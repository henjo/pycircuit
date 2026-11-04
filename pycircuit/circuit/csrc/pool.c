#define _GNU_SOURCE
#include <stdint.h>
#include <fenv.h>
#include <dlfcn.h>
#include <pthread.h>
#include <stdatomic.h>
#include <unistd.h>
#include "pcx_abi.h"
#include "pool.h"
/* THE POOL: a job's items on `nthreads` threads -- the caller and workers
   kept from call to call (spinning a while for the next call, then asleep),
   each taking chunks of the items from one counter.  The items are the
   caller's (the core's fused-group elements, the walk's limiters): each a
   function of its own inputs, writing its own outputs, so any thread may
   make it.  The workers take the caller's FP environment each call and the
   flags they raise are raised in the caller's after (the union, as one
   thread leaves them).  One call at a time (a concurrent caller makes its
   items alone: 0); a forked child starts workers of its own; the library
   pins itself (RTLD_NODELETE) before a worker runs in it. */
/* THE CHUNK: a thread takes POOL_CHUNK items at a time, POOL_CHUNK_BIG where
   every thread gets two chunks that size -- a chunk's cost is its take from
   the one counter and the output rows its neighbours write beside its own
   (measured on the gold decks, four threads: 1 or 2 items a take slower
   everywhere, jge under Radau 9.2 -> 10.7 s at 1; 8 jge 1.03x, the rest the
   same; 16 gm19 and voter25 1-2 % slower; a longer spin no change); a job of
   fewer items than a chunk one at a time (radau_tc.c's two factors) */
enum { POOL_MAX = 31, POOL_CHUNK = 4, POOL_CHUNK_BIG = 8, POOL_SPIN = 1 << 16 };

static pthread_mutex_t pool_mu = PTHREAD_MUTEX_INITIALIZER;
static pthread_cond_t pool_wake = PTHREAD_COND_INITIALIZER;
static long pool_started;
static pid_t pool_pid;
static atomic_long pool_gen, pool_next, pool_left;
static atomic_int pool_flags;
static atomic_flag pool_busy = ATOMIC_FLAG_INIT;
static pool_job_t pool_job;
static fenv_t pool_env;

static void pool_relax(void)
{
#if defined(__x86_64__) || defined(__i386__)
    __builtin_ia32_pause();
#endif
}

/* the job's items, a chunk at a time from the shared counter, until none
   is left */
static void pool_take(const pool_job_t *j)
{
    pool_range_t range = (pool_range_t) j->range;
    long i, end;
    for (;;) {
        i = atomic_fetch_add_explicit(&pool_next, j->chunk, memory_order_relaxed);
        if (i >= j->n) return;
        end = i + j->chunk < j->n ? i + j->chunk : j->n;
        range(j->ctx, i, end);
    }
}

/* a worker's wait for the call after `seen`: spinning, then asleep */
static long pool_wait(long seen)
{
    long spin;
    for (spin = 0; spin < POOL_SPIN; spin++) {
        if (atomic_load_explicit(&pool_gen, memory_order_acquire) != seen) goto out;
        pool_relax();
    }
    pthread_mutex_lock(&pool_mu);
    while (atomic_load_explicit(&pool_gen, memory_order_acquire) == seen)
        pthread_cond_wait(&pool_wake, &pool_mu);
    pthread_mutex_unlock(&pool_mu);
out:
    return atomic_load_explicit(&pool_gen, memory_order_acquire);
}

/* a worker: each call (the caller waits for every worker before the
   next, so none is missed), in the caller's FP environment, its flags
   handed back */
static void *pool_worker(void *arg)
{
    long seen = (long) (intptr_t) arg;
    int fl;
    for (;;) {
        seen = pool_wait(seen);
        fesetenv(&pool_env);
        feclearexcept(FE_ALL_EXCEPT);
        pool_take(&pool_job);
        fl = fetestexcept(FE_ALL_EXCEPT);
        if (fl) atomic_fetch_or_explicit(&pool_flags, fl, memory_order_relaxed);
        atomic_fetch_sub_explicit(&pool_left, 1, memory_order_release);
    }
    return NULL;
}

/* this library kept loaded while workers may run in it */
static int pool_pin(void)
{
    Dl_info info;
    if (!dladdr((void *) pool_worker, &info) || !info.dli_fname) return 0;
    return dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD | RTLD_NODELETE) != NULL;
}

/* `n - 1` workers (at most POOL_MAX, named "pcx-core") running in this
   process: 1 where any is.  In a forked child the parent's are gone: its
   own, from none */
static int pool_ready(long n)
{
    long want = n - 1 < POOL_MAX ? n - 1 : POOL_MAX;
    pthread_attr_t at;
    pthread_t t;
    if (pool_pid != getpid()) {
        pool_pid = getpid();
        pool_started = 0;
        pthread_mutex_init(&pool_mu, NULL);
        pthread_cond_init(&pool_wake, NULL);
        if (!pool_pin()) pool_pid = 0;
    }
    if (!pool_pid) return 0;
    while (pool_started < want) {
        pthread_attr_init(&at);
        pthread_attr_setdetachstate(&at, PTHREAD_CREATE_DETACHED);
        if (pthread_create(&t, &at, pool_worker,
                           (void *) (intptr_t) atomic_load(&pool_gen)) != 0) {
            pthread_attr_destroy(&at);
            break;
        }
        pthread_attr_destroy(&at);
        pthread_setname_np(t, "pcx-core");
        pool_started++;
    }
    return pool_started > 0;
}

/* one call on the pool: the job published, the caller's share taken, every
   worker's waited for, their flags raised here */
static void pool_run(const pool_job_t *j)
{
    int fl;
    fegetenv(&pool_env);
    atomic_store_explicit(&pool_next, 0, memory_order_relaxed);
    atomic_store_explicit(&pool_flags, 0, memory_order_relaxed);
    atomic_store_explicit(&pool_left, pool_started, memory_order_relaxed);
    pthread_mutex_lock(&pool_mu);
    atomic_fetch_add_explicit(&pool_gen, 1, memory_order_release);
    pthread_cond_broadcast(&pool_wake);
    pthread_mutex_unlock(&pool_mu);
    pool_take(j);
    while (atomic_load_explicit(&pool_left, memory_order_acquire) > 0) pool_relax();
    fl = atomic_load_explicit(&pool_flags, memory_order_relaxed);
    if (fl) feraiseexcept(fl);
}

/* THE ENTRY: `range` (pool_range_t) over the items [0, n) of `ctx` on the
   pool -- 1 -- or 0 where the caller makes them alone (one thread asked
   for, the pool another caller's or not running) */
long pcx_pool(void *range, const void *ctx, long n, long nthreads)
{
    if (nthreads < 2 || n < 1) return 0;
    if (atomic_flag_test_and_set_explicit(&pool_busy, memory_order_acquire)) return 0;
    if (!pool_ready(nthreads)) {
        atomic_flag_clear_explicit(&pool_busy, memory_order_release);
        return 0;
    }
    pool_job.range = range;
    pool_job.ctx = ctx;
    pool_job.n = n;
    pool_job.chunk = n < POOL_CHUNK ? 1
                     : n >= 2 * POOL_CHUNK_BIG * (pool_started + 1) ? POOL_CHUNK_BIG : POOL_CHUNK;
    pool_run(&pool_job);
    atomic_flag_clear_explicit(&pool_busy, memory_order_release);
    return 1;
}

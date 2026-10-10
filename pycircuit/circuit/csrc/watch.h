/* watch.h -- the dict watcher's entry: the interface `_watch` loads
 * watch.c through; watch.c includes it, so the compiler checks the
 * entry against it. */
void *hdl_fn(int which);

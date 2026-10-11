
#include <stdint.h>
#include "watch.h"
static int64_t hdl_epoch = 0;
static int hdl_dict_cb(int event, void *d, void *k, void *v)
{ (void) event; (void) d; (void) k; (void) v; hdl_epoch++; return 0; }
void *pcx_watch(int which)
{
    return which == 0 ? (void *) &hdl_epoch : (void *) hdl_dict_cb;
}

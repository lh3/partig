#include "ptpriv.h"
#include "ksort.h"

#define pt128x_key(z) ((z).x)
KRADIX_SORT_INIT(pt128x, pt128_t, pt128x_key, 8)

uint32_t pt_verbose = 3;

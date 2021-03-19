#ifndef PTPRIV_H
#define PTPRIV_H

#include "partig.h"

#define PT_MALLOC(ptr, len) ((ptr) = (__typeof__(ptr))malloc((len) * sizeof(*(ptr))))
#define PT_CALLOC(ptr, len) ((ptr) = (__typeof__(ptr))calloc((len), sizeof(*(ptr))))
#define PT_REALLOC(ptr, len) ((ptr) = (__typeof__(ptr))realloc((ptr), (len) * sizeof(*(ptr))))
#define PT_EXPAND(a, m) do { \
		(m) = (m)? (m) + ((m)>>1) : 16; \
		PT_REALLOC((a), (m)); \
	} while (0)

#ifndef kroundup64
#define kroundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))
#endif

void pt_sketch(const char *str, int len, int w, int k, uint32_t rid, int is_hpc, pt_mz1_v *p); // in sketch.c
void radix_sort_pt128x(pt128_t *st, pt128_t *en);

#endif

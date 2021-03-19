#ifndef PARTIG_H
#define PARTIG_H

#include <stdint.h>
#include "gfa.h"

#define PT_VERSION "0.0.0-dirty"

typedef struct {
	uint64_t x, y;
} pt128_t;

typedef struct {
	uint64_t x;
	uint64_t rid:28, pos:27, rev:1, span:8;
} pt_mz1_t;

typedef struct { uint32_t n, m; pt128_t *a; } pt128_v;
typedef struct { uint32_t n, m; pt_mz1_t *a; } pt_mz1_v;

typedef struct {
	int32_t k, w, is_hpc;
	int32_t max_occ;
} pt_pdopt_t;

extern uint32_t pt_verbose;

void pt_opt_init(pt_pdopt_t *opt);
void pt_pdist(const pt_pdopt_t *opt, const gfa_t *g);

double pt_cputime(void);
double pt_realtime(void);
long pt_peakrss(void);

#endif

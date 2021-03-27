#ifndef PARTIG_H
#define PARTIG_H

#include <stdint.h>
#include <stdio.h>
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
	int32_t min_cnt;
	double min_sim, diff_thres;
} pt_pdopt_t;

typedef struct {
	uint32_t sid[2], n[2], uni[2], m, rev;
	double sim;
} pt_match1_t;

typedef struct {
	uint32_t *cnt, *ucnt, n_ma, n_seg;
	uint64_t *idx, *group;
	int8_t *s; // phase
	pt_match1_t *ma;
} pt_match_t;

extern uint32_t pt_verbose;

void pt_opt_init(pt_pdopt_t *opt);
pt_match_t *pt_pdist(const pt_pdopt_t *opt, const gfa_t *g);
void pt_match_print(FILE *fp, const gfa_t *g, const pt_match_t *ma);
void pt_match_free(pt_match_t *ma);
void pt_cc(pt_match_t *ma);
void pt_solve(pt_match_t *ma, uint64_t x);

double pt_cputime(void);
double pt_realtime(void);
long pt_peakrss(void);

#endif

#include "ptpriv.h"
#include "ksort.h"

#define mz_key(z) ((z).x)
KRADIX_SORT_INIT(mz, pt_mz1_t, mz_key, 8)

void pt_opt_init(pt_pdopt_t *opt)
{
	memset(opt, 0, sizeof(pt_pdopt_t));
	opt->k = 51, opt->w = 51, opt->is_hpc = 1;
	opt->max_occ = 20;
}

void pt_pdist(const pt_pdopt_t *opt, const gfa_t *g)
{
	int32_t i;
	uint32_t st, j;
	uint64_t m_a = 0, n_a = 0;
	pt_mz1_v mz;
	pt128_t *a = 0;

	// collect minimizers
	memset(&mz, 0, sizeof(mz));
	for (i = 0; i < g->n_seg; ++i)
		pt_sketch(g->seg[i].seq, g->seg[i].len, opt->w, opt->k, i, opt->is_hpc, &mz);
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s] %d minimizers\n", __func__, mz.n);
	radix_sort_mz(mz.a, mz.a + mz.n);

	// collect anchors
	for (j = 1, st = 0; j <= mz.n; ++j) {
		if (j == mz.n || mz.a[j].x != mz.a[st].x) {
			uint32_t k, l;
			if (j - st >= 2 && j - st <= opt->max_occ) goto end_anchor;
			for (k = st; k < j; ++k) {
				for (l = k + 1; l < j; ++l) {
					uint32_t span, rev = (mz.a[k].rev != mz.a[l].rev);
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz.a[l].span;
					a[n_a].x = (uint64_t)mz.a[k].rid << 33 | mz.a[l].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz.a[k].pos << 32 | (mz.a[l].pos + 1 - span);
					if (n_a == m_a) PT_EXPAND(a, m_a);
					span = mz.a[k].span;
					a[n_a].x = (uint64_t)mz.a[l].rid << 33 | mz.a[k].rid << 1 | rev;
					a[n_a++].y = (uint64_t)mz.a[l].rid << 33 | (mz.a[k].pos + 1 - span);
				}
			}
end_anchor:	st = j;
		}
	}
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s] %ld anchors\n", __func__, (long)n_a);
	radix_sort_pt128x(a, a + n_a);

	free(mz.a);
}

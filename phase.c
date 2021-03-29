#include <zlib.h>
#include "ptpriv.h"
#include "ksort.h"
#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

#define el_key(z) ((z).x)
KRADIX_SORT_INIT(el, pt_node_t, el_key, 8)

pt128_t *pt_read_links(const gfa_t *g, const char *fn, uint32_t *n_link_)
{
	uint32_t n_link = 0, m_link = 0;
	pt128_t *link = 0;
	gzFile fp;
	kstream_t *ks;
	kstring_t str = {0,0,0};
	int dret;

	*n_link_ = 0;
	fp = gzopen(fn, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	while (ks_getuntil(ks, KS_SEP_LINE, &str, &dret) >= 0) {
		char *p, *q, *u1 = 0, *u2 = 0;
		int32_t i;
		double x = -1.0;
		for (i = 0, p = q = str.s;; ++p) {
			if (*p == '\t' || *p == 0) {
				int c = *p;
				*p = 0;
				if (i == 0) u1 = q;
				else if (i == 1) u2 = q;
				else if (i == 2) x = atof(q);
				++i, q = p + 1;
				if (c == 0) break;
			}
		}
		if (i >= 3 && x > 0.0) {
			uint64_t y = (uint64_t)(x + .499);
			int32_t sid1, sid2;
			if (y == 0) continue;
			sid1 = gfa_name2id(g, u1);
			sid2 = gfa_name2id(g, u2);
			if (sid1 < 0 || sid2 < 0) continue;
			if (n_link == m_link) PT_EXPAND(link, m_link);
			link[n_link].x = (uint64_t)sid1 << 32 | sid2;
			link[n_link++].y = y;
		}
	}
	ks_destroy(ks);
	gzclose(fp);
	*n_link_ = n_link;
	if (pt_verbose >= 3)
		fprintf(stderr, "[%s::%.3f] read %d links\n", __func__, pt_realtime(), n_link);
	return link;
}

pt_node_t *pt_merge_list(uint32_t n_ma, const pt_match1_t *ma, uint32_t n_link, const pt128_t *link, uint32_t *n_node_)
{
	uint32_t i, k, n_node, st;
	pt_node_t *node;
	PT_CALLOC(node, n_ma + n_link);
	for (i = k = 0; i < n_ma; ++i) {
		pt_node_t *e = &node[k++];
		e->x = (uint64_t)ma[i].sid[0] << 32 | ma[i].sid[1];
		e->m = ma[i].m;
	}
	for (i = 0; i < n_link; ++i) {
		pt_node_t *e = &node[k++];
		e->x = link[i].x;
		e->l = link[i].y;
	}
	n_node = k;
	radix_sort_el(node, node + n_node);
	for (st = 0, i = 1, k = 0; i <= n_node; ++i) {
		if (i == n_node || node[i].x != node[st].x) {
			assert(i - st == 1 || i - st == 2);
			if (i - st == 2) {
				assert(node[st].m * node[st+1].m == 0);
				assert(node[st].l * node[st+1].l == 0);
				node[k].m = node[st].m + node[st+1].m;
				node[k].l = node[st].l + node[st+1].l;
				node[k++].x = node[st].x;
			} else node[k++] = node[i];
			st = i;
		}
	}
	*n_node_ = n_node = k;
	return node;
}

void pt_node_print(FILE *fp, const gfa_t *g, uint32_t n_node, const pt_node_t *node)
{
	uint32_t i;
	for (i = 0; i < n_node; ++i) {
		const pt_node_t *e = &node[i];
		fprintf(fp, "E\t%s\t%s\t%ld\t%ld\n", g->seg[e->x>>32].name, g->seg[(uint32_t)e->x].name,
				(long)e->m, (long)e->l);
	}
}

void pt_phase(const gfa_t *g, const pt_match_t *ma, const char *fn_link)
{
	uint32_t n_link, n_node;
	pt128_t *link;
	pt_node_t *node;

	link = pt_read_links(g, fn_link, &n_link);
	node = pt_merge_list(ma->n_ma, ma->ma, n_link, link, &n_node);
	pt_node_print(stdout, g, n_node, node);
}

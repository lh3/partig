#include <zlib.h>
#include "ptpriv.h"
#include "kseq.h"

KSTREAM_INIT(gzFile, gzread, 0x10000)

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
	while (ks_getuntil2(ks, KS_SEP_LINE, &str, &dret, 0) >= 0) {
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
	return link;
}

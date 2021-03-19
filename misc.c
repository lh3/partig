#include <sys/resource.h>
#include <sys/time.h>
#include "ptpriv.h"
#include "ksort.h"

#define pt128x_key(z) ((z).x)
KRADIX_SORT_INIT(pt128x, pt128_t, pt128x_key, 8)

double pt_realtime0 = -1.0;
uint32_t pt_verbose = 3;

double pt_cputime(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

long pt_peakrss(void)
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
#ifdef __linux__
	return r.ru_maxrss * 1024;
#else
	return r.ru_maxrss;
#endif
}

static double pt_realtime_core(void)
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

double pt_realtime(void)
{
	if (pt_realtime0 < 0) pt_realtime0 = pt_realtime_core();
	return pt_realtime_core() - pt_realtime0;
}

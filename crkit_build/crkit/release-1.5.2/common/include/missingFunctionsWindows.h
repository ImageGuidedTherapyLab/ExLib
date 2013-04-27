#include <windows.h>
#include <winsock.h> // timeval lives here
#include <time.h>

void gettimeofday(struct timeval* p, void* tz /* IGNORED */)
{
	union {
		long long ns100; /*time since 1 Jan 1601 in 100ns units */
		FILETIME ft;
	} now;

	GetSystemTimeAsFileTime( &(now.ft) );
	p->tv_usec=(long)((now.ns100 / 10LL) % 1000000LL );
	p->tv_sec= (long)((now.ns100-(116444736000000000LL))/10000000LL);
}


double fmax(double value1, double value2)
{
	if (value1 > value2)
		return value1;
	else
		return value2;
}

double fmin(double value1, double value2)
{
	if (value1 < value2)
		return value1;
	else
		return value2;
}

//(C) Microsoft Corporation. All rights reserved.
#include "stdafx.h"
#include "Utils.h"

namespace paraF
{

LongType abs(LongType value)
{
	return value > 0 ? value : -value;
}

LongType gcd(LongType a, LongType b)
{
	while (a != 0 && b != 0)
	{
		if (a > b)
			a %= b;
		else
			b %= a;
	}
	return a + b;
}

int GetTime()
{
#ifdef __GNUC__
	rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec * 1000 + r.ru_utime.tv_usec / 1000;
#else
	FILETIME sysTime;
	FILETIME userTime;
	FILETIME dummyTime;
	GetProcessTimes(GetCurrentProcess(), &dummyTime, &dummyTime, &sysTime, &userTime);
	return (sysTime.dwLowDateTime + userTime.dwLowDateTime) / 10000;
#endif
}

}

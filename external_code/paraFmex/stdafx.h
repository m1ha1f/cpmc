//(C) Microsoft corporation. All rights reserved
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#ifndef stdafx_h_SEEN
#define stdafx_h_SEEN

#define _CRT_SECURE_NO_DEPRECATE
#define _HAS_ITERATOR_DEBUGGING 0

#ifdef __GNUC__
	typedef long long LongType;
	#define	LONGTYPE_SPEC "lld"
#else
	typedef __int64 LongType;
	#define	LONGTYPE_SPEC "I64d"
#endif

#ifdef __GNUC__
	#include <sys/time.h>
	#include <sys/resource.h>
#else
	#define WIN32_LEAN_AND_MEAN
	#define NOMINMAX
	#include <windows.h>
#endif

#include <stdio.h>
#include <algorithm>
#include <vector>
#include <list>
#include <limits>
#include <queue>
#include <map>
#include <set>
#include <stack>
#include <cmath>


#endif

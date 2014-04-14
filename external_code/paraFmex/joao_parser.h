//(C) Microsoft corporation. All rights reserved
#ifndef Parser_h_SEEN
#define Parser_h_SEEN

#include "Graph.h"

namespace paraF
{

struct Problem
{
	enum ProblemType
	{
		MAX_FLOW_PROBLEM,
		PARAMETRIC_FLOW_PROBLEM
	};

	ProblemType type;

	bool hasLambdaRange;
	int minLambda;
	int maxLambda;
};

class InputParser
{
public:
	void Parse(Graph& graph, Problem& problem, double* arg1, LongType s, LongType t, LongType m, LongType n);
};

}
#endif


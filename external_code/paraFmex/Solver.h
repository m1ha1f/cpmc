//(C) Microsoft corporation. All rights reserved
#ifndef Solver_h_SEEN
#define Solver_h_SEEN

#include "Graph.h"
#include "joao_parser.h"

namespace paraF
{

class Solver
{
public:
	Solver();

	void Solve(Graph& graph, Problem& problem);

	enum SolverAlgorithm
	{
		SIMPLE_ALGORITHM,
		GGT_ST_ONLY_ALGORITHM,
		GGT_TS_ONLY_ALGORITHM,
		GGT_BOTH_ALGORITHM,
		GGT_ALGORITHM
	};
	
	enum Verbosity
	{
		DEBUG_VERBOSITY,
		NORMAL_VERBOSITY,
		QUIET_VERBOSITY
	};

	SolverAlgorithm mode;
	Verbosity verbosity;

	bool printCuts;
	bool maxBreakpointOnly;
  int numBreakpoints;
	LongType maxBreakpoint;
	std::set<int> breakpointCutPointers;
  double lambdas[3000]; // a safe maximum number of breakpoints (should still insert some error verification code)
	int multiplier;
private:
	void Slice(Graph* stGraph, LongType stLambda, Graph* tsGraph, LongType tsLambda);
	void CopyF2(Graph* goodGraph, Graph* badGraph, MinCut& minCut);
	void CopyF2(Node* node, Graph* goodGraph, Graph* badGraph, MinCut& minCut);
	void RestoreF13(Graph& graph, MinCut& minCut, Node::MinCutType type);
	void RestoreF13(Node* node, Graph& graph, MinCut& minCut, Node::MinCutType type);
	void PrintBreakpoint(LongType nom, LongType denom);

	void SimpleSolve(Graph& graph, Problem& problem);
	void SimpleSlice(Graph* graph, LongType stLambda, LongType tsLambda);

	void MaxBreakpointSolve(Graph& graph, Problem& problem);

	MinCut minCut;
	int totalNodeCounter;

#ifdef PARANOID
	MinCut backupCut;
#endif

};

}
#endif

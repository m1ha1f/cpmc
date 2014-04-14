//(C) Microsoft Corporation. All rights reserved.
#include "stdafx.h"
#include "Solver.h"
#include "Utils.h"
#include "mex.h"

namespace paraF
{

const int PARALLEL_WORK_AMOUNT = 1000;

Solver::Solver()
{
	multiplier = 1;
	mode = GGT_ALGORITHM;
	verbosity = NORMAL_VERBOSITY;
  //verbosity = DEBUG_VERBOSITY;
	printCuts = false;
	maxBreakpointOnly = false;
	maxBreakpoint = std::numeric_limits<LongType>::min();
}

void Solver::Solve(Graph& stGraph, Problem& problem)
{
	totalNodeCounter = 0;
	
	if (maxBreakpointOnly)
	{
		MaxBreakpointSolve(stGraph, problem);
		return;
	}

	if (mode == SIMPLE_ALGORITHM)
	{
		SimpleSolve(stGraph, problem);
		return;
	}

	Graph tsGraph;
	stGraph.Reverse(tsGraph);

	LongType stLambda;
	LongType tsLambda;

	//mexPrintf("Problem has lambda range: %d", problem.hasLambdaRange);
	if (problem.hasLambdaRange)
	{    
		stLambda = problem.minLambda;
		tsLambda = problem.maxLambda;
	}
	else
	{
		stLambda = stGraph.GetLowerLambdaBound();
		tsLambda = -tsGraph.GetLowerLambdaBound();
	}

	stGraph.InitCapDelta(stLambda);
	stGraph.InitCapDelta(tsLambda);

	tsGraph.InitCapDelta(-stLambda);
	tsGraph.InitCapDelta(-tsLambda);

	multiplier = stGraph.ChooseMultiplier(stLambda, tsLambda);

	//mexPrintf("m %d\n", multiplier);

	stGraph.SetMultiplier(multiplier);
	tsGraph.SetMultiplier(multiplier);

	stLambda *= multiplier;
	tsLambda *= multiplier;

	int time0 = GetTime();

	MinCut stCut;
	MinCut tsCut;

	//mexPrintf("c Starting initial s-t maxflow computation\n");
	stGraph.InitPreflowPush(stLambda);
	stGraph.FindMaxPreflow();
	stGraph.ConvertPreflowToFlow();
#ifdef PARANOID
	stGraph.CheckMaxFlow();
#endif
	stGraph.FindMinCut(stCut);

	//mexPrintf("c Starting initial t-s maxflow computation\n");
	tsGraph.InitPreflowPush(-tsLambda);
	tsGraph.FindMaxPreflow();
	tsGraph.ConvertPreflowToFlow();
#ifdef PARANOID
	tsGraph.CheckMaxFlow();
#endif
	tsGraph.FindMinCut(tsCut);

	int time1 = GetTime();

	//mexPrintf("c Initial solve time: %d\n", time1 - time0);

	for (int i = 0; i < (int) stGraph.nodes.size(); i++)
	{
		Node* stNode = &stGraph.nodes[i];
		Node* tsNode = &tsGraph.nodes[i];
		if (!stGraph.IsInnerNode(stNode))
			continue;

		if (stCut[i] == Node::SOURCE_NODE)
		{
			stGraph.MakeSource(stNode);
			tsGraph.MakeSink(tsNode);
		}
		if (tsCut[i] == Node::SOURCE_NODE)
		{
			stGraph.MakeSink(stNode);
			tsGraph.MakeSource(tsNode);
		}
	}
	
	breakpointCutPointers.clear();
	breakpointCutPointers.insert(stGraph.state.firstSinkIndex);
	breakpointCutPointers.insert((int) stGraph.nodes.size());
	numBreakpoints = 0;

	Slice(&stGraph, stLambda, &tsGraph, tsLambda);

	if (numBreakpoints + 2 != breakpointCutPointers.size())
		throw "wrong number of cuts found";

	//mexPrintf("c Breakpoints: %d\n", numBreakpoints);

	if (printCuts)
	{
		int color = 1;
		int j = 0;
		for (std::set<int>::iterator i = breakpointCutPointers.begin();
			i != breakpointCutPointers.end(); i++)
		{
			while (j < *i)
			{
				mexPrintf("z %d %d\n", stGraph.typedNodeList[j]->index + 1, color);
				j++;
			}
			color++;
		}
	}

	//mexPrintf("c\n");
	//mexPrintf("c s-t graph statistics:\n");
	//stGraph.PrintStats("st");

	//mexPrintf("c\n");
	//mexPrintf("c t-s graph statistics:\n");
	//tsGraph.PrintStats("ts");

	//mexPrintf("c\n");
	//stGraph.AccumulateStats(&tsGraph);
	//stGraph.PrintStats("total");

	//mexPrintf("c\n");
	//mexPrintf("c Total nodes: %d\n", totalNodeCounter);
	//mexPrintf("c\n");
}

void Solver::SimpleSolve(Graph& graph, Problem& problem)
{
	LongType stLambda;
	LongType tsLambda;

	if (problem.hasLambdaRange)
	{
		stLambda = problem.minLambda;
		tsLambda = problem.maxLambda;
	}
	else
	{
		stLambda = graph.GetLowerLambdaBound();
		tsLambda = graph.GetUpperLambdaBound();
	}

	graph.InitCapDelta(stLambda);
	graph.InitCapDelta(tsLambda);

	multiplier = graph.ChooseMultiplier(stLambda, tsLambda);

	//mexPrintf("m %d\n", multiplier);

	graph.SetMultiplier(multiplier);

	stLambda *= multiplier;
	tsLambda *= multiplier;

	int time0 = GetTime();

	MinCut stCut;
	MinCut tsCut;

	//mexPrintf("c Starting first initial s-t maxflow computation\n");
	graph.InitPreflowPush(stLambda);
	graph.FindMaxPreflow();
	graph.ConvertPreflowToFlow();
	graph.FindMinCut(stCut);

	//mexPrintf("c Starting second initial s-t maxflow computation\n");
	graph.IncreaseLambda(tsLambda);
	graph.FindMaxPreflow();
	graph.ConvertPreflowToFlow();
	graph.FindMinCut(tsCut);

	int time1 = GetTime();

	//mexPrintf("c Initial solve time: %d\n", time1 - time0);

	for (int i = 0; i < (int) graph.nodes.size(); i++)
	{
		Node* node = &graph.nodes[i];
		if (!graph.IsInnerNode(node))
			continue;

		if (stCut[i] == Node::SOURCE_NODE)
		{
			graph.MakeSource(node);
		}
		if (tsCut[i] == Node::SINK_NODE)
		{
			graph.MakeSink(node);
		}
	}

	numBreakpoints = 0;
	
	SimpleSlice(&graph, stLambda, tsLambda);

	//mexPrintf("c Breakpoints: %d\n", numBreakpoints);

	//mexPrintf("c\n");
	//mexPrintf("c Graph statistics:\n");
	//graph.PrintStats("total");

	//mexPrintf("c\n");
	//mexPrintf("c Total nodes: %d\n", totalNodeCounter);
	//mexPrintf("c\n");
}

void Solver::PrintBreakpoint(LongType nom, LongType denom)
{
	maxBreakpoint = std::max(maxBreakpoint, nom / denom);
	LongType g = gcd(abs(nom), multiplier);
	nom /= g;
	denom *= (multiplier / g);
	double lambda = (double) nom / denom;
  
  lambdas[numBreakpoints] = lambda;
  
	/*if (verbosity <= DEBUG_VERBOSITY)
	{
		mexPrintf("c Found breakpoint %" LONGTYPE_SPEC "/%" LONGTYPE_SPEC " ~ %.6lf\n", nom, denom, lambda);
	}
	if (verbosity <= NORMAL_VERBOSITY)
	{
		mexPrintf("l %" LONGTYPE_SPEC " %" LONGTYPE_SPEC "\n", nom, denom);
	}*/
	numBreakpoints++;
}

void Solver::Slice(Graph* stGraph, LongType stLambda, Graph* tsGraph, LongType tsLambda)
{
	if (stGraph->GetNodeCount() == 2)
		return;

	totalNodeCounter += stGraph->GetNodeCount();
	
	LongType sourceCap = stGraph->state.sourceCap;
	LongType sourceSlope = stGraph->state.sourceSlope;

	LongType sinkCap = tsGraph->state.sourceCap;
	LongType sinkSlope = -tsGraph->state.sourceSlope;

	LongType denom = sourceSlope - sinkSlope;
	LongType nom = sinkCap - sourceCap;
	LongType g = gcd(abs(nom), abs(denom));
	denom /= g;
	nom /= g;
	LongType midLambda = nom / denom;

	if (verbosity <= DEBUG_VERBOSITY)
	{
		mexPrintf("c Slice [%.6lf, %.6lf, %.6lf]\n", (double) stLambda/multiplier, (double) midLambda/multiplier, (double) tsLambda/multiplier);
	}
	
	if (midLambda == stLambda || midLambda == tsLambda)
	{
		PrintBreakpoint(nom, denom);
		breakpointCutPointers.insert(stGraph->state.firstInnerIndex);
		return;
	}

	stGraph->InitArcLists();
	tsGraph->InitArcLists();

	stGraph->BackupFlows();
	tsGraph->BackupFlows();

 	stGraph->IncreaseLambda(midLambda);
	tsGraph->IncreaseLambda(-midLambda);

#ifdef PARANOID
	if (mode != GGT_TS_ONLY_ALGORITHM)
		stGraph->CheckPreflow();
	if (mode != GGT_ST_ONLY_ALGORITHM)
		tsGraph->CheckPreflow();

	if (stGraph->state.firstInnerIndex != (int) tsGraph->nodes.size() - tsGraph->state.firstSinkIndex)
		throw "nodes are out of sync";
	if ((int) stGraph->nodes.size() - stGraph->state.firstSinkIndex != tsGraph->state.firstInnerIndex)
		throw "nodes are out of sync";

	for (int i = stGraph->state.firstInnerIndex; i < stGraph->state.firstSinkIndex; i++)
	{
		Node* stNode = stGraph->typedNodeList[i];
		Node* tsNode = &tsGraph->nodes[stNode->index];

		if (stNode->svArc->cap != tsNode->vtArc->cap)
			throw "capacities are out of sync";

		if (stNode->vtArc->cap != tsNode->svArc->cap)
			throw "capacities are out of sync";

		for (Arc* stArc = stNode->firstOutArc; stArc != NULL; stArc = stArc->nextOutArc)
		{
			if (stGraph->IsInnerNode(stArc->headNode))
			{
				Arc* tsArc = tsGraph->arcs[stArc - &*stGraph->arcs.begin()].revArc;
				if (stArc->GetHeadNode()->index != tsArc->GetTailNode()->index ||
					stArc->GetTailNode()->index != tsArc->GetHeadNode()->index)
					throw "arcs are not of sync";
				if (stArc->cap != tsArc->cap)
					throw "capacities are out of sync";
			}
		}
	}
#endif

	Graph* goodGraph;
	Graph* badGraph;
	Node::MinCutType goodSourceType;

	for (;;)
	{
		if (mode != GGT_TS_ONLY_ALGORITHM && stGraph->FindMaxPreflow(PARALLEL_WORK_AMOUNT))
		{
			goodGraph = stGraph;
			badGraph = tsGraph;
			goodSourceType = Node::SOURCE_NODE;
			break;
		}
		if (mode != GGT_ST_ONLY_ALGORITHM && tsGraph->FindMaxPreflow(PARALLEL_WORK_AMOUNT))
		{
			goodGraph = tsGraph;
			badGraph = stGraph;
			goodSourceType = Node::SINK_NODE;
			break;
		}
	}

	goodGraph->ConvertPreflowToFlow();
#ifdef PARANOID
	goodGraph->CheckMaxFlow();
#endif
	goodGraph->FindMinCut(minCut);

#ifdef PARANOID
	backupCut.resize(stGraph->nodes.size());
	std::copy(minCut.begin(), minCut.end(), backupCut.begin());
#endif

	int sourceSideSize = 1;
	for (int i = goodGraph->state.firstInnerIndex; i < goodGraph->state.firstSinkIndex; i++)
	{
		Node* node = goodGraph->typedNodeList[i];
		if (minCut[node->index] == goodSourceType)
			sourceSideSize++;
	}

	if (verbosity <= DEBUG_VERBOSITY)
	{
		mexPrintf("c Cut size: %d of %d\n", sourceSideSize, goodGraph->GetNodeCount());		
	}
	
	bool flag = sourceSideSize > (int) goodGraph->GetNodeCount() / 2;
	
	if (mode == GGT_BOTH_ALGORITHM)
		flag = true;
	if (mode == GGT_ST_ONLY_ALGORITHM || mode == GGT_TS_ONLY_ALGORITHM)
		flag = false;
		
	if (sourceSideSize == goodGraph->GetNodeCount() - 1 ||
		sourceSideSize == 1)
	{
		if (verbosity <= DEBUG_VERBOSITY)
		{
			mexPrintf("c Cut is trivial, no amortization is needed\n");
		}
		flag = false;
	}

	if (flag)
	{
		if (verbosity <= DEBUG_VERBOSITY)
		{
			mexPrintf("c Computing reverse preflow\n");
		}
		badGraph->FindMaxPreflow();
		badGraph->ConvertPreflowToFlow();
#ifdef PARANOID
		badGraph->CheckMaxFlow();
#endif
		badGraph->FindMinCut(minCut);

#ifdef PARANOID
		for (int i = stGraph->state.firstInnerIndex; i < stGraph->state.firstSinkIndex; i++)
		{
			int j = stGraph->typedNodeList[i]->index;
			if (minCut[j] + backupCut[j] != 0)
				throw "cut check failed";
		}
#endif

		goodSourceType = (Node::MinCutType) -goodSourceType;
		std::swap(goodGraph, badGraph);
	}

	CopyF2(goodGraph, badGraph, minCut);

	RestoreF13(*stGraph, minCut, goodSourceType);
	RestoreF13(*tsGraph, minCut, (Node::MinCutType) -goodSourceType);

	int stFirstInnerIndex, stFirstSinkIndex;
	stGraph->RearrangeNodes(minCut, goodSourceType, &stFirstInnerIndex, &stFirstSinkIndex);

	int tsFirstInnerIndex, tsFirstSinkIndex;
	tsGraph->RearrangeNodes(minCut, (Node::MinCutType) -goodSourceType, &tsFirstInnerIndex, &tsFirstSinkIndex);

	GraphState stState, tsState;

	stGraph->SaveState(stState);
	tsGraph->SaveState(tsState);

	while (stGraph->state.firstSinkIndex > stFirstInnerIndex)
	{
		Node* node = stGraph->typedNodeList[stGraph->state.firstSinkIndex - 1];
		stGraph->MakeSink(node);
	}

	while (tsGraph->state.firstInnerIndex < tsFirstSinkIndex)
	{
		Node* node = tsGraph->typedNodeList[tsGraph->state.firstInnerIndex];
		tsGraph->MakeSource(node);
	}

	LongType sSlope = -tsGraph->state.sourceSlope;

	Slice(stGraph, stLambda, tsGraph, midLambda);

	stGraph->RestoreState(stState);
	tsGraph->RestoreState(tsState);

	while (stGraph->state.firstInnerIndex < stFirstSinkIndex)
	{
		Node* node = stGraph->typedNodeList[stGraph->state.firstInnerIndex];
		stGraph->MakeSource(node);
	}

	while (tsGraph->state.firstSinkIndex > tsFirstInnerIndex)
	{
		Node* node = tsGraph->typedNodeList[tsGraph->state.firstSinkIndex - 1];
		tsGraph->MakeSink(node);
	}

	LongType tSlope = stGraph->state.sourceSlope;

	Slice(stGraph, midLambda, tsGraph, tsLambda);

	stGraph->RestoreState(stState);
	tsGraph->RestoreState(tsState);

	if (sSlope != tSlope)
	{
		PrintBreakpoint(nom, denom);
		breakpointCutPointers.insert(stFirstInnerIndex);
	}
}

void Solver::SimpleSlice(Graph* graph, LongType stLambda, LongType tsLambda)
{
	if (graph->GetNodeCount() == 2)
		return;

	totalNodeCounter += graph->GetNodeCount();
	LongType sourceCap = graph->state.sourceCap;;
	LongType sourceSlope = graph->state.sourceSlope;

	LongType sinkCap = graph->state.sinkCap;
	LongType sinkSlope = graph->state.sinkSlope;

	LongType denom = sourceSlope - sinkSlope;
	LongType nom = sinkCap - sourceCap;
	LongType g = gcd(abs(nom), abs(denom));
	denom /= g;
	nom /= g;
	LongType midLambda = nom / denom;

	if (verbosity <= DEBUG_VERBOSITY)
	{
		mexPrintf("c Slice [%.6lf, %.6lf, %.6lf]\n", (double) stLambda/multiplier, (double) midLambda/multiplier, (double) tsLambda/multiplier);
	}
	
	if (midLambda == stLambda || midLambda == tsLambda)
	{
		PrintBreakpoint(nom, denom);
		return;
	}

 	graph->InitPreflowPush(midLambda);
	graph->FindMaxPreflow();
	graph->ConvertPreflowToFlow();
	graph->FindMinCut(minCut);

	int stFirstInnerIndex, stFirstSinkIndex;
	graph->RearrangeNodes(minCut, Node::SOURCE_NODE, &stFirstInnerIndex, &stFirstSinkIndex);

	GraphState state;

	graph->SaveState(state);

	while (graph->state.firstSinkIndex > stFirstInnerIndex)
	{
		Node* node = graph->typedNodeList[graph->state.firstSinkIndex - 1];
		graph->MakeSink(node);
	}

	LongType sSlope = graph->state.sinkSlope;

	SimpleSlice(graph, stLambda, midLambda);

	graph->RestoreState(state);

	while (graph->state.firstInnerIndex < stFirstSinkIndex)
	{
		Node* node = graph->typedNodeList[graph->state.firstInnerIndex];
		graph->MakeSource(node);
	}

	LongType tSlope = graph->state.sourceSlope;

	SimpleSlice(graph, midLambda, tsLambda);

	if (sSlope != tSlope)
	{
		PrintBreakpoint(nom, denom);
	}

	graph->RestoreState(state);
}

void Solver::CopyF2(Graph* goodGraph, Graph* badGraph, MinCut& minCut)
{
	CopyF2(goodGraph->GetSourceNode(), goodGraph, badGraph, minCut);
	for (int i = goodGraph->state.firstInnerIndex; i < goodGraph->state.firstSinkIndex; i++)
	{
		Node* node = goodGraph->typedNodeList[i];
		if (minCut[node->index] == Node::SOURCE_NODE)
			CopyF2(node, goodGraph, badGraph, minCut);
	}
}

void Solver::CopyF2(Node* node, Graph* goodGraph, Graph* badGraph, MinCut& minCut)
{
	for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
	{
		if (minCut[arc->headNode->index] == Node::SOURCE_NODE)
			badGraph->arcs[arc - &*goodGraph->arcs.begin()].revArc->flow = arc->flow;
	}
}

void Solver::RestoreF13(Graph& graph, MinCut& minCut, Node::MinCutType type)
{
	RestoreF13(graph.GetSourceNode(), graph, minCut, type);
	RestoreF13(graph.GetSinkNode(), graph, minCut, type);

	for (int i = graph.state.firstInnerIndex; i < graph.state.firstSinkIndex; i++)
	{
		Node* node = graph.typedNodeList[i];
		RestoreF13(node, graph, minCut, type);
	}
}

void Solver::RestoreF13(Node* node, Graph& graph, MinCut& minCut, Node::MinCutType type)
{
	if (minCut[node->index] == type)
	{
		for (Arc* arc = node->firstOutArc; arc != NULL; arc = arc->nextOutArc)
		{
			arc->flow = arc->baseFlow;
			arc->revArc->flow = arc->revArc->baseFlow;
		}
	}
}

void Solver::MaxBreakpointSolve(paraF::Graph &graph, paraF::Problem &problem)
{
	LongType loLambda = graph.GetLowerLambdaBound();
	LongType hiLambda = graph.GetUpperLambdaBound();

	graph.InitCapDelta(loLambda);
	graph.InitCapDelta(hiLambda);
	multiplier = graph.ChooseMultiplier(loLambda, hiLambda);

	mexPrintf("m %d\n", multiplier);

	graph.SetMultiplier(multiplier);

	loLambda *= multiplier;
	hiLambda *= multiplier;

	LongType tCap = graph.state.sinkCap;
	LongType tSlope = graph.state.sinkSlope;

	GraphState state;
	graph.SaveState(state);

	graph.InitPreflowPush(loLambda);
	graph.FindMaxPreflow();
	graph.ConvertPreflowToFlow();

	LongType nom = -1;
	LongType denom = -1;

	for (;;)
	{
		LongType sCap = graph.state.sourceCap;
		LongType sSlope = graph.state.sourceSlope;

		if (sSlope == tSlope)
			break;

		denom = sSlope - tSlope;
		nom = tCap - sCap;
		LongType g = gcd(abs(nom), abs(denom));
		denom /= g;
		nom /= g;
		LongType midLambda = nom / denom;

		if (verbosity <= DEBUG_VERBOSITY)
		{
			mexPrintf("c Trying lambda = %.6lf\n", (double) midLambda/multiplier);
		}

		if (mode == SIMPLE_ALGORITHM)
		{
			graph.InitPreflowPush(midLambda);
		}
		else
		{
			graph.InitArcLists();
 			graph.IncreaseLambda(midLambda);
		}

		graph.FindMaxPreflow();
		graph.ConvertPreflowToFlow();
		graph.FindMinCut(minCut);

		for (int i = graph.state.firstInnerIndex; i < graph.state.firstSinkIndex; i++)
		{
			Node* node = graph.typedNodeList[i];
			if (minCut[node->index] == Node::SOURCE_NODE)
			{
				graph.MakeSource(node);
			}
		}

		if (loLambda == midLambda)
			break;

		loLambda = midLambda;
	}

	if (nom != -1 && denom != -1)
		PrintBreakpoint(nom, denom);

	//mexPrintf("c\n");
	//mexPrintf("c Graph statistics:\n");
	//graph.PrintStats("total");

	graph.RestoreState(state);
}


}

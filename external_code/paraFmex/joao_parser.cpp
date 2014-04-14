//(C) Microsoft Corporation. All rights reserved.
#include "stdafx.h"

#include "joao_parser.h"
#include "Graph.h"
#include <string.h>
#include "mex.h"

namespace paraF
{

void InputParser::Parse(Graph& graph, Problem& problem, double* arg1, LongType s, LongType t, LongType m, LongType n)
{  
  if(n < s || n < t) {
    throw("invalid source or sink!");
  }
  
  problem.type = Problem::PARAMETRIC_FLOW_PROBLEM;
  
	graph.Reserve(m * 2 + n * 4, n);
	for (int i = 0; i < n; i++)
		graph.AddNode();

  int nArcsRead = 0;
	int sourceNodeIndex = s - 1;
	int sinkNodeIndex = t - 1;
  LongType cap = 0;
  int slope = 0;

  for(int i=0; i<m; i++)
	{    
    int tailNodeIndex;
    int headNodeIndex;

    tailNodeIndex = (int) arg1[i];
    headNodeIndex = (int) arg1[i + m]; // same row, 2nd column
    cap = arg1[i + 2*m];
    slope = arg1[i + 3*m];   

    if(n < tailNodeIndex || n < headNodeIndex || 0 >= tailNodeIndex || 0 >= headNodeIndex) {
      throw "Invalid node!";
    }
      
    //mexPrintf("tailNodeIndex: %d headNodeIndex: %d cap: %ld slope: %d\n", tailNodeIndex, headNodeIndex, cap,slope);    
    
    headNodeIndex--;
    tailNodeIndex--;

    Node* headNode = &graph.nodes[headNodeIndex];
    Node* tailNode = &graph.nodes[tailNodeIndex];

    Arc* arc;
    if (tailNodeIndex == sourceNodeIndex || headNodeIndex == sinkNodeIndex)
    {
      if (slope < 0)
        throw "slope should be non-negative";

      if (tailNodeIndex == sourceNodeIndex && headNodeIndex == sinkNodeIndex) {
			  /*mexPrintf("Found Source to Sink edge\n");*/
        continue;
			}

      if (tailNodeIndex == sourceNodeIndex) {
        arc = graph.AddArcFromSource(tailNode, headNode, cap, slope);
      } else {
        arc = graph.AddArcToSink(tailNode, headNode, cap, -slope);
      }
    }
    else // not source, not sink
    {
      if (headNodeIndex == sourceNodeIndex) {
			  /*mexPrintf("Found edge directed at the Source\n");*/
        continue;
			}
      if (tailNodeIndex == sinkNodeIndex) {
			  /*mexPrintf("Found edge pointing out of the Sink\n");*/
        continue;      
		  }
      if (tailNodeIndex == sourceNodeIndex)
        arc = graph.AddArcFromSource(tailNode, headNode, cap, 0);
      else if (headNodeIndex == sinkNodeIndex)
        arc = graph.AddArcToSink(tailNode, headNode, cap, 0);
      else
        arc = graph.AddArc(tailNode, headNode, cap);
    }
    
    nArcsRead++;
    
    Arc* revArc = graph.AddArc(headNode, tailNode, 0);
    arc->revArc = revArc;
    revArc->revArc = arc;    
  }

  /*if (nArcsRead < m) {
		mexPrintf("nArcsRead: %d , m: %d\n", nArcsRead, m);
    throw "too few arcs in the input";	
	}*/
  if (sourceNodeIndex < 0)
          throw "source node is not chosen";
  if (sinkNodeIndex < 0)
          throw "sink node is not chosen";

  graph.AddAuxArcs(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);
	graph.PrepareNodes(&graph.nodes[sourceNodeIndex], &graph.nodes[sinkNodeIndex]);
	graph.PrepareArcs();
    
//   for (ArcIterator an_arc = graph.arcs.begin(); an_arc != graph.arcs.end(); an_arc++) {
//     mexPrintf("First if: %ld\n", an_arc->cap);
//   }
}

}

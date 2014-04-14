//(C) Microsoft Corporation. All rights reserved.
// adapted by Joao Carreira, University of Bonn, January 2009

#include "stdafx.h"
#include "joao_parser.h"
#include "Solver.h"
#include "Utils.h"
#include <string.h>
#include "mex.h"

using namespace paraF;

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	Graph graph;
	InputParser parser;
	Problem problem;
	Solver solver;
		
  // INPUTS -- number of inputs should be at least the first 5:
  // - a matrix with 4 columns: id_from id_to capacity slope (the capacity is the offset, slope multiplies by lambda)
  // - source_id
  // - sink_id
  // - m edges
  // - n nodes
  // optional:
  // - lambda_min
  // - lambda_max (need to provide both)
  //
  // OUTPUTS -- number of outputs should be 2:
  // - a matrix with 2 columns and n_nodes rows
  // - a vector of lambda values for the breakpoints
  // - the precision, 1/m. it determines how small a lambda it can search (it can miss breakpoints but they will be 2/m away from some it gets).
  
  if(nrhs<5) {
    mexErrMsgTxt("Number of inputs should be at least 5!");
  }
  
  int ncols, nrows, nbreakpoints;
  LongType s, t, m, n;
  double *arg1;
  unsigned int *out_cuts; 
  double *out_lambdas;
  double *out_precision;
  
  ncols = mxGetN(prhs[0]);
  nrows = mxGetM(prhs[0]);

  //mexPrintf("ncols: %d \n nrows: %d \n", ncols, nrows);
  
  if(ncols!=4) {
    mexErrMsgTxt("First argument should have 4 columns: from to cap slope!");
  }
  
  arg1 = mxGetPr(prhs[0]);

  s = mxGetScalar(prhs[1]);
  t = mxGetScalar(prhs[2]);
  m = mxGetScalar(prhs[3]);
  n = mxGetScalar(prhs[4]);
  
  if(0 > s || 0 > t) {
    mexErrMsgTxt("s and t must be bigger than or equal to 1!");
  }
  
  if(m != nrows) {
    mexErrMsgTxt("The first argument and the number of edges don't add up!");
  }

  //
  // Now the outputs
  //
  if(nlhs > 3) {
    mexErrMsgTxt("Expecting three outputs!");
  }     
  
  //
  // Now the computation!
  //    
	try
	{
		int minLambda;
		int maxLambda;
    
    if(nrhs==7) {	
      problem.hasLambdaRange = true;
    	problem.minLambda = mxGetScalar(prhs[5]);
  		problem.maxLambda = mxGetScalar(prhs[6]);
			//mexPrintf("min_lambda: %d\n max_lambda: %d\n", problem.minLambda, problem.maxLambda);
    }
 
    //solver.mode = Solver::GGT_ST_ONLY_ALGORITHM;
    solver.mode = Solver::GGT_TS_ONLY_ALGORITHM;
    //solver.mode = Solver::GGT_BOTH_ALGORITHM;    
    //solver.mode = Solver::SIMPLE_ALGORITHM;
    solver.printCuts = false;
    solver.verbosity = Solver::QUIET_VERBOSITY;    

		//int time0 = GetTime();
		
    parser.Parse(graph, problem, arg1, s, t, m, n);        
		
    //int time1 = GetTime();
		//printf("c Parse time: %d\n", time1 - time0);
		//printf("c\n");
        
    solver.Solve(graph, problem);
    
    //int time2 = GetTime();
    //printf("c Solve time: %d\n", time2 - time1);
    
    nbreakpoints = solver.numBreakpoints;
    plhs[0] = mxCreateNumericMatrix(n, 2, mxUINT32_CLASS,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(nbreakpoints, 1, mxREAL);
    out_cuts = (unsigned int*) mxGetData(plhs[0]);
    out_lambdas = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    out_precision = mxGetPr(plhs[2]);
    
    //
    // fill in the cuts
    //
		int color = 1;
		int j = 0;
    
    std::set<int>::iterator i = (solver.breakpointCutPointers.begin());
    i++;
    std::set<int>::iterator i_end = (solver.breakpointCutPointers.end());
    i_end--;
    
		for (;i != i_end; i++){
			while (j < *i)
			{
				out_cuts[j] = graph.typedNodeList[j]->index + 1;
        out_cuts[j+n] = color;        
				j++;
			}
			color++;
		}
    
    // fill in the lambdas

    for(int i = 0; i != nbreakpoints; i++) {      
      //mexPrintf("the bps: %lf\n", solver.lambdas[i]);    
      out_lambdas[i] = solver.lambdas[i];
    }
    
    // fill in the minimum possible lambda
    out_precision[0] = (double) 1/solver.multiplier;
	}
	catch (const char* msg)
	{
		printf("ERROR: %s\n", msg);
		//return 1;
	}
}


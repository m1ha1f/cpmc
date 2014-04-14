(C) Microsoft Corporation. All rights reserved.
*** See License.txt for the use license *** 

-----
paraF version 1.0, 12/22/06
-----

paraF implements the Gallo-Grigoriadis-Tarjan [SIAM J. Comp. 18(1), 1989]
(GGT) parametric flow algorithm and some of its variants.
Developed by Maxim Babenko and Andrew Goldberg
The implementation is described in [Babenko & Goldberg,
Microsoft Tech. Report MSR-TR-2006-77, 2006].

-----

The program is designed to be compiled under
linux, cygwin, or Windows Services for Unix (SFU).
To compile run "make". Modify makefile if needed.

-----

The input/output formats are described in the Doc directory

-----

Usage: paraF [options]

By default, the program runs GGT on the standard input.

Valid options are:
  --input <filename>
  Read from <filename> instead of standard input stream.
  
  --range <lo> <hi>
  Only look for breakpoints in the range [lo, hi].
  
  --simple
  Simple mode: no flow or height amortization is performed.
  Maximum flow is always recomputed from scratch.
  Benefits: uses a single copy of graph. As a consequence,
  runs faster on most testcases.
  
  --st-only
  Uses flow amortization but always computes s-t flows.
  
  --ts-only
  Uses flow amortization but always computes t-s flows.

  --max-bp-only
  Computes the maximum breakpoint value only.
  
  --solve-max-flow-for-max-bp
  Compute maximum flow for maximum breakpoint.
  
  --print-cuts  
  Outputs a collection of minimum cuts corresponding to
  the discovered breakpoints.
  
  --normal-verbosity
  Outputs only breakpoint values (default).
  
  --debug-verbosity
  Outputs breakpoint values and some useful debugging information.
  
  --quiet-verbosity
  Does not output breakpoint values or debugging information.
    
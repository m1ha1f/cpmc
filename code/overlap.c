/* Copyright (C) 2010 Joao Carreira

 This code is part of the extended implementation of the paper:
 
 J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
 */

#include <omp.h>
#include <mex.h>

void overlap(unsigned int *intersections, unsigned int *reunions, mxLogical *segms, int nc, int nr) {
    int i, j, index_ij, index_ik, index_jk, k;
    
	{
#pragma omp parallel for private(index_ik, index_jk, index_ij, i, j, k)  
        for(i=0; i<nc; i++) { /* for each segment */
            for(j=i+1; j<nc; j++) { /* go through the others with j>i */
                index_ij = i*nc + j;
                intersections[index_ij] = 0;
                reunions[index_ij] = 0;            

                for( k=0; k<nr; k++) { /* compute intersections and reunions */
                    index_ik = i*nr + k;
                    index_jk = j*nr + k;
                    intersections[index_ij] = intersections[index_ij] + (segms[index_ik] * segms[index_jk]);
                    reunions[index_ij] = reunions[index_ij] + (segms[index_ik] | segms[index_jk]);
                }
            }
		}
	}
  return;
}

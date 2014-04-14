#ifndef __clique_h__
#define __clique_h__

#include "matrix.h"

void findCliques(const Matrix* pb, double threshold, Matrix* cliqueMembership, int* numCliques);

void getCliqueSize(const Matrix* cliqueMembership, const int numCliques, int* cliqueSize);

void findNbrs(const Matrix* cliqueMembership, const int numCliques, int** nbrs, int* nnbrs);

void lift(const double* cliques, const Matrix* cliqueMembership, Matrix* pixels);
void project(const Matrix* pixels, const Matrix* cliqueMembership, const int numCliques, double* cliques);

void lift(const int* cliques, const Matrix* cliqueMembership, Matrix* pixels);
void project(const Matrix* pixels, const Matrix* cliqueMembership, const int numCliques, int* cliques);

#endif

#ifndef BUTCHER_H
#define BUTCHER_H

void butcher(double*** a, double** b, double** c, char filename[], int* n);
void butcher_adaptive(double*** a, double** b, double** c1, double** c2, char filename[], int* n);

#endif
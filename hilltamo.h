#ifndef HILLTAMO
#define HILLTAMO

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lnutils.h"
#include "randnums.h"


void Get_Hillman_Tamori_Parameters (FILE *fp);
void Hillman_Tamori_Terminate (int tree_type, double diameter);
void Hillman_Tamori_Stem (int algorithm, int tree_type, double diameter);

#endif

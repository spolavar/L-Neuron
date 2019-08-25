#ifndef BURKE
#define BURKE

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lnutils.h"
#include "randnums.h"

void Get_Burke_Parameters (FILE *fp);
void Burke_Branch (int tree_type, double diameter);
void Burke_Stem (int tree_type, double diameter);

#endif

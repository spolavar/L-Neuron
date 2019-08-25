/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@	                                                           FILE: lnutils.h @
//@                           CHANGE HISTORY                                   @
//@                                                                            @
//@ DATE      AUTH     DESCRIPTION                                             @
//@ ----      ----     -----------                                             @
//@ 12-30-99  JLK      Added processing for mixed distributions                @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

#ifndef LNUTILS
#define LNUTILS

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "randnums.h"

#define ZERO		0.00000001
#define PI		3.141592654

#define SINGLE		0
#define MULTIPLE	1

#define GAUSSIAN	0
#define UNIFORM		1
#define CONSTANT	2
#define MIXED		3
#define GAMMA		4
#define EXP			5

#define X_TROPISM		0
#define Y_TROPISM		1
#define Z_TROPISM		2
#define SOMA_TROPISM	3

#define SOMA	1
#define AXON	2
#define APICAL	3
#define BASAL	4

#define TREETYPES	9

// Max # of different distribution for each parameter
#define MAX_MIXED	10

#define bool	int
#define false	0
#define true	1

typedef struct {
	double min;
	double max;
} Uniform_T;

typedef struct {
	double mean;
	double stdev;
	double min_range;	/* optional: min value JLK 12-29-99 */
	double max_range;	/* optional: max value JLK 12-29-99 */
} Gaussian_T ;

typedef struct {
	double alpha;
	double beta;
	double offset;		/* SCR: 05-15-01 */
	double min_range;	/* optional: min value JLK 12-29-99 */
	double max_range;	/* optional: max value JLK 12-29-99 */
} Gamma_T ;

typedef struct {
	int density_type;
	union Density_Union {
		Uniform_T uni;
		Gaussian_T gauss;
		Gamma_T gamma;
		double constant;
	} parms;
} Density_T ;

typedef struct {
	double prob;
	Density_T d;
} Density_Array_T;


typedef struct {
	int num_dists;
	Density_Array_T density [MAX_MIXED];
} Distribution_T;

bool Valid_Tropism (Density_T d);
char * Add_Tropism (int tree_type, char *tr_str);
void Add_Global_Tropism (int tree_type, int trop_type);
char * Do_Wrinkle (int tree_type, double len, char *wr_str);
void Global_Post_Processing (int tree_type);
void Convert_Str_To_Upper (char *s1, char *s2);
double Sample_Distribution (Distribution_T d);
Density_T Get_Gaussian_Distribution (int num_fields, char *parm, double val1, double val2, double val3, double val4);
Density_T Get_Uniform_Distribution (int num_fields, char *parm, double val1, double val2);
Density_T Get_Constant_Distribution (int num_fields, char *parm, double val1);
void Get_Mixed_Distribution (FILE *fp, char *line, int num_fields, char *mixparm, int num_dists, Distribution_T *d);
void Get_Distribution (FILE *fp, char *line, Distribution_T *d);

#endif

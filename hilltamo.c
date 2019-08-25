/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@	                                                          File: hilltamo.c @
//@                           CHANGE HISTORY                                   @
//@                                                                            @
//@ DATE      AUTH     DESCRIPTION                                             @
//@ ----      ----     -----------                                             @
//@ 12-29-99  JLK      Fixed end of file check and comment check               @
//@ 12-30-99  JLK      Added processing for mixed distributions                @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "genneuron.h"
#include "hilltamo.h"
#include "randnums.h"

extern char *object_s; /* production string for lparser drawer */
extern bool x_trop [];
extern bool y_trop [];
extern bool z_trop [];
extern bool soma_trop [];
extern bool wrinkle [];

/* Common Distributions */
extern Distribution_T dist_soma_diameter;
extern Distribution_T dist_initial_diameter [];
extern Distribution_T dist_number_of_trees [];
extern Distribution_T dist_tree_azimuth [];
extern Distribution_T dist_tree_elevation [];
extern Distribution_T dist_x_tropism [];
extern Distribution_T dist_y_tropism [];
extern Distribution_T dist_z_tropism [];
extern Distribution_T dist_soma_tropism [];

/* Hillman Morphometric Distributions */
extern Distribution_T dist_bifurcating_amplitude_angle [];
extern Distribution_T dist_bifurcating_orientation_angle [];
extern Distribution_T dist_bifurcating_power [];
extern Distribution_T dist_bifurcating_ratio [];
extern Distribution_T dist_consolidated_length [];
extern Distribution_T dist_terminal_length [];
extern Distribution_T dist_threshold [];
extern Distribution_T dist_taper [];
extern Distribution_T dist_contraction [];
extern Distribution_T dist_fragmentation [];
extern Distribution_T dist_power_linearization [];

char line [80];

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Get_Hillman_Tamori_Parameters (FILE *fp)                  @
//@                                                                            @
//@     Action: Gets the distribution parameters from the input file.          @
//@                                                                            @
//@     Input:  fp - input file with Hillman or Tamori parameters.             @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Get_Hillman_Tamori_Parameters (FILE *fp) {

	char parm [80];
	char PARM [80];
	int i, type;

	/* default value for power linearization constant */
	for (i = 0; i < TREETYPES; i++) {
		dist_power_linearization [i].num_dists = 1;
		dist_power_linearization [i].density [0].prob = 1.0;
		dist_power_linearization [i].density [0].d.density_type = CONSTANT;
		dist_power_linearization [i].density [0].d.parms.constant = 1.0;
	}

	fgets (line, sizeof (line), fp);

	while (!feof (fp)) {

         /* skip any invalid lines.  valid line must have a numeric in the first
		    column                                                  JLK 12-29-99 */
		if ((line [0] >= '0') && (line [0] <= '9')) { 
		
			sscanf (line, "%i %s\n", &type, parm);
			Convert_Str_To_Upper (PARM, parm);

			/* common parameters */
			if (strcmp (PARM, "SOMADIAM") == 0)
				Get_Distribution (fp, line, &dist_soma_diameter);
			else if (strcmp (PARM, "STEM_DIAMETER") == 0)
				Get_Distribution (fp, line, &dist_initial_diameter [type]);
			else if (strcmp (PARM, "N_STEMS") == 0)
				Get_Distribution (fp, line, &dist_number_of_trees [type]);
			else if (strcmp (PARM, "TREEAZIM") == 0)
				Get_Distribution (fp, line, &dist_tree_azimuth [type]);
			else if (strcmp (PARM, "TREEELEV") == 0)
				Get_Distribution (fp, line, &dist_tree_elevation [type]);

			else if (strcmp (PARM, "XTROPISM") == 0) {

				Get_Distribution (fp, line, &dist_x_tropism [type]);
				if (Valid_Tropism (dist_x_tropism [type].density [0].d))
					x_trop [type] = true;
				else {
					fprintf (stderr, "ERROR: Invalid value for X tropism: %s\n", line);
					exit (1);
				}
			}
			else if (strcmp (PARM, "YTROPISM") == 0) {

				Get_Distribution (fp, line, &dist_y_tropism [type]);
				if (Valid_Tropism (dist_y_tropism [type].density [0].d))
					y_trop [type] = true;
				else {
					fprintf (stderr, "ERROR: Invalid value for Y tropism: %s\n", line);
					exit (1);
				}
			}
			else if (strcmp (PARM, "ZTROPISM") == 0) {

				Get_Distribution (fp, line, &dist_z_tropism [type]);
				if (Valid_Tropism (dist_z_tropism [type].density [0].d))
					z_trop [type] = true;
				else {
					fprintf (stderr, "ERROR: Invalid value for Z tropism: %s\n", line);
					exit (1);
				}
			}
			else if (strcmp (PARM, "SOMATROPISM") == 0) {

				Get_Distribution (fp, line, &dist_soma_tropism [type]);
				if (Valid_Tropism (dist_soma_tropism [type].density [0].d))
					soma_trop [type] = true;
				else {
					fprintf (stderr, "ERROR: Invalid value for SOMA tropism: %s\n", line);
					exit (1);
				}
			}
			
			/* hillman specific parameters */
			else if (strcmp (PARM, "BIFAMPLITUDE") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_amplitude_angle [type]);
			else if (strcmp (PARM, "BIFORIENT") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_orientation_angle [type]);
			else if (strcmp (PARM, "RALL_POWER") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_power [type]);
			else if (strcmp (PARM, "DAUGHTER_RATIO") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_ratio [type]);
			else if (strcmp (PARM, "IBF_BRANCH_PATHLENGTH") == 0)
				Get_Distribution (fp, line, &dist_consolidated_length [type]);
			else if (strcmp (PARM, "TERM_BRANCH_PATHLENGTH") == 0)
				Get_Distribution (fp, line, &dist_terminal_length [type]);
			else if (strcmp (PARM, "DIAM_THRESHOLD") == 0)
				Get_Distribution (fp, line, &dist_threshold [type]);
			else if (strcmp (PARM, "TAPER_2") == 0)
				Get_Distribution (fp, line, &dist_taper [type]);
			else if (strcmp (PARM, "CONTRACTION") == 0) {
				Get_Distribution (fp, line, &dist_contraction [type]);
				wrinkle [type] = true;
			}
			else if (strcmp (PARM, "FRAGMENTATION") == 0)
				Get_Distribution (fp, line, &dist_fragmentation [type]);
			else if (strcmp (PARM, "PK") == 0)
				Get_Distribution (fp, line, &dist_power_linearization [type]);
			else {
				fprintf (stderr, "ERROR: unknown parameter %s\n", parm);
				exit (1);
			}
		} /* end if check for comment */
	
		fgets (line, sizeof (line), fp);
		
	} /* end while */
	
} /* end Get_Hillman_Tamori_Parameters */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Hillman_Tamori_Terminate (int tree_type, double diameter) @
//@                                                                            @
//@     Action: Writes the production string for a terminating segment.        @
//@                                                                            @
//@     Input:  tree_type - type of tree for this stem.                        @
//@             diameter - the diameter of the stem segment                    @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Hillman_Tamori_Terminate (int tree_type, double diameter) {

	double length = Sample_Distribution (dist_terminal_length [tree_type]);
	double taper = Sample_Distribution (dist_taper [tree_type]);
	char term_str [2048] = "";
	char s1 [1024], s2 [1024];

	sprintf (term_str, "?(%3.4lf)%s(%3.2lf)?(%3.2lf)%s(%3.2lf)]",
		 2*diameter/length,
		 Do_Wrinkle (tree_type, length / 2.0, s1),
		 length / 2.0,
		 1 - taper,
		 Do_Wrinkle (tree_type, length / 2.0, s2),
		 length / 2.0);

	strcat (object_s, term_str);

} /* end Hillman_Tamori_Terminate */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Hillman_Tamori_Stem (int algorithm,                       @ 
//@                                       int tree_type,                       @
//@                                       double diameter)                     @
//@                                                                            @
//@     Action: Writes the production string for a stem segment.  This is      @
//@             a recursive procedure.                                         @
//@                                                                            @
//@     Input:  algorithm - Hillman or Tamori algorithm.                       @        
//@             tree_type - type of tree for this stem.                        @
//@             diameter - the diameter of the stem segment                    @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Hillman_Tamori_Stem (int algorithm, int tree_type, double diameter) {
  double length, taper;
  double bifang1, bifang2, diam1, diam2, rall, ratio, power, totbifang, dimension;
  double eq5_4_den, eq5_4_n1, eq5_4_n2;  
  char stem_str [4096] ;
  char s1 [1024], s2 [1024];

  length = Sample_Distribution (dist_consolidated_length [tree_type]);
  taper = Sample_Distribution (dist_taper [tree_type]);

		//	 2*diameter/length,

//	fprintf(stderr,"+");
	
	sprintf (stem_str, "?(%3.4lf)%s(%3.2lf)?(%3.2lf)%s(%3.2lf)?(%3.2lf)",
			  2*diameter/length,
		 	 Do_Wrinkle (tree_type, length / 2.0, s1),
			 length / 2.0,
			 1 - taper,
		 	 Do_Wrinkle (tree_type, length / 2.0, s2),
			 length / 2.0,
			 length / (2*diameter*(1-taper)));
//	fprintf(stderr,"-");
	strcat (object_s, stem_str);
	diameter = diameter * (1-taper);
	/* If the diameter is less than a given threshold the branch ends 	*/
	/* with a terminating segment.  Otherwise, make recursive calls to 	*/
	/* create two branch stems.						*/

	if (diameter < Sample_Distribution (dist_threshold [tree_type])) 
		Hillman_Tamori_Terminate (tree_type, diameter);

	else {

		sprintf (stem_str, ">(%3.2lf)", Sample_Distribution (dist_bifurcating_orientation_angle [tree_type]));
		strcat (object_s, stem_str);
		ratio = Sample_Distribution (dist_bifurcating_ratio [tree_type]);
		power = Sample_Distribution (dist_bifurcating_power [tree_type]);
		rall = pow ((1 + pow (ratio, power)), -1/power);
		diam2 = diameter * rall * pow (Sample_Distribution (dist_power_linearization [tree_type]), 1/power);
		diam1 = ratio * diam2;

		if (algorithm == HILLMAN) {
			totbifang = Sample_Distribution (dist_bifurcating_amplitude_angle [tree_type]);
			bifang1 = rnd01 () * totbifang;
			bifang2 = bifang1 - totbifang;

			/* for display purposes, reverse the angles half the time */
			if (rnd01 () < 0.5) {
				bifang1 *= -1;
				bifang2 *= -1;
			}
		} /* end if Hillman Algorithm */

		else if (algorithm == TAMORI) { /* Tamori Algorithm */
				
			dimension = 1.0 + rnd01 () * (power-1);

			/* Equation 5.4 from Tamori, 1993 */
			eq5_4_n1 = pow (1+pow(ratio,power), 2*dimension/power);
			eq5_4_n2 = pow(ratio, 2*dimension);
			eq5_4_den = 2*pow (1+pow(ratio,power), dimension/power);
			bifang1 = acos ((eq5_4_n1 + eq5_4_n2 - 1) / (eq5_4_den * pow (ratio, dimension)));
			bifang1 = 180/PI * bifang1;
			bifang2 = acos ((eq5_4_n1 - eq5_4_n2 + 1) / (eq5_4_den));
			bifang2 = 180/PI * bifang2;

			/* for display purposes, reverse the angles half the time */
			if (rnd01 () < 0.5)
				bifang1 *= -1;
			else
				bifang2 *= -1;

		} // end else Tamori Algorithm

		/* first branch */
		sprintf (stem_str, "[+(%3.2lf)", bifang1);
		strcat (object_s, stem_str);
		Hillman_Tamori_Stem (algorithm, tree_type, diam1);

		/* second branch */
		sprintf (stem_str, "+(%3.2lf)", bifang2);
		strcat (object_s, stem_str);
		Hillman_Tamori_Stem (algorithm, tree_type, diam2);

	} /* end else diam > thr */

} /* end Hillman_Tamori_Stem */


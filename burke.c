/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@	                                                             File: burke.c @
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
#include "burke.h"
#include "randnums.h"

extern char *object_s; /* production string for lparser drawer */
extern bool x_trop [];
extern bool y_trop [];
extern bool z_trop [];
extern bool soma_trop [];

/* Common Distributions */
extern Distribution_T dist_soma_diameter;
extern Distribution_T dist_initial_diameter[];
extern Distribution_T dist_number_of_trees [];
extern Distribution_T dist_tree_azimuth [];
extern Distribution_T dist_tree_elevation [];
extern Distribution_T dist_x_tropism [];
extern Distribution_T dist_y_tropism [];
extern Distribution_T dist_z_tropism [];
extern Distribution_T dist_soma_tropism [];

/* Burke Morphometric Distributions */
extern Distribution_T dist_bifurcating_amplitude_angle [];
extern Distribution_T dist_bifurcating_orientation_angle [];
extern Distribution_T dist_bin_length [];
extern Distribution_T dist_gaussian_branch [];
extern Distribution_T dist_k1_overlap [];
extern Distribution_T dist_k2_overlap [];
extern Distribution_T dist_k1_nonoverlap [];
extern Distribution_T dist_k2_nonoverlap [];
extern Distribution_T dist_k1_terminate [];
extern Distribution_T dist_k2_terminate [];
extern Distribution_T dist_linear_branch [];
extern Distribution_T dist_taper [];
extern Distribution_T dist_extending_azimuth [];
extern Distribution_T dist_extending_elevation [];

	char line [80];

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Get_Burke_Parameters (FILE *fp)                           @
//@                                                                            @
//@     Action: Gets the distribution parameters from the input file.          @
//@                                                                            @
//@     Input:  fp - input file with Hillman parameters.                       @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Get_Burke_Parameters (FILE *fp) {

	char parm [80];
	char PARM [80];
	int type;

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
								
			/* burke specific parameters */
			else if (strcmp (PARM, "BIFAMPLITUDE") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_amplitude_angle [type]);
			else if (strcmp (PARM, "BIFORIENT") == 0)
				Get_Distribution (fp, line, &dist_bifurcating_orientation_angle [type]);
			else if (strcmp (PARM, "BINLENGTH") == 0)
				Get_Distribution (fp, line, &dist_bin_length [type]);
			else if (strcmp (PARM, "GAUSSBRANCH") == 0)
				Get_Distribution (fp, line, &dist_gaussian_branch [type]);
			else if (strcmp (PARM, "K1OVERLAP") == 0)
				Get_Distribution (fp, line, &dist_k1_overlap [type]);
			else if (strcmp (PARM, "K2OVERLAP") == 0)
				Get_Distribution (fp, line, &dist_k2_overlap [type]);
			else if (strcmp (PARM, "K1NONOVERLAP") == 0)
				Get_Distribution (fp, line, &dist_k1_nonoverlap [type]);
			else if (strcmp (PARM, "K2NONOVERLAP") == 0)
				Get_Distribution (fp, line, &dist_k2_nonoverlap [type]);
			else if (strcmp (PARM, "K1TERMINATE") == 0)
				Get_Distribution (fp, line, &dist_k1_terminate [type]);
			else if (strcmp (PARM, "K2TERMINATE") == 0)
				Get_Distribution (fp, line, &dist_k2_terminate [type]);
			else if (strcmp (PARM, "LINBRANCH") == 0)
				Get_Distribution (fp, line, &dist_linear_branch [type]);
			else if (strcmp (PARM, "TAPER_2") == 0)
				Get_Distribution (fp, line, &dist_taper [type]);
			else if (strcmp (PARM, "EXAZIM") == 0)
				Get_Distribution (fp, line, &dist_extending_azimuth [type]);
			else if (strcmp (PARM, "EXELEV") == 0)
				Get_Distribution (fp, line, &dist_extending_elevation [type]);
			else {
				fprintf (stderr, "ERROR: unknown parameter %s\n", parm);
				exit (1);
			}
		} /* end if comment */
		
		fgets (line, sizeof (line), fp);
	} /* end while */
} /* end Get_Burke_Parameters */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Burke_Branch (int tree_type, double diameter)             @
//@                                                                            @
//@     Action: Writes the production string for a branching segment.          @
//@                                                                            @
//@     Input:  tree_type - type of tree for this branch.                      @
//@             diameter - the diameter of the branch segment                  @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Burke_Branch (int tree_type, double diameter) {

	double r1 = Sample_Distribution (dist_gaussian_branch [tree_type]);
	double r2 = Sample_Distribution (dist_gaussian_branch [tree_type]);
	double burkea = Sample_Distribution (dist_linear_branch [tree_type]);
	double diam1 = diameter * (r1 + r2*burkea);
	double diam2 = diameter * (r2 + r1*burkea);
	double totbifang = Sample_Distribution (dist_bifurcating_amplitude_angle [tree_type]);
	double bifang1 = rnd01 () * totbifang;
	double bifang2 = bifang1 - totbifang;
	char branch_str [80] = "";

	sprintf (branch_str, ">(%3.2lf)", Sample_Distribution (dist_bifurcating_orientation_angle [tree_type]));
	strcat (object_s, branch_str);

	/* first branch */
	sprintf (branch_str, "[+(%3.2lf)", bifang1);
	strcat (object_s, branch_str);
	Burke_Stem (tree_type, diam1);
	
	/* second branch */
	sprintf (branch_str, "+(%3.2lf)", bifang2);
	strcat (object_s, branch_str);
	Burke_Stem (tree_type, diam2);
	
} /* end Burke_Branch */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Burke_Stem (int tree_type, double diameter)               @
//@                                                                            @
//@     Action: Writes the production string for a stem segment.  This is      @
//@             a recursive procedure.                                         @
//@                                                                            @
//@     Input:  tree_type - type of tree for this stem.                        @
//@             diameter - the diameter of the stem segment                    @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Burke_Stem (int tree_type, double diameter) {

	double length = Sample_Distribution (dist_bin_length [tree_type]);
	double taper = Sample_Distribution (dist_taper [tree_type]);
	double k1_overlap = Sample_Distribution (dist_k1_overlap [tree_type]);
	double k2_overlap = Sample_Distribution (dist_k2_overlap [tree_type]);
	double k1_nonoverlap = Sample_Distribution (dist_k1_nonoverlap [tree_type]);
	double k2_nonoverlap = Sample_Distribution (dist_k2_nonoverlap [tree_type]);
	double k1_terminate = Sample_Distribution (dist_k1_terminate [tree_type]);
	double k2_terminate = Sample_Distribution (dist_k2_terminate [tree_type]);
	double pbr, pnonoverlap, ptr;
	char   stem_str [128] = "";
	char   s [40];

	sprintf (stem_str, "?(%3.2lf)%s(%3.2lf)?(%3.2lf)",
			 diameter/length,
			 Add_Tropism (tree_type, s),
			 length,
			 length/diameter);
	strcat (object_s, stem_str);

	diameter += taper*length;
	
	pbr = k1_overlap * exp (k2_overlap * diameter);
	ptr = k1_terminate * exp (k2_terminate * diameter);
	pnonoverlap = k1_nonoverlap * exp (k2_nonoverlap * diameter);
	
	if (pnonoverlap < pbr)
		pbr = pnonoverlap;
		
	/* decide whether to branch, terminate or grow */

	if (rnd01 () < pbr*length) /* branch */
		Burke_Branch (tree_type, diameter);
		
	else if (rnd01 () < ptr*length) { /* terminate */
		sprintf (stem_str, "]");
		strcat (object_s, stem_str);
	} /* end terminate */
	
	else { /* grow a stem */
		sprintf (stem_str, ">(%3.2lf)+(%3.2lf)",
				 Sample_Distribution (dist_extending_azimuth [tree_type]),
				 Sample_Distribution (dist_extending_elevation [tree_type]));
		strcat (object_s, stem_str);
		Burke_Stem (tree_type, diameter);
	}	
} /* end Burke_Stem */


#include "genneuron.h"
#include "lnutils.h"
#include "burke.h"
#include "hilltamo.h"

/* External Variables ----------------------------------------------------- */
extern char *object_s; /* production string for lparser drawer */
extern x_trop [];
extern y_trop [];
extern z_trop [];
extern soma_trop [];
extern wrinkle [];

const int treecolor [TREETYPES] = {	DRED, SOMACOLOR, AXONCOLOR, APICALCOLOR, 
										BASALCOLOR, DBLUE, MAGENTA, DGREEN, AQUA };

/* Common Distributions */
Distribution_T dist_soma_diameter;
Distribution_T dist_initial_diameter [TREETYPES];
Distribution_T dist_number_of_trees [TREETYPES];
Distribution_T dist_tree_azimuth [TREETYPES];
Distribution_T dist_tree_elevation [TREETYPES];
Distribution_T dist_x_tropism [TREETYPES];
Distribution_T dist_y_tropism [TREETYPES];
Distribution_T dist_z_tropism [TREETYPES];
Distribution_T dist_soma_tropism [TREETYPES];
Distribution_T dist_bifurcating_amplitude_angle [TREETYPES];
Distribution_T dist_bifurcating_orientation_angle [TREETYPES];
Distribution_T dist_bifurcating_power [TREETYPES];
Distribution_T dist_bifurcating_ratio [TREETYPES];
Distribution_T dist_taper [TREETYPES];
Distribution_T dist_consolidated_length [TREETYPES];
Distribution_T dist_terminal_length [TREETYPES];
Distribution_T dist_threshold [TREETYPES];
Distribution_T dist_contraction [TREETYPES];
Distribution_T dist_fragmentation [TREETYPES];
Distribution_T dist_power_linearization [TREETYPES];

/* Burke Morphometric Distributions */
Distribution_T dist_bin_length [TREETYPES];
Distribution_T dist_gaussian_branch [TREETYPES];
Distribution_T dist_k1_overlap [TREETYPES];
Distribution_T dist_k2_overlap [TREETYPES];
Distribution_T dist_k1_nonoverlap [TREETYPES];
Distribution_T dist_k2_nonoverlap [TREETYPES];
Distribution_T dist_k1_terminate [TREETYPES];
Distribution_T dist_k2_terminate [TREETYPES];
Distribution_T dist_linear_branch [TREETYPES];
Distribution_T dist_extending_azimuth [TREETYPES];
Distribution_T dist_extending_elevation [TREETYPES];

int num_tree_types = 0;

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Grow_Soma (double d)                  			   	       @
//@                                                                            @
//@     Action: Generates the production string for a soma of size d.          @
//@                                                                            @
//@     Input:  d - diameter of the soma.                                      @
//@     Output: none                                                           @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Grow_Soma (double d) {

	char soma_str [128];
	sprintf (soma_str, "c(%i)[?(0.1)'(%3.2lf)", SOMACOLOR, d*OCTAGON_FACTOR); 
	strcpy (object_s, soma_str);
	sprintf (soma_str, "[+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(45)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(135)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(90)+(45)>(-90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(90)+(90)>(-90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[>(90)+(135)>(-90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[+(135)>(90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	sprintf (soma_str, "[+(45)>(90)+(22.5)z(%3.2lf)+(112.5)F+F+F+F+F+F+F+F]", d/2.0);
	strcat (object_s, soma_str);
	strcat (object_s, "]");

} // end Grow_Soma


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Initialize_Neuron ()                  			           @
//@                                                                            @
//@     Action: Initializes any state data used by the neuron.                 @
//@                                                                            @
//@     Input:  none.                                                          @
//@     Output: none.                                                          @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Initialize_Neuron (int seed) {

	int i;
	
	rndseed (seed);
	
	for (i = 0; i < TREETYPES; i++) {
		x_trop [i] = false;
		y_trop [i] = false;
		z_trop [i] = false;
		soma_trop [i] = false;
		wrinkle [i] = false;
	}
	
} // end Initialize_Neuron

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  double Grow_Neuron (int ntype)                  			   @
//@                                                                            @
//@     Action: Generates a neuron using the specified neuronal algorithm      @
//@                                                                            @
//@     Input:  ntype - type of neuronal growth algorithm.                     @
//@     Output: returns the soma diameter of the neuron                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

double Grow_Neuron (int ntype) {

	int numtrees;
	double azim, elev;
	char neuron_str [80];
	double soma_diameter = Sample_Distribution (dist_soma_diameter);
	int i, j;
	
	Grow_Soma (soma_diameter);
	
	for (i = APICAL; i <= num_tree_types; i++) {
	
		numtrees = (int) Sample_Distribution (dist_number_of_trees [i]);
		if(numtrees<1){
					fprintf (stderr, "ERROR:NrTree==0!!!\n");
				numtrees=1;
		}
		for (j = 0; j < numtrees; j++) {

			sprintf (neuron_str, "c(%i)", treecolor [i]); 
			strcat (object_s, neuron_str);

			elev = Sample_Distribution (dist_tree_elevation [i]);
			azim = Sample_Distribution (dist_tree_azimuth [i]);
			sprintf (neuron_str, "[>(%3.2lf)+(%3.2lf)z(%3.2lf)", azim, elev, soma_diameter*0.924/2); /* somad*cos (22.5)/2 */
			strcat (object_s, neuron_str);
	
			switch (ntype) {
			
				case BURKEMARKS:
					Burke_Stem (i, Sample_Distribution (dist_initial_diameter [i]));
					break;
				case HILLMAN:
					Hillman_Tamori_Stem (HILLMAN, i, Sample_Distribution (dist_initial_diameter [i]));
					break;
				case TAMORI:
					Hillman_Tamori_Stem (TAMORI, i, Sample_Distribution (dist_initial_diameter [i]));
					break;

			} /* end switch */
		} /* end for j */

		/* global post processing not currently used */
		// Global_Post_Processing (i);	/* alter string for tropism, metabolic contraints, etc. */	

	} /* end for i */

	return soma_diameter;
	
} /* end Grow_Neuron */

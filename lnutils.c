/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@	                                                           File: lnutils.c @
//@                           CHANGE HISTORY                                   @
//@                                                                            @
//@ DATE      AUTH     DESCRIPTION                                             @
//@ ----      ----     -----------                                             @
//@ 12-29-99  JLK      Added max value for Gaussian distribution               @
//@ 12-30-99  JLK      Added processing for mixed distributions                @
//@ 12-30-99  JLK      Added error checking for zero tropism                   @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

#include "lnutils.h"

extern Distribution_T dist_x_tropism [];
extern Distribution_T dist_y_tropism [];
extern Distribution_T dist_z_tropism [];
extern Distribution_T dist_soma_tropism [];
extern Distribution_T dist_contraction [];
extern Distribution_T dist_fragmentation [];
extern char *object_s; /* production string for lparser drawer */
extern int num_tree_types;

#define LINESIZE 80

bool x_trop [TREETYPES];
bool y_trop [TREETYPES];
bool z_trop [TREETYPES];
bool soma_trop [TREETYPES];
bool wrinkle [TREETYPES];

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  bool Valid_Tropism (Distribution_T d)                          @
//@                                                                            @
//@     Action: Checks that zero tropism is not requested.                     @
//@                                                                            @
//@     Input:  d - distribution to validate.                                  @
//@     Output: returns true if constant value for tropism is non-zero         @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

bool Valid_Tropism (Density_T d) {

	return ((d.density_type != CONSTANT) || (d.parms.constant > ZERO));

} /* end Valid_Tropism

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Add_Tropism (int trop_type, char *tr_str)                 @
//@                                                                            @
//@     Action: Creates tropism for a segment of the tree.                     @
//@                                                                            @
//@     Input:  tree_type - apply tropism to this tree type                    @
//@             tr_str - string to add tropism to.                             @
//@     Output: returns tr_str                 								   @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

char * Add_Tropism (int tree_type, char *tr_str) {

	char x_str [20], y_str [20], z_str [20], s_str [20];
	tr_str [0] = '\0';
		
	if (x_trop [tree_type]) {	
		sprintf (x_str, "x(%3.3lf)", Sample_Distribution (dist_x_tropism [tree_type]));
		strcpy (tr_str, x_str);
	}
	if (y_trop [tree_type]) {	
		sprintf (y_str, "y(%3.3lf)", Sample_Distribution (dist_y_tropism [tree_type]));
		if (!x_trop [tree_type])
			strcpy (tr_str, y_str);
		else
			strcat (tr_str, y_str);
	}
	if (z_trop [tree_type]) {	
		sprintf (z_str, "w(%3.3lf)", Sample_Distribution (dist_z_tropism [tree_type]));
		if (!x_trop [tree_type] && !y_trop [tree_type]) 
			strcpy (tr_str, z_str);
		else
			strcat (tr_str, z_str);
	}
	if (soma_trop [tree_type]) {	
		sprintf (s_str, "s(%3.3lf)", Sample_Distribution (dist_soma_tropism [tree_type]));
		if (!x_trop [tree_type] && !y_trop [tree_type] && !z_trop [tree_type] )
			strcpy (tr_str, s_str);
		else
			strcat (tr_str, s_str);
	}

	strcat (tr_str, "F");
	return tr_str;
		
} /* end Add_Tropism */


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Add_Global_Tropism (int trop_type, int trop_type)         @
//@                                                                            @
//@     Action: Adds tropism globally to the production string.                @
//@                                                                            @
//@     Input:  tree_type - apply tropism to this tree type                    @
//@             trop_type - tropism type (x, y, z, or soma)                    @
//@     Output: none                         								   @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Add_Global_Tropism (int tree_type, int trop_type) {

	Distribution_T d;
	char t_sym;
	unsigned int ui, uj;
	int i;
	int inx;
	char *tmp_str;
	char trop_str [80] = "";
	int x = strlen (object_s);
	tmp_str = (char *) malloc (4*strlen (object_s));
	
	switch (trop_type) {
	
		case X_TROPISM:
			d = dist_x_tropism [tree_type];
			t_sym = 'x';
			break;
			
		case Y_TROPISM:
			d = dist_y_tropism [tree_type];
			t_sym = 'y';
			break;
			
		case Z_TROPISM:
			d = dist_z_tropism [tree_type];
			t_sym = 'w';
			break;
			
		case SOMA_TROPISM:
			d = dist_soma_tropism [tree_type];
			t_sym = 's';
			break;
			
	} /* end switch */
	
	/* skip over the soma - the end of the soma production string 	*/
	/* is marked by a second c										*/
	inx = 0;
	while (object_s [inx] != 'c') {
		tmp_str [inx] = object_s [inx];
		++inx;
	}
	tmp_str [inx] = object_s [inx];
	++inx;

	/* skip over the value "(nn)" */
	if (object_s [inx] == '(' ) {
		while (object_s [inx] != ')' ) {
			tmp_str [inx] = object_s [inx];
			++inx;
		}
		tmp_str [inx] = object_s [inx];
		++inx;
	}

	/* for the given tree type skip over the other trees 		*/
	/* trees are marked by their color. it is assumed that the  */
	/* order is soma then apical tree then basal etc followed   */
	/* by the axon if there is one.								*/
	for (i = 0; i < tree_type - AXON; i++) {
		while (object_s [inx] != 'c') {
			tmp_str [inx] = object_s [inx];
			++inx;
		}
		tmp_str [inx] = object_s [inx];
		++inx;

		/* skip over the value "(nn)" */
		if (object_s [inx] == '(' ) {
			while (object_s [inx] != ')' ) {
				tmp_str [inx] = object_s [inx];
				++inx;
			}
			tmp_str [inx] = object_s [inx];
			++inx;
		}
	} /* end for i */
	
	/*	next parse through the tree string and replace any 'F' command with 
		a tropism x(nn)F command.											*/
		
	for (ui = inx; ui < strlen (object_s); ui++) {
	
		if (object_s [ui] == 'F') {
			sprintf (trop_str, "%c(%3.3lf)F", t_sym, Sample_Distribution (d));
			
			for (uj = 0; uj < strlen (trop_str); uj++) {
				tmp_str [inx] = trop_str [uj];
				++inx;
			}
		}
		else {
			tmp_str [inx] = object_s [ui];
			++inx;
		}
	} /* end for i */
		
	tmp_str [inx] = '\0';
	strcpy (object_s, tmp_str);
	
} /* end Add_Global_Tropism */


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  char * Do_Wrinkle (int trop_type, int len, char *tr_str)       @
//@                                                                            @
//@     Action: Puts contraction and fragmentation into a segment.             @
//@                                                                            @
//@     Input:  tree_type - apply tropism to this tree type                    @
//@             len - length of segment to split up.                           @
//@             wr_str - string with wrinkled segment and tropism.             @
//@     Output: returns tr_str                 								   @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

char * Do_Wrinkle (int tree_type, double len, char *wr_str) {

	double contract;
	double fragments;
	double alpha, plusminus;
	double evenodd = -1;
	char   tr_str [40], s [40];
	int i;
	static int l=0;

/*	//fprintf(stderr,"%d ",l);
	strcpy(tr_str,"");
	strcpy(s,"");

	/*
l1=strlen(w);
if(l1>l){
		l=l1;
		fprintf(stderr,"%d\n",l);
}
*/

	//SCR remove wrinklefor debug
	if (!wrinkle [tree_type]) 
		Add_Tropism (tree_type, wr_str);
	else {
		contract = Sample_Distribution (dist_contraction [tree_type]);
		fragments = Sample_Distribution (dist_fragmentation [tree_type]) - 1;
		
		if (fragments < 1)
			fragments = 1;
	//SCR remove wrinklefor debug
			//fragments = 1;

		if (rnd01 () < 0.5)
			plusminus = 1;
		else
			plusminus = -1;

		alpha = acos (contract)*180.0*plusminus/PI;
		sprintf (wr_str, "+(%3.3lf)?(%3.3lf)%s(%3.3lf)?(0.5)",
			alpha,
			2*fragments,
			Add_Tropism (tree_type, tr_str),
			len/(2*fragments));

		for (i = 0; i < (int) (fragments-1); i++) {
			sprintf (s, "+(%3.3lf)%s(%3.3lf)",
				2*evenodd*alpha,
				Add_Tropism (tree_type, tr_str),
				len/fragments);
			strcat (wr_str, s);
			evenodd *= -1;
			
		}
		sprintf (s, "?(2)+(%3.3lf)%s(%3.3lf)?(%3.3lf)+(%3.3lf)",
				2*evenodd*alpha,
				Add_Tropism (tree_type, tr_str),
				len/(2*fragments),
				1/(2*fragments),
				-evenodd*alpha);
		strcat (wr_str, s);

	} /* end else */

	return wr_str;

} // end Do_Wrinkle

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  void Global_Post_Processing (int tree_type)	                   @
//@                                                                            @
//@     Action: Applies global changes after the production string is complete.@
//@                                                                            @
//@     Input:  tree_type - indicates the tree to apply post processing        @                      									   @
//@     Output: none                         								   @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Global_Post_Processing (int tree_type) {

	if (x_trop [tree_type])
		Add_Global_Tropism (tree_type, X_TROPISM);
		
	if (y_trop [tree_type])
		Add_Global_Tropism (tree_type, Y_TROPISM);
		
	if (z_trop [tree_type])
		Add_Global_Tropism (tree_type, Z_TROPISM);
		
	if (soma_trop [tree_type])
		Add_Global_Tropism (tree_type, SOMA_TROPISM);
		
	
} /* end Global_Post_Processing

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  double Sample_Distribution (Distribution_T d)                  @
//@                                                                            @
//@     Action: Gets a sample value from the given distribution.               @
//@                                                                            @
//@     Input:  d - parameter and distribution type                            @
//@     Output: returns a value from the distribution.                         @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

double Sample_Distribution (Distribution_T d) {

	double x, x_d;
	int i,constant=1;

	//check distribution is compused of only one constant value 
	//then don't generate a random number!!!

	if(d.num_dists==1 && d.density [0].d.density_type==CONSTANT)
		return d.density [0].d.parms.constant;

	x_d = rnd01 ();

	for (i = 0; i < d.num_dists; i++) {

		if (x_d < d.density [i].prob) {
			if (d.density [i].d.density_type == UNIFORM) 
				x = rnduniform (d.density [i].d.parms.uni.min, d.density [i].d.parms.uni.max);

			else if (d.density [i].d.density_type == GAUSSIAN) {  /* Gaussian */	
				do {
					x = rndgauss (d.density [i].d.parms.gauss.mean, d.density [i].d.parms.gauss.stdev);
				} while ((x < d.density [i].d.parms.gauss.min_range) || (x > d.density [i].d.parms.gauss.max_range));
			}
			else if (d.density [i].d.density_type == GAMMA) {  /* Gamma */	
				do {
					//x = rndgamma (d.density [i].d.parms.gauss.mean, d.density [i].d.parms.gauss.stdev);
					x = rndgamma (d.density [i].d.parms.gamma.alpha, d.density [i].d.parms.gamma.beta, d.density [i].d.parms.gamma.offset);  /* SCR 05-15-01 */
				} while ((x < d.density [i].d.parms.gamma.min_range) || (x > d.density [i].d.parms.gamma.max_range));
			}
			else /* constant value */
				x = d.density [i].d.parms.constant;

			i = d.num_dists; /* exit the loop */
		}
	} /* end for i */

	return x;

} /* end double Sample_Distribution */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title:  Convert_Str_To_Upper (char *s1, char *s2)                      @
//@                                                                            @
//@     Action: Converts an input string to all upper case.                    @
//@                                                                            @
//@     Input:  s1 - the converted string                                      @
//@             s2 - string to be converted                                    @
//@     Output: converted string                                               @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Convert_Str_To_Upper (char *s1, char *s2) {

	unsigned int i;
	for (i = 0; i <= strlen (s2); i++)
		s1 [i] = toupper (s2 [i]);

} /* end Convert_Str_To_Upper */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Density_T Get_Gaussian_Distribution (                           @
//@                               int num_fields, double val1, double val2,    @
//@                               double val3, double val4)                    @
//@                                                                            @
//@     Action: Sets the parameters for the Gaussian distribution and returns  @
//@             the completed distribution.                                    @
//@                                                                            @
//@     Input:  num_fields - the number of fields read from the input line     @
//@             parm - name of the parameter getting the distribution          @
//@             val1, val2, val3, val4 - fields read from the input line       @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

Density_T Get_Gaussian_Distribution (int num_fields, char *parm, double val1, double val2, double val3, double val4) {

	Density_T d;

	if ((num_fields < 5) || (num_fields > 7)) {
		fprintf (stderr, "ERROR: incorrect fields for GAUSSIAN %s\n", parm);
		exit (1);
	}

	/* no min or max for gaussian */
	if (num_fields == 5) {
		val3 = val1 - 1000 * val2;	/* make the min extremely small */
		val4 = val1 + 1000 * val2;	/* make the max extremely large */
	}
	/* min only for gaussian */
	if (num_fields == 6) {
		val4 = val1 + 1000 * val2;	/* make the max extremely large */
	}
	d.density_type = GAUSSIAN;
	d.parms.gauss.mean = val1;
	d.parms.gauss.stdev = val2;
	d.parms.gauss.min_range = val3;
	d.parms.gauss.max_range = val4;

	return d;

} /* end Get_Gaussian_Distribution */
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Density_T Get_Gamma_Distribution (                           @
//@                               int num_fields, double val1, double val2,    @
//@                               double val3, double val4)                    @
//@                                                                            @
//@     Action: Sets the parameters for the Gaussian distribution and returns  @
//@             the completed distribution.                                    @
//@                                                                            @
//@     Input:  num_fields - the number of fields read from the input line     @
//@             parm - name of the parameter getting the distribution          @
//@             val1, val2, val3, val4 - fields read from the input line       @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

Density_T Get_Gamma_Distribution (int num_fields, char *parm, double val1, double val2, double val3, double val4, double val5) {

	Density_T d;

	if ((num_fields < 6) || (num_fields > 8)) {
		fprintf (stderr, "ERROR: incorrect fields for Gamma %s\n", parm);
		exit (1);
	}

	/* no min or max for Gamma */
	if (num_fields == 6) {
		val4 = val1 - 1000 * val2;	/* make the min extremely small */
		val5 = val1 + 1000 * val2;	/* make the max extremely large */
	}
	/* min only for Gamma */
	if (num_fields == 7) {
		val4 = val1 + 1000 * val2;	/* make the max extremely large */
	}
	d.density_type = GAMMA;
	d.parms.gamma.alpha = val1;
	d.parms.gamma.beta = val2;
	d.parms.gamma.offset = val3;
	d.parms.gamma.min_range = val4;
	d.parms.gamma.max_range = val5;

	return d;

} /* end Get_Gamma_Distribution */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Density_T Get_Uniform_Distribution (int num_fields,             @
//@                               char *parm, double val1, double val2)        @
//@                                                                            @
//@     Action: Sets the parameters for the Uniform distribution and returns   @
//@             the completed distribution.                                    @
//@                                                                            @
//@     Input:  num_fields - the number of fields read from the input line     @
//@             parm - name of the parameter getting the distribution          @
//@             val1, val2,  - fields read from the input line                 @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

Density_T Get_Uniform_Distribution (int num_fields, char *parm, double val1, double val2) {

	Density_T d;

	if (num_fields != 5) {
		fprintf (stderr, "ERROR: incorrect fields for UNIFORM %s\n", parm);
		exit (1);
	}
	d.density_type = UNIFORM;
	d.parms.uni.min = val1;
	d.parms.uni.max = val2;

	return d;

} /* end Get_Uniform_Distribution */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Density_T Get_Constant_Distribution (                           @
//@                               int num_fields, char *parm, double val1)     @
//@                                                                            @
//@     Action: Sets the parameters for the Constant distribution and returns  @
//@             the completed distribution.                                    @
//@                                                                            @
//@     Input:  num_fields - the number of fields read from the input line     @
//@             parm - name of the parameter getting the distribution          @
//@             val1  - fields read from the input line                        @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

Density_T Get_Constant_Distribution (int num_fields, char *parm, double val1) {

	Density_T d;
	
	if (num_fields != 4) {
		fprintf (stderr, "ERROR: incorrect fields for %s\n", parm);
		exit (1);
	}
	d.density_type = CONSTANT;
	d.parms.constant = val1;

	return d;

} /* Get_Constant_Distribution */

/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Distribution_T Get_Mixed_Distribution (FILE *fp, char *line     @
//@                               int num_fields, char *parm, int num_dists)   @
//@                                                                            @
//@     Action: Sets the parameters for the Mixed distribution and returns     @
//@             the completed distribution.                                    @
//@                                                                            @
//@     Input:  fp - input file of distributions                               @
//@             line - current input line                                      @                          
//@             num_fields - the number of fields read from the input line     @
//@             parm - name of the parameter getting the distribution          @
//@             num_dists - the number of distributions needed to fill         @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Get_Mixed_Distribution (FILE *fp, char *line, int num_fields, char *mixparm, int num_dists, Distribution_T *d) {

	char parm [80];
	double prob, val1, val2, val3, val4, val5;
	int type;
	char dist;
	double total_prob = 0.0;
	int i;

	if (num_fields != 4) {
		fprintf (stderr, "ERROR: incorrect fields for %s\n", parm);
		exit (1);
	}
	d->num_dists = num_dists;
	
	for (i = 0; i < num_dists; i++) {
	
		fgets (line, LINESIZE, fp);
		sscanf (line, "%i %s %c\n", &type, parm, &dist);

		if (strcmp (mixparm, parm) != 0) {
			fprintf (stderr, "ERROR: invalid mixed distribution %s\n", mixparm);
			exit (1);
		}

		if (toupper (dist) == 'G') {
			num_fields = sscanf (line, "%i %s %c %lf %lf %lf %lf %lf\n", &type, parm, &dist, &prob, &val1, &val2, &val3, &val4);
			total_prob += prob;
			d->density [i].prob = total_prob;
			d->density [i].d = Get_Gaussian_Distribution (num_fields-1, parm, val1, val2, val3, val4);
		}
		else if (toupper (dist) == 'Y') {
			num_fields = sscanf (line, "%i %s %c %lf %lf %lf %lf %lf %lf\n", &type, parm, &dist, &prob, &val1, &val2, &val3, &val4, &val5);
			total_prob += prob;
			d->density [i].prob = total_prob;
			d->density [i].d = Get_Gamma_Distribution (num_fields-1, parm, val1, val2, val3, val4,val5);
		}
		else if (toupper (dist) == 'U') {
			num_fields = sscanf (line, "%i %s %c %lf %lf %lf\n", &type, parm, &dist, &prob, &val1, &val2);
			total_prob += prob;
			d->density [i].prob = total_prob;
			d->density [i].d = Get_Uniform_Distribution (num_fields-1, parm, val1, val2);
		}
		else if (toupper (dist) == 'K') {
			num_fields = sscanf (line, "%i %s %c %lf %lf\n", &type, parm, &dist, &prob, &val1);
			total_prob += prob;
			d->density [i].prob = total_prob;
			d->density [i].d = Get_Constant_Distribution (num_fields-1, parm, val1);
		}
		else {
			fprintf (stderr, "ERROR: invalid distribution type %c\n", dist);
			exit (1);
		}

	} /* end for i */

	if (fabs (total_prob - 1.0) > ZERO) {
		fprintf (stderr, "ERROR: probabilities must add to 1 in mixed %s distribution\n", parm);
		exit (1);
	}
	
} /* end Get_Mixed_Distribution */


/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@                                                                            @
//@     Title: Distribution_T Get_Distribution (char *line)                    @
//@                                                                            @
//@     Action: Gets the distribution parameters from the file.                @
//@                                                                            @
//@     Input:  line - string with the distribution information.               @
//@     Output: the distribution filled in with values.                        @
//@                                                                            @
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/

void Get_Distribution (FILE *fp, char *line, Distribution_T *d) {

	int num_fields, num_dists;
	int type;
	char parm [80];
	char dist;
	double val1, val2, val3, val4,val5;

	sscanf (line, "%i %s %c\n", &type, parm, &dist);
	
	if (toupper (dist) == 'G') {
		num_fields = sscanf (line, "%i %s %c %lf %lf %lf %lf\n", &type, parm, &dist, &val1, &val2, &val3, &val4);
		d->num_dists = 1;
		d->density [0].prob = 1.0;
		d->density [0].d = Get_Gaussian_Distribution (num_fields, parm, val1, val2, val3, val4);		
	}
	else if (toupper (dist) == 'U') {
		num_fields = sscanf (line, "%i %s %c %lf %lf\n", &type, parm, &dist, &val1, &val2);
		d->num_dists = 1;
		d->density [0].prob = 1.0;
		d->density [0].d = Get_Uniform_Distribution (num_fields, parm, val1, val2);
	}
	else if (toupper (dist) == 'Y') {
		num_fields = sscanf (line, "%i %s %c %lf %lf %lf %lf %lf\n", &type, parm, &dist, &val1, &val2, &val3, &val4,&val5);
		d->num_dists = 1;
		d->density [0].prob = 1.0;
		d->density [0].d = Get_Gamma_Distribution (num_fields, parm, val1, val2, val3, val4, val5);		
	}
	else if (toupper (dist) == 'K') {
		num_fields = sscanf (line, "%i %s %c %lf\n", &type, parm, &dist, &val1);
		d->num_dists = 1;
		d->density [0].prob = 1.0;
		d->density [0].d = Get_Constant_Distribution (num_fields, parm, val1);
	}
	else if (toupper (dist) == 'M') {
		//Mixed distribution	
		num_fields = sscanf (line,"%i %s %c %i\n", &type, parm, &dist, &num_dists);
		Get_Mixed_Distribution (fp, line, num_fields, parm, num_dists, d);
	}
	else {
		fprintf (stderr, "ERROR: invalid distribution type %c\n", dist);
		exit (1);
	}
	
	if (type > num_tree_types)
		num_tree_types = type;
					
} /* end Get_Distribution */


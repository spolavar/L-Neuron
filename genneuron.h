#ifndef GENNEURON
#define GENNEURON

#define BURKEMARKS	0
#define HILLMAN		1
#define TAMORI		2

/* colors used by the lparser routine */
#define DRED		1	/* dark red */
#define LRED		2	/* light red */
#define YELLOW		3	/* yellow */
#define LGREEN		4	/* light green */
#define LBLUE		5	/* light blue */
#define DBLUE		6	/* dark blue */
#define MAGENTA		7	/* magenta */
#define DGREEN		8	/* dark green */
#define AQUA		9	/* aqua */
#define NBLUE		10	/* night blue */

#define	SOMACOLOR	LRED
#define AXONCOLOR	YELLOW
#define APICALCOLOR	LGREEN
#define BASALCOLOR	LBLUE

#define TREETYPES	9

// the soma is a serie of 45 degree lines used to draw the soma (with octagnola shape)
#define SOMA_ANGLE	45.0
//Factor to obtain the diameter of the lines of the octagon to draw a soma of the right diameter
#define OCTAGON_FACTOR	0.3827

void Initialize_Neuron ();
double Grow_Neuron (int ntype);
void Grow_Soma (double d);

#endif
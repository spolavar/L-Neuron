/* -------------------------------------------------------------------------
   LNEURON, L-System Neuronal Generator
   -------------------------------------------------------------------------
   Jeff Krichmar
   jkrichma@gmu.edu
   
   LNEURON is based on LPARSER code
   -------------------------------------------------------------------------
*/
/* -------------------------------------------------------------------------
   LPARSER, L-System Parser/Mutator
   -------------------------------------------------------------------------
   Laurens Lapre
   ljlapre@xs4all.nl
   http://www.xs4all.nl/~ljlapre/
   -------------------------------------------------------------------------
*/


/* Includes --------------------------------------------------------------- */


#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <malloc.h>
#include <time.h>
#include "burke.h"
#include "genneuron.h"
#include "hilltamo.h"

//comment about version!
#include "notes.h"
/* Basic types ------------------------------------------------------------ */


/* Simple types */
#define u8               unsigned char
#define u16              unsigned short int
#define u32              unsigned long int
#define s8               signed char
#define s16              signed short int
#define s32              signed long int
#define r32              float
#define r64              double

/* My own boolean type */
#define boolean          s16
#ifndef TRUE
#define TRUE             (s16) 1
#endif
#ifndef FALSE
#define FALSE            (s16) 0
#endif


/* Constants -------------------------------------------------------------- */


/* Max char size of filename and large string */
#define max_file         512

/* Max vectors per polygon */
#define vectors_p_poly   15

/* Max polygons per object */
#define max_p_object     400

/* Max size of the [] and {} stacks during drawing */
#define max_stack        1024L

/* Max number of the rules in the ls file */
#define rule_n           200

/* Max size of the rules in the ls file */
#define rule_s           500

/* Version id for the VOL file format */
#define VersionID        11

/* Vector indices */
#define _x               0
#define _y               1
#define _z               2

/* Some often used consts */
#define zero             (r32) 0.0
#define one              (r32) 1.0
#define half_pi          (r32) 1.570796
#define pi               (r32) 3.141592
#define two_pi           (r32) 6.283185
#define LF               ((char) 10)
#define CR               ((char) 13)
#define min_bar          "---------------------------------------------------------"
#define r32_epsilon		 (r32) 0.001


/* Bounding space for floats */
#define float_min        (r32) -1e30
#define float_max        (r32)  1e30


/* Array and vector types ------------------------------------------------- */
char	*object_s; 


/* A large string */
typedef char            string_file[max_file];

/* Vector arrays */
typedef r32             vector[3];
typedef vector          vectors_4[4];
typedef vector          vectors_max[vectors_p_poly];
typedef r32             float_array[max_p_object + 1];

/* Polygon arrays */
typedef s16             polygon_type[4];
typedef polygon_type    polygon_array[max_p_object + 1];


/* Intrinsics ------------------------------------------------------------- */


#define Abs_s16(A)       ((s16) abs((A)))
#define Abs_s32(A)       ((s32) abs((A)))
#define Abs_r32(A)       ((r32) fabs((A)))
#define Abs_r64(A)       ((r64) fabs((A)))

#define MIN(A,B)         (((A) < (B)) ? (A) : (B))
#define MAX(A,B)         (((A) > (B)) ? (A) : (B))

#define Vector_copy_max_r32(n,A,B) memcpy(((void *) (B)), ((void *) (A)), n*12L)
#define Vector_copy_r32(A,B)       memcpy(((void *) (B)), ((void *) (A)),   12L)
#define Mem_copy(A,B,C)            memcpy(((void *) (B)), ((void *) (A)), (C))
#define Mem_clear(A,B,C)           memset(((void *) (A)), (C), (B))


/* Vector utils and inlines ----------------------------------------------- */


#define Clip_low(a,b)   a = ((a) < (b)) ? (b) : (a)
#define Clip_high(a,b)  a = ((a) > (b)) ? (b) : (a)
#define Clamp(a,b,c)    a = ((a) < (b)) ? (b) : ((a) > (c)) ? (c) : (a)
#define Wrap(a,b,c)     a = ((a) < (b)) ? (a) + (c) : ((a) > (c)) ? (a) - (c) : (a)
#define Lerp(a,b,c)     ((b) + (((c) - (b)) * (a)))
#define Swap(a,b)       {(a) ^= (b); (b) ^= (a); (a) ^= (b);}

#define Vector_length(A)    ((r32) sqrt(   (r32)(A[_x] * A[_x])\
                                         + (r32)(A[_y] * A[_y])\
                                         + (r32)(A[_z] * A[_z]) ))

#define Vector_equal(A,B)   ((A[_x] == B[_x]) && (A[_y] == B[_y]) && (A[_z] == B[_z]))
#define Scalar_product(A,B) (A[_x] * B[_x] + A[_y] * B[_y] + A[_z] * B[_z])

#define Vector_make(A,a,b,c) {\
  A[_x] = a;\
  A[_y] = b;\
  A[_z] = c;\
}

#define Vector_break(A,a,b,c) {\
  a = A[_x];\
  b = A[_y];\
  c = A[_z];\
}

#define Vector_lerp(a,A,B) {\
  A[_x] = Lerp(a, A[_x], B[_x]);\
  A[_y] = Lerp(a, A[_y], B[_y]);\
  A[_z] = Lerp(a, A[_z], B[_z]);\
}

#define Vector_normalize(A)\
{ r32 Dist = (r32) 1.0 / Vector_length(A);\
\
  A[_x] *= Dist;\
  A[_y] *= Dist;\
  A[_z] *= Dist;\
}

#define Vector_copy(A,B) {\
  B[_x] = A[_x];\
  B[_y] = A[_y];\
  B[_z] = A[_z];\
}

#define Vector_product(A,B,C) {\
  C[_x] = A[_y] * B[_z] - A[_z] * B[_y];\
  C[_y] = A[_z] * B[_x] - A[_x] * B[_z];\
  C[_z] = A[_x] * B[_y] - A[_y] * B[_x];\
}

#define Vector_min(A,B,C) {\
  C[_x] = A[_x] - B[_x];\
  C[_y] = A[_y] - B[_y];\
  C[_z] = A[_z] - B[_z];\
}

#define Vector_plus(A,B,C) {\
  C[_x] = A[_x] + B[_x];\
  C[_y] = A[_y] + B[_y];\
  C[_z] = A[_z] + B[_z];\
}

#define Vector_dec(A,B) {\
  A[_x] -= B[_x];\
  A[_y] -= B[_y];\
  A[_z] -= B[_z];\
}

#define Vector_neg(A) {\
  A[_x] = (-A[_x]);\
  A[_y] = (-A[_y]);\
  A[_z] = (-A[_z]);\
}

#define Vector_inc(A,B) {\
  A[_x] += B[_x];\
  A[_y] += B[_y];\
  A[_z] += B[_z];\
}

#define Vector_plus_fac(A,B,t,C) {\
  C[_x] = A[_x] + (t) * B[_x];\
  C[_y] = A[_y] + (t) * B[_y];\
  C[_z] = A[_z] + (t) * B[_z];\
}

#define Vector_plus_fac2(A,B,b,C,c,D) {\
  D[_x] = A[_x] + (b) * B[_x] + (c) * C[_x];\
  D[_y] = A[_y] + (b) * B[_y] + (c) * C[_y];\
  D[_z] = A[_z] + (b) * B[_z] + (c) * C[_z];\
}

#define Vector_combine(A,a,B,b,C) {\
  C[_x] = (a) * A[_x] + (b) * B[_x];\
  C[_y] = (a) * A[_y] + (b) * B[_y];\
  C[_z] = (a) * A[_z] + (b) * B[_z];\
}

#define Vector_add(A,d) {\
  A[_x] += d;\
  A[_y] += d;\
  A[_z] += d;\
}

#define Vector_sub(A,d) {\
  A[_x] -= d;\
  A[_y] -= d;\
  A[_z] -= d;\
}

#define Vector_div(A,d) {\
  A[_x] /= d;\
  A[_y] /= d;\
  A[_z] /= d;\
}

#define Vector_mul(A,d) {\
  A[_x] *= d;\
  A[_y] *= d;\
  A[_z] *= d;\
}


/* Vector procs ----------------------------------------------------------- */


static vector           M1, M2, M3;         /* The current movetransform matrix */


r32
Do_angle(r32 x1, r32 y1, r32 x2, r32 y2)
{                                           /* Calculate the angle between
                                             * x-axis and x1,y1 -> x2,y2. It
                                             * can handle all kinds of weird
                                             * exceptions */
    r32                     temp, x, y;

    x = x2 - x1;
    y = y2 - y1;

    if (x == zero) {
        if (y < zero)
            temp = -half_pi;
        else
            temp = half_pi;
    } else {
		//SCR cast
        temp = (float) atan(y / x);
        if (x < zero) {
            if (y < zero)
                temp = -pi + temp;
            else
                temp = pi + temp;
        }
    }

    if (Abs_r32(temp) < (r32) 0.0001)
        temp = zero;

    if (temp < zero)
        temp += two_pi;
    else if (temp > two_pi)
        temp -= two_pi;

    return temp;
}


void
Move_transform(vector v)
{                                           /* Transform the vector according
                                             * to the current movetransform
                                             * matrix */
    vector                  t;

    Vector_copy_r32(v, t);
    v[_x] = Scalar_product(M1, t);
    v[_y] = Scalar_product(M2, t);
    v[_z] = Scalar_product(M3, t);
}


void
Set_move_transform(r32 a, vector no)
{                                           /* Set a movetransformation matrix
                                             * based on an angle rotation of
                                             * 'a' around the vector 'no' */
    r32                     n11, n22, n33, nxy, nxz, nyz, sina, cosa;

    cosa = (r32) cos(a);
    sina = (r32) sin(a);

    n11 = no[_x] * no[_x];
    n22 = no[_y] * no[_y];
    n33 = no[_z] * no[_z];

    nxy = no[_x] * no[_y];
    nxz = no[_x] * no[_z];
    nyz = no[_y] * no[_z];

    M1[_x] = n11 + (one - n11) * cosa;
    M1[_y] = nxy * (one - cosa) - no[_z] * sina;
    M1[_z] = nxz * (one - cosa) + no[_y] * sina;

    M2[_x] = nxy * (one - cosa) + no[_z] * sina;
    M2[_y] = n22 + (one - n22) * cosa;
    M2[_z] = nyz * (one - cosa) - no[_x] * sina;

    M3[_x] = nxz * (one - cosa) - no[_y] * sina;
    M3[_y] = nyz * (one - cosa) + no[_x] * sina;
    M3[_z] = n33 + (one - n33) * cosa;
}


/* File and conio procs --------------------------------------------------- */


static boolean          native_mode = TRUE;


void
Set_lowhigh(boolean b)
{                                           /* TRUE native Intel mode Low-High,
                                             * FALSE High-Low */
    native_mode = b;
}


void
User_error(char *s,...)
{                                           /* Displays and error messages and
                                             * exits the program */
    string_file             buf;
    va_list                 args;

    va_start(args, s);
    vsprintf(buf, s, args);
    va_end(args);

    fprintf(stdout, "\n\nError: %s\n\n", buf);
    fflush(stdout);

    exit(EXIT_FAILURE);
}


void
Message(char *s,...)
{                                           /* Sends a message to the output
                                             * stream */
    string_file             buf;
    va_list                 args;

    va_start(args, s);
    vsprintf(buf, s, args);
    va_end(args);

    fprintf(stdout, "%s", buf);
    fflush(stdout);
    fflush(stdout);
}


void
Fget_string(FILE * f, char *s)
{                                           /* Get a string from a file. Used
                                             * for parsing the LS files */
    s16                     i = 0;
    char                    c;

    s[0] = '\0';
    for (;;) {
        c = (char) getc(f);
        if (feof(f))
            return;
        if (c == '\r')
            continue;
        if (c == LF)
            break;
        s[i++] = c;
    }
    s[i] = '\0';
}


void
Fget_bin_r32(FILE * f, r32 *val)
{                                           /* Get a r32 value, check for order */
    s32                     temp, ta, tb, tc, td;
    r32                    *tempr = ((r32 *) ((void *) &temp));

    ta = (s32) getc(f);
    tb = (s32) getc(f);
    tc = (s32) getc(f);
    td = (s32) getc(f);

    if (native_mode) {
        temp = td << 8;
        temp += tc;
        temp = temp << 8;
        temp += tb;
        temp = temp << 8;
        temp += ta;
    } else {
        temp = ta << 8;
        temp += tb;
        temp = temp << 8;
        temp += tc;
        temp = temp << 8;
        temp += td;
    }

    *val = *tempr;
}


void
Fget_bin_s16(FILE * f, s16 *val)
{                                           /* Get a s16 value, check for order */
    s32                     ta, tb;

    ta = (s32) getc(f);
    tb = (s32) getc(f);

    if (native_mode) {
        *val = (s16) tb << 8;
        *val += (s16) ta;
    } else {
        *val = (s16) ta << 8;
        *val += (s16) tb;
    }
}


void
Fget_bin_s32(FILE * f, s32 *val)
{                                           /* Get a s32 value, check for order */
    s32                     ta, tb, tc, td;

    ta = (s32) getc(f);
    tb = (s32) getc(f);
    tc = (s32) getc(f);
    td = (s32) getc(f);

    if (native_mode) {
        *val = td << 8;
        *val += tc;
        *val = *val << 8;
        *val += tb;
        *val = *val << 8;
        *val += ta;
    } else {
        *val = ta << 8;
        *val += tb;
        *val = *val << 8;
        *val += tc;
        *val = *val << 8;
        *val += td;
    }
}


void
Fget_bin_u16(FILE * f, u16 *val)
{                                           /* Get a u16 value, check for order */
    s32                     ta, tb;

    ta = (s32) getc(f);
    tb = (s32) getc(f);

    if (native_mode) {
        *val = (u16) tb << 8;
        *val += (u16) ta;
    } else {
        *val = (u16) ta << 8;
        *val += (u16) tb;
    }

}


void
Fget_bin_s8(FILE * f, s8 *val)
{
    *val = (s8) getc(f);
}


void
Fget_bin_u8(FILE * f, u8 *val)
{
    *val = (u8) getc(f);
}


void
Fget_bin_u32(FILE * f, u32 *val)
{                                           /* Get a u32 value, check for order */
    s32                     ta, tb, tc, td;

    ta = (s32) getc(f);
    tb = (s32) getc(f);
    tc = (s32) getc(f);
    td = (s32) getc(f);

    if (native_mode) {
        *val = td << 8;
        *val += tc;
        *val = *val << 8;
        *val += tb;
        *val = *val << 8;
        *val += ta;
    } else {
        *val = ta << 8;
        *val += tb;
        *val = *val << 8;
        *val += tc;
        *val = *val << 8;
        *val += td;
    }
}


void
Fput_bin_u8(FILE * f, u8 r)
{
    if (fwrite(&r, sizeof(u8), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_s8(FILE * f, s8 r)
{
    if (fwrite(&r, sizeof(s8), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_r32(FILE * f, r32 r)
{
    if (fwrite(&r, sizeof(r32), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_s16(FILE * f, s16 r)
{
    if (fwrite(&r, sizeof(s16), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_u16(FILE * f, u16 r)
{
    if (fwrite(&r, sizeof(u16), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_s32(FILE * f, s32 r)
{
    if (fwrite(&r, sizeof(s32), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


void
Fput_bin_u32(FILE * f, u32 r)
{
    if (fwrite(&r, sizeof(u32), 1, f) != 1)
        User_error("Can't continue writing outputfile");
}


/* File buffer procs ------------------------------------------------------ */


/* Sinces there is a lot of fileio we use large buffers */
#define f_buffer_size  30L * 1024L
static char            *f_buffer[3] = {NULL, NULL, NULL};


void
Buffer_IO(FILE * f, s16 i)
{                                           /* Attach filebuffer i to open file */
    setvbuf(f, f_buffer[i], _IOFBF, f_buffer_size);
}


void
Init_file_buf(s16 i)
{                                           /* Init filebuffer i */
    if (f_buffer[i] != NULL)
        return;
    f_buffer[i] = (char *) malloc(f_buffer_size);
    if (f_buffer[i] == NULL)
        User_error("Not enough memory to allocate file buffer");
}


/* Feedback percentage counter -------------------------------------------- */


/* These vars are used to calculate the percentages counters for feedback */
static r32              bar_fac2 = zero;
static s16              old_bar2 = 0;
static u32              bar_max2 = 0;


void
Process_start2(u32 max)
{                                           /* Start bar 2 with the maximum
                                             * value it's going to get */
	//SCR cast
    bar_fac2 = (float) 100.0 / (r32) max;
    bar_max2 = max;
}


void
Process_update2(u32 now)
{                                           /* Update the percentage counter
                                             * when needed */
    s16                     bar = (s16) (bar_fac2 * (r32) now);

    if (bar != old_bar2) {
        old_bar2 = bar;
        Message("\r%3d%%\r", bar);
    }
}


void
Process_end2(void)
{                                           /* Close bar */
    Message("\r    \r");
}


/* Comline procs ---------------------------------------------------------- */


#define OPTCHAR '-'

static char             empty[] = "";
static char            *optarg = empty;
static int              optind = 1, opterr = 1;
static char             opts[150] = "";     /* option string */
static s16              s_argc = 0;         /* pointers to comline */
static char           **s_argv = NULL;


static int
lgetopt(int argc, char *argv[], const char *optstring)
{                                           /* Taken from a source lib
                                             * somewhere */
    static char            *in_pointer = empty;
    char                   *find;

    if (!*in_pointer) {
        if ((!argv[optind]) || (optind >= argc) || (argv[optind][0] != OPTCHAR))
            return -1;
        in_pointer = argv[optind];
        in_pointer++;
        optind++;
        if (*in_pointer == OPTCHAR) {
            return -1;
        }
    };

    if (*in_pointer == '\0') {
        return 0;
    };

    find = strchr(optstring, *in_pointer);
    if (find == NULL) {
        if (opterr)
            User_error("Option -%c not known", *in_pointer);
        in_pointer = empty;
        return '?';
    };

    if (*(find + 1) == ':') {
        if (!*(in_pointer + 1)) {
            if (optind >= argc) {
                if (opterr)
                    User_error("No argument for option -%c", *in_pointer);
                optarg = empty;
            } else {
                if (argv[optind][0] == OPTCHAR) {
                    if (opterr)
                        User_error("No argument for option -%c but found %s instead", *in_pointer, argv[optind]);
                }
                optarg = argv[optind++];
            }
        } else {
            optarg = ++in_pointer;
        }
        in_pointer = empty;
    } else {
        optarg = empty;
        in_pointer++;
    }

    return *find;
}


void
Get_comline_opt(char *c, boolean * found, char *result)
{                                           /* Check if comline option 'c' has
                                             * been used and return parameter
                                             * if any */
    int                     f, argc;
    char                  **argv;

    argc = s_argc;
    argv = s_argv;
    optind = 1;

    while ((f = lgetopt(argc, argv, opts)) != -1) {

        if (f == *c) {
            strcpy(result, optarg);
            *found = TRUE;
            return;
        }
    };

    *found = FALSE;
    strcpy(result, "");
}


void
Get_comline_filename(char *c)
{                                           /* Get the filename argument from
                                             * the comline */
    int                     argc;
    char                  **argv;

    argc = s_argc;
    argv = s_argv;
    optind = 1;

    while (lgetopt(argc, argv, opts) != -1);

    if (optind == argc)
        User_error("Ran out of arguments before finding file name");

    strcpy(c, argv[optind]);
}


void
Get_comline_progname(char *c)
{                                           /* Get the program name from the
                                             * comline */
    strcpy(c, s_argv[0]);
}


/* Main lparser vars ------------------------------------------------------ */


/* Settings stack used for solving [] references */
typedef struct s_rec {
    vector                  pos;            /* position in 3space of turtle
                                             * origin */
    vector                  fow;            /* forward direction */

    vector                  lef;            /* left direction */
    vector                  upp;            /* up direction */
    vector                  last;           /* last position used for
                                             * connecting cylinders */
    vector                  last_v[9];      /* last vertices of object used for
                                             * connecting cylinders */
    r32                     dis;            /* value of F distance */
    r32                     ang;            /* value of basic angle */
    r32                     thick;          /* value of thickness */
    r32                     dis2;           /* value of Z distance */
    r32                     tr;             /* trope value */
    s16                     col;            /* current color */
    s16                     last_col;       /* color of last object */
    
    /* changes needed for GENESIS output */
    /* JLK 2/10/99                       */
    s16						parent;			/* parent compartment number */
} s_rec;

/* Polygon stack used for solving {} references */
typedef struct p_rec {
    s16                     count;          /* number of vertices */
    vector                 *ver;            /* vertex store */
} p_rec;

/* Flags */
static boolean          trope_set = FALSE;  /* see at comline scannign */
static boolean          rand_set = FALSE;
static boolean          user_form = FALSE;
static boolean          closed_form = FALSE;
static boolean          pov_form = FALSE;
static boolean          pov_form2 = FALSE;
static boolean          pov_form3 = FALSE;
static boolean          blb_form = FALSE;
static boolean          inc_out = FALSE;
static boolean          dxf1 = FALSE;
static boolean          dxf2 = FALSE;
static boolean          dxf3 = FALSE;
static boolean          vrml = FALSE;
static boolean          turner = FALSE;
static boolean          turner_tree = FALSE;
static boolean          growing = FALSE;    /* real is used for recursion level */
static boolean          last_recur = FALSE; /* processing the last recursion
                                             * step */
static boolean          seed_set = FALSE;  /* JLK: random seed set switch -s */

/* Rule and strings storage */
static char             rule[rule_n][rule_s];   /* the rule store */
static char				*otemp;   /* the two strings used for
                                             * building the production string */
static s16              size[rule_n];       /* size of the rule in the rule
                                             * store */
static boolean          mark[rule_n];       /* marked rules need special
                                             * processing when growing shapes
                                             * are active */

/* Init vars */
//SCR cast
static r32              zmin = (float) 1e30, thick, min_thick =(float) zero, rand_amount =(float) zero;
static r32              trope_amount = zero;
static u32              polcount = 0;
static u32              poly_limit = 500000L, max_string;
static s16              num = 0, col = 2, lev, last_col = 0;
static r32              dis, ang, dis2, tr = (float) 0.2;
static vector           sky = {0.0, 0.0, 1.0}, trope;
static vector           last = {1.0, 1.0, 1.0}, last_v[9];
static r32              recursion, fraction;
static vector           axis_x, axis_y, axis_z;
static r32				soma_d; /* JLK - diameter of the neuron's soma */
static int				seed_amount; /* JLK - diameter of the neuron's soma */

/* Stacks [] and {} */
static s_rec           *stack, org, save;
static s16              scount = 0;
static p_rec           *pstack;
static s16              pscount = 0;

/* Current active transform matrix for drawing */
static vector           C1, C2, C3;

/* Var for rnd */
static r32              rm = (r32) 1.0 / (r32) RAND_MAX;

/* Ouput files */
static FILE            *volume_file = NULL; /* the basic open geometry file */
static FILE            *vf[8];              /* the 8 files when writing
                                             * multiple povray blob files */

/* Object stores */
static polygon_array    poly_store;         /* the store where polygons
                                             * accumulate before saved to disc */
static vector           ver[max_p_object];  /* the store where vertices
                                             * accumulate */

/* Storage of a loaded shape for the -X option */
static char             x_name[max_file];   /* filename of VOL file */
static vector           form_c[max_p_object];   /* vertices */
static polygon_array    form_s;             /* polygons */
static s16              form_ver, form_pol; /* vertices and polygon counts */


/* Check for weird polygons ----------------------------------------------- */


 /*
  * Sometimes polygons and/or vertices end up containing floating pint
  * exception like NAN etc. These routines find these problems. They can create
  * havoc on input parsers which expect the geometry to be flawless. Typical,
  * normalization routines blow up on NAN in vectors.
  */


/* IEEE coded floating point exceptions */
static u32              fp_exp1 = 0x7f800000L;
static u32              fp_exp2 = 0xff800000L;
static u32              fp_exp3 = 0xffc00000L;
static u32              fp_exp4 = 0x7fc00000L;


static                  boolean
Bad_vertex(r32 x, r32 y, r32 z)
{                                           /* Does this vertex contain a
                                             * floation point exception ? */
    union {
        r32                     f;
        u32                     i;
    }                       u;

    u.f = x;
    if ((u.i == fp_exp1) ||
            (u.i == fp_exp2) ||
            (u.i == fp_exp3) ||
            (u.i == fp_exp4)) {
        return TRUE;
    };

    u.f = y;
    if ((u.i == fp_exp1) ||
            (u.i == fp_exp2) ||
            (u.i == fp_exp3) ||
            (u.i == fp_exp4)) {
        return TRUE;
    };

    u.f = z;
    if ((u.i == fp_exp1) ||
            (u.i == fp_exp2) ||
            (u.i == fp_exp3) ||
            (u.i == fp_exp4)) {
        return TRUE;
    };

    return FALSE;
}


static                  boolean
Invalid_polygon(s16 p)
{                                           /* Can a normal be created on this
                                             * polygon ? */
    vector                  X, Y, Z, N;
    s16                     i;
    r32                     D, x, y, z;

    for (i = 0; i < 4; i++) {
        x = ver[poly_store[p][i]][_x];
        y = ver[poly_store[p][i]][_y];
        z = ver[poly_store[p][i]][_z];
        if (Bad_vertex(x, y, z))
            return TRUE;
    }

    for (i = 0; i < 3; i++) {
        X[i] = ver[poly_store[p][i]][_x];
        Y[i] = ver[poly_store[p][i]][_y];
        Z[i] = ver[poly_store[p][i]][_z];
    }

    N[_x] = Y[_x] * (Z[_y] - Z[_z]) + Y[_y] * (Z[_z] - Z[_x]) + Y[_z] * (Z[_x] - Z[_y]);
    N[_y] = ((-X[_x])) * (Z[_y] - Z[_z]) - X[_y] * (Z[_z] - Z[_x]) - X[_z] * (Z[_x] - Z[_y]);
    N[_z] = X[_x] * (Y[_y] - Y[_z]) + X[_y] * (Y[_z] - Y[_x]) + X[_z] * (Y[_x] - Y[_y]);

    D = N[_x] * N[_x] + N[_y] * N[_y] + N[_z] * N[_z];

    if (D <= (r32) 0.0001)
        return TRUE;
    else
        return FALSE;
}


/* Output data file procs ------------------------------------------------- */


static void
Open_datafile(void)
{                                           /* Open and setup the different
                                             * output files depending on the
                                             * flags */

    char                    S[max_file] = "";
    s16                     i;

    if (pov_form || pov_form2) {
        strcpy(S, inc_out ? "output.inc" : "output.pov");
        Message("Pov file       : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);
        return;

    } else if (pov_form3) {
        for (i = 0; i < 8; i++) {
            sprintf(S, inc_out ? "output%d.inc" : "output%d.pov", i);
            Message("Pov file       : %s\n", S);
            vf[i] = fopen(S, "wt");
            if (!vf[i])
                User_error("Cannot open file [%s]", S);
            fprintf(vf[i], "component 1.0 1.0 <0, 0, 0>\n");
        }
        return;

    } else if (blb_form) {
        strcpy(S, "output.blb");
        Message("Blob file      : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);
        fprintf(volume_file, "[blob]\nThreshold = 0.5\n");
        return;

    } else if (dxf1) {
        strcpy(S, "output.dxf");
        Message("Dxf file       : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

        fprintf(volume_file, "999\nL-System Parser/Mutator\n");
        fprintf(volume_file, "999\nPolyline Polyface Meshes\n");

        if (user_form) {                   /* in this case build a block
                                            * section, include the loaded shape
                                            * and use only block inserts in the
                                            * entities section of the dxf file */

            fprintf(volume_file, "0\nSECTION\n2\nTABLES\n0\nTABLE\n2\nAPPID\n70\n4\n0\nAPPID\n2\nLPARSER\n70\n0\n0\nENDTAB\n0\nENDSEC\n");

            fprintf(volume_file, "0\nSECTION\n2\nBLOCKS\n0\nBLOCK\n8\n0\n2\nBLOCK\n70\n0\n");
            fprintf(volume_file, "10\n0.0\n20\n0.0\n30\n0.0\n3\nBLOCK\n");

            fprintf(volume_file, "0\nPOLYLINE\n66\n1\n8\n0\n62\n0\n70\n64\n");
            fprintf(volume_file, "1001\nLPARSER\n1071\n18500\n1070\n11003\n1000\n<byblock>\n1070\n10999\n");

            for (i = 1; i <= form_ver; i++) {
                fprintf(volume_file, "0\nVERTEX\n8\n0\n62\n0\n");
                fprintf(volume_file, "10\n%g\n", form_c[i][_x]);
                fprintf(volume_file, "20\n%g\n", form_c[i][_y]);
                fprintf(volume_file, "30\n%g\n", form_c[i][_z]);
                fprintf(volume_file, "70\n192\n");
            }

            for (i = 1; i <= form_pol; i++) {
                fprintf(volume_file, "0\nVERTEX\n8\n0\n62\n0\n");
                fprintf(volume_file, "10\n0\n20\n0\n30\n0\n70\n128\n");
                fprintf(volume_file, "71\n%d\n", form_s[i][0]);
                fprintf(volume_file, "72\n%d\n", form_s[i][1]);
                fprintf(volume_file, "73\n%d\n", form_s[i][2]);
                fprintf(volume_file, "74\n%d\n", form_s[i][3]);
            };

            fprintf(volume_file, "0\nSEQEND\n8\n0\n0\nENDBLK\n8\n0\n0\nENDSEC\n");
        };

        fprintf(volume_file, "0\nSECTION\n2\nENTITIES\n");

    } else if (dxf2) {
        strcpy(S, "output.dxf");
        Message("Dxf file       : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

        fprintf(volume_file, "999\nL-System Parser/Mutator\n");
        fprintf(volume_file, "999\n3d Faces List\n");

        fprintf(volume_file, "0\nSECTION\n2\nENTITIES\n");

    } else if (dxf3) {
        strcpy(S, "output.raw");
        Message("Raw file       : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

    } else if (vrml) {
        strcpy(S, "output.wrl");
        Message("VRML file      : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

        fprintf(volume_file, "#VRML V1.0 ascii\n");     /* vrml header */
        fprintf(volume_file, "\nSeparator {\n");
        fprintf(volume_file, "\tShapeHints {\n");
        fprintf(volume_file, "\t\tvertexOrdering UNKNOWN_ORDERING\n");
        fprintf(volume_file, "\t\tshapeType UNKNOWN_SHAPE_TYPE\n");
        fprintf(volume_file, "\t\tfaceType CONVEX\n");
        fprintf(volume_file, "\t\tcreaseAngle 0.5\n");
        fprintf(volume_file, "\t}\n");
        fprintf(volume_file, "\tDirectionalLight {\n");
        fprintf(volume_file, "\t\tdirection -0.3 -0.6 -0.9\n");
        fprintf(volume_file, "\t}\n");

    } else if (turner) {
        strcpy(S, "output.swc");
        Message("Turner Southampton Archive file      : %s\n", S);
        volume_file = fopen(S, "wt");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

    } else {                               /* default output in Lviewer VOL
                                            * format */
        strcpy(S, "output.vol");
        Message("Datafile       : %s\n", S);
        volume_file = fopen(S, "wb");
        if (!volume_file)
            User_error("Cannot open file [%s]", S);
        Buffer_IO(volume_file, 0);

        Fput_bin_u8(volume_file, (u8) VersionID);       /* VOL header */
        Fput_bin_r32(volume_file, (r32) 45.0);
        Fput_bin_r32(volume_file, (r32) 45.0);
        Fput_bin_r32(volume_file, (r32) 90.0);
        Fput_bin_r32(volume_file, (r32) 45.0);
        Fput_bin_r32(volume_file, (r32) 0.0);
        Fput_bin_s16(volume_file, 0);
        Fput_bin_s16(volume_file, 0);
        Fput_bin_s16(volume_file, 100);
        Fput_bin_s16(volume_file, 3000);
        Fput_bin_u16(volume_file, (u16) 0);/* stream format */
    }
}


static void
Close_datafile(void)
{                                           /* Close the different open output
                                             * files */
    s16                     i;

    if (pov_form3) {
        Message("Objects        : %ld\n", polcount);
        for (i = 0; i < 8; i++)
            fclose(vf[i]);
        return;
    };

    if (dxf1 || dxf2)
        fprintf(volume_file, "0\nENDSEC\n0\nEOF\n");

    if (vrml)
        fprintf(volume_file, "}\n");

    Message("Objects        : %ld\n", polcount);
    fclose(volume_file);
}


static void
Save_object(s16 vertices, s16 polygons, s16 color)
{                                           /* Save an object from store to
                                             * disc */
    s32                     t, i, max;

    if (pov_form2 || pov_form3 || blb_form)/* in these case no 'real' geometry
                                            * is saved */
        return;

    polcount += polygons;

    if (pov_form) {
        for (t = 1; t <= polygons; t++) {

			//SCR cast
            if (Invalid_polygon((short)t))
                continue;

            if (poly_store[t][2] == poly_store[t][3]) { /* 3 vertex triangle */
                fprintf(volume_file, "object{triangle{");
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][0]][_x],
                        ver[poly_store[t][0]][_z], ver[poly_store[t][0]][_y]);
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][1]][_x],
                        ver[poly_store[t][1]][_z], ver[poly_store[t][1]][_y]);
                fprintf(volume_file, "<%g, %g, %g>}", ver[poly_store[t][2]][_x],
                        ver[poly_store[t][2]][_z], ver[poly_store[t][2]][_y]);
                fprintf(volume_file, "finish{t_leaf} pigment{color col_%d}}\n", color % 16);

            } else {                       /* 4 vertex polygon == 2x triangle */
                fprintf(volume_file, "object{triangle{");
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][0]][_x],
                        ver[poly_store[t][0]][_z], ver[poly_store[t][0]][_y]);
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][1]][_x],
                        ver[poly_store[t][1]][_z], ver[poly_store[t][1]][_y]);
                fprintf(volume_file, "<%g, %g, %g>}", ver[poly_store[t][2]][_x],
                        ver[poly_store[t][2]][_z], ver[poly_store[t][2]][_y]);
                fprintf(volume_file, "finish{t_leaf} pigment{color col_%d}}\n", color % 16);
                fprintf(volume_file, "object{triangle{");
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][2]][_x],
                        ver[poly_store[t][2]][_z], ver[poly_store[t][2]][_y]);
                fprintf(volume_file, "<%g, %g, %g>", ver[poly_store[t][3]][_x],
                        ver[poly_store[t][3]][_z], ver[poly_store[t][3]][_y]);
                fprintf(volume_file, "<%g, %g, %g>}", ver[poly_store[t][0]][_x],
                        ver[poly_store[t][0]][_z], ver[poly_store[t][0]][_y]);
                fprintf(volume_file, "finish{t_leaf} pigment{color col_%d}}\n", color % 16);
            }
        }

    } else if (dxf1) {                     /* dxf 3d mesh object */
        fprintf(volume_file, "0\nPOLYLINE\n66\n1\n8\n%d\n70\n64\n", color);

        for (i = 1; i <= vertices; i++) {
            fprintf(volume_file, "0\nVERTEX\n8\n%d\n", color);
            fprintf(volume_file, "10\n%g\n", ver[i][_x]);
            fprintf(volume_file, "20\n%g\n", ver[i][_y]);
            fprintf(volume_file, "30\n%g\n", ver[i][_z]);
            fprintf(volume_file, "70\n192\n");
        }

        for (i = 1; i <= polygons; i++) {
			//SCR cast
            if (Invalid_polygon((short)i))
                continue;
            fprintf(volume_file, "0\nVERTEX\n8\n%d\n", color);
            fprintf(volume_file, "10\n0\n20\n0\n30\n0\n70\n128\n");
            fprintf(volume_file, "71\n%d\n", poly_store[i][0]);
            fprintf(volume_file, "72\n%d\n", poly_store[i][1]);
            fprintf(volume_file, "73\n%d\n", poly_store[i][2]);
            fprintf(volume_file, "74\n%d\n", poly_store[i][3]);
        }

        fprintf(volume_file, "0\nSEQEND\n8\n%d\n", color);

    } else if (dxf2) {                     /* dxf 3dface object */
        for (i = 1; i <= polygons; i++) {
			//SCR cast
            if (Invalid_polygon((short)i))
                continue;

            if (poly_store[i][2] == poly_store[i][3]) { /* 3 vertex face */
                fprintf(volume_file, "0\n3DFACE\n8\n%d\n", color);
                fprintf(volume_file, "10\n%g\n", ver[poly_store[i][0]][_x]);
                fprintf(volume_file, "20\n%g\n", ver[poly_store[i][0]][_y]);
                fprintf(volume_file, "30\n%g\n", ver[poly_store[i][0]][_z]);
                fprintf(volume_file, "11\n%g\n", ver[poly_store[i][1]][_x]);
                fprintf(volume_file, "21\n%g\n", ver[poly_store[i][1]][_y]);
                fprintf(volume_file, "31\n%g\n", ver[poly_store[i][1]][_z]);
                fprintf(volume_file, "12\n%g\n", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "22\n%g\n", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "32\n%g\n", ver[poly_store[i][2]][_z]);
                fprintf(volume_file, "13\n%g\n", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "23\n%g\n", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "33\n%g\n", ver[poly_store[i][2]][_z]);

            } else {                       /* 4 vertex face */
                fprintf(volume_file, "0\n3DFACE\n8\n%d\n", color);
                fprintf(volume_file, "10\n%g\n", ver[poly_store[i][0]][_x]);
                fprintf(volume_file, "20\n%g\n", ver[poly_store[i][0]][_y]);
                fprintf(volume_file, "30\n%g\n", ver[poly_store[i][0]][_z]);
                fprintf(volume_file, "11\n%g\n", ver[poly_store[i][1]][_x]);
                fprintf(volume_file, "21\n%g\n", ver[poly_store[i][1]][_y]);
                fprintf(volume_file, "31\n%g\n", ver[poly_store[i][1]][_z]);
                fprintf(volume_file, "12\n%g\n", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "22\n%g\n", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "32\n%g\n", ver[poly_store[i][2]][_z]);
                fprintf(volume_file, "13\n%g\n", ver[poly_store[i][3]][_x]);
                fprintf(volume_file, "23\n%g\n", ver[poly_store[i][3]][_y]);
                fprintf(volume_file, "33\n%g\n", ver[poly_store[i][3]][_z]);
            }
        }

    } else if (vrml) {
        fprintf(volume_file, "\tSeparator {\n");
        fprintf(volume_file, "\t\tMaterial {\n");
        fprintf(volume_file, "\t\t\tdiffuseColor ");

        switch (color) {                   /* translate colors from lparser to
                                            * RGB values */
          case 1:
            fprintf(volume_file, "%f %f %f", 0.3, 0.3, 0.3);
            break;
          case 2:
            fprintf(volume_file, "%f %f %f", 0.8, 0.4, 0.4);
            break;
          case 3:
            fprintf(volume_file, "%f %f %f", 0.8, 0.8, 0.4);
            break;
          case 4:
            fprintf(volume_file, "%f %f %f", 0.4, 0.8, 0.4);
            break;
          case 5:
            fprintf(volume_file, "%f %f %f", 0.4, 0.8, 0.8);
            break;
          case 6:
            fprintf(volume_file, "%f %f %f", 0.4, 0.4, 0.8);
            break;
          case 7:
            fprintf(volume_file, "%f %f %f", 0.8, 0.4, 0.8);
            break;
          case 8:
            fprintf(volume_file, "%f %f %f", 0.2, 0.5, 0.2);
            break;
          case 9:
            fprintf(volume_file, "%f %f %f", 0.2, 0.5, 0.5);
            break;
          case 10:
            fprintf(volume_file, "%f %f %f", 0.2, 0.2, 0.5);
            break;
          case 11:
            fprintf(volume_file, "%f %f %f", 0.5, 0.2, 0.5);
            break;
          case 12:
            fprintf(volume_file, "%f %f %f", 0.6, 0.2, 0.2);
            break;
          case 13:
            fprintf(volume_file, "%f %f %f", 0.5, 0.5, 0.5);
            break;
          case 14:
            fprintf(volume_file, "%f %f %f", 0.7, 0.7, 0.7);
            break;
          case 15:
            fprintf(volume_file, "%f %f %f", 0.9, 0.9, 0.9);
            break;
          default:
            fprintf(volume_file, "%f %f %f", 0.5, 0.5, 0.5);
            break;
        };

        fprintf(volume_file, "\n\t\t}\n");

    /* Write vertices */
        fprintf(volume_file, "\t\tCoordinate3 {\n");
        fprintf(volume_file, "\t\t\tpoint [\n");
        for (t = 1; t < vertices; t++)
            fprintf(volume_file, "\t\t\t\t%f %f %f,\n", ver[t][_x], ver[t][_z], ver[t][_y]);
        for (t = vertices; t <= vertices; t++)
            fprintf(volume_file, "\t\t\t\t%f %f %f\n", ver[t][_x], ver[t][_z], ver[t][_y]);
        fprintf(volume_file, "\t\t\t]\n");
        fprintf(volume_file, "\t\t}\n");

        fprintf(volume_file, "\t\tMaterialBinding {\n\t\t\tvalue OVERALL\n\t\t}\n");

    /* Write polygons */
        fprintf(volume_file, "\t\tIndexedFaceSet {\n");
        fprintf(volume_file, "\t\t\tcoordIndex [\n");
        for (t = 1; t < polygons; t++) {
            fprintf(volume_file, "\t\t\t\t");
            if (poly_store[t][2] == poly_store[t][3])
                max = 2;
            else
                max = 3;
            for (i = 0; i <= max; i++)
                fprintf(volume_file, "%d, ", poly_store[t][i] - 1);
            fprintf(volume_file, "-1,\n");
        }
        for (t = polygons; t <= polygons; t++) {
            fprintf(volume_file, "\t\t\t\t");
            if (poly_store[t][2] == poly_store[t][3])
                max = 2;
            else
                max = 3;
            for (i = 0; i <= max; i++)
                fprintf(volume_file, "%d, ", poly_store[t][i] - 1);
            fprintf(volume_file, "-1\n");
        }
        fprintf(volume_file, "\t\t\t]\n");
        fprintf(volume_file, "\t\t}\n");

        fprintf(volume_file, "\t}\n");

    } else if (dxf3) {                     /* simple raw format */
        for (i = 1; i <= polygons; i++) {
			//SCR cast
            if (Invalid_polygon((short)i))
                continue;

            if (poly_store[i][2] == poly_store[i][3]) { /* 3 vertex triangle */
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "%g\n", ver[poly_store[i][2]][_z]);

            } else {                       /* 4 vertex polygon = 2x triangle */
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][1]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "%g\n", ver[poly_store[i][2]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][0]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_y]);
                fprintf(volume_file, "%g ", ver[poly_store[i][2]][_z]);
                fprintf(volume_file, "%g ", ver[poly_store[i][3]][_x]);
                fprintf(volume_file, "%g ", ver[poly_store[i][3]][_y]);
                fprintf(volume_file, "%g\n", ver[poly_store[i][3]][_z]);
            }
        }

    } else {                               /* standard VOL output */
        Fput_bin_s16(volume_file, 20);
        Fput_bin_s16(volume_file, vertices);
        Fput_bin_s16(volume_file, polygons);
        Fput_bin_s16(volume_file, color);

        for (t = 1; t <= vertices; t++) {
            Fput_bin_r32(volume_file, ver[t][_x]);
            Fput_bin_r32(volume_file, ver[t][_y]);
            Fput_bin_r32(volume_file, ver[t][_z]);
        }

        for (t = 1; t <= polygons; t++) {
            for (i = 0; i <= 3; i++)
                Fput_bin_s16(volume_file, poly_store[t][i]);
        }
    }
}


/* Add object ------------------------------------------------------------- */


static void
Inverse(vector t, vector v)
{                                           /* Inverse vector transform of a
                                             * matrix built in C123 */
    v[_x] = C1[_x] * t[_x] + C2[_x] * t[_y] + C3[_x] * t[_z];
    v[_y] = C1[_y] * t[_x] + C2[_y] * t[_y] + C3[_y] * t[_z];
    v[_z] = C1[_z] * t[_x] + C2[_z] * t[_y] + C3[_z] * t[_z];
}


static void
Set_ECS(vector n)
{                                           /* Build an ECS transform in the
                                             * axis_xyz vars, used for dxf1
                                             * output */
    vector                  Wy, Wz;
    r32                     fac = (r32) 0.015625;

    Wy[_x] = (r32) 0.0;
    Wy[_y] = (r32) 1.0;
    Wy[_z] = (r32) 0.0;

    Wz[_x] = (r32) 0.0;
    Wz[_y] = (r32) 0.0;
    Wz[_z] = (r32) 1.0;

    Vector_copy_r32(n, axis_z);
    Vector_normalize(axis_z);

    if ((Abs_r32(n[_x]) < fac) && (Abs_r32(n[_y]) < fac)) {
        Vector_product(Wy, axis_z, axis_x);
    } else {
        Vector_product(Wz, axis_z, axis_x);
    }

    Vector_normalize(axis_x);
    Vector_product(axis_z, axis_x, axis_y);
    Vector_normalize(axis_y);
}


static void
Inverse_ECS(vector p1, vector p2)
{                                           /* Used for dxf1 output */
    p2[_x] = axis_x[_x] * p1[_x] + axis_x[_y] * p1[_y] + axis_x[_z] * p1[_z];
    p2[_y] = axis_y[_x] * p1[_x] + axis_y[_y] * p1[_y] + axis_y[_z] * p1[_z];
    p2[_z] = axis_z[_x] * p1[_x] + axis_z[_y] * p1[_y] + axis_z[_z] * p1[_z];
}


static void
Read_form(void)
{                                           /* Read in and store an external
                                             * object */
    s32                     i, j;
    s16                     ver, pol, k;
    r32                     LR, LE, LT;
    r32                     Dum32;
    s16                     Dum16;

    if (pov_form || pov_form2 || pov_form3 || blb_form)
        return;

    strcat(x_name, ".vol");
    volume_file = fopen(x_name, "rb");
    if (!volume_file)
        User_error("Cannot open file [%s]", x_name);

    Buffer_IO(volume_file, 0);

    Dum16 = getc(volume_file);
	//SCR cast
    Set_lowhigh((short)(Dum16 == 11));            /* get storage mode */

    Fget_bin_r32(volume_file, &Dum32);     /* header */
    Fget_bin_r32(volume_file, &Dum32);
    Fget_bin_r32(volume_file, &LR);
    Fget_bin_r32(volume_file, &LE);
    Fget_bin_r32(volume_file, &LT);
    fread(&Dum16, sizeof(s16), 1, volume_file);
    fread(&Dum16, sizeof(s16), 1, volume_file);
    fread(&Dum16, sizeof(s16), 1, volume_file);
    fread(&Dum16, sizeof(s16), 1, volume_file);
    fread(&Dum16, sizeof(u16), 1, volume_file);

    form_ver = 0;
    form_pol = 0;
    k = 0;

    for (;;) {
        if (feof(volume_file))
            break;

        Fget_bin_s16(volume_file, &Dum16);
        Fget_bin_s16(volume_file, &ver);   /* vertex count */
        Fget_bin_s16(volume_file, &pol);   /* polygon count */
        Fget_bin_s16(volume_file, &Dum16);

        for (i = 1; i <= ver; i++) {       /* vertices */
            form_ver++;
            Fget_bin_r32(volume_file, &form_c[form_ver][_x]);
            Fget_bin_r32(volume_file, &form_c[form_ver][_y]);
            Fget_bin_r32(volume_file, &form_c[form_ver][_z]);
        }

        for (i = 1; i <= pol; i++) {       /* polygons */
            form_pol++;
            for (j = 0; j <= 3; j++)
                Fget_bin_s16(volume_file, &form_s[form_pol][j]);
            for (j = 0; j <= 3; j++)
                form_s[form_pol][j] += k;
        }

        k = form_ver;
    }

    fclose(volume_file);
}


static void
Define_form(vector p1, vector p2, vector up, s16 c)
{                                           /* Insert external object. The
                                             * basis is to create the resulting
                                             * normalized direction vectors
                                             * from the turtle vectors and use
                                             * this a matrix for inverse
                                             * transforming the object from its
                                             * location round the origin to its
                                             * location at the turtle origin */
    vector                  dis, d1, d2, d3, in, ext, Q, P;
    s16                     i;
    r32                     s, d, r1, r2;
    char                    layer[10] = "";

 /* p1 = location, p2 = forward, up = up */
    zmin = MIN(zmin, p1[_z]);
    zmin = MIN(zmin, p2[_z]);

 /* setup */
    Vector_min(p2, p1, dis);
    d = Vector_length(dis);
    if (d == (r32) 0.0)
        return;
    s = d * thick;
    s = (s < (r32) min_thick) ? min_thick : s;

 /* d1 */
    Vector_copy_r32(dis, d1);
    Vector_normalize(d1);

    if (dxf1) {                            /* setup ECS and insert the block
                                            * reference */
        Set_ECS(d1);
        Inverse_ECS(p1, d2);
        Vector_copy(d1, ext);
        Vector_copy(d2, in);
        sprintf(layer, "%d", c);
        fprintf(volume_file, "0\nINSERT\n8\n%s\n2\nBLOCK\n", layer);
        fprintf(volume_file, "10\n%g\n20\n%g\n30\n%g\n", in[_x], in[_y], in[_z]);
        fprintf(volume_file, "41\n%g\n42\n%g\n43\n%g\n", s, s, d);
        fprintf(volume_file, "210\n%g\n220\n%g\n230\n%g\n", ext[_x], ext[_y], ext[_z]);
        polcount++;
        return;
    };

 /* d2 */
    Vector_copy_r32(up, d2);
    Vector_normalize(d2);

 /* d3 */
    Vector_product(d1, d2, d3);
    Vector_normalize(d3);

 /* setup transform */
    Vector_copy_r32(d3, C1);               /* new x-axis */
    Vector_copy_r32(d2, C2);               /* new y-axis */
    Vector_copy_r32(d1, C3);               /* new z-axis */

    if (pov_form) {                        /* insert a reference to l_base */
        //SCR cast
		d =(float) (d* 0.7);
        s =(float) (s* 0.7);
        r1 =(float) (57.0 * Do_angle(0.0, 0.0, d1[_z], (float)sqrt(d1[_x] * d1[_x] + d1[_y] * d1[_y])));
        r2 =(float) (57.0 * Do_angle(0.0, 0.0, d1[_x], d1[_y]));
        fprintf(volume_file, "object{l_base ");
        fprintf(volume_file, "finish{t_base} pigment{color col_%d}", c % 16);
        fprintf(volume_file, "scale<%g, %g, %g>", s, d, s);
        fprintf(volume_file, "rotate<%g, %g, %g>", 0.0, 0.0, -r1);
        fprintf(volume_file, "rotate<%g, %g, %g>", 0.0, -r2, 0.0);
        fprintf(volume_file, "translate<%g, %g, %g>}\n", p1[_x], p1[_z], p1[_y]);
        polcount++;
        return;

    } else if (pov_form2) {                /* write out a blob origin */
        fprintf(volume_file, "component 1.0 %g <%g, %g, %g>\n", d, p1[_x], p1[_z], p1[_y]);
        polcount++;
        return;

    } else if (pov_form3) {                /* write out a blob origin to the
                                            * file based on the color */
        fprintf(vf[c % 8], "component 1.0 %g <%g, %g, %g>\n", d, p1[_x], p1[_z], p1[_y]);
        polcount++;
        return;

    } else if (blb_form) {                 /* write out a blob in BLB format */
        fprintf(volume_file, "Sphere = %g %g %g 1.0 %g\n", p1[_x], p1[_z], p1[_y], d);
        polcount++;
        return;
    };

 /* Else use the good old Lviewer VOL format */
    for (i = 1; i <= form_ver; i++) {      /* vertices */
        Q[_x] = form_c[i][_x] * s;
        Q[_y] = form_c[i][_y] * s;
        Q[_z] = form_c[i][_z] * d;
        Inverse(Q, P);
        Vector_plus(P, p1, ver[i]);
    }

    for (i = 1; i <= form_pol; i++) {      /* polygons */
        poly_store[i][0] = form_s[i][0];
        poly_store[i][1] = form_s[i][1];
        poly_store[i][2] = form_s[i][2];
        poly_store[i][3] = form_s[i][3];
    }

    Save_object(form_ver, form_pol, c);    /* save the stored object */
}


static void
Define_block(vector p1, vector p2, vector up, s16 c)
{                                           /* Insert basic block. Here we
                                             * build a cube shape directly on
                                             * the input vectors. */
    vector                  dis, d1, d2, d3;
    s16                     i;
    r32                     s, d;

    if (pov_form || pov_form2 || pov_form3 || blb_form)
        return;

    zmin = MIN(zmin, p1[_z]);
    zmin = MIN(zmin, p2[_z]);

 /* setup */
    Vector_min(p2, p1, dis);
    d = Vector_length(dis);
    if (d == (r32) 0.0)
        return;
    s = d * thick;
    s = (s < (r32) min_thick) ? min_thick : s;
    s *= 0.5;

 /* d1 */
    Vector_copy_r32(dis, d1);
    Vector_normalize(d1);

 /* d2 */
    Vector_copy_r32(up, d2);
    Vector_normalize(d2);

 /* d3 */
    Vector_product(d1, d2, d3);
    Vector_normalize(d3);

 /* base 1, 3 */
    Vector_plus(d2, d3, d1);
    Vector_normalize(d1);
    Vector_plus_fac(p1, d1, s, ver[1]);
    Vector_plus_fac(p1, d1, -s, ver[3]);

 /* base 2, 4 */
    Vector_min(d2, d3, d1);
    Vector_normalize(d1);
    Vector_plus_fac(p1, d1, s, ver[2]);
    Vector_plus_fac(p1, d1, -s, ver[4]);

 /* end */
    for (i = 1; i <= 4; i++)
        Vector_plus(ver[i], dis, ver[i + 4]);

 /* polygons */
    poly_store[1][0] = 1;
    poly_store[1][1] = 5;
    poly_store[1][2] = 6;
    poly_store[1][3] = 2;

    poly_store[2][0] = 2;
    poly_store[2][1] = 6;
    poly_store[2][2] = 7;
    poly_store[2][3] = 3;

    poly_store[3][0] = 3;
    poly_store[3][1] = 7;
    poly_store[3][2] = 8;
    poly_store[3][3] = 4;

    poly_store[4][0] = 4;
    poly_store[4][1] = 8;
    poly_store[4][2] = 5;
    poly_store[4][3] = 1;

    poly_store[5][0] = 1;
    poly_store[5][1] = 2;
    poly_store[5][2] = 3;
    poly_store[5][3] = 4;

    poly_store[6][0] = 8;
    poly_store[6][1] = 7;
    poly_store[6][2] = 6;
    poly_store[6][3] = 5;

    Save_object(8, 6, c);
}


static void
Define_closed(vector p1, vector p2, vector up, s16 c)
{                                           /* Insert connected cylinder shape.
                                             * The lastxxx vars are used to
                                             * store the previous top of the
                                             * cylinder for conntecing a next
                                             * one. Since the vars are stacked
                                             * for [] we can connect correctly
                                             * according to current nesting
                                             * level. */
    vector                  dis, d1, d2, d3, t1, t2;
    s16                     i, ii;
    r32                     s, d, dd = float_max;

    if (pov_form || pov_form2 || pov_form3 || blb_form)
        return;

    zmin = MIN(zmin, p1[_z]);
    zmin = MIN(zmin, p2[_z]);

 /* setup */
    Vector_min(p2, p1, dis);
    d = Vector_length(dis);
    if (d == (r32) 0.0)
        return;
    s = d * thick;
    s = (s < (r32) min_thick) ? min_thick : s;
    s *= 0.5;

 /* d1 */
    Vector_copy_r32(dis, d1);
    Vector_normalize(d1);

 /* d2 */
    Vector_copy_r32(up, d2);
    Vector_normalize(d2);

 /* d3 */
    Vector_product(d1, d2, d3);
    Vector_normalize(d3);

    Vector_plus(d2, d3, t1);
    Vector_normalize(t1);
    Vector_min(d2, d3, t2);
    Vector_normalize(t2);

    Vector_plus_fac(p1, t1, s, ver[1]);
    Vector_plus_fac(p1, t1, -s, ver[5]);
    Vector_plus_fac(p1, t2, s, ver[3]);
    Vector_plus_fac(p1, t2, -s, ver[7]);

    s *= (r32) 0.7071;
    Vector_plus_fac2(p1, t1, s, t2, s, ver[2]);
    Vector_plus_fac2(p1, t1, -s, t2, s, ver[4]);
    Vector_plus_fac2(p1, t1, -s, t2, -s, ver[6]);
    Vector_plus_fac2(p1, t1, s, t2, -s, ver[8]);

 /* end */
    for (i = 1; i <= 8; i++)
        Vector_plus(ver[i], dis, ver[i + 8]);

    if (last_col == c) {
        Vector_min(p1, last, dis);
        d = Vector_length(dis);

        if (d < (r32) 1.0) {
            for (i = 1; i <= 8; i++) {
                Vector_min(ver[1], last_v[i], dis);
                d = Vector_length(dis);
                if (d < dd) {
                    dd = d;
                    ii = i;
                }
            }
            for (i = 1; i <= 8; i++) {
                Vector_copy_r32(last_v[ii], ver[i]);
                ii = (ii + 1) % 9;
                if (ii == 0)
                    ii = 1;
            }
        }
    };

 /* polygons */
    poly_store[1][0] = 1;
    poly_store[1][1] = 9;
    poly_store[1][2] = 10;
    poly_store[1][3] = 2;

    poly_store[2][0] = 2;
    poly_store[2][1] = 10;
    poly_store[2][2] = 11;
    poly_store[2][3] = 3;

    poly_store[3][0] = 3;
    poly_store[3][1] = 11;
    poly_store[3][2] = 12;
    poly_store[3][3] = 4;

    poly_store[4][0] = 4;
    poly_store[4][1] = 12;
    poly_store[4][2] = 13;
    poly_store[4][3] = 5;

    poly_store[5][0] = 5;
    poly_store[5][1] = 13;
    poly_store[5][2] = 14;
    poly_store[5][3] = 6;

    poly_store[6][0] = 6;
    poly_store[6][1] = 14;
    poly_store[6][2] = 15;
    poly_store[6][3] = 7;

    poly_store[7][0] = 7;
    poly_store[7][1] = 15;
    poly_store[7][2] = 16;
    poly_store[7][3] = 8;

    poly_store[8][0] = 8;
    poly_store[8][1] = 16;
    poly_store[8][2] = 9;
    poly_store[8][3] = 1;

    Save_object(16, 8, c);

    last_col = c;
    Vector_copy_r32(p2, last);
    for (i = 1; i <= 8; i++)
        Vector_copy_r32(ver[i + 8], last_v[i]);
}


static void
Ground_plane(void)
{                                           /* Add a simple large groundplane */
    r32                     l = (r32) 1e5;

    ver[1][_x] = -l;
    ver[1][_y] = l;
    ver[1][_z] = zmin;

    ver[2][_x] = l;
    ver[2][_y] = l;
    ver[2][_z] = zmin;

    ver[3][_x] = l;
    ver[3][_y] = -l;
    ver[3][_z] = zmin;

    ver[4][_x] = -l;
    ver[4][_y] = -l;
    ver[4][_z] = zmin;

    poly_store[1][0] = 4;
    poly_store[1][1] = 3;
    poly_store[1][2] = 2;
    poly_store[1][3] = 1;

    Save_object(4, 1, 1);
}


/* L-system routines ------------------------------------------------------ */


static r32
Rnd(void)
{                                           /* Get a random number */
    return (r32) rand() * rm;
}



#define Util_t(In,C1,C2,C3,Out) {\
  Out[_x] = Scalar_product(C1,In);\
  Out[_y] = Scalar_product(C2,In);\
  Out[_z] = Scalar_product(C3,In);\
}


static void
Set_rot(r32 a, vector n)
{                                           /* Set up a rotation matrix */
    r32                     n11, n22, n33, nxy, nxz, nyz, sina, cosa;
	//SCR cast
    cosa = (float)cos(a);
    sina = (float)sin(a);

    n11 = n[_x] * n[_x];
    n22 = n[_y] * n[_y];
    n33 = n[_z] * n[_z];

    nxy = n[_x] * n[_y];
    nxz = n[_x] * n[_z];
    nyz = n[_y] * n[_z];

    C1[_x] = n11 + (one - n11) * cosa;
    C1[_y] = nxy * (one - cosa) - n[_z] * sina;
    C1[_z] = nxz * (one - cosa) + n[_y] * sina;

    C2[_x] = nxy * (one - cosa) + n[_z] * sina;
    C2[_y] = n22 + (one - n22) * cosa;
    C2[_z] = nyz * (one - cosa) - n[_x] * sina;

    C3[_x] = nxz * (one - cosa) - n[_y] * sina;
    C3[_y] = nyz * (one - cosa) + n[_x] * sina;
    C3[_z] = n33 + (one - n33) * cosa;
}


static r32
Get_value(u32 *j)
{                                           /* Read a (xx) value from a
                                             * production string at location j
                                             * and make it into a real */
    s16                     i = 0;
    r32                     r = 0.0;
    char                    val[40] = "";

    (*j)++;
    (*j)++;

    for (;;) {
        if (object_s[*j] == ')')
            break;
        val[i] = object_s[*j];
        i++;
        (*j)++;
    }

    val[i] = '\0';
    sscanf(val, "%f", &r);

    if (last_recur)
        r *= fraction;

    return r;
}


static
L_init(void)
{                                           /* Process a ls file and setup */
    FILE                   *f;
    char                    name[max_file], temp[rule_s];
    boolean                 found;
    s16                     i;

 /* Init mem */
    object_s = (char *) malloc(max_string);
    otemp = (char *) malloc(max_string);
    stack = (s_rec *) malloc(sizeof(s_rec) * max_stack);
    pstack = (p_rec *) malloc(sizeof(p_rec) * max_stack);
    if ((object_s == NULL) || (otemp == NULL) || (stack == NULL) || (pstack == NULL))
        User_error("Not enough memory to startup");

 /* Get file name */
    Get_comline_filename(name);
    strcat(name, ".ls");
    f = fopen(name, "rt");
    if (!f)
        User_error("Cannot find file [%s]", name);
    Message("L-Neuron file  : %s\n", name);

 /* Recursion level */
    do {
        Fget_string(f, temp);
    } while (temp[0] == '#');
    sscanf(temp, "%f", &recursion);
    Get_comline_opt("r", &found, temp);    /* Overrule ? */
    if (found)
        sscanf(temp, "%f", &recursion);
    Message("Recursion depth: %g\n", recursion);
    lev = (s16) recursion;
    fraction = recursion - (r32) lev;
    if (fraction > zero) {                 /* Check for fraction */
        lev++;
        growing = TRUE;
    } else {
        growing = FALSE;
    }

 /* Basic angle */
    do {
        Fget_string(f, temp);
    } while (temp[0] == '#');
    sscanf(temp, "%f", &ang);
    Get_comline_opt("a", &found, temp);    /* Overrule ? */
    if (found)
        sscanf(temp, "%f", &ang);
    Message("Basic angle    : %g\n", ang);
	//SCR cast
    ang = (float)((ang / 180.0) * 3.141592654);

 /* Thickness */
    do {
        Fget_string(f, temp);
    } while (temp[0] == '#');
    sscanf(temp, "%f", &thick);
    Message("Thickness      : %g\n", thick);
    thick /= 100.0;

 /* Axiom */
    do {
        Fget_string(f, temp);
    } while (temp[0] == '#');
    strcpy(object_s, strtok(temp, " \r\n\t#"));
    Message("Axiom          : %s\n", object_s);

 /* Get rules */
    num = 0;
    for (i = 0; i < 150; i++) {
        do {
            Fget_string(f, temp);
        } while (temp[0] == '#');
        strcpy(rule[num], strtok(temp, " \r\n\t#"));
        if (rule[num][0] == '\0')
            continue;
        if (rule[num][0] == '@')
            break;
        num++;
    };

    fclose(f);

 /* Add default rules */
    strcpy(rule[num++], "+=+");
    strcpy(rule[num++], "-=-");
    strcpy(rule[num++], "&=&");
    strcpy(rule[num++], "^=^");
    strcpy(rule[num++], "<=<");
    strcpy(rule[num++], ">=>");

    strcpy(rule[num++], "%=%");
    strcpy(rule[num++], "|=|");
    strcpy(rule[num++], "!=!");
    strcpy(rule[num++], "?=?");
    strcpy(rule[num++], ":=:");
    strcpy(rule[num++], ";=;");
    strcpy(rule[num++], "\'=\'");
    strcpy(rule[num++], "\"=\"");
    strcpy(rule[num++], "c=c");

    strcpy(rule[num++], "[=[");
    strcpy(rule[num++], "]=]");
    strcpy(rule[num++], "{={");
    strcpy(rule[num++], "}=}");

    strcpy(rule[num++], "F=F");
    strcpy(rule[num++], "f=f");
    strcpy(rule[num++], "t=t");
    strcpy(rule[num++], "x=x");
    strcpy(rule[num++], "y=y");
    strcpy(rule[num++], "w=w");
    strcpy(rule[num++], "s=s");
    strcpy(rule[num++], "g=g");
    strcpy(rule[num++], "Z=Z");
    strcpy(rule[num++], "z=z");
    strcpy(rule[num++], "*=*");
    strcpy(rule[num++], "$=$");
    strcpy(rule[num++], "~=~");

    strcpy(rule[num++], ".=.");
    strcpy(rule[num++], "1=1");
    strcpy(rule[num++], "2=2");
    strcpy(rule[num++], "3=3");
    strcpy(rule[num++], "4=4");
    strcpy(rule[num++], "5=5");
    strcpy(rule[num++], "6=6");
    strcpy(rule[num++], "7=7");
    strcpy(rule[num++], "8=8");
    strcpy(rule[num++], "9=9");
    strcpy(rule[num++], "0=0");
    strcpy(rule[num++], "(=(");
    strcpy(rule[num++], ")=)");

    strcpy(rule[num++], "_=_");            /* closer default */

 /* Set start values fir F and Z distances */
    dis = 100.0;
	//SCR cast
    dis2 = (float) (dis * 0.5);

 /* Get sizes and marks */
    for (i = 0; i < num; i++) {
        size[i] = strlen(rule[i]) - 2;
        mark[i] = FALSE;
    };

 /* Check which rules need to be marked for last recursion when growing */
    for (i = 0; i < num; i++) {
        if (rule[i][0] == '+')
            break;
        mark[i] = TRUE;

    /* All rules with basic move/block before '=' mark false */
        if (rule[i][0] == 'F')
            mark[i] = FALSE;
        if (rule[i][0] == 'f')
            mark[i] = FALSE;
        if (rule[i][0] == 'Z')
            mark[i] = FALSE;
        if (rule[i][0] == 'z')
            mark[i] = FALSE;

        Message("Rule           : %s\n", rule[i]);
    }
}


static void
L_mutate(void)
{                                           /* Apply mutations to the rules */
    s16                     i, j, k, rules, ii, max = 1000;
    char                    T, R, S[10] = "";
    char                    rulet[100] = "";

    for (i = 0; i < num; i++) {
        if (rule[i][0] == '+')
            break;
    }
    rules = i;

    switch ((s16) (Rnd() * 6.0)) {

      default:
        return;

      case 1:                              /* Insert */
        i = (s16) (Rnd() * (r32) rules);
        T = rule[i][0];
        j = (s16) (Rnd() * (r32) rules);
        k = (s16) (Rnd() * (r32) strlen(rule[j]));
        k = (k < 2) ? 2 : k;
        strcpy(rulet, &(rule[j][k]));
        rule[j][k] = '[';
        rule[j][k + 1] = T;
        rule[j][k + 2] = ']';
        rule[j][k + 3] = '\0';
        strcat(rule[j], rulet);
        size[j] = strlen(rule[j]) - 2;
        break;

      case 0:
      case 2:                              /* Replace */
        do {
            i = (s16) (Rnd() * (r32) rules);
            j = (s16) (Rnd() * (r32) rules);
            T = rule[i][0];
            R = rule[j][0];
        } while (T == R);
        for (ii = 0; ii < max; ii++) {
            i = (s16) (Rnd() * (r32) rules);
            for (j = 2; j < (int)strlen(rule[i]); j++) {
                if (rule[i][j] == T) {
                    rule[i][j] = R;
                    return;
                }
            }
        }
        break;

      case 3:                              /* Append */
        i = (s16) (Rnd() * (r32) rules);
        S[0] = rule[i][0];
        i = (s16) (Rnd() * (r32) rules);
        strcat(rule[i], S);
        size[i] = strlen(rule[i]) - 2;
        break;

      case 4:                              /* Swap directions */
        for (ii = 0; ii < max; ii++) {
            i = (s16) (Rnd() * (r32) rules);
            for (j = 2; j < (int)strlen(rule[i]); j++) {
                switch ((s16) (Rnd() * 12.0)) {
                  default:
                    return;
                  case 0:
                    if (rule[i][j] == '+') {
                        rule[i][j] = '-';
                        return;
                    }
                    break;
                  case 1:
                    if (rule[i][j] == '-') {
                        rule[i][j] = '+';
                        return;
                    }
                    break;
                  case 2:
                    if (rule[i][j] == '&') {
                        rule[i][j] = '^';
                        return;
                    }
                    break;
                  case 3:
                    if (rule[i][j] == '^') {
                        rule[i][j] = '&';
                        return;
                    }
                    break;
                  case 4:
                    if (rule[i][j] == '>') {
                        rule[i][j] = '<';
                        return;
                    }
                    break;
                  case 5:
                    if (rule[i][j] == '<') {
                        rule[i][j] = '>';
                        return;
                    }
                    break;
                  case 6:
                    if (rule[i][j] == '|') {
                        rule[i][j] = '%';
                        return;
                    }
                    break;
                  case 7:
                    if (rule[i][j] == '%') {
                        rule[i][j] = '|';
                        return;
                    }
                    break;
                  case 8:
                    if (rule[i][j] == ':') {
                        rule[i][j] = ';';
                        return;
                    }
                    break;
                  case 9:
                    if (rule[i][j] == ';') {
                        rule[i][j] = ':';
                        return;
                    }
                    break;
                  case 10:
                    if (rule[i][j] == '\'') {
                        rule[i][j] = '\"';
                        return;
                    }
                    break;
                  case 11:
                    if (rule[i][j] == '\"') {
                        rule[i][j] = '\'';
                        return;
                    }
                    break;
                }
            }
        }
        break;

      case 5:                              /* Swap sizes */
        for (ii = 0; ii < max; ii++) {
            i = (s16) (Rnd() * (r32) rules);
            for (j = 2; j < (int)strlen(rule[i]); j++) {
                switch ((s16) (Rnd() * 6.0)) {
                  default:
                    return;
                  case 0:
                    if (rule[i][j] == 'F') {
                        rule[i][j] = 'Z';
                        return;
                    }
                    break;
                  case 1:
                    if (rule[i][j] == 'Z') {
                        rule[i][j] = 'F';
                        return;
                    }
                    break;
                  case 2:
                    if (rule[i][j] == 'f') {
                        rule[i][j] = 'z';
                        return;
                    }
                    break;
                  case 3:
                    if (rule[i][j] == 'z') {
                        rule[i][j] = 'f';
                        return;
                    }
                    break;
                  case 4:
                    if (rule[i][j] == '!') {
                        rule[i][j] = '?';
                        return;
                    }
                    break;
                  case 5:
                    if (rule[i][j] == '?') {
                        rule[i][j] = '!';
                        return;
                    }
                    break;
                }
            }
        }
        break;

    }
}


static
L_system(void)
{                                           /* Expand l-system into production
                                             * string. Object_s is read with
                                             * the k counter and the next
                                             * generation is build up in otemp
                                             * with the ot counter. */

    u32                     k, st, s, ss, max = max_string - 10L;
    char                   *ot;
    char                    m_1 = '@';
    s16                     S[256], i, l;
    boolean                 incomplete = FALSE, marker = FALSE;

 /*
  * Clear fast access array. This array is to quickly find the rule asociated
  * with a char.
  */
    for (i = 0; i < 256; i++)
        S[i] = num - 1;                    /* Num -1 contains the default rule
                                            * which does nothing */

 /* Each char gets a rule number */
    for (i = num - 1; i >= 0; i--)
        S[(s16) rule[i][0]] = i;

    Process_start2(lev - 1);

    for (l = 0; l < lev; l++) {            /* For each recursion */
        Process_update2(l);

        marker = ((l == (lev - 1)) && growing); /* Need markers ? */

        st = strlen(object_s);
        ot = otemp;
        ss = 0;

        for (k = (u32) 0; k < st; k++) {   /* For each char in the string */
            i = S[object_s[k]];            /* i = rule number attached to
                                            * current char */
            s = size[i];                   /* s = size of current rule */
            ss += s;

            if (ss >= max) {               /* Overflow */
                l = lev;
                incomplete = TRUE;
                break;

            } else {
                if (marker && mark[i]) {   /* Add mark char */
                    Mem_copy(&m_1, ot, 1);
                    ot += 1;
                    Mem_copy(&(rule[i][2]), ot, s);
                    ot += s;
                    Mem_copy(&m_1, ot, 1);
                    ot += 1;
                } else {
                    Mem_copy(&(rule[i][2]), ot, s);     /* Copy */
                    ot += s;
                }
            }
        };

        *ot = '\0';
        strcpy(object_s, otemp);           /* Copy the temp string to object_s
                                            * and repeat cycle */
    };

    Process_end2();
    Message("Size of string : %ld chars %s\n", strlen(object_s), incomplete ? "(incomplete)" : "");
}


static void
L_save(void)
{                                           /* Save mutated ls-system for
                                             * re-run */
    FILE                   *f;
    s16                     i;

    remove("mutation.ls");
    f = fopen("mutation.ls", "wt");
    if (!f)
        User_error("Cannot open file [mutation.ls]");

    Message("Saving ls file : mutation.ls\n");

    fprintf(f, "%d\n", lev);
    fprintf(f, "%g\n", (ang / 3.141592654) * 180.0);
    fprintf(f, "%g\n", thick * 100.0);

    fprintf(f, "%s\n", object_s);

    for (i = 0; i < num; i++) {
        if (rule[i][0] == '+')
            break;
        fprintf(f, "%s\n", rule[i]);
    }
    fprintf(f, "@\n");

    fclose(f);
}



static
L_draw(void)
{   
	r32					dia;
                                        /* Process a production string and
                                             * generate form */
    vector                  pos, end, v, fow, upp, lef;
    u32                     i, max = strlen(object_s);
    r32                     r, a, thick_l, ang_l, dis_l, dis2_l, trope_l;
    s16                     vcount, pcount, j;
    char                    temp[max_file], next;
    boolean                 found, poly_on = FALSE;
    //u32                     one_star, star_num;
    s16						cmpt, parent, type;  	/* compartment numbers JLK 2/10/99 */
	cmpt = 2;
    
 /* Echo production string */
    Get_comline_opt("l", &found, temp);
    if (found) {
        Message("Production     : ");
        for (i = 0; object_s[i] != '\0'; i++) {
            Message("%c", object_s[i]);
        }
        Message("\n");
    };

 /* Get user form if needed */
    if (user_form)
        Read_form();

 /* Setup vectors */
    pos[_x] = 0.0;
    pos[_y] = 0.0;
    pos[_z] = 0.0;
    fow[_x] = 0.0;
    fow[_y] = 0.0;
    fow[_z] = 1.0;
    lef[_x] = 0.0;
    lef[_y] = -1.0; /* JLK puts the turtle in standard right hand coordinates */
    lef[_z] = 0.0;
    upp[_x] = 1.0;
    upp[_y] = 0.0;
    upp[_z] = 0.0;
    Vector_normalize(trope);

 /* Do it */
    Open_datafile();

 /* Start values */
    org.col = col;
    org.dis = dis;
    org.dis2 = dis2;
    org.ang = ang;
    org.thick = thick;
    org.tr = tr;
		
 /* Feedback */
    Process_start2(strlen(object_s));
      //fprintf (volume_file, "%s\n",object_s);

    for (i = 0; object_s[i] != '\0'; i++) {
        Process_update2(i);

        if (polcount > poly_limit)         /* overflow */
            break;

        next = object_s[i + 1];            /* the next char in de string */

        switch (object_s[i]) {             /* the current char in the string */
          default:
            break;

          case '@':                        /* Marks last recursion level during
                                            * growing fase */
            last_recur = !last_recur;
            if (last_recur) {              /* Store all vars and do fraction */
                thick_l = thick;
                ang_l = ang;
                dis_l = dis;
                dis2_l = dis2;
                trope_l = trope_amount;
                dis *= fraction;
                dis2 *= fraction;
                thick *= fraction;
                ang *= fraction;
                trope_amount *= fraction;
            } else {                       /* Restore */
                thick = thick_l;
                ang = ang_l;
                dis = dis_l;
                dis2 = dis2_l;
                trope_amount = trope_l;
            }
            break;

          case '+':
            save.ang = ang;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            Set_rot(-ang, upp);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Vector_normalize(fow);
            Vector_normalize(lef);
            ang = save.ang;
            break;

          case '-':
            save.ang = ang;
            /* check for error condition. 
               -n is a negative number and not a turn */
            if (isdigit (next))
            	break;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            
            Set_rot(ang, upp);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Vector_normalize(fow);
            Vector_normalize(lef);
            ang = save.ang;
            break;

          case '~':
            if (next == '(')
                r = ((r32) 0.017453) * Get_value(&i);
            else if (rand_set)
                r = ((r32) 0.017453) * rand_amount;
            else
                r = (r32) 6.0;
            a = (Rnd() * r * (r32) 2.0) - r;
            Set_rot(a, upp);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Vector_normalize(fow);
            Vector_normalize(lef);
            a = (Rnd() * r * (r32) 2.0) - r;
            Set_rot(a, lef);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(fow);
            Vector_normalize(upp);
            a = (Rnd() * r * (r32) 2.0) - r;
            Set_rot(a, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(lef);
            Vector_normalize(upp);
            break;

          case 't':
            if (next == '(') {
                tr = Get_value(&i);
                if (last_recur)
                    tr *= fraction;
            }
            if ((fow[_x] == (r32) 0.0) && (fow[_y] == (r32) 0.0))
                break;
            save.tr = tr;
            if (trope_set)
                tr = trope_amount;
            Vector_copy_r32(fow, trope);
            trope[_x] = -trope[_x];
            trope[_y] = -trope[_y];
            trope[_z] = (r32) 0.0;
            Vector_normalize(trope);
            r = tr * Scalar_product(fow, trope);
            Set_rot(-r, lef);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(fow);
            Vector_normalize(upp);
            tr = save.tr;
            break;

/* tropism along the "x", "y", "z" axes, and "somatofugal" tropism -- GAA 2/14/99 */
/* corresponding symbols are x, y, z, s, respectvely. */

/* NB I think that the coordinate system of LParser is xzy i.o. xyz !! */
/* We could and probably should fix it by making the initial left position aligned */
/* with -y i.o. y, I think, i.e. by setting lef[_y] = -1 at page 45 of the code... */

		/*  tropism along the "x" axis  */
		case 'x':
   			save.tr = tr; 
   			if (next == '(') 
      			tr = Get_value (&i);
   			if ((fabs (fow[_y]) < r32_epsilon) && (fabs (fow[_z]) < r32_epsilon))
      			break;
   			if (trope_set)
      			tr = trope_amount;
			/* we don't need last_recur and fraction because we don't use fraction recursion */
   			
   			Vector_copy_r32 (fow, trope);
   			fow[_x] = fow[_x]+tr; /* tropism is given in units of intrinsic growth */
		   	Vector_normalize (fow); /* fow is done; now we need a rotation for upp and lef */
   			Vector_product(trope, fow, v); /* rotation axis for upp and lef */
   			//SCR cast
			r = (float)acos (Scalar_product (trope, fow)); /* rotation angle for upp and lef */
   			Set_rot (r, v);
   			Util_t (upp, C1, C2, C3, v);
   			Vector_copy_r32 (v, upp);
   			Util_t (lef, C1, C2, C3, v);
   			Vector_copy_r32 (v, lef);
   			Vector_normalize (lef);
   			Vector_normalize (upp);
   			tr = save.tr;
   			break;

		/*  tropism along the "y" axis  */
		case 'y':
   			save.tr = tr; 
   			if (next == '(') {
      			tr = Get_value (&i);
			}
   			if ((fabs (fow[_x]) < r32_epsilon) && (fabs (fow[_z]) < r32_epsilon))
      			break;
   			if (trope_set)
      			tr = trope_amount;
   			Vector_copy_r32 (fow, trope);
   			fow[_y] = fow[_y]+tr; /* tropism is given in units of intrinsic growth */
   			Vector_normalize (fow); /* fow is done; now we need a rotation for upp and lef */
   			Vector_product(trope, fow, v); /* rotation axis for upp and lef */
			//SCR cast
			r = (float) acos (Scalar_product (trope, fow)); /* rotation angle for upp and lef */
   			Set_rot (r, v);
   			Util_t (upp, C1, C2, C3, v);
   			Vector_copy_r32 (v, upp);
   			Util_t (lef, C1, C2, C3, v);
   			Vector_copy_r32 (v, lef);
   			Vector_normalize (lef);
   			Vector_normalize (upp);
   			tr = save.tr;
   			break;

		/*  tropism along the "z" axis: we call it "w" because its taken */
		case 'w':
   			save.tr = tr; 
   			if (next == '(') 
      			tr = Get_value (&i);
   			if ((fabs (fow[_y]) < r32_epsilon) && (fabs (fow[_x]) < r32_epsilon))
      			break;
   			if (trope_set)
      			tr = trope_amount;
   			
   			Vector_copy_r32 (fow, trope);
   			fow[_z] = fow[_z]+tr; /* tropism is given in units of intrinsic growth */
   			Vector_normalize (fow); /* fow is done; now we need a rotation for upp and lef */
   			Vector_product(trope, fow, v); /* rotation axis for upp and lef */
			//SCR cast
			r = (float) acos (Scalar_product (trope, fow)); /* rotation angle for upp and lef */
   			Set_rot (r, v);
   			Util_t (upp, C1, C2, C3, v);
   			Vector_copy_r32 (v, upp);
   			Util_t (lef, C1, C2, C3, v);
   			Vector_copy_r32 (v, lef);
   			Vector_normalize (lef);
   			Vector_normalize (upp);
   			tr = save.tr;
   			break;

		/*  somatofugal tropism */
		case 's':
   			if (next == '(') {
      			tr = Get_value (&i);
			/* we don't need last_recur and fraction because we don't use fractionary recursion */
   			}
   			Vector_copy_r32 (pos, v);
   			
   			/* check for zero vector */
   			if ((fabs (v[_x]) < r32_epsilon) && (fabs (v[_y]) < r32_epsilon) && (fabs (v[_z]) < r32_epsilon))
   				break;
   				
   			Vector_normalize (v);
   			save.tr = tr; 
   			
   			/* check for equal vectors */
   			if (((fabs (fow[_x]-v[_x]) < r32_epsilon) &&
   				 (fabs (fow[_y]-v[_y]) < r32_epsilon) &&
   				 (fabs (fow[_z]-v[_z]) < r32_epsilon)) ||
   			    ((fabs (fow[_x]+v[_x]) < r32_epsilon) &&
   				 (fabs (fow[_y]+v[_y]) < r32_epsilon) &&
   				 (fabs (fow[_z]+v[_z]) < r32_epsilon)))
      			break;
      			
   			if (trope_set)
      			tr = trope_amount;
      			
   			Vector_copy_r32 (fow, trope);
   			fow[_x] = fow[_x]+tr * v[_x]; 
   			fow[_y] = fow[_y]+tr * v[_y];
   			fow[_z] = fow[_z]+tr * v[_z]; /* tropism is given in units of intrinsic growth */
   			Vector_normalize (fow); /* fow is done; now we need a rotation for upp and lef */
   			Vector_product(trope, fow, v); /* rotation axis for upp and lef */

			//SCR 8-1-01
   			/* check for zero vector */
   			if ((fabs (v[_x]) < r32_epsilon) && (fabs (v[_y]) < r32_epsilon) && (fabs (v[_z]) < r32_epsilon))
   				break;

   			//SCR cast
			r = (float)acos (Scalar_product (trope, fow)); /* rotation angle for upp and lef */
   			Set_rot (r, v);
   			Util_t (upp, C1, C2, C3, v);
   			Vector_copy_r32 (v, upp);
   			Util_t (lef, C1, C2, C3, v);
   			Vector_copy_r32 (v, lef);
   			Vector_normalize (lef);
   			Vector_normalize (upp);
   			tr = save.tr;
   			break;

/* END of tropism along the "x', "y", "z" axes, and somatofugal tropism -- GAA */
   
          case '$':
            Vector_min(fow, sky, v);
            if (Vector_length(v) == (r32) 0.0)
                break;
            Vector_product(fow, sky, lef);
            Vector_product(fow, lef, upp);
            if (upp[_z] < (r32) 0.0) {
                upp[_x] = -upp[_x];
                upp[_y] = -upp[_y];
                upp[_z] = -upp[_z];
                lef[_x] = -lef[_x];
                lef[_y] = -lef[_y];
                lef[_z] = -lef[_z];
            }
            break;

          case '&':
            save.ang = ang;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            Set_rot(ang, lef);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(fow);
            Vector_normalize(upp);
            ang = save.ang;
            break;

          case '^':
            save.ang = ang;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            Set_rot(-ang, lef);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(fow);
            Vector_normalize(upp);
            ang = save.ang;
            break;

          case '<':
            save.ang = ang;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            Set_rot(-ang, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(lef);
            Vector_normalize(upp);
            ang = save.ang;
            break;

          case '>':
            save.ang = ang;
            if (next == '(') {
                ang = ((r32) 0.017453) * Get_value(&i);
                if (last_recur)
                    ang *= fraction;
            }
            Set_rot(ang, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(lef);
            Vector_normalize(upp);
            ang = save.ang;
            break;

          case '%':
			//SCR cast
            Set_rot((float)3.141592654, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Util_t(upp, C1, C2, C3, v);
            Vector_copy_r32(v, upp);
            Vector_normalize(lef);
            Vector_normalize(upp);
            break;

          case '|':
			//SCR cast
			Set_rot((float)3.141592654, upp);
            Util_t(fow, C1, C2, C3, v);
            Vector_copy_r32(v, fow);
            Util_t(lef, C1, C2, C3, v);
            Vector_copy_r32(v, lef);
            Vector_normalize(fow);
            Vector_normalize(lef);
            break;

          case '!':
            if (next == '(') {
                if (last_recur)
                    thick *= one + fraction * (Get_value(&i) - one);
                else
                    thick *= Get_value(&i);
            } else {
                if (last_recur)
                    thick *= one + fraction * ((r32) 0.7 - one);
                else
                    thick *= (r32) 0.7;
            }
            break;

          case '?':
            if (next == '(') {
                if (last_recur)
                    thick *= one + fraction * (Get_value(&i) - one);
                else
                    thick *= Get_value(&i);
            } else {
                if (last_recur)
                    thick /= one + fraction * ((r32) 0.7 - one);
                else
                    thick /= (r32) 0.7;
            }
            break;

          case ':':
            if (next == '(') {
                if (last_recur)
                    ang *= one + fraction * (Get_value(&i) - one);
                else
                    ang *= Get_value(&i);
            } else {
                if (last_recur)
                    ang *= one + fraction * ((r32) 0.9 - one);
                else
                    ang *= (r32) 0.9;
            }
            break;

          case ';':
            if (next == '(') {
                if (last_recur)
                    ang *= one + fraction * (Get_value(&i) - one);
                else
                    ang *= Get_value(&i);
            } else {
                if (last_recur)
                    ang /= one + fraction * ((r32) 0.9 - one);
                else
                    ang /= (r32) 0.9;
            }
            break;

          case '\'':
            if (next == '(') {
                r = Get_value(&i);
                if (last_recur) {
                    dis *= one + fraction * (r - one);
                    dis2 *= one + fraction * (r - one);
                } else {
                    dis *= r;
                    dis2 *= r;
                }
            } else {
                if (last_recur) {
                    dis *= one + fraction * ((r32) 0.9 - one);
                    dis2 *= one + fraction * ((r32) 0.9 - one);
                } else {
                    dis *= (r32) 0.9;
                    dis2 *= (r32) 0.9;
                }
            }
            break;

          case '\"':
            if (next == '(') {
                r = Get_value(&i);
                if (last_recur) {
                    dis *= one + fraction * (r - one);
                    dis2 *= one + fraction * (r - one);
                } else {
                    dis *= r;
                    dis2 *= r;
                }
            } else {
                if (last_recur) {
                    dis /= one + fraction * ((r32) 0.9 - one);
                    dis2 /= one + fraction * ((r32) 0.9 - one);
                } else {
                    dis /= (r32) 0.9;
                    dis2 /= (r32) 0.9;
                }
            }
            break;

          case 'Z':
            save.dis2 = dis2;
            if (next == '(') {
                dis2 = Get_value(&i);
                if (last_recur)
                    dis2 *= fraction;
            }
            Vector_plus_fac(pos, fow, dis2, end);
            if (user_form)
                Define_form(pos, end, upp, col);
            else if (closed_form)
                Define_closed(pos, end, upp, col);
       		else if (turner) {
       			if (turner_tree) {
       				fprintf (volume_file, "%i %i %f %f %f %f %i\n", cmpt, type, end [_x], end [_y], end [_z], thick/2.0, parent);
       				parent = cmpt;
       				++cmpt;
       			}
       		}
            else
                Define_block(pos, end, upp, col);
            Vector_copy_r32(end, pos);
            dis2 = save.dis2;
            break;

          case 'F':
            save.dis = dis;
            if (next == '(') {
                dis = Get_value(&i);
                if (last_recur)
                    dis *= fraction;
            }
            Vector_plus_fac(pos, fow, dis, end);
            if (user_form)
                Define_form(pos, end, upp, col);
            else if (closed_form)
                Define_closed(pos, end, upp, col);
       		else if (turner) {
       			if (turner_tree) {
					if(thick>100){
						int ertrt=cmpt;
					}
					if(cmpt==41){
						int hhh=0;
					}
					dia= dis * thick /2.0;
					if(cmpt==8){
						int ii=0;
					}
       				fprintf (volume_file, "%i %i %f %f %f %f %i\n", cmpt, type, end [_x], end [_y], end [_z], dia, parent);
					parent = cmpt;
       				++cmpt;
       			}
       		}
       		else
                Define_block(pos, end, upp, col);
            Vector_copy_r32(end, pos);
            dis = save.dis;
            break;

          case '[':
            if (scount > max_stack)
                User_error("Ran out of stack");
            Vector_copy_r32(pos, stack[scount].pos);
            Vector_copy_r32(fow, stack[scount].fow);
            Vector_copy_r32(lef, stack[scount].lef);
            Vector_copy_r32(upp, stack[scount].upp);
            stack[scount].col = col;
            stack[scount].dis = dis;
            stack[scount].dis2 = dis2;
            stack[scount].ang = ang;
            stack[scount].thick = thick;
            stack[scount].tr = tr;
            stack[scount].parent = parent;
            if (closed_form) {
                Vector_copy_r32(last, stack[scount].last);
                stack[scount].last_col = last_col;
                for (j = 1; j <= 8; j++)
                    Vector_copy_r32(last_v[j], stack[scount].last_v[j]);
            }
            scount++;
            break;

          case ']':
            scount--;
            Vector_copy_r32(stack[scount].pos, pos);
            Vector_copy_r32(stack[scount].fow, fow);
            Vector_copy_r32(stack[scount].lef, lef);
            Vector_copy_r32(stack[scount].upp, upp);
            col = stack[scount].col;
            dis = stack[scount].dis;
            dis2 = stack[scount].dis2;
            ang = stack[scount].ang;
            thick = stack[scount].thick;
            parent = stack[scount].parent;
            tr = stack[scount].tr;
            if (closed_form) {
                Vector_copy_r32(stack[scount].last, last);
                last_col = stack[scount].last_col;
                for (j = 1; j <= 8; j++)
                    Vector_copy_r32(stack[scount].last_v[j], last_v[j]);
            }
            break;

          case '{':
            if (poly_on) {
                pstack[pscount].count = vcount;
                pstack[pscount].ver = (vector *) malloc(vcount * 12L);
                if (pstack[pscount].ver == NULL)
                    User_error("Ran out of memory");
                Vector_copy_max_r32(vcount, ver, pstack[pscount].ver);
                pscount++;
                if (pscount > max_stack)
                    User_error("Ran out of stack");
            }
            poly_on = TRUE;
            vcount = (s16) 1;
            pcount = (s16) 1;
            Vector_copy_r32(pos, ver[vcount++]);
            break;

          case 'f':
            save.dis = dis;
            if (next == '(') {
                dis = Get_value(&i);
                if (last_recur)
                    dis *= fraction;
            }
            Vector_plus_fac(pos, fow, dis, pos);
            if (poly_on)
                Vector_copy_r32(pos, ver[vcount++]);
            dis = save.dis;
            break;

          case '.':
            if (poly_on)
                Vector_copy_r32(pos, ver[vcount++]);
            break;

          case 'g':
            save.dis = dis;
            if (next == '(') {
                dis = Get_value(&i);
                if (last_recur)
                    dis *= fraction;
            }
            Vector_plus_fac(pos, fow, dis, pos);
            dis = save.dis;
            break;

          case 'z':
            save.dis2 = dis2;
            if (next == '(') {
                dis2 = Get_value(&i);
                if (last_recur)
                    dis2 *= fraction;
            }
            Vector_plus_fac(pos, fow, dis2, pos);
            if (poly_on)
                Vector_copy_r32(pos, ver[vcount++]);
            dis2 = save.dis2;
            break;

          case '}':
            if (vcount > (s16) 3) {
                for (j = 1; j < vcount - 2; j++) {
                    poly_store[pcount][0] = 1;
                    poly_store[pcount][1] = j + 1;
                    poly_store[pcount][2] = j + 2;
                    poly_store[pcount++][3] = j + 2;
                }
                Save_object((short) (vcount - 1), (short) (pcount - 1),(short) col);
            }
            poly_on = FALSE;
            if (pscount > 0) {
                pscount--;
                Vector_copy_max_r32(pstack[pscount].count, pstack[pscount].ver, ver);
                vcount = pstack[pscount].count;
                poly_on = TRUE;
            }
            break;

          case 'c':
            if (next == '(') {
                col = (s16) Get_value(&i);
            } else {
                col++;
            }

            /*	JLK 2/10/99 */
            /* 		Start writing out the Turner swc file after the soma
            		is drawn.  The tree starts with a "c"					*/
			if (turner) {
				parent = 1;
				type = col-1;
				if ((!turner_tree) && (col != SOMACOLOR)) {
					turner_tree = TRUE;
					fprintf (volume_file, "1 1 %f %f %f %f -1\n", pos [_x], pos [_y], pos [_z], soma_d/2.0);
				}
			}
            break;
          
        }
    };
     	
    Process_end2();
}


/* Lparser main ----------------------------------------------------------- */


static void
Help(void)
{
    char                    s[max_file];

    Get_comline_progname(s);
    Message("%s [options] ls-filename\n\n", s);

    Message("%s \tnone      (default) output Lviewer VOL file\n", "vol");
    Message("%s \t-v        output POV object file\n", "pov");
    Message("%s \t-b        output POV blob file\n", "pov");
    Message("%s \t-B        output multiple POV blob files\n", "pov");
    Message("%s \t-c        output inc files instead of pov files\n", "pov");
    Message("%s \t-d        output polyface meshes DXF file\n", "dxf");
    Message("%s \t-3        output 3dfaces DXF file\n", "dxf");
    Message("%s \t-R        output triangles in RAW file\n", "raw");
    Message("%s \t-O        output Blob Sculptor BLB file\n", "blb");
    Message("%s \t-V        output VRML world file\n", "wrl");

    Message("%s \t-X [name] use name.vol as base element\n", "base");
    Message("%s \t-i        link base elements together\n", "base");

    Message("%s \t-S [num]  set string size to num Kbytes\n", "set");
    Message("%s \t-t [num]  set minimum thickness\n", "set");
    Message("%s \t-r [num]  overrule recursion depth\n", "set");
    Message("%s \t-a [num]  overrule angle\n", "set");

    Message("%s \t-u [num]  mutate [num] times\n", "ls");
    Message("%s \t-l        show final L-string\n", "ls");
    Message("%s \t-g        add ground plane\n", "ls");
    Message("%s \t-L [num]  set amount for ~ command\n", "ls");
    Message("%s \t-P [num]  set amount for t command\n", "ls");

    Message("%s \t-p [num]  limit polygons to [num]\n", "limit");
    Message("%s \t-M        use Burke-Marks algorithm\n", "bur");
    Message("%s \t-H        use Hillman algorithm\n", "hil");
    Message("%s \t-Y        use Tamori algorithm\n", "tam");
    Message("%s \t-e        extend also with Hillman and Tamori algorithm\n", "ext");
    Message("%s \t-T        output in Southampton Archive format\n", "swc");
    Message("%s \t-G        output in GENESIS format\n", "p");
    Message("%s \t-N        output in NEURON format\n", "???");
    Message("%s \t-m [name] morphometric measurement of a Southampton Archive file\n", "swc");
    Message("%s \t-2        use Burke's model 2 algorithm\n", "???");
    Message("%s \t-D        use generic developer algorithm\n", "???");
    Message("%s \t-s [num]  set random number generator seed to num\n", "???");
}


void
main(int argc, char *argv[])
{
    char                    temp[max_file];
    boolean                 found, burke_found, hillman_found, tamori_found;
    s16                     tim, i;
	char                    name[max_file];
	FILE					*f;
	int ii=0;

	/*
	//gamma random generator
	if(argc==5){
			//rndgamma(argv[0],1,0);
		double alpha,lambda,offset,num;
		alpha=atof(argv[1]);
		lambda=atof(argv[2]);
		offset=atof(argv[3]);
		num=atof(argv[4]);
		printf("alpha:%g \nlambda:%g\noffset:%g\n#points:%g\n",alpha,lambda,offset,num);
		
		for(ii=0;ii<num;ii++){
			printf("%g\n",rndgamma(alpha,lambda,multip));
		}	
		
	}
	exit(0);

  */
 /* Store the pointers to the comline */
    s_argc = argc;
    s_argv = argv;

 /* Display header and help file if needed */
 	//char *ver=VERSION;
    Message("%s\n", min_bar);
    Message("%s v%s (%s)\n", "L-Neuron Generator", VERSION, __DATE__);
    Message("%s v%s\n", "\tbased on L-System Parser/Mutator", "5.0");
    Message("%s\n", min_bar);

    if (argc < 2) {
        Help();
        User_error("Need arguments\n");
    };

 /* Set the option string */
    strcpy(opts, "VP:OhL:S:iecBbRd3vX:p:t:u:r:a:glHTYMs:");

 /* Init files */
    Init_file_buf(0);
    Init_file_buf(1);
    Init_file_buf(2);

 /* Check for all the comline options */
   	Get_comline_opt("t", &found, temp);    /* set minimum thickness */
   	if (found)
       	sscanf(temp, "%f", &min_thick);

   	Get_comline_opt("p", &found, temp);    /* limit total generated polygons */
   	if (found)
       	sscanf(temp, "%ld", &poly_limit);

   	Get_comline_opt("L", &rand_set, temp); /* set random amount */
   	if (rand_set) {
       	sscanf(temp, "%f", &rand_amount);
       	srand(0);                          /* when using changing amounts of
                                            * randomness for animations we want
                                            * the seed to be the same so we set
                                            * it here */
   	} else {
       	srand(time(NULL));
   	};

   	Get_comline_opt("P", &trope_set, temp);/* set amount of trope */
   	if (trope_set) {
       	sscanf(temp, "%f", &trope_amount);
       	srand(0);                          /* same for trope animations */
   	};

   	Get_comline_opt("S", &found, temp);    /* set maximum production string
                                            * size in kbytes */
   	if (found) {
       	sscanf(temp, "%ld", &max_string);
       	max_string *= 1024L;               /* we will need two string of this
            	                                * size ! */
   	} else {
       	max_string = 2L * 1024L * 1024L;   /* default is 2 mbytes */
   	};

   	Get_comline_opt("X", &user_form, x_name);   /* use a vol file as user shape */
   	Get_comline_opt("i", &closed_form, temp);   /* create closed connected
                                                 * cylinders */
   	Get_comline_opt("d", &dxf1, temp);     /* create a dxf file with inserts */
   	Get_comline_opt("3", &dxf2, temp);     /* create a dxf with only polygons */
   	Get_comline_opt("R", &dxf3, temp);     /* create a RAW triangle file */
   	Get_comline_opt("V", &vrml, temp);     /* create a WRL VRML v1.0 file */
   	Get_comline_opt("T", &turner, temp);   /* create a Turner Southampton Archive file */
   	Get_comline_opt("B", &pov_form3, temp);/* create multiple povray blob files
                                            * for multi colord blobs */

   	Get_comline_opt("O", &blb_form, temp); /* create a blb blob file */
   	Get_comline_opt("b", &pov_form2, temp);/* create a povray blob file */
   	Get_comline_opt("v", &pov_form, temp); /* create a povray file */
   	Get_comline_opt("c", &inc_out, temp);  /* create a povray inc file */

 	/*
  	* In these cases an user defined form is used. In the other cases a simple
  	* block or connected cylinder shape will be used.
  	*/
    if (pov_form || pov_form2 || pov_form3 || blb_form)
        user_form = TRUE;

	/* check for a non l-system algorithm */
    Get_comline_opt("M", &burke_found, temp);
    Get_comline_opt("H", &hillman_found, temp);
   	Get_comline_opt("Y", &tamori_found, temp);
   	Get_comline_opt("s", &seed_set, temp);
   	if (seed_set)
	//changing %f to %d to fix the random seed error- sridevi 1/20/2012
       	sscanf(temp, "%d", &seed_amount);
    else 
        seed_amount = 0;
   	
   	if (burke_found || hillman_found || tamori_found) {

 		/* Init mem */
 		srand (0); /* make output deterministic for now */
    	object_s = (char *) malloc(max_string);
    	otemp = (char *) malloc(max_string);
    	stack = (s_rec *) malloc(sizeof(s_rec) * max_stack);
    	pstack = (p_rec *) malloc(sizeof(p_rec) * max_stack);
    	if ((object_s == NULL) || (otemp == NULL) || (stack == NULL) || (pstack == NULL))
        	User_error("Not enough memory to startup");

		 /* Get file name */
    	Get_comline_filename(name);
     	f = fopen(name, "rt");
    	if (!f)
        	User_error("Cannot find file [%s]", name);

		lev = 2;
		//SCR cast
		ang = (float) ((SOMA_ANGLE / 180.0) * 3.141592654);
		dis = 1;
		thick = 1;
		
		Initialize_Neuron (seed_amount);
		
   	}  /* end non L-System algorithm */
   	
    if (burke_found) { /* Burke's algorithm */
 
    	Message("Burke-system file  : %s\n", name);
    	Get_Burke_Parameters (f);
		soma_d = (float) Grow_Neuron (BURKEMARKS);
	}
    else if (hillman_found) { /* Hillman */
 
      Message("Hillman-system file  : %s\n", name);
   
      Get_Hillman_Tamori_Parameters (f);
      
      soma_d = (float)Grow_Neuron (HILLMAN);
   
    }
    else if (tamori_found) {	/* check for Yoshihide Tamori's algorithm */
      
      Message("Tamori-system file  : %s\n", name);
		Get_Hillman_Tamori_Parameters (f);
    	soma_d = (float)Grow_Neuron (TAMORI);
	} /* if found Tamori */
		
	else { /* l-system */
		/* Read ls file and setup rules */
    	L_init();

		/* Execute mutations */
    	Get_comline_opt("u", &found, temp);
    	if (found) {
       		sscanf(temp, "%hd", &tim);
       		Message("Mutating       : ..");
       		for (i = 0; i < tim; i++) {
           		Message("\b\b%-2d", i + 1);
           		L_mutate();                    /* perform mutations on stored rule
    	           	                            * set */
        	}
        	Message("\n");
        	L_save();                          /* save the mutated ls file */
    	};

 		/* Create L-system production string */
    	L_system();
	} /* end else l-system */
	
 /* Parse production string and create geometry */
    L_draw();

 /* Add groundplane */
    Get_comline_opt("g", &found, temp);
    if (found)
        Ground_plane();

    Close_datafile();

    Message("\n");
}


/* -------------------------------------------------------------------------
   End of file.
   -------------------------------------------------------------------------
*/

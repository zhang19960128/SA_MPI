#include 	"string.h"
#include	"stddef.h"
#include	"stdlib.h"
#include  <stdio.h>
#include  <math.h>
#include  "simann.h"
struct {
    float u[97], c, cd, cm;
    int i97, j97;
} raset;

/* Subroutine */ int rmarin(int *ij, int *kl)
{
    /* Local variables */
    static int i, j, k, l, m;
    static float s, t;
    static int ii, jj;

/*  This subroutine and the next function generate random numbers. See */
/*  the comments for SA for more information. The only changes from the */
/*  orginal code is that (1) the test to make sure that RMARIN runs first */
/*  was taken out since SA assures that this is done (this test didn't */
/*  compile under IBM's VS Fortran) and (2) typing ivec as int was */
/*  taken out since ivec isn't used. With these exceptions, all following */
/*  lines are original. */
/* This is the initialization routine for the random number generator */
/*     RANMAR() */
/* NOTE: The seed variables can have values between:    0 <= IJ <= 31328 */
/*                                                      0 <= KL <= 30081 */
    if (*ij < 0 || *ij > 31328 || *kl < 0 || *kl > 30081) {
	printf(" The first random number seed must have a value between 0 and 31328");
	printf(" The second seed must have a value between 0 and 30081");
	return(0);
    }
    i = *ij / 177 % 177 + 2;
    j = *ij % 177 + 2;
    k = *kl / 169 % 178 + 1;
    l = *kl % 169;
    for (ii = 1; ii <= 97; ++ii) {
	s = 0.f;
	t = .5f;
	for (jj = 1; jj <= 24; ++jj) {
	    m = i * j % 179 * k % 179;
	    i = j;
	    j = k;
	    k = m;
	    l = (l * 53 + 1) % 169;
	    if (l * m % 64 >= 32) {
		s += t;
	    }
	    t *= .5f;
	}
	raset.u[ii - 1] = s;
    }
    raset.c = .021602869033813477f;
    raset.cd = .45623308420181274f;
    raset.cm = .99999982118606567f;
    raset.i97 = 97;
    raset.j97 = 33;
    return 0;
} /* rmarin */

double ranmar(void)
{
    /* System generated locals */
    float ret_val;

    /* Local variables */
    static float uni;

    uni = raset.u[raset.i97 - 1] - raset.u[raset.j97 - 1];
    if (uni < 0.f) {
	uni += 1.f;
    }
    raset.u[raset.i97 - 1] = uni;
    --raset.i97;
    if (raset.i97 == 0) {
	raset.i97 = 97;
    }
    --raset.j97;
    if (raset.j97 == 0) {
	raset.j97 = 97;
    }
    raset.c -= raset.cd;
    if (raset.c < 0.f) {
	raset.c += raset.cm;
    }
    uni -= raset.c;
    if (uni < 0.f) {
	uni += 1.f;
    }
    ret_val = uni;
    return ret_val;
} /* ranmar */

double exprep(double rdum)
{
    /* System generated locals */
    double ret_val;

    /* Builtin functions */
    double exp(double);

/*  This function replaces exp to avoid under- and overflows and is */
/*  designed for IBM 370 type machines. It may be necessary to modify */
/*  it for other machines. Note that the maximum and minimum values of */
/*  EXPREP are such that they has no effect on the algorithm. */
    if (rdum > 174.f) {
	ret_val = 3.69e75;
    } else if (rdum < -180.f) {
	ret_val = 0.f;
    } else {
	ret_val = exp(rdum);
    }
    return ret_val;
} /* exprep */

double maxval(double x, double y)
{
   if (x >= y)
   {
      return x;

   }
   return y;
} /* maxval */

double minval(double x, double y)
{
   if (x <= y)
   {
      return x;

   }
   return y;
} /* minval */

//char    *get_line(line, len, file)
//char    *line;
//int     len;
//FILE    *file;
//{
//   char *s;
//   int  i;
//   do
//   {
//      s = fgets(line, len, file);               /* Read one line of input     */
//      if(s == NULL) break;                      /* exit if end of file        */
//      i = strlen(s) - 1;
//      while(i >= 0 && (s[i] == ' ' || s[i] == '\t' || s[i] == '\n'))
//         s[i--] = '\0';                         /* Strip trailing white space */
//   }
//   while(*s == '\0' || *s == '#');              /* Repeat if blank or comment */
//   if(s == NULL)
//      *line = '\0';                             /* Return null at eof         */
//   return(line);
//}

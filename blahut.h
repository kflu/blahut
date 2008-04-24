#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifndef __BLAHUT_H__
#define __BLAHUT_H__
struct blahut_constrained_capacity {
    /* These two members must be specified */
    gsl_matrix *	Q;	/* the foward transition matrix */
    gsl_vector *	e;	/* the expense vector */

    gsl_vector *	p;	/* the source probability distribution */
    gsl_matrix *	P;	/* the backward transition matrix */
    gsl_vector *	c;	/* the capacity vector c_j */

    /* a M x 2 matrix storing the Expense and corresponding Capacity.
       M = floor((s_U - s_L)/s_d), only used when iterating over
       s, otherwise this field is always set to NULL */
    gsl_matrix *	ce_curve;	

    int		numIn;	/* This is the # rows in Q matrix
				   and the # col. in P matrix */
    int 	numOut;	/* This is the # col. in Q matrix 
				   and the # rows in P matrix */

    double 	I_L;	/* the lower bound of I */
    double 	I_U;	/* the upper bound of I */
    double 	E;	/* the Expense */
    double 	C;	/* the Capacity */

    double 	s;	/* the parameter s */
    /* These three are only used when calculate the C(E) curve 
     * (iterates on s) */
    double	s_L;	/* the lowest value of s */
    double 	s_U;	/* the highest value of s */
    double 	s_d;	/* the iteration step on s */


    double 	epsilon;/* the error toleration */

    unsigned int	it;	/* the interation index */
    unsigned int	maxNumIt; /* the max # of iterations allowed.
				     The default value set in 
				     blahut_cap_init() is 
				     the limit of unsigned int: UINT_MAX */
    int		exceedsMaxNumIt;  /* this flag is set to true if after maxNumIt 
				     iterations, the terminating creteria still 
				     haven't been met. */
};
typedef struct blahut_constrained_capacity blahut_cap; /* a shorthand */

blahut_cap * 
blahut_cap_init( const gsl_matrix* Q, 
		   const gsl_vector* e );

void 
blahut_cap_free( blahut_cap* cap);

/*
blahut_cap * 
blahut_cap_set_p( blahut_cap * cap,
		    const gsl_vector* p );
*/

blahut_cap *
blahut_cap_calc( blahut_cap * cap );

blahut_cap * 
blahut_cap_set_p_uniform( blahut_cap * cap );

gsl_matrix * 
blahut_cap_iterate_over_s( blahut_cap * cap, const char* filename);

blahut_cap * 
blahut_cap_setSRange(blahut_cap * cap, double s_L, double s_U, double step);
#endif /*  __BLAHUT_H__ */


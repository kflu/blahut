/*
 * -------------------------------------------------------------------------
 *
 * Blahut Algorithm for Estimating Channel Capacity
 *
 * Copyright (C) 2008 - Kefei Lu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifndef __BLAHUT_H__
#define __BLAHUT_H__

/*
 * This type is used to determine the unit used to measure information.
 */
typedef enum {
    BITS, 
    NATS
} blahut_unit;

/* This structure stores the estimated capacity-exapense curve
 * and the optimizing input probability distribution.
 * The length of the three vectors are the same. */
struct blahut_ce_curve_str {
    unsigned int len;
    gsl_matrix * p;
    gsl_vector * E;
    gsl_vector * C;
};
typedef struct blahut_ce_curve_str blahut_ce_curve;

struct blahut_constrained_capacity {
    /* Unit used to express the channel capacity */
    blahut_unit	unit;
    /* These two members must be specified */
    gsl_matrix *	Q;	/* the foward transition matrix */
    gsl_vector *	e;	/* the expense vector */

    gsl_vector *	p;	/* the source probability distribution */
    gsl_matrix *	P;	/* the backward transition matrix */
    gsl_vector *	c;	/* the capacity vector c_j */

    /* A M x (2+numIn) matrix storing the Expense and corresponding Capacity, 
     *   and the maximizing p vector. M = floor((s_U - s_L)/s_d). 
     * Only used when iterating over s, DO NOT use this varaible otherwise. */
    /*
    gsl_matrix *	ce_curve;	
    */
    blahut_ce_curve 	ce_curve;

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

blahut_ce_curve 
blahut_cap_iterate_over_s( blahut_cap * cap, const char* filename);

blahut_cap * 
blahut_cap_setSRange(blahut_cap * cap, double s_L, double s_U, double step);
#endif /*  __BLAHUT_H__ */


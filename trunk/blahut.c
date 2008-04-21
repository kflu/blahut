#include <assert.h>
#include "blahut.h"

#define DEBUG 
// #define DEBUG_PRINT_WARNING /* print warning message */
#define DEFAULT_MAX_IT UINT_MAX
#define DOUBLE_COMP_LIMIT 1e8

#ifdef XXXXXXXXXXXXXX
/* 
 * Many systems already have log2 function. If not, remove the 
 * ifdef, endif preprocessor to enable this function.
 */
inline static double log2(const double value)
{
    return log10(value)/log10(2);
}
#endif

static int 
vector_isnonneg(const gsl_vector* vec)
{
    unsigned int i;
    for (i=0;i<vec->size;i++) {
	if (gsl_vector_get(vec,i) < 0) {
	    return 0;
	}
    }
    return 1;
}

static int 
matrix_isnonneg(const gsl_matrix* mat)
{
    unsigned int j,k;
    for (j=0;j<mat->size1;j++) {
	for (k=0;k<mat->size2;k++) {
	    if (gsl_matrix_get(mat,j,k) < 0.0) {
		return 0;
	    }
	}
    }
    return 1;
}

blahut_cap * 
blahut_cap_init( const gsl_matrix* Q, 
		   const gsl_vector* e )
{
    unsigned int i=0, k=0;
    blahut_cap * cap = (blahut_cap*) malloc (sizeof(blahut_cap));
    if (!cap) {
	fprintf(stderr, "Not enough memory when allocating a blahut_cap.\n");
	exit(1);
    }
    memset(cap, 0, sizeof(blahut_cap));

    /* Check validity of Q and e */
    if (Q == NULL || e == NULL) {
	fprintf(stderr, "Q or e is NULL pointer.\n");
	exit(1);
    }
    if (Q->size1 != e->size) {
	fprintf(stderr, "Q's # rows is not equal to e's size.\n");
	exit(1);
    } else if (Q->size2 <= 0) {
	fprintf(stderr, "Q's # columns should not be negative.\n");
	exit(1);
    } else if (!vector_isnonneg(e) || !matrix_isnonneg(Q)) {
	fprintf(stderr, "Q or e contains negative elements.\n");
	exit(1);
    }

    /* Check the validity of Q and e */
    double sum_e=0, sum_Q=0;

#ifdef XXXXXXXXXXXXXXXX
    for (i=0;i<e->size;i++) {
	sum_e += gsl_vector_get(e,i);
    }
    if (fabs(sum_e - 1.0) > DOUBLE_COMP_LIMIT) {
	fprintf(stderr, "Sum of e seems not to be 1.\n");
	exit(2);
    }
#endif

    /* each row of Q should sum to 1 */
    for (i=0;i<Q->size1;i++) {
	sum_Q=0;
	for (k=0;k<Q->size2;k++) {
	    sum_Q += gsl_matrix_get(Q, i, k);
	}
	if (fabs(sum_Q - 1.0) > DOUBLE_COMP_LIMIT) {
	    fprintf(stderr, "Sum over row %d of Q seems not to be 1.\n", k);
	    exit(2);
	}
    }

    cap->Q = (gsl_matrix*)Q;
    cap->e = (gsl_vector*)e;
    
    /* Process Q and e */
    cap->numIn = Q->size1;
    cap->numOut = Q->size2;

    /* Initialize p, P, c */
    cap->p = gsl_vector_alloc(cap->numIn);
    gsl_vector_set_all(cap->p, 1.0/cap->numIn); /* init. the input distribution
					      as uniform */
    cap->P = gsl_matrix_calloc(cap->numOut, cap->numIn); /* numOut x numIn */
    cap->c = gsl_vector_calloc(cap->numIn); /* numIn x 1 */

    /* The default values */
    cap->s_L = 0.0;
    cap->s_U = 1e4;
    cap->s_d = 0.001;

    cap->epsilon = 1e-5;
    cap->maxNumIt = DEFAULT_MAX_IT;

    return cap;
}

void 
blahut_cap_free( blahut_cap* cap)
{
    gsl_matrix_free(cap->P);
    gsl_vector_free(cap->p);
    gsl_vector_free(cap->c);
    gsl_matrix_free(cap->ce_curve);

    free(cap);
}

/* Calculate:
 * sum_j {p_j * Q_{k|j}} */
inline static double
sum_p_Q (const blahut_cap * cap, int k)
{
    register int j = 0;
    register double sum = 0;
    for (j; j<cap->numIn; j++) {
	sum += gsl_vector_get(cap->p, j) * gsl_matrix_get(cap->Q, j, k);
    }
    return sum;
}

/* Calculate the sum of the first term of exp(...) in c_j
 * expression */
inline static double 
sum_Q_log (const blahut_cap * cap, int j)
{
    register int k = 0;
    register double Q_kj;
    register double sum = 0;

    for (k; k<cap->numOut; k++) {
	Q_kj = gsl_matrix_get(cap->Q, j, k);
	/* use the convention that 0log0 = 0 */
	sum += (Q_kj == 0 ? 0 : Q_kj * log2 (Q_kj/sum_p_Q(cap,k)));
    }

    return sum;
}

/* Calculate c_j over all j */
static blahut_cap * calc_c_j ( blahut_cap * cap )
{
    int j;
    for (j=0; j<cap->numIn; j++) {
	gsl_vector_set(cap->c, j, 
		exp(sum_Q_log(cap, j)
		    - cap->s * gsl_vector_get(cap->e,j)));
    }
    return cap;
}

inline static double calc_I_L (blahut_cap * cap)
{
    int j;
    double sum=0;
    for (j=0; j<cap->numIn; j++) {
	sum += gsl_vector_get(cap->p,j) * gsl_vector_get(cap->c,j);
    }
    cap->I_L = log2(sum);
    return cap->I_L;
}

inline static double calc_I_U (blahut_cap * cap)
{
    cap->I_U = log2(gsl_vector_max ( cap->c ));
    return cap->I_U;
}

static blahut_cap * update_p (blahut_cap * cap)
{
    double sum=0;
    int j;
    for(j=0;j<cap->numIn;j++) {
	sum += gsl_vector_get(cap->p,j) * gsl_vector_get(cap->c,j);
    }

    for (j=0;j<cap->numIn; j++) {
	gsl_vector_set(cap->p,j, 
		gsl_vector_get(cap->p,j)
		* (gsl_vector_get(cap->c,j) / sum));
    }

    return cap;
}

static double calc_E(blahut_cap * cap)
{
    int j;
    double sum = 0;
    for(j=0;j<cap->numIn; j++) {
	sum += gsl_vector_get(cap->p,j) * gsl_vector_get(cap->e,j);
    }
    cap->E = sum;
    return sum;
}

inline static double calc_C(blahut_cap * cap)
{
    cap->C = cap->s * cap->E + cap->I_L;
    return cap->C;
}
    
blahut_cap *
blahut_cap_calc( blahut_cap * cap )
{
    for(cap->it=0 ; cap->it < cap->maxNumIt ; cap->it++ ) {
	calc_c_j(cap);
	if (calc_I_U(cap) - calc_I_L(cap) < cap->epsilon) {
	    calc_E(cap);
	    calc_C(cap);
	    break;
	} else if ( isnan(cap->I_U) || isnan(cap->I_L) ) {
#ifdef DEBUG_PRINT_WARNING
	    fprintf(stdout, "I_U or I_L is NaN, terminate loop.\n");
#endif
	    break;
	}
	update_p(cap);
    }

    /* to check if for() is terminated by a break or by
     * cap->it >= cap->maxNumIt */
    if (cap->it >= cap->maxNumIt) {
	/* for() is terminated by exceeding the max # iterations */
	calc_E(cap);
	calc_C(cap);
	cap->exceedsMaxNumIt = 1;
    }

    return cap;
}

inline blahut_cap * 
blahut_cap_set_p_uniform( blahut_cap * cap )
{
    gsl_vector_set_all(cap->p, 1.0/cap->p->size);
    return cap;
}

gsl_matrix * 
blahut_cap_iterate_over_s( blahut_cap * cap, const char* filename)
{
    unsigned int NumS = floor((cap->s_U - cap->s_L)/cap->s_d);
    cap->ce_curve = gsl_matrix_calloc(NumS,2);

#ifdef DEBUG
    fprintf(stdout, "Calculating C(E) curve ...\n");
#endif

    unsigned int i;
    double step = cap->s_d;
    for (i=0,cap->s = cap->s_L; i<NumS; i++, cap->s+=step) {
	//printf("%d, %g\n",i,cap->s);
	blahut_cap_set_p_uniform(cap);
	blahut_cap_calc(cap);

	/* stop iterating */
	if (cap->C == 0 || cap->E == 0) {
#ifdef DEBUG_PRINT_WARNING
	    fprintf(stdout, "C = 0 or E = 0 encountered, break the iteration.\n");
#endif
	    break;
	}
	gsl_matrix_set(cap->ce_curve, i, 0, cap->E);
	gsl_matrix_set(cap->ce_curve, i, 1, cap->C);
    }
#ifdef DEBUG
    fprintf(stdout, "Finished.\n");
#endif

    /* Write to file */
    if (filename != NULL) {
#ifdef DEBUG
	fprintf(stdout, "Writing data to file ...\n");
#endif
	FILE * file = fopen(filename,"w");
	for (i=0;i<cap->ce_curve->size1;i++) {
	    fprintf(file, "%g ", gsl_matrix_get(cap->ce_curve, i, 0));
	    fprintf(file, "%g ", gsl_matrix_get(cap->ce_curve, i, 1));
	    fprintf(file, "\n");
	}
	fclose(file);
#ifdef DEBUG
	fprintf(stdout, "Finished.\n");
#endif
    }

    return cap->ce_curve;
}

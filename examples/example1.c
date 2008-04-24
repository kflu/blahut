#include <stdio.h>
#include <stdlib.h>

#include "blahut.h"

int main(int argc, char * argv[])
{
    /* At most one argument from command line is allowd to specify
     * the cross over probability */
    char ** end;
    double q;	/* Setting q */
    if (argc > 2) {
	fprintf(stderr, "(E) At most one argument is allowd to specify"
		"the cross over probability p\n");
	exit(1);
    } else if (argc == 2) {
	q = strtod(argv[1], end);
	if ( q==HUGE_VAL||q==-HUGE_VAL ) {
	    fprintf(stderr, "(E) Double overflows in %s.\n", argv[1]);
	    exit(1);
	} else if (*end == argv[1]) {
	    fprintf(stderr, "(E) Not a valid double: %s.\n", argv[1]);
	    exit(1);
	}
	fprintf(stdout, "(I) Cross over probability set to %g.\n", q);
    } else if (argc == 1) {
	fprintf(stdout, "(W) Default value of q is used.\n");
	q = 0.2;	/* Take the default value */
    } else {
	/* no else */
	exit(10);
    }

    /* Specify the forward transition matrix
     * and the expense schedule vector.
     *
     * Note that in C and gsl, matrix are row-major,
     * meaning that eacho row forms a contiguous memory
     * block*/
    double QQ[] = {1, 0,  0, 0,
		   0, q,1-q, 0,
		   0, 0,  0, 1};
    double ee[] = {0, 1, 0};

    /* Prepare a gsl_matrix representation of Q,
     * and a gsl_vector representation of e.
     * For more information on gsl vector/matrix,
     * refer to gsl manual */
    gsl_matrix_view Q_view = gsl_matrix_view_array(QQ,3,4);
    gsl_matrix * Q = &(Q_view.matrix);

    gsl_vector_view e_view = gsl_vector_view_array(ee,3);
    gsl_vector * e = &(e_view.vector);

    /* Initialize a muanipulating object of Blahut algorithm */
    blahut_cap * cap = blahut_cap_init(Q,e);
    blahut_cap_setSRange(cap, 0, 1000, 0.01);

    /* Begin the estimation of C(E) curve.
     * The curve is stored in gsl_matrix * cap->ce_curve 
     * Meanwhile the curve is also write into file "ce.txt" */
    blahut_cap_iterate_over_s(cap, "ce.txt"); 

    /* Free the memory allocated for cap */
    blahut_cap_free(cap);

    exit(0);
}


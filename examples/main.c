#include "blahut.h"
#define NDEBUG

int main()
{
    double p = 0; /* cross over prob */
    // double QQ[4] = {1-p,p,p,1-p};
    // gsl_matrix_view Q = gsl_matrix_view_array(QQ,2,2);

    // double ee[2] = {1,0};
    // gsl_vector_view e =gsl_vector_view_array(ee, 2);

    // blahut_cap * cap = blahut_cap_init(&(Q.matrix),&(e.vector));
    
    double QQ[6]={.5,.3,.2,.3,.3,.4};
    gsl_matrix_view Q = gsl_matrix_view_array(QQ,2,3);

    double ee[2] = {1,0};
    gsl_vector_view e =gsl_vector_view_array(ee, 2);

    blahut_cap * cap = blahut_cap_init(&(Q.matrix),&(e.vector));
    cap->s_d = 0.001;




    blahut_cap_iterate_over_s(cap,"expense-capacity.dat");

    blahut_cap_free(cap);

    return 0;
}

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/**
 * The measurement function
 */
void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht,gsl_vector *y){

    gsl_matrix_set(Ht,0,0,1.0);
    gsl_matrix_set(Ht,1,0,param[0]);
    gsl_matrix_set(Ht,2,0,param[1]);
    gsl_matrix_set(Ht,3,1,1.0);
    gsl_matrix_set(Ht,4,1,param[2]);
    gsl_matrix_set(Ht,5,1,param[3]);
	
	gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);

}

void function_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *param, size_t n_gparam,const gsl_vector *co_variate,
        void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *),
        gsl_vector *x_tend){
			
		/*abs->fabs, exp=exp*/	
		gsl_vector_set(x_tend,0,param[4]*gsl_vector_get(xstart,0)+param[6]*(exp(fabs(gsl_vector_get(xstart,1)))/(1+exp(fabs(gsl_vector_get(xstart,1)))))*gsl_vector_get(xstart,1));
		gsl_vector_set(x_tend,1,param[5]*gsl_vector_get(xstart,1)+param[7]*(exp(fabs(gsl_vector_get(xstart,0)))/(1+exp(fabs(gsl_vector_get(xstart,0)))))*gsl_vector_get(xstart,0));	

}


void function_jacob_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart,
        double *param, size_t num_func_param, const gsl_vector *co_variate,
        void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *),
		gsl_matrix *Jx){
        
		/*the sign function, ^2->pow(double, double)*/
		double sign_x0=(gsl_vector_get(xstart,0)>0?1.0:0.0)-(gsl_vector_get(xstart,0)<0?1.0:0.0);
		double sign_x1=(gsl_vector_get(xstart,1)>0?1.0:0.0)-(gsl_vector_get(xstart,1)<0?1.0:0.0);
		gsl_matrix_set(Jx,0,0,param[4]);
		gsl_matrix_set(Jx,0,1,(param[6]*fabs(gsl_vector_get(xstart,1)))/(fabs(gsl_vector_get(xstart,1)) + 1.0) + (param[6]*gsl_vector_get(xstart,1)*sign_x1)/(fabs(gsl_vector_get(xstart,1)) + 1.0)- (param[6]*gsl_vector_get(xstart,1)*fabs(gsl_vector_get(xstart,1))*sign_x1)/pow((fabs(gsl_vector_get(xstart,1))+1.0),2.0));
		gsl_matrix_set(Jx,1,0,(param[7]*fabs(gsl_vector_get(xstart,0)))/(fabs(gsl_vector_get(xstart,0)) + 1.0) + (param[7]*gsl_vector_get(xstart,0)*sign_x0)/(fabs(gsl_vector_get(xstart,0)) + 1.0) - (param[7]*gsl_vector_get(xstart,0)*fabs(gsl_vector_get(xstart,0))*sign_x0)/pow((fabs(gsl_vector_get(xstart,0))+1.0),2.0));
		gsl_matrix_set(Jx,1,1,param[5]);

}

/**
 * Set the initial condition
 */

void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector *pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0){

    gsl_vector_set(pr_0,0,1);

    size_t num_regime=pr_0->size;
    size_t dim_latent_var=error_cov_0[0]->size1;
    size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
    /*printf("%lu %lu %lu\n",num_regime,dim_latent_var,num_sbj);*/
    size_t i,j;
    for(j=0;j<num_regime;j++){
        for(i=0;i<num_sbj;i++){
            /*printf("%lu %lu\n",i,j);*/
            gsl_vector_set((eta_0)[j],i*dim_latent_var,0.0);
            gsl_vector_set((eta_0)[j],i*dim_latent_var+1,0.0);
        }/*statevar_1_p1 statevar_2_p1 statevar_1_p2 statevar_2_p2 ..., eta_0[] with a length of num_sbj*dim_latent_var*/

        gsl_matrix_set((error_cov_0)[j],0,0,1.0);
        gsl_matrix_set((error_cov_0)[j],1,1,1.0);
    }
}

/**
 * Set the regime-switch transition probability matrix
 */

void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){
     /*if there is only one regime*/
	gsl_matrix_set_identity(regime_switch_mat);
}

/**
 * Set the noise covariance matrix
 *They are to be LDL' transformed.
 *e.g., [a b
         b c]
 *-->LDL', L=[1 0;b 1], D=diag(a,c)
 */

void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){
    
    gsl_matrix_set(y_noise_cov,0,0, param[8]);
    gsl_matrix_set(y_noise_cov,1,1, param[9]);
    gsl_matrix_set(y_noise_cov,2,2, param[10]);
    gsl_matrix_set(y_noise_cov,3,3, param[11]);
    gsl_matrix_set(y_noise_cov,4,4, param[12]);
    gsl_matrix_set(y_noise_cov,5,5, param[13]);
	
	gsl_matrix_set(eta_noise_cov,0,0,param[14]);
	gsl_matrix_set(eta_noise_cov,1,1,param[15]);

}

/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 */
void function_transform(double *param){

}


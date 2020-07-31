/*This code will run the dynamics for active Brownian particles (ABPs) and calculate the bound in equation 7 from the paper "Entropy production fluctuations encode collective behavior in active matter" by Trevor GrandPre, Katie Klymko, Kranthi K. Mandadapu, and David T. Limmer.

This code was written by Katie Klymko and Trevor GrandPre. 

The inputs are lambda, the self-propulsion v, and time in that order. In order to compile, you will need GNU scientific library. */


#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>	
#include <algorithm>	
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

using namespace std;

/*within this code I will run ABPS with interactions and measure the coefficient froma running average and get the lambda dependence as well.*/

int main(int argc, char *argv[]) {

/*define variables and vectors*/
    double time1 = atof(argv[3]);//Time 
	double Pe = atof(argv[2]);//self-propulsion
    double lams = atof(argv[1]);//lambda bias in the form of lambda*N where N is the number of particles.
  	double h = 0.0001;//this is the timestep
  	int Nstep =int(time1/h); //number of total steps.
  	int i,j,m;
    int N=10;//number of partilcles
    double L=5;// box size is 2*L
    double den=0.1;//the density
	double beta = 1.0;//invere temperature
	double gamma = 1.0;//friction
	double T = 1.0;//Temperature.
    double rcut=pow(2,1./6.);//cutoff for WCA
    double ecut = 1; //this is epsilon in the WCA formula
    //setting the observables equal to zero
    double fluxp=0;
    double flux=0;
    
  	vector< vector<double> > f(N, vector<double>(2));
    vector<double> x(N);
    vector<double> y(N);
    vector<double> theta(N);
    //defining the non-interacting CGF for ABPs
    double pf=Pe*Pe*(lams/N)+Pe*Pe*(lams/N)*(lams/N);
    //create a file to write into called "density.txt". 
   FILE * den1;
    den1 = fopen("density.txt","w");
//read in the initial conditions from file named "equi"
FILE *init;//restart file
 	if ( ( init = fopen( "equi", "r" ) ) == NULL){
   		printf ("initial data could not be opened\n");}
 	else {
      for(i=0;i<N;i++){
				fscanf(init, "%lf %lf %lf", &x[i], &y[i], &theta[i]);
     	}
		
   	 }
   	rewind(init);
   	fclose(init);
//generate a random number to start the stochastic dynamics. 
   gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL));

	srand(time(0));



	for(i=0;i<Nstep;i++){//time loop
 
       //Calculate forces

		for (m=0; m<N; m++) {
			f[m][0] = 0;
			f[m][1] = 0;
		}

		for (m=0; m<N-1; m++) {
			for (j=m+1; j<N; j++) {
				double xr = x[m]-x[j];
				if(xr>2*L*0.5) {
					 xr-=2*L;}
				else if(xr<(-2*L*0.5)) {
					 xr=2*L+xr;}
				double yr = y[m]-y[j];
				if(yr>2*L*0.5) {
					 yr-=2*L;}
				else if(yr<(-2*L*0.5)) {
					 yr=2*L+yr;}

				double r2 = xr*xr+yr*yr,r = sqrt(r2);

				if (r < rcut) {

						double r2i = 1/r2;
						double r6i = pow(r2i,3);
						double ff = 48*r2i*r6i*(r6i-0.5); // L-J force

						f[m][0] = f[m][0]+ff*xr; // update pair force in components
						f[m][1] = f[m][1]+ff*yr;

						f[j][0] = f[j][0]-ff*xr;
						f[j][1] = f[j][1]-ff*yr;

				}
		}
	}
             
        
    
		//update positions
		for(j=0;j<N;j++){
       
            
            double fd_term = sqrt(2 * h *T/ (gamma ));
			double fd_term_theta = sqrt(3. * 2. * h *T/ (gamma ));
            
			double noise_theta = fd_term_theta * gsl_ran_gaussian(rng,1);
      		double noise_0 = fd_term * gsl_ran_gaussian(rng,1);
      		double noise_1 = fd_term * gsl_ran_gaussian(rng,1);

            
            
      		x[j] +=h / gamma * f[j][0]+ noise_0+Pe*cos(theta[j])*h+2*Pe*cos(theta[j])*h*(lams/N);
            y[j] +=h / gamma * f[j][1]+ noise_1+Pe*sin(theta[j])*h+2*Pe*sin(theta[j])*h*(lams/N);
      		theta[j] += noise_theta;
            
        
			//periodic boundaries o f 2L%
            if(x[j]>L){ x[j]=-L+(x[j]-L); }
            if(x[j]<-L){ x[j]=L+(L+x[j]); }
            if(y[j]>L){ y[j]=-L+(y[j]-L); }
            if(y[j]<-L){ y[j]=L+(L+y[j]); }

		
/*calculate the observable*/

            flux+=(Pe*cos(theta[j]))*(h / gamma * f[j][0])+Pe*sin(theta[j])*(h / gamma * f[j][1]);
            fluxp+=(Pe*cos(theta[j]))*(h / gamma * f[j][0]+ noise_0+Pe*cos(theta[j])*h+2*Pe*cos(theta[j])*h*(lams/N))+(Pe*sin(theta[j]))*(h / gamma * f[j][1]+ noise_1+Pe*sin(theta[j])*h+2*Pe*sin(theta[j])*h*(lams/N));
          
            
            
	                       }   
        printf("%lf %lf %lf %lf %lf\n",i*h,flux/(Pe*Pe)/(i*h)/(den*N),fluxp/N/(i*h),pf,pf+(lams/N)*flux/(i*h)/N);
         // time, effective drage, observable-EP, free CGF, bound
                           
                        }
    //save the running average of all the quantities to a file, the density.txt.
       fprintf(den1,"%lf %lf %lf %lf\n",i*h,flux/(Pe*Pe)/(i*h)/(den*N),fluxp/N/(i*h),pf+(lams/N)*flux/(i*h)/N);



}

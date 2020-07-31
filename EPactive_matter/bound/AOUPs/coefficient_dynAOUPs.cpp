/*This code will run the dynamics for active Ornstein Uhlenbeck particles (AOUPs) and calculate the bound in equation 7 from the paper "Entropy production fluctuations encode collective behavior in active matter" by Trevor GrandPre, Katie Klymko, Kranthi K. Mandadapu, and David T. Limmer (2020).

This code was written by Katie Klymko and Trevor GrandPre. 

The inputs are lambda, active diffusion constant, and time of system in that order. In order to compile, you will need GNU scientific library. */


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


int main(int argc, char *argv[]) {

/*define variables and vectors*/
    double lams = atof(argv[1]);//lambda bias in the form of lambda*N where N is the number of particles.
    double da = atof(argv[2]);//active diffusion 
    double time1 = atof(argv[3]);//Time 


  	double h = 0.0001;//this is the timestep

  	int Nstep =int(time1/h); //total number of steps
  	int i,j,m;
    int N=20;//number of particles
    double rho=0.1;//density of the syste 
    double  L=7.07106781187;//length of system [-L L]
	double T = 1.0;
    double rcut=pow(2,1./6.);
    double ecut = 1; //this is epsilon in the WCA formula
    
    //initialize quantities
    double fluxp=0;
    double flux=0;
   double fluxint=0;
  	vector< vector<double> > f(N, vector<double>(2));
    vector<double> x(N);
    vector<double> y(N);
    vector<double> theta(N);
    vector<double> theta1(N);
    vector<double> theta2(N);
     
   FILE * den1;
    den1 = fopen("density.txt","w");
//read in the initial conditions
FILE *init;//restart file
 	if ( ( init = fopen( "pe_3.3_den0.5.txt", "r" ) ) == NULL){
   		printf ("initial data could not be opened\n");}
 	else {
      for(i=0;i<N;i++){
				fscanf(init, "%lf %lf %lf", &x[i], &y[i], &theta[i]);
          theta1[i]=0.0;
          theta2[i]=0.0;
          //set both Fps initially to zero.
     	}
		
   	 }
   	rewind(init);
   	fclose(init);
//generate a random number
   gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL));

	srand(time(0));



	for(i=0;i<Nstep;i++){
        /*calculate the forces every time step*/
       

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
             
   //set the rotational and translational diffusion constants to 1
        int dr=1;
        int dt=1;
        //the non-interacting CGF
        double p=dr*(1-sqrt(1-4*da*(lams/N)-4*da*(lams/N)*(lams/N)));
		//update positions
		for(j=0;j<N;j++){
            
            double fd_term = sqrt(2 * h *T);
			double fd_term_theta = sqrt(2.*da*h*dr*dr);
			double noise_theta1 = fd_term_theta * gsl_ran_gaussian(rng,1);
            double noise_theta2 = fd_term_theta * gsl_ran_gaussian(rng,1);
      		double noise_0 = fd_term * gsl_ran_gaussian(rng,1);
      		double noise_1 = fd_term * gsl_ran_gaussian(rng,1);

            

      		x[j] +=noise_0+theta1[j]*h+h* f[j][0]+2*theta1[j]*(lams/N)*h;
            y[j] +=noise_1+theta2[j]*h+h* f[j][1]+2*theta2[j]*(lams/N)*h;
            theta1[j] += -dr*theta1[j]*h+noise_theta1+theta1[j]*p*h;
      		theta2[j] += -dr*theta2[j]*h+noise_theta2+theta2[j]*p*h;

        
			//periodic boundaries o f 2L%
            if(x[j]>L){ x[j]=(x[j]-2*L); }
            if(x[j]<-L){ x[j]=(2*L+x[j]); }
            if(y[j]>L){ y[j]=(y[j]-2*L); }
            if(y[j]<-L){ y[j]=(2*L+y[j]); }

		
/*calculate the observable*/
            //non-interacting CGF in terms of a time average
            flux+=((((lams/N)*theta1[j]*theta1[j])/(dt))+(((lams/N)*(lams/N)*theta1[j]*theta1[j])/(dt))+((p*p*theta1[j]*theta1[j])/(4*da*dr*dr))-((p*theta1[j]*theta1[j])/(2*da*dr))+p)*h/2+((((lams/N)*theta2[j]*theta2[j])/(dt))+(((lams/N)*(lams/N)*theta2[j]*theta2[j])/(dt))+((p*p*theta2[j]*theta2[j])/(4*da*dr*dr))-((p*theta2[j]*theta2[j])/(2*da*dr))+p)*h/2;
      //the EP
            fluxp+=theta1[j]*(noise_0+theta1[j]*h+h* f[j][0]+2*theta1[j]*(lams/N)*h)+theta2[j]*(noise_1+theta2[j]*h+h* f[j][1]+2*theta2[j]*(lams/N)*h);
       
            //the interactions that will be used for the bound 
            fluxint+=theta1[j]*(lams/N)*(h* f[j][0])/dt+theta2[j]*(lams/N)*(h* f[j][1])/dt;
        
	                       }
    //save to a file time , the observable-EP, the free CGF time average, free CGF prediction, difference between prediction and time average for free CGF, the interaction correction, the bound.      
     printf("%lf %lf %lf %lf %lf %lf %lf\n",i*h,fluxp/(i*h)/N,flux/(i*h)/(N),p,flux/(i*h)/(N)-p, fluxint/(i*h)/(N),(fluxint+flux)/(i*h)/(N));
     fprintf(den1,"%lf %lf %lf %lf %lf %lf %lf\n",i*h,fluxp/(i*h)/N,flux/(i*h)/(N),p,flux/(i*h)/(N)-p, fluxint/(i*h)/(N),(fluxint+flux)/(i*h)/(N));

 
                        }
    }







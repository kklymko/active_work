#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>	
#include <algorithm>	//std::copy
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <mpi.h>

using namespace std;

/*within this code I will run ABPS with interactions and measure the coefficient froma running average and get the lambda dependence as well.*/

int main(int argc, char *argv[]) {

/*define variables and vectors*/
    double lams = atof(argv[1]);//Temperature
    double da = atof(argv[2]);//bonding energy for 2grb2 to 
    double time1 = atof(argv[3]);//Time in reduced units
     double realizations = atof(argv[4]);

  	double h = 0.0005;//think this is the timestep
    //double T=10;
  	int Nstep =int(time1/h); //number of steps in the interval
  	int i,j,m;
    int N=10;
    double rho=0.1;
    double L=5;
	double beta = 1.0;
	double gamma = 1.0;
	double T = 1.0;
    double rcut=pow(2,1./6.);
    double ecut = 1; //this is epsilon in the WCA formula
    double fluxp=0;
    double flux=0;
   
    
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
 	if ( ( init = fopen( "equi", "r" ) ) == NULL){
   		printf ("initial data could not be opened\n");}
 	else {
      for(i=0;i<N;i++){
				fscanf(init, "%lf %lf %lf", &x[i], &y[i], &theta[i]);
          theta1[i]=theta[i];
          theta2[i]=theta[i];
          //set both Fps initially at theta of restart.
     	}
		
   	 }
   	rewind(init);
   	fclose(init);

   gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL));

	srand(time(0));

//for(int k=0 k,walkers;k++){
    
for(int k=1;k<=realizations;k++){

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
             
        
        //zero every time step for the instantaneous value
      //  double flux = 0.0;
        double p=-(-1+sqrt(1-4*da*(lams/N)-4*da*(lams/N)*(lams/N)));
        double ps2=-(-1+sqrt(1-4*da*(lams/N)-4*da*(lams/N)*(lams/N)))+(da*0.86*0.86*(lams/N)*(lams/N)*rho*rho)/((1-4*da*(lams/N)*(1+(lams/N))));
		//update positions
		for(j=0;j<N;j++){
            int dr=1;
            double fd_term = sqrt(2 * h *T/ (gamma ));
			double fd_term_theta = sqrt(2.*da*h);
			double noise_theta1 = fd_term_theta * gsl_ran_gaussian(rng,1);
            double noise_theta2 = fd_term_theta * gsl_ran_gaussian(rng,1);
      		double noise_0 = fd_term * gsl_ran_gaussian(rng,1);
      		double noise_1 = fd_term * gsl_ran_gaussian(rng,1);

            

      		x[j] +=noise_0+theta1[j]*h+h / gamma * f[j][0]+2*theta1[j]*(lams/N)*h;
            y[j] +=noise_1+theta2[j]*h+h / gamma * f[j][1]+2*theta2[j]*(lams/N)*h;
            theta1[j] += -dr*theta1[j]*h+dr*noise_theta1+theta1[j]*p*h;
      		theta2[j] += -dr*theta2[j]*h+dr*noise_theta2+theta2[j]*p*h;

        
			//periodic boundaries o f 2L%
            if(x[j]>L){ x[j]=-L+(x[j]-L); }
            if(x[j]<-L){ x[j]=L+(L+x[j]); }
            if(y[j]>L){ y[j]=-L+(y[j]-L); }
            if(y[j]<-L){ y[j]=L+(L+y[j]); }

		
/*calculate the observable*/
            //friction
            flux+=(theta1[j])*(h / gamma * f[j][0])+theta2[j]*(h / gamma * f[j][1]);
            //average
            fluxp+=theta1[j]*(noise_0+theta1[j]*h+h / gamma * f[j][0]+2*theta1[j]*(lams/N)*h)+theta2[j]*(theta2[j]*h+h / gamma * f[j][1]+2*theta2[j]*(lams/N)*h);
            //fluxp[j]+=0;
            
            
	                       }
      
        double dels=flux*(lams/N)*(1/(i*h))/N;
     
    printf("%lf %lf %lf %lf %lf %lf %lf\n",i*h,fluxp/N/(i*h)/k,flux/(2*da)/(i*h)/(rho*N)/k,p,dels/k,p+dels/k, ps2);
    fprintf(den1,"%lf %lf %lf %lf %lf %lf %lf\n",i*h,fluxp/N/(i*h)/k,flux/(2*da)/(i*h)/(rho*N)/k,p,dels/k,p+dels/k, ps2);
 
                        }
    }

//print the running average over walkers.

//printf("%d %lf\n",k,);

//}	
                     // fluxp+=flux/(i*h);//this is the running average.
                           
                           //print the instantaneous and running average of friction
                          // printf("%lf %lf %lf\n",i*h,flux/(Pe*Pe),fluxp/(Pe*Pe));





}

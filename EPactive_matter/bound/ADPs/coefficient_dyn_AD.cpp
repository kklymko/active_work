/*This code will run the dynamics for active dumbell particles (ADPs) and calculate the bound in equation 7 from the paper "Entropy production fluctuations encode collective behavior in active matter" by Trevor GrandPre, Katie Klymko, Kranthi K. Mandadapu, and David T. Limmer (2020).

This code was written by Katie Klymko and Trevor GrandPre. 

The inputs are lambda, the self-propulsion v, time, and Temperature of system in that order. In order to compile, you will need GNU scientific library. */


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
    double time1 = atof(argv[3]);//Time 
	double Pe = atof(argv[2]);//self-propulsion magnitude
    double T = atof(argv[4]);//temperature
    double lams = atof(argv[1]);//lambda bias in the form of lambda*N where N is the number of particles.
  	double h = 0.0001;//this is the timestep
  	int Nstep =int(time1/h); //total number of steps
    

	double D = 1.0;//set the diffusion coefficient which is set by kT/(m*gamma) 
	double kappa =100;// the spring constant
	double lr = 1.5; //the rest length of the spring. 
    double den=0.2;
    
  	int i,j,m;
    int N=20;//number of particles
    double L = 14.1421356237;//box length in L
	double beta = 1.0;//inverse temp
	double gamma = 1.0;//friction
    double rcut=pow(2,1./6.);//cutoff for WCA
    double ecut = 1; //this is epsilon in the WCA formula
    //setting observables to zero
    double fluxp=0;
    double flux=0;
    //initialize  forces, positions, and orientation vector d.
  	vector< vector<double> > F(N, vector<double>(2));
    vector<double> x(N);
    vector<double> y(N);
    vector<vector<double> > d(N, vector<double>(2));
    //non-interacting CGF for ADPs.
      double pfad=2*Pe*Pe*(lams/N)+2*Pe*Pe*(lams/N)*(lams/N);
    //create a file to write into
   FILE * den1;
    den1 = fopen("density.txt","w");
//read in the initial conditions
FILE *init;
 	if ( ( init = fopen( "restart_AD10M.txt", "r" ) ) == NULL){
	  printf ("initial data could not be opened\n");}
 	else {
	  for(i=0;i<N;i++){
	    fscanf(init, "%lf %lf", &x[i], &y[i]);
	  }
	  
	}
   	rewind(init);
   	fclose(init);
//generate a random seed
   gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL));

	srand(time(0));



	for(m=0;m<Nstep;m++){
        /*calculate the forces every time step*/
      
	  /*first calculate the spring forces and WCA forces between the two dumbells.*/
	  for (int i=0; i<N; i++) {//initialize the force fields to zero.
	    F[i][0] = 0;
	    F[i][1] = 0;
	    d[i][0]=0;
	    d[i][1]=0;
	  }
	  
	  for (int i=0; i<N-1; i++) {
	    for (int j=i+1; j<N; j++) {//Periodic boundaries of the separation distance between the particles of the same dumbell
       
            
	      double xr = x[i]-x[j];
	      if(xr>L*0.5) {
	        xr-=L;}
	      else if(xr<(-L*0.5)) {
	        xr=L+xr;}
	      double yr = y[i]-y[j];
	      if(yr>L*0.5) {
	        yr-=L;}
	      else if(yr<(-L*0.5)) {
	        yr=L+yr;}      
	      
	      double r2 = xr*xr+yr*yr,r = sqrt(r2);

	      if(j==i+1 & i%2==0){
	      double pre=(2.*kappa*(lr-r))/r;//harmonic bond when particles bond is nonzero. 
	      F[i][0] = F[i][0]+pre*xr; // update pair force in components
	      F[i][1] = F[i][1]+pre*yr;
	      
	      F[j][0] =F[j][0]-pre*xr;//Newtons third law.
	      F[j][1] =F[j][1]-pre*yr;
	      //double ecut = 0.0;
	      
	      //calculation of the direction for both of the particles of the dumbell
	      d[i][0]=xr/r;
	      d[i][1]=yr/r;
	      
	      d[j][0]=xr/r;
	      d[j][1]=yr/r;
	 
	      }
	     
	      if (r < rcut) {//WCA force
	        
	        double r2i = 1/r2;
	        double r6i = pow(r2i,3);
	        double ff = 48*ecut*r2i*r6i*(r6i-0.5); // L-J force
	        
	        F[i][0] = F[i][0]+ff*xr; // update pair force in components
	        F[i][1] = F[i][1]+ff*yr;
	        
	        F[j][0] =F[j][0]-ff*xr;//Newtons third law.
	        F[j][1] =F[j][1]-ff*yr;
	      }
	      
	      
	     
	    }
	  }
		for(j=0;j<N;j++){
       
            
           double fd_term = sqrt(2 *D* h);//noise of the brownian motion
	    double noise_0 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);
	    double noise_1 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);
	    
        x[j] +=(h / gamma) * F[j][0]+ Pe*d[j][0]*h+noise_0+2*Pe*d[j][0]*h*(lams/N);
	    y[j] +=(h / gamma) * F[j][1]+ Pe*d[j][1]*h+noise_1+2*Pe*d[j][1]*h*(lams/N);
	    
	    
	    //periodic boundaries for 0 to L
	    if (x[j] > L) {x[j] -= (floor(x[j]/L))*L;}
	    if (x[j] < 0) {x[j] += (1+(floor(-x[j]/L)))*L;}
	    if (y[j] > L) {y[j] -= (floor(y[j]/L))*L;}
	    if (y[j] < 0) {y[j] += (1+(floor(-y[j]/L)))*L;}

		
/*calculate the observable*/
            
         

            flux+=(Pe*d[j][0])*(h / gamma * F[j][0])+Pe*d[j][1]*(h / gamma * F[j][1]);
         
            
            fluxp+=((Pe*d[j][0])*((h / gamma) * F[j][0]+ Pe*d[j][0]*h+noise_0+2*Pe*d[j][0]*h*(lams/N))+(Pe*d[j][1])*((h / gamma) * F[j][1]+ Pe*d[j][1]*h+noise_1+2*Pe*d[j][1]*h*(lams/N)));
      
            
            
	                       }
    //save to a file time, effective drage, observable-EP, non-interacting CGF, bound
        printf("%lf %lf %lf %lf %lf\n",m*h,flux/(m*h)/(Pe*Pe)/2/(den*N),fluxp/(N/2)/(m*h),pfad,pfad+(lams/N)*flux/(m*h)/N);
         fprintf(den1,"%lf %lf %lf\n",m*h,flux/(Pe*Pe)/(m*h),lams*flux/(Pe*Pe)/(m*h));
                           
                        }





}

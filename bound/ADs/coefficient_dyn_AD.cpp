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
    double time1 = atof(argv[3]);//Time in reduced units
	double Pe = atof(argv[2]);
    double T = atof(argv[4]);
    double lams = atof(argv[1]);//Temperature
  	double h = 0.0001;//think this is the timestep
    //double T=10;
  	int Nstep =int(time1/h); //number of steps in the 
    
    //direction in the x and y
	double D = 1.0;//set the diffusion coefficient which is set by kT/(m*gamma) 
	double kappa =100;// the spring constant
	double lr = 1.5; //the rest length of the spring. 
    double den=0.2;
    
  	int i,j,m;
    int N=20;
    double L = 14.1421356237;
	double beta = 1.0;
	double gamma = 1.0;
	//double T = 1.0;
    double rcut=pow(2,1./6.);
    double ecut = 1; //this is epsilon in the WCA formula
    double fluxp=0;
    double flux=0;
    
  	vector< vector<double> > F(N, vector<double>(2));
    vector<double> x(N);
    vector<double> y(N);
    vector<vector<double> > d(N, vector<double>(2));//the 
    //vector<double> theta(N);
      double pfad=2*Pe*Pe*(lams/N)+2*Pe*Pe*(lams/N)*(lams/N);
   FILE * den1;
    den1 = fopen("density.txt","w");
//read in the initial conditions
FILE *init;//restart file
 	if ( ( init = fopen( "restart_AD10M.txt", "r" ) ) == NULL){
	  printf ("initial data could not be opened\n");}
 	else {
	  for(i=0;i<N;i++){
	    fscanf(init, "%lf %lf", &x[i], &y[i]);
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
            // double ecut = 1;
            
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
	     // printf("%d %d %lf %lf %lf %lf\n",i,j,d[i][0],d[j][0],d[i][1],d[j][1]);
	      }
	      //printf("%lf %lf %lf %lf\n",d[i][0],d[j][0],d[i][1],d[j][1]);
	      if (r < rcut) {//WCA force
	        
	        double r2i = 1/r2;
	        double r6i = pow(r2i,3);
	        double ff = 48*ecut*r2i*r6i*(r6i-0.5); // L-J force
	        
	        F[i][0] = F[i][0]+ff*xr; // update pair force in components
	        F[i][1] = F[i][1]+ff*yr;
	        
	        F[j][0] =F[j][0]-ff*xr;//Newtons third law.
	        F[j][1] =F[j][1]-ff*yr;
	      }
	      
	      
	      //printf("%d %d %lf\n",i,j,rc[i][j]);
	    }//printf("%lf\n",i,my_vec[i][2]);
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
            
           // if(j%2==0){

            flux+=(Pe*d[j][0])*(h / gamma * F[j][0])+Pe*d[j][1]*(h / gamma * F[j][1]);
            //}
            
          //  flux=flux/2;
            
            fluxp+=((Pe*d[j][0])*((h / gamma) * F[j][0]+ Pe*d[j][0]*h+noise_0+2*Pe*d[j][0]*h*(lams/N))+(Pe*d[j][1])*((h / gamma) * F[j][1]+ Pe*d[j][1]*h+noise_1+2*Pe*d[j][1]*h*(lams/N)));
            //fluxp[j]+=0;
            
            
	                       }
    //printf("%lf %lf %lf %lf %lf\n",i*h,flux/(Pe*Pe)/(i*h)/(den*N),fluxp/N/(i*h),pf,pf+(lams/N)*flux/(i*h)/N);
        printf("%lf %lf %lf %lf %lf\n",m*h,flux/(m*h)/(Pe*Pe)/2/(den*N),fluxp/(N/2)/(m*h),pfad,pfad+(lams/N)*flux/(m*h)/N);
         fprintf(den1,"%lf %lf %lf\n",m*h,flux/(Pe*Pe)/(m*h),lams*flux/(Pe*Pe)/(m*h));
                           
                        }

//print the running average over walkers.

//printf("%d %lf\n",k,);

//}	
                     // fluxp+=flux/(i*h);//this is the running average.
                           
                           //print the instantaneous and running average of friction
                          // printf("%lf %lf %lf\n",i*h,flux/(Pe*Pe),fluxp/(Pe*Pe));





}

/*This code will run the dynamics for active dumbell particles (ADPs) and run cloning calculations from the paper "Entropy production fluctuations encode collective behavior in active matter" by Trevor GrandPre, Katie Klymko, Kranthi K. Mandadapu, and David T. Limmer (2020).

This code was written by Katie Klymko and Trevor GrandPre. It was adapted from a cloning code written by Chloe Ya Gao found here:https://github.com/chloegao12/Nonlinear-transport-coefficients-from-LDF

The inputs are lambda, self-propulsion v, and time of system in that order. In order to compile, you will need GNU scientific library. */

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

// The walker structure is not necessary in this code----we define these quantities further down
struct Walker_info{
	vector<double> x;
	vector<double> theta;
	int parent;
	double Qavg;
};

double fRand(double fMin, double fMax){
	double f = (double)rand()/RAND_MAX;
	return fMin + f*(fMax-fMin);
}


double propagation(vector<double>& x, int N, double L, double tint, gsl_rng* rng,vector<double>& y)
{
    
    /*STEP 1 define the variables that will be used in the simulation*/ 
  double h=0.0001;
  int Nstep = tint/h; //number of steps in the interval
  int i,j;
	double beta = 1.0;
	double gamma = 1.0;
	double T = 1.0;
	double a,b;
	double flux = 0.0;
	double vnewstart,vnewend;
	double Pe=10; // this is the magnitude of the active force.
	vector<vector<double> > F(N, vector<double>(2));//forces in the x and y
	vector<vector<double> > d(N, vector<double>(2));//the direction in the x and y
	double D = 1.0;//set the diffusion coefficient which is set by kT/(m*gamma) 
	double kappa =100;// the spring constant
	double lr = 1.5; //the rest length of the spring. 
	double rcut=pow(2,1./6.);//WCA cutoff
	double ecut =1; //this is epsilon in the WCA formula 

	for(int f=0;f<Nstep;f++){
	  
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
	  
	  
	  /*UPDATE PARTICLE POSTIONS*/
	  for(int j=0;j<N;j++){
	    
	    double fd_term = sqrt(2 *D* h);//noise of the brownian motion
	    double noise_0 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);
	    double noise_1 = fd_term * gsl_ran_gaussian_ziggurat(rng,1);
	    
	    x[j] +=(h / gamma) * F[j][0]+ Pe*d[j][0]*h+noise_0;
	    y[j] +=(h / gamma) * F[j][1]+ Pe*d[j][1]*h+noise_1;
	    
	    //periodic boundaries for 0 to L
	    if (x[j] > L) {x[j] -= (floor(x[j]/L))*L;}
	    if (x[j] < 0) {x[j] += (1+(floor(-x[j]/L)))*L;}
	    if (y[j] > L) {y[j] -= (floor(y[j]/L))*L;}
	    if (y[j] < 0) {y[j] += (1+(floor(-y[j]/L)))*L;}
	    
	    flux+=((Pe*d[j][0])*((h / gamma) * F[j][0]+ Pe*d[j][0]*h+noise_0)+(Pe*d[j][1])*((h / gamma) * F[j][1]+ Pe*d[j][1]*h+noise_1));
	 
	  }//N
}//time
	  return flux/(N/2);//returns the scaled observable per molecule
}

int main(int argc, char *argv[])
{       
        int nprocs,rank;
        gsl_rng* rng;
        gsl_rng_env_setup();
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng,time(NULL)+rank);

	srand(time(0));
	
	//initializing the MPI environment
	MPI_Init(NULL,NULL);
	
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Status status;	

	//parameters in the ADPs model
	int N = 20; //number of particles in the model
	double L = 14.1421356237;//this is half the box of 2L
	vector<double> xinit(N), yinit(N);//allocation
	int i,j,k,s,index;//initialization
	
	
	FILE *init;//restart file
	if ( ( init = fopen( "restart_AD10M.txt", "r" ) ) == NULL){
	  printf ("initial data could not be opened\n");}
	else {
	  for(i=0;i<N;i++){
	    fscanf(init, "%lf %lf", &xinit[i], &yinit[i]);
	  }
	  
	}
	rewind(init);
	fclose(init);


 	
	//parameters for DMC
    
	int Nw =120000;
	int n = Nw/nprocs; //gives the number of walkers per processor
	double lambda = atof(argv[1]);//the strength of the bias is given by the input
	double tint = 0.005, tobs = 10.;//tried 700 try twice as long with a million walkers
	int Ntint = 2000;
//	int Ntint = tobs/tint;
	double weightsum = 0.0;
	int walkersum = 0;
	vector<double> globalq(Nw),localq(n);
	double Q,Q2,localQ,localQ2,sigma;
	double phi = 0.0;

	vector< vector<double> > x(n, vector<double>(N)); //memory allocation
	vector< vector<double> > y(n, vector<double>(N));

	vector< vector<double> > newx(n, vector<double>(N));
	vector< vector<double> > newy(n, vector<double>(N));

    
	vector<double> weight(Nw);
	vector<int> number(Nw);
	vector<int> parent(Nw);
	vector<double> Qavg(Nw);
	vector<int> oldparent(Nw);
	vector<double> oldQavg(Nw);

	vector<double> mean,var,ldf;
	vector<int> multiplicity,m;
	vector<int> table;
	table.reserve(Nw);

	//Initialization of all the walkers


	for(i=0;i<n;i++){
		x[i] = xinit; //this is copying the full vector; each
		y[i] = yinit;
	}

	for(i=0;i<Nw;i++){
		parent[i] = i;
		Qavg[i] = 0.0;
	}
	
	for(i=0;i<Ntint;i++){ //this is the time loop we are interested in
		walkersum = 0;
		weightsum = 0.0;
			printf("time = %d\n",i);


		for(j=0;j<n;j++){
			localq[j] = propagation(x[j],N,L,tint,rng,y[j]);//this returns the flux of a trajectory
		}	

		// Gather all the currents generated in tint to rank 0
		MPI_Gather(&localq[0],n,MPI_DOUBLE,&globalq[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// Reweight the walkers on rank 0 and write the lookup table 
		if(rank==0){
		for(j=0;j<Nw;j++){
			Qavg[j] += globalq[j];
			weight[j] = exp(lambda*globalq[j]);
			weightsum += weight[j];
	}
		for(j=0;j<Nw;j++){
			number[j] = floor(Nw*weight[j]/weightsum + fRand(0,1));
			walkersum += number[j];
		}
		if(walkersum < Nw){
			while(walkersum < Nw){
				number[fRand(0,Nw)] += 1;
				walkersum += 1;
			}	
		}
		if(walkersum > Nw){
			while(walkersum > Nw){
				s = floor(fRand(0,Nw));
				if(number[s]>0){
					number[s] -= 1;
					walkersum -= 1;
				}
			}
		}		 	
		for(j=0;j<Nw;j++){
			for(s=0;s<number[j];s++){
				table.push_back(j);
			}
		}
		}

					
		MPI_Bcast(&table[0],Nw,MPI_INT,0,MPI_COMM_WORLD);
		
		
		for(j=0;j<Nw;j++){
			if(table[j]!=j){
			//if the replaced walker is not on the same core
			if((table[j]/n)!=(j/n)){
			if(rank==table[j]/n){
				MPI_Send(&x[table[j]%n][0],N,MPI_DOUBLE,j/n,j,MPI_COMM_WORLD);
				MPI_Send(&y[table[j]%n][0],N,MPI_DOUBLE,j/n,j,MPI_COMM_WORLD);
			
			}
			if(rank==j/n){
				MPI_Recv(&newx[j%n][0],N,MPI_DOUBLE,table[j]/n,j,MPI_COMM_WORLD,&status);
				MPI_Recv(&newy[j%n][0],N,MPI_DOUBLE,table[j]/n,j,MPI_COMM_WORLD,&status);
				
			}
			}
			else{
				if(rank==j/n){
					newx[j%n] = x[table[j]%n];
					newy[j%n] = y[table[j]%n];
				}
			}
			}
		}

		// Search in the lookup table for replacing walkers	
		for(j=0;j<n;j++){
			index = rank*n+j;
			if(table[index]!=index){
				x[index%n] = newx[index%n];
				y[index%n] = newy[index%n];
			}
		}

		// replace the Q's at rank 0 in the inverse order
		if(rank==0){
			Q = 0.0;
			Q2 = 0.0;
			for(j=0;j<Nw;j++){
				oldQavg[j] = Qavg[j];
				oldparent[j] = parent[j];
			}
			for(j=0;j<Nw;j++){
				if(table[j]!=j){
					Qavg[j] = oldQavg[table[j]];
					parent[j] = oldparent[table[j]];
				}
				Q += Qavg[j];
				Q2 += Qavg[j]*Qavg[j]; 
			}
		}

		table.erase(table.begin(),table.end());
		
		//evaluating average observable and large deviation function
		
		if(rank==0){
		Q = Q/Nw;
		Q2 = Q2/Nw;
		sigma = (Q2-Q*Q)/(i+1)/(tint);
		Q = Q/(i+1)/(tint);//this is Q/tobs making it intensive
		phi += log(weightsum/Nw);
		mean.push_back(Q);
		var.push_back(sigma);
		ldf.push_back(phi/(i+1)/(tint));//1/tobs*vol.
		//compute the multiplicity
		for(j=0;j<Nw;j++){
			m.push_back(parent[j]);
		}
		sort(m.begin(),m.end());
		multiplicity.push_back(distance(m.begin(),unique(m.begin(),m.end())));
		m.erase(m.begin(),m.end());
/*
		if(i%(Ntint/1)==0){ //not sure what this is
                        printf("%lf completed\n",i*100.0/Ntint);
                }
*/		}
	}
	
	if(rank==0){	
	//write out the result
	FILE *result,*current;
	char buffer[800],buffer1[800];
	sprintf(buffer,"basicOutput_%.5lf",lambda);
	sprintf(buffer1,"Qcurrent_%.5lf.txt",lambda);
	result = fopen(buffer,"w");
	fprintf(result,"t	ldf	mean	var	multiplicity\n");
	for(i=0;i<Ntint;i++){
		fprintf(result,"%lf	%lf	%lf	%lf	%lf\n",(i+1)*tint,ldf[i],mean[i],var[i],multiplicity[i]*1.0);
	}
	fclose(result);

        current= fopen(buffer1,"w");
        for(i=0;i<Nw;i++){
                fprintf(current,"%lf\n",Qavg[i]/tobs);
        }
        fclose(current);
	}

	MPI_Finalize();
	return 0;
}

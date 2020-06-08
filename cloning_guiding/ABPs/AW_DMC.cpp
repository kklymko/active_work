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
//calling the force in the integration step
vector< vector<double> > getforce(const vector<double> x, const vector<double> y, int N, double L){
    int i,j;
    vector< vector<double> > f(N, vector<double>(2));
    double rcut=pow(2,1./6.);
    double ecut = 1; //this is epsilon in the WCA formula

		for (i=0; i<N; i++) {
			f[i][0] = 0;
			f[i][1] = 0;
		}

		for (i=0; i<N-1; i++) {
			for (j=i+1; j<N; j++) {
				double xr = x[i]-x[j];
				if(xr>2*L*0.5) {
					 xr-=2*L;}
				else if(xr<(-2*L*0.5)) {
					 xr=2*L+xr;}
				double yr = y[i]-y[j];
				if(yr>2*L*0.5) {
					 yr-=2*L;}
				else if(yr<(-2*L*0.5)) {
					 yr=2*L+yr;}

				double r2 = xr*xr+yr*yr,r = sqrt(r2);

				if (r < rcut) {

						double r2i = 1/r2;
						double r6i = pow(r2i,3);
						double ff = 48*r2i*r6i*(r6i-0.5); // L-J force

						f[i][0] = f[i][0]+ff*xr; // update pair force in components
						f[i][1] = f[i][1]+ff*yr;

						f[j][0] = f[j][0]-ff*xr;
						f[j][1] = f[j][1]-ff*yr;

				}
		}
	}

	return f;
}
//(NEED-we might be able to get rid of the propagationstep and have lammps dump the observable)
// this function takes the initial value x&v and modify it
// and returns the cumulative current Q within tint 
vector<double> propagation(vector<double>& x, int N, double L, double tint, gsl_rng* rng, vector<double>& theta,vector<double>& y,double lams,double Pe)
{
    
    
    //run lammps and have it output the trajectories or the observable. It might make more sense for it to just output the observable since we know that it works. in propagation have lammps dump trajectories and then we compute the local density as flux. THis is also where theentropy production would go. 
  	double h = 0.00001; //think this is the timestep
  	int Nstep = tint/h; //number of steps in the interval
  	int i,j;
	double beta = 1.0;
	double gamma = 1.0;
	double T = 1.0;
	double a,b;
	double flux = 0.0;
	double vnewstart,vnewend;
    	int g=1;
  	vector< vector<double> > f(N, vector<double>(2));
    // vector<double> flamx(N),flamy(N);
	double sum_flambda=0.0;
	double sum_theta=0.0;
double fluxp=0.0;
    
    


	for(i=0;i<Nstep;i++){
			 f = getforce(x,y,N,L);
		//update positions
		for(j=0;j<N;j++){
       
            
            double fd_term = sqrt(2 * h *T/ (gamma ));
			double fd_term_theta = sqrt(3. * 2. * h *T/ (gamma ));
			double noise_theta = fd_term_theta * gsl_ran_gaussian(rng,1);
      		double noise_0 = fd_term * gsl_ran_gaussian(rng,1);
      		double noise_1 = fd_term * gsl_ran_gaussian(rng,1);

            //flamx[j]=2*lams*Pe*cos(theta[j]);
            //flamy[j]=2*lams*Pe*sin(theta[j]);
            
      	    x[j] +=h / gamma * f[j][0]+ noise_0+Pe*cos(theta[j])*h+2*lams*Pe*cos(theta[j])*h/(N);
            y[j] +=h / gamma * f[j][1]+ noise_1+Pe*sin(theta[j])*h+2*lams*Pe*sin(theta[j])*h/(N);
      	    theta[j] += noise_theta;
            
        
			//periodic boundaries o f 2L%
            if(x[j]>L){ x[j]=-L+(x[j]-L); }
            if(x[j]<-L){ x[j]=L+(L+x[j]); }
            if(y[j]>L){ y[j]=-L+(y[j]-L); }
            if(y[j]<-L){ y[j]=L+(L+y[j]); }

		
//The observable is the active work. We take the tagged particle to be the zeroth one. and look at only the x direction. 

            fluxp+=(Pe*cos(theta[j]))*(h / gamma * f[j][0]+ noise_0+Pe*cos(theta[j])*h+2*lams*Pe*cos(theta[j])*h/(N))+Pe*sin(theta[j])*(h / gamma * f[j][1]+ noise_1+Pe*sin(theta[j])*h+2*lams*Pe*sin(theta[j])*h/(N));
            //fluxp+=0.0;
            sum_flambda+=Pe*Pe*lams*h/(N)+Pe*Pe*lams*lams*h/(N*N)+h*f[j][0]*Pe*cos(theta[j])*lams/(N)+h*f[j][1]*Pe*sin(theta[j])*lams/(N);
	                       }   
                        }

					//calculate polarization below--make sure it should be double sum 

					for(i=0;i<N;i++){

						for(j=0;j<N;j++){

							sum_theta+=cos(theta[i]-theta[j]);

						}

					}

				
	
            vector<double> weights(3);
            weights[0]= fluxp/N;
            weights[1]=sum_flambda;
            weights[2]=sqrt(sum_theta)/N;
    	//			printf("%lf\n ",weights[2]);
    		//		printf("HERE");

		        return weights;//returns the scaled observable
}

int main(int argc, char *argv[])
{       int nprocs,rank;
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

	//parameters in the ABPs model
	int N = 10.0; //number of particles in the model
	double L =5;//this is half the box of 2L
	vector<double> xinit(N), yinit(N), thetainit(N);//allocation
	int i,j,k,s,index;//initialization
	
	FILE *init;//restart file
 	if ( ( init = fopen( "equi", "r" ) ) == NULL){
   		printf ("initial data could not be opened\n");}
 	else {
      for(i=0;i<N;i++){
				fscanf(init, "%lf	%lf %lf", &xinit[i], &yinit[i], &thetainit[i]);
     	}
		
   	 }
   	rewind(init);
   	fclose(init);
 	
	//parameters for DMC

	
    
	double lambda = atof(argv[1]);//the strength of the bias is given by the input
    double peclet=atof(argv[2]);
     int Nw=atof(argv[3]);
    // double den=atof(argv[4]);
    // double zeta=atof(argv[5]); 
 	int n = Nw/nprocs; 
    
	double tint = 0.0005, tobs =5.0;//tried 700 try twice as long with a million walkers
	int Ntint =10000;
//	int Ntint = tobs/tint;
	double weightsum = 0.0;
	int walkersum = 0;
	vector<double> globalq(Nw),localq(n);
    vector<double> globalq2(Nw),localq2(n), localq3(n), globalq3(Nw);
    
	double Q,Q2,localQ,localQ2,sigma,P;
	double phi = 0.0;

	vector< vector<double> > x(n, vector<double>(N)); //memory allocation
	vector< vector<double> > y(n, vector<double>(N));
	vector< vector<double> > theta(n, vector<double>(N));
	vector< vector<double> > newx(n, vector<double>(N));
	vector< vector<double> > newy(n, vector<double>(N));
	vector< vector<double> > newtheta(n, vector<double>(N));
	vector<double> weight(Nw);
	vector<int> number(Nw);
	vector<int> parent(Nw);
	vector<double> Qavg(Nw);
	vector<double> Pavg(Nw);
	vector<int> oldparent(Nw);
	vector<double> oldQavg(Nw);
	vector<double> oldPavg(Nw);
    vector<double> values(2);

	vector<double> mean,var,ldf,meanP;
	vector<int> multiplicity,m;
	vector<int> table;
	table.reserve(Nw);

	//Initialization of all the walkers


	for(i=0;i<n;i++){
		x[i] = xinit; //this is copying the full vector; each
		y[i] = yinit;
		theta[i] = thetainit;
	}


	for(i=0;i<Nw;i++){
		parent[i] = i;
		Qavg[i] = 0.0;
		Pavg[i] = 0.0;
	}
	
	for(i=0;i<Ntint;i++){ //this is the time loop we are interested in
		walkersum = 0;
		weightsum = 0.0;
			printf("time = %d\n",i);


		for(j=0;j<n;j++){
			values = propagation(x[j],N,L,tint,rng,theta[j],y[j],lambda,peclet);//this returns the flux of a trajectory
            		localq[j] = values[0];

			localq2[j] = values[1];
			localq3[j] = values[2];
//			printf("%d %lf %lf %lf\n",j,localq[j],localq2[j],localq3[j]);
		}	


		// Gather all the currents generated in tint to rank 0
		MPI_Gather(&localq[0],n,MPI_DOUBLE,&globalq[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&localq2[0],n,MPI_DOUBLE,&globalq2[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Gather(&localq3[0],n,MPI_DOUBLE,&globalq3[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// Reweight the walkers on rank 0 and write the lookup table 
		if(rank==0){
		for(j=0;j<Nw;j++){
			Qavg[j] += globalq[j];
			weight[j] = exp(+globalq2[j]);
			weightsum += weight[j];      
			Pavg[j] += globalq3[j];
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
				MPI_Send(&theta[table[j]%n][0],N,MPI_DOUBLE,j/n,j,MPI_COMM_WORLD);
			}
			if(rank==j/n){
				MPI_Recv(&newx[j%n][0],N,MPI_DOUBLE,table[j]/n,j,MPI_COMM_WORLD,&status);
				MPI_Recv(&newy[j%n][0],N,MPI_DOUBLE,table[j]/n,j,MPI_COMM_WORLD,&status);
				MPI_Recv(&newtheta[j%n][0],N,MPI_DOUBLE,table[j]/n,j,MPI_COMM_WORLD,&status);
				
			}
			}
			else{
				if(rank==j/n){
					newx[j%n] = x[table[j]%n];
					newy[j%n] = y[table[j]%n];
					newtheta[j%n] = theta[table[j]%n];
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
				theta[index%n] = newtheta[index%n];
			}
		}

		// replace the Q's at rank 0 in the inverse order
		if(rank==0){
			P = 0.0;
			Q = 0.0;
			Q2 = 0.0;
			for(j=0;j<Nw;j++){
				oldPavg[j] = Pavg[j];
				oldQavg[j] = Qavg[j];
				oldparent[j] = parent[j];
			}
			for(j=0;j<Nw;j++){
				if(table[j]!=j){
					Pavg[j] = oldPavg[table[j]];
					Qavg[j] = oldQavg[table[j]];
					parent[j] = oldparent[table[j]];
				}
				P +=Pavg[j];
				Q += Qavg[j];
				Q2 += Qavg[j]*Qavg[j]; 
			}
		}

		table.erase(table.begin(),table.end());
		
		//evaluating average observable and large deviation function
		
		if(rank==0){
		P = P/Nw;
		Q = Q/Nw;
		Q2 = Q2/Nw;
		sigma = (Q2-Q*Q)/(i+1)/(tint);
		Q = Q/(i+1)/(tint);//this is Q/tobs making it intensive
		phi += log(weightsum/Nw);
		meanP.push_back(P);
		mean.push_back(Q);
		var.push_back(sigma);
		ldf.push_back(phi/(i+1)/(tint));//1/tobs*vol.
		//compute the multiplicit
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
	sprintf(buffer1,"Qcurrent_%.5lf",lambda);
	result = fopen(buffer,"w");
	fprintf(result,"t	ldf	mean	var	multiplicity\n");
	for(i=0;i<Ntint;i++){
		fprintf(result,"%lf	%lf %lf	%lf	%lf	%lf \n",(i+1)*tint,ldf[i],mean[i],var[i],multiplicity[i]*1.0,meanP[i]);
	}
	fclose(result);

        current= fopen(buffer1,"w");
        for(i=0;i<Nw;i++){
                fprintf(current,"%lf\n",Qavg[i]/tobs);
        }
        fclose(current);
	}
/*
	for(i=n;i<Nw;i++){
                if(rank==i/n){
                        MPI_Send(&x[i][0],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
                        MPI_Send(&v[i][0],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD);
                }
                if(rank==0){
                        MPI_Recv(&x[i][0],N,MPI_DOUBLE,i/n,i,MPI_COMM_WORLD,&status);
                        MPI_Recv(&v[i][0],N,MPI_DOUBLE,i/n,i,MPI_COMM_WORLD,&status);
                }
        }
*/
/*
        FILE *restart;
	for(i=0;i<nprocs;i++){
		if(rank==i){
        	restart = fopen("~/ActiveMatter/Trevor/cloning/restart_file","a");
        	for(j=0;j<n;j++){
                	for(k=0;k<N;k++) fprintf(restart,"%lf\n",x[j][k]);
		}
		fclose(restart);
		}
		MPI_Barrier(MPI_COMM_WORLD);	
        }
  */      
	MPI_Finalize();
	return 0;
}

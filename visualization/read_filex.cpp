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
int main(int argc, char *argv[])
{

//read in the traj file. 
//initialize variables
int Ntint=500;
int Nw=30000;
int N=10;
int Number=29652510;//size of traj raw 


int dummy,dummy1,dummy2;
int i,j,k;

vector<vector<vector<double> > > traj1 (N,vector<vector<double> >(Ntint,vector <double>(Nw)));
//vector<vector<vector<double> > > traj2 (N,vector<vector<double> >(Ntint,vector <double>(Nw)));
//vector<vector<vector<double> > > traj3 (N,vector<vector<double> >(Ntint,vector <double>(Nw)));

vector<vector<vector<double> > > x (N,vector<vector<double> >(Nw,vector <double>(Ntint)));
//vector<vector<vector<double> > > y (N,vector<vector<double> >(Nw,vector <double>(Ntint)));
//vector<vector<vector<double> > > theta (N,vector<vector<double> >(Nw,vector <double>(Ntint)));

vector< vector<int> > history(Nw, vector<int>(Ntint));
    
//what is this number supposed to be?       


FILE *data;
if((data=fopen("traj_raw_test","r"))==NULL){
    printf("data file cannot be opened\n");}
else{
    for(i=0;i<Number;i++){
        fscanf(data,"%d %d",&dummy1,&dummy2);
        for(j=0;j<N;j++){//ive changed the indices of the vectors to have N first.
 //           fscanf(data,"%lf %lf %lf",&traj1[j][dummy1][dummy2],&traj2[j][dummy1][dummy2],&traj3[j][dummy1][dummy2]);
            fscanf(data,"%lf %*lf %*lf",&traj1[j][dummy1][dummy2]);
        }
    }
}

//read in the history file 
fclose(data);
printf("finished reading raw data\n");

FILE *table;
if((table=fopen("historytable","r"))==NULL){
    printf("data file cannot be opened\n");}
else{
    for(i=0;i<Nw;i++){
				int dummm;
        for(j=0;j<Ntint;j++)    fscanf(table,"%*d %d",&history[i][j]);
     //   for(j=0;j<Ntint;j++)    fscanf(table,"%d %d",&dummm,&history[i][j]);
		//		printf("%d %d %d %d\n",i,j,dummm,history[i][j]);
    }
}

fclose(table);
//construct
for(i=0;i<Nw;i++){
    for(j=0;j<Ntint;j++){
        dummy=history[i][j];
        for(k=0;k<N;k++){
            x[k][i][j]=traj1[k][j][dummy];//particle number, Ntint, Nw(right). left is particle number NW, Ntint
  //          y[k][i][j]=traj2[k][j][dummy];
 //           theta[k][i][j]=traj3[k][j][dummy];
            //if(v[i][j][k]>10)       print("something is wrong with input data i=%d,j=%d\n",i,j);
        }
    }
}
printf("finished constructing trajectories\n");

//how many clone videos to make
int sizee = 100;

for(k=0;k<sizee;k++){
				FILE *output;
				char buffer[80];

				sprintf(buffer,"videoX%d.xyz",k);
				output = fopen(buffer,"w");

				for(j=0;j<Ntint;j++){
					fprintf(output,"%d\n\n",N);
					for(i=0;i<N;i++){
//						fprintf(output,"O %lf %lf %lf 0.0\n",x[i][k][j],y[i][k][j],theta[i][k][j]);
//						if(i==0){fprintf(output,"O %lf\n",x[i][k][j]);}
						fprintf(output,"N %lf\n",x[i][k][j]);
					}
				}

				fclose(output);
}
} //end of main

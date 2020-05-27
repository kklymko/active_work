#include <stdio.h>
#include <time.h>
#include <string.h>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <errno.h>
using namespace std;

#define N 2000
#define n_steps 8000
#define BOXX 100
#define BOXY 100

int main(int argc, char *argv[]){

	FILE *data,*data2;
	int n,m;
	int i,j;
	double* X = new double[N+1];
	double* Y = new double[N+1];
	double* theta = new double[N+1];
	int* type = new int[N+1];
//	double X[N+1],Y[N+1],theta[N+1];
//	int type[N+1];
	double sigma=60.0;
	int nhis=4000;
	int nhistheta=400;
	double delg=BOXX/(double)(2*nhis);
	double delgtheta=2*M_PI/(double)(nhistheta);
	double rho=0.5;
//	double probe[nhis][nhistheta];
	char buffer[800];
	char lines[800];
	int t;
	int ngr=0;

  double** probe = new double*[nhis];

  for(i = 0; i < nhis; i++){
      probe[i] = new double[nhistheta];
  }


//	data=fopen("THETA60.lammpstrj","r");
	data = fopen("test4.xyz","r");


	int junk1,junk2,junk;
	double work,work2;
	double xx,yy;

	for(i=0;i<nhis;i++){
		for(j=0;j<nhistheta;j++){
			probe[i][j]=0;
		}
	}

int counter=0;

for(t=0;t<n_steps;t++){

	ngr++;

	printf("%d\n",t);

  for(i=0;i<28;i++){

        fscanf(data,"%s",&lines);
  }

	for(j=1;j<=N;j++){

        fscanf(data, "%d", &junk1);
        fscanf(data, "%d", &junk2);
        fscanf(data, "%lf", &xx);
        fscanf(data, "%lf", &yy);
        fscanf(data, "%d", &junk);
        fscanf(data, "%lf", &work);
			//	if(xx>100 || yy>100){
			//		printf("%d %d %lf %lf %d %lf\n",junk1,junk2,xx,yy,junk,work);
		//		}
				type[junk1]=junk2;
				X[junk1]=xx;
				Y[junk1]=yy;
				theta[junk1]=work;
	//			printf("%d %lf %lf %lf\n",type[j],X[j],Y[j],theta[j]);

	}

	for(j=1;j<=N-1;j++){
		for(i=j+1;i<=N;i++){
					

				double distx=X[j]-X[i];
						if(distx>BOXX*.5){
								distx-=BOXX;
						}
						if(distx<-BOXX*.5){
								distx+=BOXX;
						}
				double disty=Y[j]-Y[i];
						if(disty>BOXY*.5){
								disty-=BOXY;
						}
						if(disty<-BOXY*.5){
								disty+=BOXY;
						}


					double dr=sqrt(distx*distx+disty*disty);

					if(dr<BOXX/2.){
						counter++;
						double thetax=cos(theta[i]);
						double thetay=sin(theta[i]);

						double argument=max((thetax*distx+thetay*disty)/(dr),-1.0);
						double argument2=min(argument,1.0);

						double angle=acos(argument2);
						if(disty>0){
							angle+=M_PI;
						}
						if(disty<0){
							angle=M_PI-angle;
						}
						//printf("%lf\n",angle);

			//			printf("%d %d %lf %lf %lf %lf\n",i,j,dr,angle,delg,delgtheta);
			//		printf("%lf %lf %lf %lf\n",X[i],X[j],Y[i],Y[j]);

						//think normalization should be the same!
						int ig=(int)(dr/delg);
						int ig2=(int)(fmod((fabs(angle)),2*M_PI)/delgtheta);
						if(ig2>=nhistheta){
							printf("PROBLEM 1 %d %d %lf %lf\n",nhistheta,ig2,fabs(angle),fmod(fabs(angle),M_PI)/delgtheta);
							exit(0);
						}
						if(ig>=nhis){
							printf("PROBLEM 2 %d %d\n",nhis,ig);
							exit(0);
						}
            probe[ig][ig2]=probe[ig][ig2]+1;

				}
			}
	}
}

	
fclose(data);

	for(i=0;i<nhis;i++){
		for(j=0;j<nhistheta;j++){

			printf("%d %d %lf\n",i,j,probe[i][j]);

		}
	}
	char buffer10[80];
	char buffer2[80];

	sprintf(buffer10,"corr_dat_dil_%d_%d",nhis,nhistheta);
	data=fopen(buffer10,"w");

	for(j=0;j<nhistheta;j++){
	//	sprintf(buffer10,"corr_dat_dil%d",j);
	//	data2=fopen(buffer10,"w");
		for(i=0;i<nhis;i++){
				double vb=M_PI*(pow((i*delg+delg),2)-pow((i*delg),2));
								printf("%lf %lf %lf\n",probe[i][j],(double)(nhistheta)*ngr*vb*rho*(((double)(N)-1)/2.),probe[i][j]/((double)(nhistheta)*ngr*vb*rho*(((double)(N)-1)/2.)));
                //probe[i][j]=probe[i][j]/((double)(nhistheta)*n_steps*vb*rho*(((double)(N)-1)/2.));
    //            fprintf(data2,"%lf %lf %lf %lf\n",(double)(i)*delg,(double)(j)*delgtheta,probe[i][j],probe[i][j]/(ngr*vb*rho*(((double)(N)-1)/2.)));
                fprintf(data,"%lf %lf %lf %lf\n",(double)(i)*delg,(double)(j)*delgtheta-M_PI,probe[i][j],probe[i][j]/(ngr*vb*rho*(((double)(N)-1)/2.)));
//                fprintf(data2,"%lf\n",probe[i][j]);
		}
//		fclose(data2);
	}

	fclose(data);


	delete [] X;
	delete [] Y;
	delete [] theta;
	delete [] type;

      for(i = 0; i < nhis; i++) {
            delete [] probe[i];
        }
        delete [] probe;

	return 0;

}

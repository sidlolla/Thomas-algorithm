#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define aE_T 1
#define aW_T 1
#define aN_T 1
#define aS_T 1
#define aP_T 4	

double TDMA(int Nx,int Ny,double aE[][Ny+1],double aW[][Ny+1],double aN[][Ny+1],double aS[][Ny+1],double aP[][Ny+1],double Phi[][Ny+1])
{	
	int i,j,k=0;
	double a[Nx+1][Ny+1],b[Nx+1][Ny+1],c[Nx+1][Ny+1],d[Nx+1][Ny+1],t[Nx+1][Ny+1],g,e=1,f;
	
		
	for(i=0; i<=Nx; i++)
  		for(j=0; j<=Ny; j++) 
   		  	t[i][j]=Phi[i][j];		//'t' stores old values
		
	while(e>0.001){
		e=0;	
		k++;
		for(j=1;j<Ny;j++){
			
		/* input values of coefficients
		SOLVING EQUATION "aT(i-1)+bT(i)+cT(i-1)=d"   */	
			
			d[1][j] = -aN[1][j]*Phi[1][j+1] -aS[1][j]*Phi[1][j-1] -aW[1][j]*Phi[0][j];
			d[Nx-1][j] = -aN[Nx-1][j]*Phi[Nx-1][j+1] -aS[Nx-1][j]*Phi[Nx-1][j-1] -aE[Nx-1][j]*Phi[Nx][j];
			
			for(i=2;i<Nx-1;i++)	
				d[i][j] = -aN[i][j]*Phi[i][j+1] -aS[i][j]*Phi[i][j-1];
				
			for(i=2;i<Nx;i++)
				a[i][j] = aW[i][j];
			for(i=1;i<Nx;i++)
				b[i][j] = -aP[i][j];						
			for(i=1;i<Nx-1;i++)
				c[i][j] = aE[i][j];	
		
		//elimination
			for(i=2;i<Nx;i++){
				g = (double)a[i][j]/b[i-1][j];
				b[i][j] = b[i][j] -g*c[i-1][j];
				d[i][j] = d[i][j] -g*d[i-1][j];
			}			
			
		//back substitution		
			Phi[Nx-1][j]=(double)d[Nx-1][j]/b[Nx-1][j];
			
			for(i=Nx-2;i>=1;i--){
				Phi[i][j]=(double)(d[i][j]-c[i][j]*Phi[i+1][j])/(double)b[i][j];
				f = t[i][j]-Phi[i][j];
     			e=e+fabsf(f);		//'e' stores the iteration error
     		}}
		
			for(i=1; i<=Nx-1; i++)
  				for(j=1; j<=Ny-1; j++) 
   		  			t[i][j]=Phi[i][j];	//updating old values
	}   		  			
}


int main()
{
	int Nx=100,Ny=100,j,i;
	double Phi[Nx+1][Ny+1],aE[Nx+1][Ny+1],aW[Nx+1][Ny+1],aN[Nx+1][Ny+1],aS[Nx+1][Ny+1],aP[Nx+1][Ny+1];
		
	
	for(i=1; i<=Nx-1; i++){
		for(j=1; j<=Ny-1; j++){
			aE[i][j] = aE_T;
			aW[i][j] = aW_T;
			aN[i][j] = aN_T;
			aS[i][j] = aS_T;
			aP[i][j] = aP_T;
		} 	
	}
  		
	
	//Boundary conditions	
		for(i=0;i<=Nx;i++){			
			Phi[i][Ny] = 100;		//top face
			Phi[i][0] = 0;		//bottom face
		}
		for(j=0;j<=Ny;j++){
			Phi[0][j] = 0;		//left face
			Phi[Nx][j] = 0;		//right face
		}
		
	//Initial guess
		for(i=1; i<Nx; i++){
   	 		for(j=1; j<Ny; j++){
   	  			Phi[i][j]=0;  		
   	  		}}
	
	/*calling the solver*/
	TDMA(Nx,Ny,aE,aW,aN,aS,aP,Phi);
	
	
	//Printing result	
	for(j=Ny; j>=0; j--){
			printf("i=50	j=%d	%lf\n",j,Phi[50][j]);
		}	
    return 0;
}

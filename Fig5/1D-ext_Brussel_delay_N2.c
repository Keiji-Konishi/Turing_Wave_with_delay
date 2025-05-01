/** 1D-ext_Brussel **/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f_1(double,double,double,double,double,double); 	/** Local Map **/
double f_2(double,double,double,double,double,double); 	/** Local Map **/
double f_3(double,double,double,double,double,double); 	/** Local Map **/

double 	ep, a, b, Dt, Dx, D1, D2, D3; 	/** Global Parameter **/

int main(void)
{
	static double		u_a1[1000][5000], u_a2[1000][5000], u_a3[1000][5000], uu_a1[1000], uu_a2[1000], uu_a3[1000];
	static double		u_b1[1000][5000], u_b2[1000][5000], u_b3[1000][5000], uu_b1[1000], uu_b2[1000], uu_b3[1000];
	double			ep0, tau, L, t_max, t_min, t_start;
	int 			i, it, i_delay, j, t_skip, N;  

	a=1.08; b=3.08;					/** a, b **/
	/** a=0.96; b=2.85;				/** a, b **/
    D1=0.01; D2=0.1; D3=1.0;        /** D_1,2,3 **/
	ep0=0.00;					    /** Coupling strength **/
	tau=0.50;					    /** Time delay **/
	/** ep0=0.05;					/** Coupling strength **/
    /** tau=1.50;					/** Time delay **/
    L=16.0;						    /** Length of media **/
	N=256;						    /** Number of grid **/
	t_skip=2000;					/** Time skip for print**/
	Dt=5e-4;					    /** Time Step **/
	Dx=L/((double) N);				/** Grid width **/
	i_delay=((int) (tau/Dt));		/** Time step delay **/
	t_max=1000;					    /** Time max**/
	t_min=800;      				/** Time max**/
    t_start=850;			        /** Time to start **/

	/** Dt, Dx check**/
	if ( Dt/(Dx*Dx) >= 1.0/6.0 )
	{ printf("Dt and Dx should be changed %f \n", Dt/(Dx*Dx)); return 0; } 

	/** Initial condition **/
	for (i=0; i<= N+1; i++)
	{ 
		for (j=1; j <= i_delay+1; j++)
		{
		u_a1[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+a;
	      	u_a2[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+(b/a);
	      	u_a3[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+a;
		u_b1[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+a;
	      	u_b2[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+(b/a);
	      	u_b3[i][j]=0.2*( 2.0*((double) rand())/((double) RAND_MAX + 1.0) - 1.0)+a;
		}
	}

	/** Time development **/
	for (it=0; it <= ((int) t_max/Dt); it++)
       	{
		/* Boundary condition (Periodic) */
		u_a1[0][1]=u_a1[N][1]; u_a1[N+1][1]=u_a1[1][1];
		u_a2[0][1]=u_a2[N][1]; u_a2[N+1][1]=u_a2[1][1];
		u_a3[0][1]=u_a3[N][1]; u_a3[N+1][1]=u_a3[1][1];
		u_b1[0][1]=u_b1[N][1]; u_b1[N+1][1]=u_b1[1][1];
		u_b2[0][1]=u_b2[N][1]; u_b2[N+1][1]=u_b2[1][1];
		u_b3[0][1]=u_b3[N][1]; u_b3[N+1][1]=u_b3[1][1];

		/* Coupling start */
		if ( ((double) it)*Dt >= t_start)
		{ep=ep0;} else {ep=0.0;}
						
		/* Map */
		for (i=1; i<=N; i++)
		{   
			uu_a1[i]=f_1(u_a1[i][1],u_a2[i][1],u_a3[i][1],u_a1[i-1][1],u_a1[i+1][1],u_b1[i][i_delay+1]);
			uu_a2[i]=f_2(u_a1[i][1],u_a2[i][1],u_a3[i][1],u_a2[i-1][1],u_a2[i+1][1],u_b2[i][i_delay+1]);
			uu_a3[i]=f_3(u_a1[i][1],u_a2[i][1],u_a3[i][1],u_a3[i-1][1],u_a3[i+1][1],u_b3[i][i_delay+1]);
			uu_b1[i]=f_1(u_b1[i][1],u_b2[i][1],u_b3[i][1],u_b1[i-1][1],u_b1[i+1][1],u_a1[i][i_delay+1]);
			uu_b2[i]=f_2(u_b1[i][1],u_b2[i][1],u_b3[i][1],u_b2[i-1][1],u_b2[i+1][1],u_a2[i][i_delay+1]);
			uu_b3[i]=f_3(u_b1[i][1],u_b2[i][1],u_b3[i][1],u_b3[i-1][1],u_b3[i+1][1],u_a3[i][i_delay+1]);

			/* Print */	
			if ( (it % t_skip) == 0 && ((double) it)*Dt >= t_min ) 
/*			{printf("%f %f %f %f \n", ((double) i)*Dx, ((double) it)*Dt, v_R1[i][1], v_R2[i][1]);} */
			{printf("%f %f %f %f \n", ((double) i)*Dx, ((double) it)*Dt, u_a1[i][1], u_b1[i][1]);}
		}
		/* Mapping */
		for (i=1; i<=N; i++) 
		{
			for (j=i_delay; j>=1; j=j-1)
			{
				u_a1[i][j+1]=u_a1[i][j]; u_a2[i][j+1]=u_a2[i][j]; u_a3[i][j+1]=u_a3[i][j];
				u_b1[i][j+1]=u_b1[i][j]; u_b2[i][j+1]=u_b2[i][j]; u_b3[i][j+1]=u_b3[i][j];
			}
			u_a1[i][1]=uu_a1[i]; u_a2[i][1]=uu_a2[i]; u_a3[i][1]=uu_a3[i];
			u_b1[i][1]=uu_b1[i]; u_b2[i][1]=uu_b2[i]; u_b3[i][1]=uu_b3[i];
		}
		
		/* Divergence */
		if (fabs(u_a1[1][1]) > 5.0){printf("Divergence! \n");return 0;}

	}
return 0;
}

double f_1(	double u_1, double u_2, double u_3, 
                double u_1l, double u_1r, double u_1c_other)
{
	return 
	( 	u_1 
		+ Dt*(a-(2.0+b)*u_1+u_1*u_1*u_2+u_3)
		+ Dt*ep*(u_1c_other- u_1)
		+ ( Dt*D1/(Dx*Dx))*( u_1r+u_1l-2.0*u_1 ) 
	);  
}

double f_2(	double u_1, double u_2, double u_3, 
                double u_2l, double u_2r, double u_2c_other)
{
	return 
	( 	u_2 
		+ Dt*(b*u_1-u_1*u_1*u_2)
		+ Dt*ep*(u_2c_other- u_2)
		+ ( Dt*D2/(Dx*Dx))*( u_2r+u_2l-2.0*u_2 ) 
	);  
}

double f_3(	double u_1, double u_2, double u_3, 
                double u_3l, double u_3r, double u_3c_other)
{
	return 
	( 	u_3 
		+ Dt*(u_1-u_3)
		+ Dt*ep*(u_3c_other- u_3)
		+ ( Dt*D3/(Dx*Dx))*( u_3r+u_3l-2.0*u_3 ) 
	);  
}


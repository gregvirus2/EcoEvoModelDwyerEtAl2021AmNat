// --------------------------------Begin ODE system of White model --------------------------------------------//
int fast_odes(double t, const double y[], double dydt[],void *Paramstuff)
{
//struct STRUCTURE *Params=(struct STRUCTURE *)Paramstuff;
ParamStruct* ODEParams;
ODEParams = (ParamStruct*) Paramstuff;


double k		= ODEParams->FxdPars[1];
double V 		= 1/k; //C^2
double nu		= ODEParams->FxdPars[2];		
double mu		= ODEParams->FxdPars[3];
//4 is sigma and 5 is ratio, neither is needed here		
double delta		= ODEParams->FxdPars[6];
//7 is h and 8 is correlation interval, neither is needed here
int m			= ((int) ODEParams->FxdPars[9]);  //CK//  the number of exposed classes.  adjustable and fit.  yay!
double S0 = ODEParams->FxdPars[11];

				

int i;

double hetero = y[m+1]*nu*pow( ((double) (y[0]/S0)), ((double) V));	// save time by doing once (heterogeneity term), where y[m+1] is P, y[0] is S, and V = C^2, 
// ------------------------------------------ ODEs -------------------------------------------- //


if(y[0]< 1.0e-6){
	dydt[0]=0;
}else{				
dydt[0]  = -y[0]*hetero;	 //dS/dt
}

dydt[1]  = y[0]*hetero - m*delta*y[1]; //dE1/dt

for(i=2; i <= m; i++){  //dEi/dt for i = 2,m
	dydt[i]=m*delta*(y[i-1] -y[i]);
}

dydt[m+1]  = m*delta*y[m] - mu*y[m+1];  //dP/dt, where y[m] is E_m
dydt[m+2] = m*delta*y[m];  //dR/dt

	   
/*for (i=0;i<15;i++)	{
	printf("dydt(%d)=%e\n",i,dydt[i]);
}
*/
return GSL_SUCCESS;
}

////////////////////Begin Jacobian of model///////////////////////////
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *Paramstuff)
{
	//printf("beginning of Jacobian\n");	getc(stdin);


ParamStruct* ODEParams;
ODEParams = (ParamStruct*) Paramstuff;


double k		= ODEParams->FxdPars[1];
double V 		= 1/k;
double nu		= ODEParams->FxdPars[2];		
double mu		= ODEParams->FxdPars[3];
//4 is sigma and 5 is ratio, neither is needed here		
double delta		= ODEParams->FxdPars[6];
//7 is h and 8 is corr'l'n interval, neither is needed
int m			= ODEParams->FxdPars[9];  //CK//  the number of exposed classes.  adjustable and fit.  yay!	

				

//int i;
int DIM = m + 3;

double S0 = ODEParams->FxdPars[10];


double hetero = nu*pow((y[0]/S0),(V));	// save time by doing once (heterogeneity term)
//double hetero_pow = pow((y[0]/S0),V);  //What the hell is this? It looks like garbage

int i;
int a;
int b;

/*
for (i=1;i<=18;i++)	{
	printf("parm(%d)=%e\n",i,Params->PARS[i]);
}getc(stdin);
*/



//dS/dt
*(dfdy + 0 * DIM +0) = -y[m+1]*hetero;  // (\partial/\partial S) of dS/dt (IF YOU ARE SKEPTICAL WRITE IT DOWN YOURSELF, YOU'LL SEE) y[m+1] is P
for(i=1; i <= m; i++){
	*(dfdy + 0 * DIM +i) = 0;  //(\partial/\partial E_i) of dS/dt  

}
*(dfdy + 0 * DIM + m + 1) = -y[0]*hetero;   //(\partial/\partial P) of dS/dt, y[0] is S

*(dfdy + 0 * DIM + m + 2) = 0;  // (\partial/\partial R) of dS/dt


//dE_1/dt
*(dfdy + 1 * DIM +0) = y[m+1]*hetero; // (\partial/\partial S) of dE_1/dt, y[m+1] is P
*(dfdy + 1 * DIM +1) = -m*delta;  // (\partial/\partial E_1) of dE_1/dt

for(i=2; i <= m; i++){            // (\partial/\partial E_i) of dE_1/dt for i = 2,...,m

	*(dfdy + 1 * DIM +i) = 0;
}
*(dfdy + 1 * DIM + m + 1) = y[0]*hetero;  // (\partial/\partial P) of dE_1/dt, y[0] is S

*(dfdy + 1 * DIM + m + 2) = 0;  // (\partial/\partial R) of dE_1/dt


//dE_i/dt for i = 2  //PROBABLY POINTLESS - IF YOU DELETE THIS, YOU MUST DO a = 2,...,m below
*(dfdy + 2 * DIM +0) = 0;  // (\partial/\partial S) of dE_i/dt
*(dfdy + 2 * DIM +1) = m*delta; // (\partial/\partial E_2) of dE_2/dt   
*(dfdy + 2 * DIM +2) = -m*delta; // (\partial/\partial E_{3}) of dE_2/dt  
for(i=3; i <= m+1; i++){  // (\partial/\partial E_i) of dE_2/dt for i = 3,...,m, and (\partial /\partial P) of dE_2/dt
	*(dfdy + 2 * DIM +i) = 0;
}
*(dfdy + 2 * DIM + m + 2) = 0;   // (\partial/\partial R) of dE_2/dt for  


//dE_i/dt for i = 3,...,m 
for(a=3; a <= m; a++){
	for(b=0; b <= m+2; b++){
		if(b == (a-1)){*(dfdy + a * DIM +b) =  m*delta;}   // (\partial/\partial E_i) of dE_i/dt 
		else if (b == a){*(dfdy + a * DIM +b) = -m*delta;}  // (\partial/\partial E_{i+1}) of dE_i/dt 
 		else {*(dfdy + a * DIM +b) = 0;}   // (\partial/\partial S) and (\partial/\partial P) and (\partial/\partial R) of dE_i/dt 
	}
}


//dP/dt
for(i=0; i <= m-1; i++){    // (\partial/\partial S) and (\partial/\partial E_i) for i = 1,...,m-1,  of dP/dt 

	*(dfdy + m+1 * DIM +i) = 0;
}
*(dfdy + m+1 * DIM +m) = m*delta; // (\partial/\partial E_{m}) of dP/dt
*(dfdy + m+1 * DIM +m+1) = -mu; // (\partial/\partial P) of dP/dt
*(dfdy + m+1 * DIM + m+2) = 0;  // (\partial/\partial R) of dP/dt



//dR/dt
for(i=0; i <= m-1; i++){   // (\partial/\partial S) and (\partial/\partial E_i) for i = 1,...,m-1,  of dR/dt 
	*(dfdy + (m+2) * DIM +i) = 0;
}
*(dfdy + (m+2) * DIM + m) = m*delta;  // (\partial/\partial E_{m}) of dR/dt

*(dfdy + (m+2) * DIM + m+1) = 0;
*(dfdy + (m+2) * DIM + m+2) = 0; 


for(i=0; i <= DIM-1; i++){  //there are no functions in the ODE's that are continuous functions of time
	dfdt[i] = 0;
}

return GSL_SUCCESS;
}

// ------------------------------------------  ODE Solver  ----------------------------------------------- //
double ODE_Solver(double t_ode,double t_end,void *Paramstuff,double *y_ode)
{
int i;
int status_ode;
double h_init=1.0e-5;

ParamStruct* ODEParams;
ODEParams = (ParamStruct*) Paramstuff;

int DIM = ODEParams->FxdPars[9]+3; //The total number of ODEs is the number of E classes, plus 3 because there are also the S, P and R classes


//const gsl_odeiv2_step_type *solver_ode	= gsl_odeiv2_step_rkf45; // Runge-Kutta Fehlberg (4, 5)
const gsl_odeiv_step_type *solver_ode = gsl_odeiv_step_rk4; 

// returns pointer to a newly allocated instance of a stepping function of type 'solver_ode' for a system of DIM dimensions //
//gsl_odeiv2_step *step_ode	= gsl_odeiv2_step_alloc(solver_ode, DIM);
gsl_odeiv_step *step_ode	= gsl_odeiv_step_alloc(solver_ode, DIM);


gsl_odeiv_control *tol_ode	= gsl_odeiv_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);
//gsl_odeiv2_control *tol_ode	= gsl_odeiv2_control_standard_new(1.0e-10, 1.0e-5, 1.0, 0.2);

//gsl_odeiv2_evolve *evol_ode	= gsl_odeiv2_evolve_alloc(DIM);
gsl_odeiv_evolve *evol_ode	= gsl_odeiv_evolve_alloc(DIM);


gsl_odeiv_system sys_ode;
//gsl_odeiv2_system sys_ode;
sys_ode.function  = fast_odes;
sys_ode.jacobian  = jacobian;
sys_ode.dimension = (size_t)(DIM);
sys_ode.params	  = ODEParams;

//double y_err[DIM]; double dydt_in[DIM];	double dydt_out[DIM];

// ----------------------------------- Integrate Over Time ------------------------------------ //



while (t_ode<t_end)	{


	status_ode = gsl_odeiv_evolve_apply(evol_ode, tol_ode, step_ode, &sys_ode, &t_ode, t_end, &h_init, y_ode);
	//status_ode = gsl_odeiv2_step_apply(step_ode, t_ode, h_init, y_ode, y_err, dydt_in, dydt_out, &sys_ode);



	for (i=0;i<DIM;i++)	{
	
		/*
		if (y_ode[i]>0)		{
			// keep y_ode as is
		}
		else				{
			//printf("y(%d) NEGATIVE or not a number\n",i);
			y_ode[i]=0;
		}
		*/
		if(y_ode[i]<0) y_ode[i] = 0.0; 
		//This next bit holds only because initial cadavers is a fraction in RK45SEIR 
		/*
		if (y_ode[i]>ODEParams->FxdPars[11])	{  //FxdPars[11] is S0
			//printf("y(%d) TOO LARGE!!\n",i);  
			y_ode[i]=ODEParams->FxdPars[11];
		}
		*/
		
	}
	
	//t_ode+=h_init;


}
// -------------------------------------- Clear Memory ----------------------------------------- //


//gsl_odeiv2_evolve_free(evol_ode);
//gsl_odeiv2_control_free(tol_ode);
//gsl_odeiv2_step_free(step_ode);

gsl_odeiv_evolve_free(evol_ode);
gsl_odeiv_control_free(tol_ode);
gsl_odeiv_step_free(step_ode);

return (t_end);
}



extern int DataSet;
extern int Figures;

extern float FractI2[20];


#define IntMax 600000
// tau divided by the step size


//int PrintMe = 0; //whether to print output to a file



double DistdRK45(double *Pars,float **Covars,float ***ModelFI,float ***ModelS,double ***RandNums,int *MaxDays,int MaxPlots,int DataInterval, int FractDead){


ParamStruct ODEParams;

 
	//This version assumes that Covars[2][1] is the actual virus pop, NOT the fractinf 
 double initz;
 int flag = 1;
 int Plt;
 //int MaxPlts = 5;

 double ratio;
 int Split;
 double Hch;
 double maxData;
 int plt;
 FILE *fp,*fp2,*fp3;
 double maxt;    // number of days to run for


 double z;
 int ind; //indicates where you are in storage arrays
 double *Intx;
 double *SStor, *PStor, *nuStor;
 double Stmp, Ptmp, nutmp;

 
 double s; 
 double *inf; 
 //double inf[100];
 /* *initx; */
  /* z is pibs, and s is suscept.'s */
 /* Intx is sum over the suscepts. times nu times pibs z, for storage */
 double zinit;  /* initial pibs, to simulate overwintering */
 //double y,oldy;      /* infect.'s - no need to store'em */
 double SumE; //sum of the E classes
 int i;
 double t;
 int tau, tau0;
 float WM[1000], Suscept[1000];
 //double *WM, *Suscept;
 int week;
 double FractI;
 double S0;
 double nubar, mu, k, gamma, delta;
 double nubarAdj; 
 int MaxStages;
 int CorrlnIntrvl, TimeIntrvl;
 int DataDay;
 double S0Adj;
 double TooLow = -1;
 double TooHigh = 1e10;
 int TooHighOut;
 double tempInt;
 int OutputCount;
 int ReDo, ReDoCount = 0;
 int Bang, FirstCrash;
 float Z0;
 double SExact, IExact, PExact, dummy1, dummy2, dummy3;
int Trial = 1;
 

  //if(Figures==1) fp = fopen("FigureOutSEIR9Aug14lokLoDens.dat","a"); //Det version
  //if(Figures==1) fp = fopen("StochOutSEIR9Aug14lokHiDens.dat","a");

  Figures = 0;
  if(Figures==1) fp = fopen("RK452WM.dat","w");




  k = Pars[1];  //model parameter: Squared CV of dist'n of transmission rates
  nubar = Pars[2];   //average transmission rate
  
  mu = Pars[3];

//For the control plots in the Otvos dataset, the epidemic starts so late that all larvae are already at least 4ths
  if((Pars[0]<4)&&(DataSet==16)){
	ratio = 1;
  }else{
  	ratio = Pars[5];
  }
  //printf("Pop:%f ratio:%f\n",Pars[0],ratio); getc(stdin);
	
  delta = Pars[6];
  
  CorrlnIntrvl = Pars[8];  //Correlation Interval

  double StagesTemp = Pars[9] + 0.5;
  int m = ( (int) StagesTemp);   //Number of Stages, rounded off
  int DIM = m; //DIM may be used elsewhere?  Otherwise this is a waste of time




//for(Trial=1;Trial<=MaxTrials;Trial++){


ReDo = 1; 

	TooHighOut = 0;


Plt=Pars[0]; 


	for(i=0;i<20;i++){
		WM[i] = 0.0;
		Suscept[i] = 0.0;
	}



    S0 = Covars[1][Plt]; //Initial population size
    FractI = Covars[2][Plt]; //Initial fraction infected
    Split = Covars[3][Plt];  //Whether hatch time is spread out, yes or no

   
if(FractI<0){
 FractI = pow(10,Pars[11+Plt]);
}
	 //maxData = 9;  //Maximum data in any plot

	 maxt = MaxDays[Plt];



   //z = S0*FractI*Hch*ratio;  //Initial density of virus OBSOLETE
   z = S0*FractI*ratio;  //Initial density of virus - this is the new version, where if there's a split, it is completely handled below.  Search on FractI2 to find it
   

   Z0 = z;


   if(Split)
   	s = S0*(1-FractI)*(1-FractI2[Plt]);
   else
   	s = S0*(1-FractI);  //Initial density of hosts


   S0Adj = S0;  //Initial density of hosts, saved for later use
   initz = z;  //Saving initial virus density, not sure why?

  
   Z0 = z;  //Emu
  

  week = 1;  //Initializing week number, need it for data storage

  flag = 1;  //Says whether or not additional hatch has occurred

  //y= 0.0;  //Initializing number dead per unit time
  //Dead = 0.0;  //Initializing cumulative number dead per week
  int OldTimeIntrvl = -1;  //Useful for checking the time interval action
 

  OutputCount = 1; 

double Oldnubar = -1;
Bang = 0; FirstCrash = 1;
double StageAdjust;


ind = 0;

//for(t=h;t<=maxt;t+=h){  //time loop

t = 0; TooHighOut = 0;


  double y_ode[DIM];
  for(i=1;i<=m;i++) //These are the E classes
	y_ode[i] = 0;
  y_ode[0] = s;
  y_ode[m+1] = z;
  y_ode[m+2] = 0; //Dead
 
  DataDay = 0;



ODEParams.FxdPars[1] = k;
ODEParams.FxdPars[3] = mu;
ODEParams.FxdPars[6] = delta;
ODEParams.FxdPars[9] = (double) m;
ODEParams.FxdPars[11] = S0;
//Dead = 0.0;
double oldt = 0.0;
int Day = 0;
DataDay = 0;

	gamma = Pars[6]*m;
	

SumE = 0.0;



while((t<maxt)&&(TooHighOut!=1)){

	
  if(Figures==1){
	fprintf(fp,"%e %e %e %e\n",t,s,z,SumE);
	printf("%f %e %e %e\n",t,s,z,SumE);
	
  }
  

   //TimeIntrvl = (t/CorrlnIntrvl);  //TimeIntrvl indexes the Random Numbers - you want it to truncate the floats
   //nubarAdj = (RandNums[1][Plt][TimeIntrvl]);

   nubarAdj = (RandNums[1][Plt][Day]); //TEMP
   
   ODEParams.FxdPars[2] = nubarAdj;

   ///////////////////  Here is the ODE solver at work ///////////////////    
   double t_next = t + 1;
   y_ode[m+1] = z;
   //y_ode[m+2] = Dead;
   //printf("t:%f s:%f z:%f\n",t,y_ode[0],y_ode[m+1]); //TEMP

   t=ODE_Solver(t,t_next,&ODEParams,y_ode); //Numerical integration from t to t_next, y_ode holds pop densities
   s = y_ode[0];	
   z = y_ode[m+1]; //set pop densities equal to ode output

   //Dead = y_ode[m+2];
  

   if(Figures==1){
      	SumE = 0.00;
   	for(i=1;i<=m;i++) SumE = SumE + y_ode[i];
   }
 
  
   //DataDay = DataDay + (t - oldt);
   //printf("t:%f oldt:%f t-oldt:%f \n",t,oldt,t-oldt);

   DataDay++;

   oldt = t;
   Day++;




  // printf("DataDay:%d\n",DataDay);

   if((DataDay>= DataInterval)&&(TooHighOut!=1)){


	   //SumE = 0.00;
          SumE = 1e-300; //For some horrible reason, this works better than the preceding line
   	   for(i=1;i<=m;i++) SumE = SumE + y_ode[i];

	   Suscept[week] = s;
	   if(FractDead==1){ //Data are fraction dead (Woods & Elk)

 

			//if((Dead<=0)&&((s<=0)&&(y<=0))) 
			if((y_ode[m+2]<=0)&&((s<=0)&&(SumE<=0))) 
				WM[week] = 0.0;  
			else{
				//WM[week] = Dead/(s+y+Dead);  //note y is sum of exposed classes
				WM[week] = y_ode[m+2]/(s+SumE+y_ode[m+2]);  //note y is sum of exposed classes

				//if(Plt==8) printf("week:%d WM:%e s:%f y:%f Dead:%f \n", week,WM[week],s,y_ode[m+2]);
				//getc(stdin);
				
			}
		  //if(Figures==1){
			//fprintf(fp,"%e %e %e %e %e %e\n",t,s,z,y,y_ode[m+2],WM[week]);
			//if(Plt==8) printf("week:%d WM:%f s:%f y:%f Dead:%f \n", week,WM[week],s,y,y_ode[m+2]);
  		 // }
   } else { 

		   if((s<=0)&&(SumE<=0))
				WM[week] = 0.0;  
			else
				WM[week] = SumE/(s+SumE);  

	   }
           if(WM[week]>1) WM[week] = 1.0;
           if(WM[week]<0) WM[week] = 0.0;
 
	   

    DataDay = 0.0;  //end of week
    //Dead = 0.0;  //re-start cumulative dead
 
    y_ode[m+2] = 0.0;
    week++;      //update week number
   }

   if(((week==1)&&(Split))&&(flag)){
      //z += FractI*Hch*S0*ratio; OBSOLETE
	z += FractI2[Plt]*(1-FractI)*S0*ratio;
	
      flag = 0;
    }



   //t += h;  //Not any more - now the ode solver does the updating

   


 }    // t loop 


//printf("week:%d\n",week);
for(i=0;i<week;i++){
		if((isnan(WM[i])!=1)&&(TooHighOut!=1))
			ModelFI[Trial][Plt][i] = WM[i];
		else
			ModelFI[Trial][Plt][i] = 0.0;
		
		ModelS[Trial][Plt][i] = -1; //Suscept[i];
                //ModelFI[Trial][Plt][i] = 0.2;
                //fprintf(fp3,"%d\t%d\t%d\t%f\t%f\n",Trial,Plt,i,WM[i],Suscept[i]);

		//printf("ModelFI:%f WM:%f\n",ModelFI[Trial][Plt][i],WM[i]);
		

	}



	if(Figures==1){ 
		fclose(fp);
		exit(1);
	}
	
	return 1.0 -(s/S0Adj);


}  //The End







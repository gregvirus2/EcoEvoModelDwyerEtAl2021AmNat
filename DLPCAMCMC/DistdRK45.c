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



  k = Pars[1];  //model parameter: Squared CV of dist'n of transmission rates
  nubar = Pars[2];   //average transmission rate
  
  mu = Pars[3];

  ratio = Pars[5];
	
  
  //ratio = pow(10,Pars[5]); 

  delta = Pars[6]/Pars[9];
  



  CorrlnIntrvl = Pars[8];  //Correlation Interval

  double StagesTemp = Pars[9] + 0.5;
  int m = ( (int) StagesTemp);   //Number of Stages, rounded off
  int DIM = m; //DIM may be used elsewhere?  Otherwise this is a waste of time





ReDo = 1; 

	TooHighOut = 0;


Plt=Pars[0];

	for(i=0;i<20;i++){
		WM[i] = 0.0;
		Suscept[i] = 0.0;
	}



    S0 = Covars[1][Plt]; //Initial population size
   // FractI = Covars[2][Plt]; //Initial fraction infected

    z = Covars[2][Plt]; //THIS IS THE ONLY (MAIN?) DIFF OF DDE_SSE9B from DDE_SSE9
    Split = Covars[3][Plt];  //Whether hatch time is spread out, yes or no

	
    maxt = MaxDays[Plt]; //Epidemic end should be 56 days...
	

   
     

  
   Z0 = z;



   S0Adj = S0;  //Initial density of hosts, saved for later use
   initz = z;  //Saving initial virus density, not sure why?

   s = S0;
   Z0 = z;  //Emu
  

  week = 1;  //Initializing week number, need it for data storage

  flag = 1;  //Says whether or not additional hatch has occurred


  int OldTimeIntrvl = -1;  //Useful for checking the time interval action
 

  OutputCount = 1; 

double Oldnubar = -1;
Bang = 0; FirstCrash = 1;
double StageAdjust;


ind = 0;

//for(t=h;t<=maxt;t+=h){  //time loop

t = 0; TooHighOut = 0;



  double y_ode[DIM];


//  printf("z before y_ode:%f\n",z);

  for(i=1;i<=m;i++) //These are the E classes
	y_ode[i] = 0;
  y_ode[0] = s;
  y_ode[m+1] = z;
  y_ode[m+2] = 0; //Dead

// printf("s:%f z:%f\n",s,z);
 
  DataDay = 0;

//Parameters for ODE solver...
ODEParams.FxdPars[1] = k;
ODEParams.FxdPars[3] = mu;
ODEParams.FxdPars[6] = delta;
ODEParams.FxdPars[9] = (double) m;
ODEParams.FxdPars[11] = S0;
//Dead = 0.0;
double oldt = 0.0;
int Day = 0;
DataDay = 0;
	

SumE = 0.0;


//printf("s:%f z:%f k:%f nubar:%f mu:%f delta:%f ratio:%f m:%d maxt:%f\n",s,z,k,nubar,mu,delta,ratio,m,maxt); //getc(stdin);

//fp = fopen("junkRK45.dat","w");
 
while((t<maxt)&&(TooHighOut!=1)){

	
  // printf("%f %f %f\n",t,s,z);
   //fprintf(fp,"%f %f %f\n",t,s,z); //getc(stdin);
   //exit(1);

   nubarAdj = (RandNums[1][Plt][Day]); //TEMP  - NOT RELEVANT FOR DP
   //printf("Day:%d nubarAdj:%f\n",Day,nubarAdj); getc(stdin); 

      
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
  
 
  
   //DataDay = DataDay + (t - oldt);
   //printf("t:%f oldt:%f t-oldt:%f \n",t,oldt,t-oldt);

   DataDay++;

   oldt = t;
   Day++;




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

				//printf("week:%d WM:%f s:%f y:%f Dead:%f \n", week,WM[week],s,y,y_ode[m+2]);
				//getc(stdin);
				
			}
		  //if(Figures==1){
			//fprintf(fp,"%e %e %e %e %e %e\n",t,s,z,y,y_ode[m+2],WM[week]);
			//printf("week:%d WM:%f s:%f y:%f Dead:%f \n", week,WM[week],s,y,y_ode[m+2]);
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


	/*

   if(((week==1)&&(Split))&&(flag)){
      //z += FractI*Hch*S0*ratio; OBSOLETE
	z += FractI2[Plt]*(1-FractI)*S0*ratio;
	
      flag = 0;
    }

	*/



   //t += h;  //Not any more - now the ode solver does the updating

   


 }    // t loop 


//fclose(fp);
//exit(1);


//printf("week:%d\n",week);
/*
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
*/



	
	return 1.0 -(s/S0Adj);


}  //The End







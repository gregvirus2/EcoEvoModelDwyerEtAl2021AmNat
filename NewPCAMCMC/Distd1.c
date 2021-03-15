extern int DataSet;
extern int Figures;

extern float FractI2[20];


#define IntMax 600000
// tau divided by the step size


int PrintMe = 0; //whether to print output to a file


double F(double sus,double nubar, double path, double sus0, double k)
 {
	
    double temp, temp2;
    if(k<1e5){
	temp = pow(sus/sus0,1/k);
    }else{
	temp = 1.0;
    }
    //return( - sus * nubar * path * pow(sus/sus0,1/k));
    temp2 = (- sus * nubar * path * temp);
    if((isnan(temp2)==0)&&(isinf(temp2)==0)){
	 return(temp2);
    }else{
	return(0.0);
    }

  }

double G(double sus,double nubar, double path, double sus0, double k, double inf, double gamma)
  {

     double temp, temp2;
    if(k<1e5){
	temp = pow(sus/sus0,1/k);
    }else{
	temp = 1.0;
    }
   //return(sus * nubar * path * pow(sus/sus0,1/k) - gamma*inf);
   temp2 = (sus * nubar * path * temp - gamma*inf);
     if((isnan(temp2)==0)&&(isinf(temp)==0)){
	 return(temp2);
    }else{
	return(0.0);
    }


  }

double G2(double inf1, double gamma, double inf2)
{
    double temp;

    temp = (gamma*inf1 - gamma*inf2);
     if((isnan(temp)==0)&&(isinf(temp)==0)){
	 return(temp);
    }else{
	return(0.0);
    }

}




double H(double inf, double gamma, double path, double mu)
  {
   double temp;
   temp = (gamma*inf - mu*path);
    if((isnan(temp)==0)&&(isinf(temp)==0)){
	 return(temp);
    }else{
	return(0.0);
    }

  }

double Death(double inf, double gamma)
{

	double temp;

    	temp = (gamma*inf);
    if((isnan(temp)==0)&&(isinf(temp)==0)){
	 return(temp);
    }else{
	return(0.0);
    }


}


double Distd1(double *Pars,float **Covars,float ***ModelFI,float ***ModelS,double ***RandNums,int *MaxDays,int MaxPlots,int DataInterval, int FractDead){
 
	//This version assumes that Covars[2][1] is the actual virus pop, NOT the fractinf 

 double DataIntervald = DataInterval;
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
 //double l,l0,l1,l2,l3;
 double n,n0,n1,n2,n3; //Dead
 //double p,p0,p1,p2,p3;
 double q,q0,q1,q2,q3;
 double h;  /* time step */

 double z;
 int ind; //indicates where you are in storage arrays
 double *Intx;
 double *SStor, *PStor, *nuStor;
 double Stmp, Ptmp, nutmp;
 double l, p;
 //double m[100];
 double l0, l1, l2, l3;
 //double m0[100], m1[100], m2[100], m3[100];
double *m, *m0, *m1, *m2, *m3;
 double p0, p1, p2, p3;
double dd,dd0, dd1, dd2, dd3;
 //double l0tmp, l1tmp, l2tmp, l3tmp;
 //double p0tmp, p1tmp, p2tmp, p3tmp;
 
 double s; 
 double *inf; 
 //double inf[100];
 /* *initx; */
  /* z is pibs, and s is suscept.'s */
 /* Intx is sum over the suscepts. times nu times pibs z, for storage */
 double zinit;  /* initial pibs, to simulate overwintering */
 double Dead,y,oldy;      /* infect.'s - no need to store'em */
 int i;
 double t;
 int tau, tau0;
 float WM[1000], Suscept[1000];
 //double *WM, *Suscept;
 int week;
 double FractI;
 double S0;
 double nubar, mu, k, gamma;
 double nubarAdj; 
 int MaxStages;
 int CorrlnIntrvl, TimeIntrvl;
 double DataDay;
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

  //Figures = 1;
  //if(Figures==1) fp = fopen("Reg2WM.dat","w");


    

 //Intx = dvector(0,IntMax-1);  //allocates memory, more than needed
 //SStor = dvector(0,IntMax-1); 
 //PStor = dvector(0,IntMax-1); 
 //nuStor = dvector(0,IntMax-1);


  k = Pars[1];  //model parameter: Squared CV of dist'n of transmission rates
  nubar = Pars[2];   //average transmission rate
  
  mu = Pars[3];

  h = Pars[7];   //DDE time step

  CorrlnIntrvl = Pars[8];  //Correlation Interval

  double StagesTemp = Pars[9] + 0.5;
  MaxStages = ( (int) StagesTemp);   //Number of Stages
 
 

  m = dvector(0,MaxStages);
  m0 = dvector(0,MaxStages);
  m1 = dvector(0,MaxStages);
  m2 = dvector(0,MaxStages);
  m3 = dvector(0,MaxStages);
  inf = dvector(0,MaxStages);


  //ratio = pow(10,Pars[5]); //size of early stage larva to size of later stage larva, in # of virus particles OLD
  
  
     //if(DataSet>=15){
	//ratio = pow(10,-2.5);
     //}else{
	ratio = Pars[5];
     //}

//WM = dvector(0,50);
//Suscept = dvector(0,50);

/*
fp3 = fopen("WMOut.dat","wt");
fclose(fp3);
*/
/*
fp2 = fopen("WMOut12.dat","wt");
fclose(fp2);
*/


//for(Trial=1;Trial<=MaxTrials;Trial++){


ReDo = 1; h = Pars[7];


	TooHighOut = 0;


Plt=Pars[0]; 


	for(i=0;i<20;i++){
		WM[i] = 0.0;
		Suscept[i] = 0.0;
	}

	//tau0 = Pars[6];

	//tau = Pars[6]/h;  //time delay is tau, divide by h to get how much memory you need
	gamma = Pars[6]*MaxStages;
	//printf("delta:%f MaxStages:%d gamma:%f\n",Pars[6],MaxStages,gamma);
	//getc(stdin);


    S0 = Covars[1][Plt]; //Initial population size
	FractI = Covars[2][Plt]; //Initial fraction infected
	Split = Covars[3][Plt];  //Whether hatch time is spread out, yes or no

if(FractI<0){
 FractI = pow(10,Pars[11+Plt]);
}
	 //maxData = 9;  //Maximum data in any plot

	 maxt = MaxDays[Plt]; 
	 

/* OBSOLETE
  if(Split)
    Hch = 0.5;  //If there is some spread, half hatches in first week
   else
    Hch = 1.0;  //No spread, everyone hatches at once.
*/

   //z = S0*FractI*Hch*ratio;  //Initial density of virus, EMU only, I think
	//s = S0;
	//z = FractI; //FractI for Emma's code is actually the virus density...
	

/* OBSOLETE
   if(Split)
    Hch = 0.5;  //If there is some spread, half hatches in first week
   else
    Hch = 1.0;  //No spread, everyone hatches at once.
*/

  //ratio = pow(10,-1.0*Pars[5]); //size of early stage larva to size of later stage larva, in # of virus particles OLD
   ratio = Pars[5];

   //z = S0*FractI*Hch*ratio;  //Initial density of virus OBSOLETE
   z = S0*FractI*ratio;  //Initial density of virus - this is the new version, where if there's a split, it is completely handled below.  Search on FractI2 to find it
   
   /*
   Z0 = z;
   if(Split)
   	s = S0*(1-FractI)*(1-FractI2[Plt]);
   else
   	s = S0*(1-FractI);  //Initial density of hosts
   */


   S0Adj = S0;  //Initial density of hosts, saved for later use
   initz = z;  //Saving initial virus density, not sure why?

   //z = Covars[2][Plt]; //EMU
   Z0 = z;  //Emu
   //s = S0*(1-FractI);  //Initial density of hosts
   //s = S0; //this is for graphs for Emma's paper, or Karl when sprayed bugs start off in E class

   //z = S0*FractI*Hch*ratio;  //Initial density of virus
   //z = 0;
   //Z0 = z;
   for(i=0;i<MaxStages;i++)
	inf[i] = 0;
   //inf[0] = S0*FractI;  //Karl only?  I think...look out for ratio, below...



   /*
   for(i=0;i<tau;i++){
    //Intx[i] = 0.0;  //Initializing storage variable for delay
	SStor[i] = 0.0; PStor[i] = 0.0; nuStor[i] = 0.0;
	//p0[i] = 0.0; p1[i] = 0.0; p2[i] = 0.0; p3[i] = 0.0;
	//l0[i] = 0.0; l1[i] = 0.0; l2[i] = 0.0; l3[i] = 0.0; 
   }
   */

  week = 1;  //Initializing week no, need it for data storage

  flag = 1;  //Says whether or not additional hatch has occurred

  y= 0.0;  //Initializing number dead per unit time
  Dead = 0.0;  //Initializing cumulative number dead per week
  int OldTimeIntrvl = -1;  //Useful for checking the time interval action
  DataDay = 0.0; 

  OutputCount = 1; 

double Oldnubar = -1;
Bang = 0; FirstCrash = 1;
double StageAdjust;


ind = 0;

//for(t=h;t<=maxt;t+=h){  //time loop

//t = h; 
t = 0;
TooHighOut = 0;

//TEMP
//s = 25.000000; z=0.100000; Split = 0;
//k=20.836960; nubar=0.100000; mu=0.283986; double delta=0.283251; MaxStages=5; S0Adj=25.000000; initz=0.100000; maxt=56.000000; //TEMP
//gamma = MaxStages*delta;

while((t<=maxt)&&(TooHighOut!=1)){

	
	if(Figures==1){
		fprintf(fp,"%e %e %e %e\n ",t,s,z,y);
		printf("%f %e %e %e\n",t,s,z,y);
	}
	
	
		TimeIntrvl = (t/CorrlnIntrvl);  //TimeIntrvl indexes the Random Numbers - you want it to truncate the floats

		//nubarAdj = exp(RandNums[Trial][Plt][TimeIntrvl]);
		//nubarAdj = (RandNums[1][Plt][TimeIntrvl]);
		nubarAdj = nubar; //TEMP

		
		//nutmp = nuStor[ind]; //old garbage, I think as of 24 May 2012
		//if(t<=14) nubarAdj*=ratio;
		
	//nuStor[ind] = nubarAdj; //old garbage, I think as of 24 May 2012


    //DDE solver: method of steps using 4th order R-K
		
	//l0tmp = l0[ind]; l1tmp = l1[ind]; l2tmp = l2[ind]; l3tmp = l3[ind]; 
	//p0tmp = p0[ind]; p1tmp = p1[ind]; p2tmp = p2[ind]; p3tmp = p3[ind];
	//Stmp = SStor[ind]; Ptmp = PStor[ind];

    l0 = F(s,nubarAdj,z,S0Adj,k);
	m0[0] = G(s,nubarAdj,z,S0Adj,k,gamma,inf[0]);
	for(i=1;i<MaxStages;i++){
		m0[i] = G2(inf[i-1],gamma,inf[i]);
       }

    p0 = H(inf[MaxStages-1],gamma,mu,z);
    dd0 = Death(inf[MaxStages-1],gamma);
		
    l1 = F(s+0.5*l0*h,nubarAdj,z+0.5*p0*h,S0Adj,k);
	m1[0] = G(s+0.5*l0*h,nubarAdj,z+0.5*p0*h,S0Adj,k,inf[0]+0.5*m0[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m1[i] = G2(inf[i-1]+0.5*m0[i-1]*h,gamma,inf[i]+0.5*m0[i]*h);
	p1 = H(inf[MaxStages-1]+0.5*m0[MaxStages-1]*h,gamma,z+0.5*p0*h,mu);
	dd1 = Death(inf[MaxStages-1]+0.5*m0[MaxStages-1]*h,gamma);
	
    l2 = F(s+0.5*l1*h,nubarAdj,z+0.5*p1*h,S0Adj,k);
	m2[0] = G(s+0.5*l1*h,nubarAdj,z+0.5*p1*h,S0Adj,k,inf[0]+0.5*m1[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m2[i] = G2(inf[i-1]+0.5*m1[i-1]*h,gamma,inf[i]+0.5*m1[i]*h);
    p2 = H(inf[MaxStages-1]+0.5*m1[MaxStages-1]*h,gamma,z+0.5*p1*h,mu);
    dd2 = Death(inf[MaxStages-1]+0.5*m1[MaxStages-1]*h,gamma);
	
    l3 = F(s+l2*h,nubarAdj,z+p2*h,S0Adj,k);
	m3[0] = G(s+l2*h,nubarAdj,z+p2*h,S0Adj,k,inf[0]+m2[0]*h,gamma);
	for(i=1;i<MaxStages;i++)
		m3[i] = G2(inf[i-1]+m2[i-1]*h,gamma,inf[i]+m2[i]*h);
	p3 = H(inf[MaxStages-1]+m2[MaxStages-1]*h,gamma,z+p2*h,mu);
	dd3 = Death(inf[MaxStages-1]+m2[MaxStages-1]*h,gamma);


   
    
    l = (l0 + 2*l1 + 2*l2 + l3)/6.0;
	for(i=0;i<MaxStages;i++)
		m[i] = (m0[i] + 2*m1[i] + 2*m2[i] + m3[i])/6.0;
    p = (p0 + 2*p1 + 2*p2 + p3)/6.0;

    dd = (dd0 + 2*dd1 + 2*dd2 + dd3)/6.0;

	int Stop = 0;
/*
		s += l*h;
		z += p*h;
	Dead += dd*h;
*/


	if((isnan(l)==0)&&(isinf(l)==0)){
		s += l*h;
	}else{
		
       	TooHighOut = 1;
		s = 0.0;
	}
	
	if((isnan(p)==0)&&(isinf(p)==0)){
	
		z += p*h;
	}else{

		

		TooHighOut = 1;
		z = 0.0;
	}
	

	
	if((isnan(dd)==0)&&(isinf(dd)==0)){
		Dead += dd*h;
	}else{
		

		Dead = 0.0;
		TooHighOut = 1;
	}
	
	
       y = 0;
	for(i=0;i<MaxStages;i++){
			
			//if(m[i]==m[i]){
			
			if( (isnan(m[i])==0) && (isinf(m[i])==0)){
				inf[i] += m[i]*h;
			}else{
				inf[i] = 0.0;
				TooHighOut = 1;
			}
			if(inf[i]>0){
				y += inf[i];
			}else{
				inf[i] = 0.0;
			}
			//if(Figures==1) fprintf(fp,"%f ",inf[i]);	
         } //for loop

		
	

   
  // printf("t:%f DataDay - DataIntervald:%e TooHighOut:%d\n",t,DataDay-DataIntervald,TooHighOut);
   if((DataDay>= DataIntervald)&&(TooHighOut!=1)){

	   Suscept[week] = s;
	   if(FractDead==1){ //Data are fraction dead (Woods & Elk)

 

			if((Dead<=0)&&((s<=0)&&(y<=0))) 
				WM[week] = 0.0;  
			else{
				WM[week] = Dead/(s+y+Dead);  //note y is sum of exposed classes
				
			}
			//if(Figures==1){
			//	fprintf(fp,"%e %e %e %e %e %e\n",t,s,z,y,Dead,WM[week]);
				//printf("week:%d WM:%f s:%f y:%f Dead:%f \n", week,WM[week],s,y,Dead);
				//getc(stdin);
  		  	//}
   } else { 

		   if((s<=0)&&(y<=0))
				WM[week] = 0.0;  
			else
				WM[week] = y/(s+y);  

	   }
           if(WM[week]>1) WM[week] = 1.0;
           if(WM[week]<0) WM[week] = 0.0;
/*
      fp2 = fopen("WMOut12.dat","at");
      //fprintf(fp2,"%e\t%d\t%d\t%f\t%f\t%f\t%f\n",t,Plt,week,WM[week],Dead,s,y);
      fprintf(fp2,"%e %e %e \n",t,s,z);
      fclose(fp2);
*/
      
	   

    DataDay = 0.0;  //end of week
    Dead = 0.0;  //re-start cumulative dead
    week++;      //update week number
   }

    /*
   if(((week==1)&&(Split))&&(flag)){
      //z += FractI*Hch*S0*ratio; OBSOLETE
	z += FractI2[Plt]*(1-FractI)*S0*ratio;
	
      flag = 0;
    }
    */



   t += h;
   DataDay +=h;  


 }    // t loop 


	for(i=0;i<week;i++){
		if((isnan(WM[i])!=1)&&(TooHighOut!=1))
			ModelFI[Trial][Plt][i] = WM[i];
		else
			ModelFI[Trial][Plt][i] = 0.0;
		
		ModelS[Trial][Plt][i] = -1; //Suscept[i];
                //ModelFI[Trial][Plt][i] = 0.2;
                //fprintf(fp3,"%d\t%d\t%d\t%f\t%f\n",Trial,Plt,i,WM[i],Suscept[i]);
		

	}



	if(Figures==1){
		fclose(fp); 
		exit(1);
	}
	
	return 1.0 -(s/S0Adj);


}  //The End







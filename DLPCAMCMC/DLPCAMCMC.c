int Verbose;
int Verbose2; 
int CrashOut=0; 
int ExtRlzns = 0;
int Count = 0;
int Evoln = 0;   //Search on .txt
int gammaFix = 0; //gammaFix = 1 means set gamma to 0
int sigmaFix = 1; //sigmaFix = 1 means set sigma to 0 (0.01 for no-evolution model)
int WAIC = 1; //EPSILON!!!
int PrdAmpStor = 1; //Have to check the file names; search on PrdAmpStor
int LineSearch = 0;
int NiceVerbose = 0;  //use this to make time series figs.  Everything else must be turned off, including LineSearch, and WAIC, but gammaFix and sigmaFix still matter
int sFix = 1;
int aFix = 1;
int bFix = 0;
int DIC = 0; //As of 19 August 2020, this is obsolete; it's no longer DIC, instead its NiceVerbose
int ModelVsData = 0; //19 August 2020, this is obsolet, it used to be ture tha tDIC must also equal 1, and NiceVerbose, but now it's just NiceVerbose
int sampSize = 6;  //5 is for return times, 6 is for visual inspection
double GlobalAvgAmp, GlobalAvgPrd;
int PopWiseNumSamp;
int DataSet = 2; //16;  //2 is woods & elk, 16 is Otvos et al //Not relevant here, right?

int GM = 1; 
int MAXN = 300; //1500 for time series plots, 300 for fitting
int ReturnTimeStat = 1;

double PopWiseMean[100], PopWiseM2Log[100], PopWiseVarLog[100], PopWiseMeanLog[100];



typedef struct {

	double FxdPars[30];

} ParamStruct;

#include <time.h>
#include "stdio.h"
#include "string.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_odeiv.h>	// ODE solver


#include "math.h"
#include "stdio.h"

#include "fast_odesGD.h"


//#include "DataFuncs.c"

#include "Nrutil.c"
#include "NRUTIL.H"
//#include "LHoodFunc.c"

long seed;
#include "Distd1.c"

#include "DistdRK45.c"

#include "Zbrentd2.c"
#include "FX.c"


     #include <stdlib.h>
     #include <gsl/gsl_math.h>
     #include <gsl/gsl_monte.h>
     #include <gsl/gsl_monte_plain.h>
     #include <gsl/gsl_monte_miser.h>
     #include <gsl/gsl_monte_vegas.h>


gsl_rng *r2;



double mean(int length, double *x){

	int i;
	double answer = 0.0;
	for(i=0;i<length;i++) answer += x[i];
	answer /= length;
	
	return(answer);
}

int ReturnTime(int End, double *TimeSeries, double *Amp, double *Prd){

	double Periods[1000];
	int FirstTime = 0; 
	int UpFirst = 0; 
	int DownFirst = 0;
	int Flag1 = 0;
	int Flag2 = 0;
	int Above = 0;
	int Below = 0;
	int i, Last;
	int MaxTimeStart, MaxTimeStop, MinTimeStart, MinTimeStop;

	int FirstTrans = 0;

	double UpTime[1000];
	double DownTime[1000];
	double AmpDiff[1000];
	double PrdDiff[1000];
	double Biggest[1000], Smallest[1000];
	double MaxPop[1000], MinPop[1000];
	
	int j = 0; 
	int k = 0;
	double AvgAmp = 0;

	MaxPop[0] = -1e10;
	MinPop[0] = 1e10;
		
	double TestStat = mean(End,TimeSeries);
	//printf("TestStat:%f\n",TestStat);

	for(i=1;i<(End-1);i++){
		
		//printf("FirstTrans:%d TimeSeries:%f\n",FirstTrans,TimeSeries[i]);

		
	
		//Until a crossing occurs, just ignore the numbers
	if(FirstTrans==0){
			if((TimeSeries[i-1]<TestStat)&&(TimeSeries[i]>TestStat)){ FirstTrans = 1; Below = 1;}  //You were above, now you're below
			if((TimeSeries[i-1]>TestStat)&&(TimeSeries[i]<TestStat)){FirstTrans = 1; Above = 1;} //You were below,  now you're above
	}

	if(FirstTrans==1){
		if(TimeSeries[i]>TestStat){
							
				
				
			

				if(Below==1){ //if you just switched from below to above

						UpTime[j] = i;
						MaxPop[j] = TimeSeries[i]; //Set MaxPop equal to the current pop size
						Below = 0; //switch to below mode
						Above = 1; //switch to above mode
						if(Flag2 == 1) k = k + 1; //next time you go below, you want to have a higher "below number"
						//printf("1st Above:%d i:%d TimeSeries[i]:%f j:%d UpTime[j]:%f MaxPop[j]:%f\n",Below,i,TimeSeries[i],j,UpTime[j],MaxPop[j]);



				}else{ 
						
						if(TimeSeries[i]>MaxPop[j])  //if you didn't just switch, only save the pop size if it's bigger
							MaxPop[j] = TimeSeries[i];


						//printf("2nd Above:%d i:%d TimeSeries[i]:%f j:%d UpTime[j]:%f MaxPop[j]:%f\n",Below,i,TimeSeries[i],j,UpTime[j],MaxPop[j]);
						//getc(stdin);
				}


				Flag1 = 1;
				if(FirstTime==0){

					FirstTime = 1;
					UpFirst = 1;
					
					
				}
				//getc(stdin);


			
		} //TimeSeries[i]>TestStat
		
		if(TimeSeries[i]<TestStat){

				//printf("i:%d\n",i);

				
				

				//printf("Below\n");
				if(Above==1){ //if you just switched from below to above

					DownTime[k] = i;
					MinPop[k]= TimeSeries[i]; //Set MaxPop equal to the current pop size
					Above = 0;
					Below = 1;
					if(Flag1==1) j = j + 1;
					//printf("1st Below:%d i:%d TimeSeries[i]:%f k:%d DownTime[k]:%f MinPop[k]:%f\n",Below,i,TimeSeries[i],k,DownTime[k],MinPop[k]);

				}else{ 

					if(TimeSeries[i]<MinPop[k])
						MinPop[k] = TimeSeries[i];

					//printf("2nd Below:%d i:%d TimeSeries[i]:%f k:%d DownTime[k]:%f MinPop[k]:%f\n",Below,i,TimeSeries[i],k,DownTime[k],MinPop[k]);

				}




				Flag2 = 1;		
				
				if(FirstTime==0){
					FirstTime = 1;
					UpFirst = 0;
				}

		}

		
          }//FirstTrans
	} //For loop?

//exit(1);


//printf("j:%d k:%d \n",j,k);

	if(k>j){
		Last = j;
	}else{
		Last = k;
	}

	//printf("Last:%d\n",Last);
	//for(j=0;j<Last;j++)
		//printf("j:%d UpTime[j]:%e\n",j,UpTime[j]);

	for(i=0;i<Last;i++){

		Amp[i] = MaxPop[i] - MinPop[i];
		if(i==0){
			Prd[i] = -1;
		} else{
			Prd[i] = ((int) fabs((double) (UpTime[i] - UpTime[i-1])));
		}

		//printf("i:%d Max:%f Min:%f Amp:%f Prd:%f UpTime[i]:%f UptTime[i-1]:%f\n",i,MaxPop[i],MinPop[i],Amp[i],Prd[i],UpTime[i],UpTime[i-1]);

	}
	//getc(stdin);

	//double MaxPop[100], MinPop[100];
	

	//MaxPop = numeric(); MinPop = numeric();	
	/*
	k = 1;
	//printf("TimeSeries:%f \n",TimeSeries[0]);
	int time;
	for(j=0;j<Last;j++){
	
		//Up = UpTime[j]; Down = DownTime[j];
		if(j==0){
			PrdDiff[j] = -1;
			if(UpFirst==1){

				//MaxPopTime = c(1:UpTime[j]);
				//MinPopTime = c(UpTime[j]:DownTime[j]);
				MaxTimeStart = UpTime[j]; MaxTimeStop = UpTime[j];
				MinTimeStart = UpTime[j]; MinTimeStop = DownTime[j];
							

			} else{ //UpFirst == 0
				
				//MinPopTime = c(1:DownTime[j]);
				//MaxPopTime = c(DownTime[j]:UpTime[j]);

				MinTimeStart = 0; MinTimeStop = DownTime[j];
				MaxTimeStart = DownTime[j]; MaxTimeStop = UpTime[j];
				
			} //else UpFirst == 0
		}else{ //j !=0
			PrdDiff[j] = UpTime[j] - UpTime[j-1];
			
			//MaxPopTime = c(DownTime[j-1]:UpTime[j]);
			//MinPopTime = c(UpTime[j-1]:UpTime[j]);
			MaxTimeStart = DownTime[j-1]; MaxTimeStop = UpTime[j];
			MinTimeStart = UpTime[j-1]; MinTimeStop = DownTime[j];
		}

		Biggest[j] = TimeSeries[MaxTimeStart];
		//printf("j:%d MaxTimeStart:%d TimeSeries:%f Biggest:%f\n",j,MaxTimeStart,TimeSeries[j],Biggest[j]);
		for(time=(MaxTimeStart+1);time<MaxTimeStop;time++)
			if(TimeSeries[time]>Biggest[j]) Biggest[j] = TimeSeries[time];

		Smallest[j] = TimeSeries[MinTimeStart];
		for(time=(MinTimeStart+1);time<MinTimeStop;time++)
			if(TimeSeries[time]<Smallest[j]) Smallest[j] = TimeSeries[time];

		//MaxPop[j] = max(TimeSeries[MaxPopTime]);
		//MinPop[j] = min(TimeSeries[MinPopTime]);
				
		
			
	} //for j

	for(j=0;j<Last;j++) {
		AmpDiff[j] = Biggest[j]-Smallest[j];
		//printf("j:%d Biggest:%f Smallest:%f AmpDiff:%f PrdDiff:%f \n",j,Biggest[j],Smallest[j],AmpDiff[j],PrdDiff[j]);
		
	}
	*/

	//for(j in 1:length(MaxPop))
		//AmpDiff = MaxPop - MinPop;
	//printf("Last:%d\n",Last);
	if((Flag1==1)&&(Flag2==1)){
		return(Last);
	}else{
		//cat("Flag2:",Flag2,"\n");
		//printf("returning 0...\n"); //getc(stdin);
		return(0);
	}
						
} //#End of ReturnTime

int AmpCalc(int End, double *TimeSeries, double *Amp, double *Prd){

	double Periods[1000];
	int FirstTime = 0; 
	int UpFirst = 0; 
	int DownFirst = 0;
	int Flag1 = 0;
	int Flag2 = 0;
	int i, Last;
	int MaxTimeStart, MaxTimeStop, MinTimeStart, MinTimeStop;

	double UpTime[1000];
	double DownTime[1000];
	double AmpDiff[1000];
	double PrdDiff[1000];
	double Biggest[1000], Smallest[1000];
	double MaxPop[1000], MinPop[1000];
	
	int j = 0; 
	int k = 0;
	double AvgAmp = 0;
		
	
	for(i=1;i<(End-1);i++){


		//if((TimeSeries[i-1]<mean(TimeSeries))&&(TimeSeries[i]>mean(TimeSeries))){
		if((TimeSeries[i-1]<TimeSeries[i])&&(TimeSeries[i+1]<TimeSeries[i])){
			
				
				UpTime[j] = i;
				//printf("Inside if TimeSeries[i-1] - TimeSeries[i]:%e TimeSeries[i]:%f TimeSeries[i+1]-TimeSeries[i]:%e UpTime:%f\n",TimeSeries[i-1]-TimeSeries[i],TimeSeries[i],TimeSeries[i+1]-TimeSeries[i],UpTime[j]);
				//printf("Inside if j:%d UpTime[j]:%f\n",j,UpTime[j]);
				MaxPop[j] = TimeSeries[i];
			 	Flag1 = 1;
				j = j + 1;
				if(FirstTime==0){
					FirstTime = 1;
					UpFirst = 1;
				}
			
		}
		//if((TimeSeries[i-1]>mean(End,TimeSeries))&&(TimeSeries[i]<mean(End,TimeSeries))){
		if((TimeSeries[i-1]>TimeSeries[i])&&(TimeSeries[i+1]>TimeSeries[i])){

				//printf("i:%d\n",i);

				DownTime[k] = i;
				//printf("Inside if k:%d DownTime[k]:%f\n",k,DownTime[k]);

				
				
				MinPop[k] = TimeSeries[i];
				Flag2 = 1;		
				k = k + 1;
				if(FirstTime==0){
					FirstTime = 1;
					UpFirst = 0;
				}

		}

		//printf("j:%d UpTime[j]:%f\n",j,UpTime[j]);

	}


//printf("j:%d k:%d \n",j,k);

	if(k>j){
		Last = j;
	}else{
		Last = k;
	}

	//printf("Last:%d\n",Last);
	//for(j=0;j<Last;j++)
		//printf("j:%d UpTime[j]:%e\n",j,UpTime[j]);

	for(i=0;i<Last;i++){

		Amp[i] = MaxPop[i] - MinPop[i];
		if(i==0){
			Prd[i] = -1;
		} else{
			Prd[i] = ((int) fabs((double) (UpTime[i] - UpTime[i-1])));
		}

		//printf("i:%d Max:%f Min:%f Amp:%f Prd:%f UpTime[i]:%f UptTime[i-1]:%f\n",i,MaxPop[i],MinPop[i],Amp[i],Prd[i],UpTime[i],UpTime[i-1]);

	}
	//getc(stdin);

	//double MaxPop[100], MinPop[100];
	

	//MaxPop = numeric(); MinPop = numeric();	
	/*
	k = 1;
	//printf("TimeSeries:%f \n",TimeSeries[0]);
	int time;
	for(j=0;j<Last;j++){
	
		//Up = UpTime[j]; Down = DownTime[j];
		if(j==0){
			PrdDiff[j] = -1;
			if(UpFirst==1){

				//MaxPopTime = c(1:UpTime[j]);
				//MinPopTime = c(UpTime[j]:DownTime[j]);
				MaxTimeStart = UpTime[j]; MaxTimeStop = UpTime[j];
				MinTimeStart = UpTime[j]; MinTimeStop = DownTime[j];
							

			} else{ //UpFirst == 0
				
				//MinPopTime = c(1:DownTime[j]);
				//MaxPopTime = c(DownTime[j]:UpTime[j]);

				MinTimeStart = 0; MinTimeStop = DownTime[j];
				MaxTimeStart = DownTime[j]; MaxTimeStop = UpTime[j];
				
			} //else UpFirst == 0
		}else{ //j !=0
			PrdDiff[j] = UpTime[j] - UpTime[j-1];
			
			//MaxPopTime = c(DownTime[j-1]:UpTime[j]);
			//MinPopTime = c(UpTime[j-1]:UpTime[j]);
			MaxTimeStart = DownTime[j-1]; MaxTimeStop = UpTime[j];
			MinTimeStart = UpTime[j-1]; MinTimeStop = DownTime[j];
		}

		Biggest[j] = TimeSeries[MaxTimeStart];
		//printf("j:%d MaxTimeStart:%d TimeSeries:%f Biggest:%f\n",j,MaxTimeStart,TimeSeries[j],Biggest[j]);
		for(time=(MaxTimeStart+1);time<MaxTimeStop;time++)
			if(TimeSeries[time]>Biggest[j]) Biggest[j] = TimeSeries[time];

		Smallest[j] = TimeSeries[MinTimeStart];
		for(time=(MinTimeStart+1);time<MinTimeStop;time++)
			if(TimeSeries[time]<Smallest[j]) Smallest[j] = TimeSeries[time];

		//MaxPop[j] = max(TimeSeries[MaxPopTime]);
		//MinPop[j] = min(TimeSeries[MinPopTime]);
				
		
			
	} //for j

	for(j=0;j<Last;j++) {
		AmpDiff[j] = Biggest[j]-Smallest[j];
		//printf("j:%d Biggest:%f Smallest:%f AmpDiff:%f PrdDiff:%f \n",j,Biggest[j],Smallest[j],AmpDiff[j],PrdDiff[j]);
		
	}
	*/

	//for(j in 1:length(MaxPop))
		//AmpDiff = MaxPop - MinPop;
	//printf("Last:%d\n",Last);
	if((Flag1==1)&&(Flag2==1)){
		return(Last);
	}else{
		//cat("Flag2:",Flag2,"\n");
		printf("returning 0...\n");
		return(0);
	}
						
} //#End of AmpCalc

double LHoodFunc(double ModelPrd, double ModelAmp){

	int i;
	
	//double DataPrds[] = {10,9,9}; //Not sure of origin of theses
	int NumPrds;
	//DataPrds[] = {10,9,9,9,9,9,9,9,9,10}; //What are these?
	//double DataAmps[] = {2.5742, 2.468, 2.4043, 4.879798728, 4.006998728, 4.238607315}; //These use every sequential peak-trough combination
	//GM first
 	//double DataAmps[] = {3.48, 3.67, 3.54, 3.3, 4.56}; // first 4 are Williams et al 91, last is Ostfeld
	//double DataPrds[] = {11,11,11,11,10};
	double DataAmps[10], DataPrds[10];

	if(GM==1){ //GM
		     
		//First two lines seem to have errors in the amps, and in 2 of the Williams et al prds, so...
		//DataAmps[0] = 3.48; DataAmps[1] = 3.67; DataAmps[2] = 3.54; DataAmps[3] = 3.30; DataAmps[4] = 4.56; DataAmps[5] = 2.86;  //First 4 are Williams et al. 91, 5th is Ostfeld, 6th is Skaller
		//DataPrds[0] = 11; DataPrds[1] = 11; DataPrds[2] = 11; DataPrds[3] = 10; DataPrds[4] = 10;  DataPrds[5] = 10;
		
		//Here are the corrected versions - August 2020 - based on visual inspection
		//DataAmps[0] = 3.48; DataAmps[1] = 3.67; DataAmps[2] = 3.54; DataAmps[3] = 3.30; DataAmps[4] = 4.56; DataAmps[5] = 2.86;  //First 4 are Williams et al. 91, 5th is Ostfeld, 6th is Skaller
		//DataPrds[0] = 11; DataPrds[1] = 11; DataPrds[2] = 10; DataPrds[3] = 10; DataPrds[4] = 10;  DataPrds[5] = 10;
	    
		//This set based on return times, but including Skaller from visual inspection - August 2020
		DataAmps[0] = 3.97; DataAmps[1] = 3.97; DataAmps[2] = 3.67; DataAmps[3] = 3.37; DataAmps[4] = 4.24; DataAmps[5] = 2.86;  //First 4 Williams et al. 91, 5th Ostfeld, 6th Skaller
		DataPrds[0] = 11; DataPrds[1] = 12; DataPrds[2] = 10; DataPrds[3] = 9; DataPrds[4] = 9; DataPrds[5] = 10;
	     
	} else{ //DFTM

		DataAmps[0] = log10(exp(7.53)); DataAmps[1] = log10(exp(7.05)); DataAmps[2] = log10(exp(6.06)); DataAmps[3] = log10(exp(5.61));
		DataPrds[0] = 9; DataPrds[1] = 9; DataPrds[2] = 6; DataPrds[3] = 8;		

	}	

	int NumAmps;
	double LHPrd, sigPrd;
	double LHAmp, sigAmp;
	double LH; //Overall LH
	double rho; //Correlation coefficient
	int samp;


	NumPrds = NumAmps = sampSize;

	sigPrd = 0.0;
	for(i=0;i<NumPrds;i++){
		sigPrd += (ModelPrd-DataPrds[i])*(ModelPrd-DataPrds[i]);
	}
	sigPrd /= NumPrds;
	sigPrd = pow(sigPrd,0.5);

	/*
	LHPrd = 0.0; 
	for(i=0;i<NumPrds;i++){
		printf("ModelPrds:%f\n",ModelPrd);
		LHPrd += log(gsl_ran_gaussian_pdf((ModelPrd-DataPrds[i])*(ModelPrd-DataPrds[i]),sigPrd));
	}
	*/

	sigAmp = 0.0;
	for(i=0;i<NumAmps;i++){
		sigAmp += (ModelAmp-DataAmps[i])*(ModelAmp-DataAmps[i]);
	}
	sigAmp /= NumAmps;
	sigAmp = pow(sigAmp,0.5);

	rho = 0;
	for(samp=0;samp<sampSize;samp++){

		rho += ((DataAmps[samp]-ModelAmp)*(DataPrds[samp]-ModelPrd))/(sigAmp*sigPrd);

	}
	rho /= sampSize;

	LH = 0.0;
	
	for(samp=0;samp<sampSize;samp++){

		double diff1 = DataAmps[samp]-ModelAmp;
		double diff2 = DataPrds[samp]-ModelPrd;
		
		double dummy = (gsl_ran_bivariate_gaussian_pdf(diff1,diff2,sigAmp,sigPrd,rho));
		LH += log(gsl_ran_bivariate_gaussian_pdf(diff1,diff2,sigAmp,sigPrd,rho));

		if(WAIC==1){
			if((isnan(dummy)!=1)&&(isinf(dummy)!=1)){
				PopWiseMean[samp] += dummy;
				PopWiseMeanLog[samp] += log(dummy);
				PopWiseM2Log[samp] += log(dummy)*log(dummy);
				//printf("PopWiseNumSamp:%d dummy:%e PopWiseMean:%e PopWiseMeanLog:%e PopWiseM2Log:%e\n",PopWiseNumSamp,dummy,PopWiseMean[samp],PopWiseMeanLog[samp],PopWiseM2Log[samp]); //getc(stdin);
			}
		}
		

		
		
		//if(WAIC==1)
			//printf("ModelAmp:%f ModelPrd:%f diff1:%f diff2:%f sigAmp:%f sigPrd:%f rho:%f dummy:%e LH:%f \n",ModelAmp,ModelPrd,diff1,diff2,sigAmp,sigPrd,rho,dummy,LH);

	}
	//getc(stdin);
	//printf("ModelAmp:%f ModelPrd:%f \n",ModelAmp,ModelPrd);
	
	if(isnan(ModelAmp)!=1) GlobalAvgAmp += ModelAmp;
	if(isnan(ModelPrd)!=1) GlobalAvgPrd += ModelPrd;


	
	
	//return(LHPrd+LHAmp);

	if(isnan(LH)==1) LH = -1e10;

	return(LH);

	

}

double LHood(double *FxdPars, double kIn, double nuIn, double muIn, double sigmaIn, double ratioIn, double deltaIn, double mIn){

//double kIn[5000],nuIn[5000],muIn[5000],ratioIn[5000],deltaIn[5000],mIn[5000],
			
	long n, maxn = MAXN;
	long nMin = MAXN-200;
	if(nMin>=maxn) nMin = 0;
	double NStor[2000];
	int Dump = 0;

	double Pars[10]; //Plt, k, nubar, mu, ?, log10(ratio), tau, h, CorrlnIntrvl, MaxTrials
	float **Covars;
	float ***ModelFI, ***ModelS;
	double ***RandNums;
	int *Split;
	int i, j, k;
	int *MaxDays;
	int MaxPlots = 1;
	int DataInterval = 7;
	int FractDead = 0;
	int MaxRlzns = 1;
	FILE *fp, *fp2, *fp3, *fpModel;
	double FIOut;
	double FI;
	
	ModelFI = f3tensor(1,MaxRlzns,1,25,0,25); //indices:Rlzns,plots,weeks
       ModelS = f3tensor(1,MaxRlzns,1,25,0,25);
	RandNums = d3tensor(1,MaxRlzns,1,25,0,1000);
	Split = ivector(1,100);
	MaxDays = ivector(1,100);
	Covars = matrix(1,25,1,25);

		double cx = 0.0, dx = 1.0, tol = 1e-60;
	double kParam = 1.35;
	Pars[1] = kParam;
	double numu = Pars[2]/Pars[3];
	double N, oldN; 
	double Z, oldZ;
	double nu, oldnu;
	double denom, numer;
	double InfN, oldInfN;
	double InfZ, oldInfZ;
	double Infnu, oldInfnu;
	double Amp[1000], Prd[1000];	
		 
	double FractInf;
	double lambda = 5.5;
	double r = 0.2;
        double V, predation;
	double b = 1.0;
	double a = 0.14;
	double w = 0.967;
        double gamma = 0.2;
	double phi = 10;  //lambda 2 and phi 17.5 - cycles, apparently, but phi 17.8 crashes (k = 1.1)
 	double sigma = 0.0;
	double sigmaLongTerm;
	double Neqm, Zeqm, nuEqm, iEqm, lamV;
	double epsilon = 1e-5; //1e-3; //5e-6; //1e-4; //5e-5; // 1e-7; //5e-6; //5e-5; //1e-4; //1e-6; //Default: 1e-5;

	Pars[0] = 1.0; //Plt
	Pars[1] = 0.5; //k
	Pars[2] = 0.45;
	Pars[3] = 0.412361; 
	Pars[4] =  -1; //sigma
	Pars[5] =  -2.8; //-1.7; Colin //-2.8;  //best-fit GM //ratio
	Pars[7] = 2e-2; //1e-4; //h
	Pars[8] = 1;  //CorrlnIntrvl
	Pars[9] = 20; //number of stages, let's say
	

	Pars[6] = Pars[9]/12.0;
	MaxDays[1] = 56;

	double gauss;
	double s;
	
		
	/* //Debug Test 1
	V = 4;
	r = 0.1;
	s = 2.0;
	gamma = 0.0;
	phi = 7; //TEMP
	*/
	
	/* //Debug Test 2	
	Evoln = 0;
	lambda = 9; //TEMP
	nuIn = 2;
	muIn = 1.0/5.0;
	phi = 25.0;
	kIn = kParam = 1.48721; //TEMP
	gamma = 0.0;
	double lamV = pow(lambda,V);
	double Neqm = (lamV - 1)*lambda/(V*(1+phi-gamma)*(lambda-1));
	Neqm *= (muIn/nuIn);
	double Zeqm = phi*(lamV-1)/(V*(1+phi-gamma));
	Zeqm *= (muIn/nuIn);
	*/

	//Debug Test 3
	/*
	Evoln = 1;
	lambda = 9; //TEMP
	nuIn = 2;
	muIn = 1.0;
	kIn = kParam = 1.0/5.0; //TEMP
	ratioIn = 0.02;
	deltaIn = 0.064;
	mIn = 9;
	V = 1/kIn;
	gamma = 0.1;
	r = FxdPars[1] = 0.2;
	s = FxdPars[2] = 45.0;
	phi = FxdPars[3] = 2.0;
	//maxn = 300;
	*/
	
	//b = 1;
	/*
	double cx2 = 0.0;
	double dx2 = 0.9999;
	double tol2 = 1e-60;
	double iEqm = zbrentd3(nuEqmfx,s,b,V,r,cx2,dx2,tol2);
	double nuEqm = (1 - r*(1-iEqm))/(r*s*pow(1-iEqm,V+1));
	double N1 = 1/pow(1-iEqm,V);
	double N2 = 1 - pow(1-iEqm,V);
	double N3 = (muIn/(nuEqm*V*iEqm));
	double N4 = (1-gamma)/(1 - gamma + phi);
	Neqm = N1*N2*N3*N4;
	Zeqm = phi*Neqm*iEqm/(1-gamma);
	*/

        Pars[1] = kIn; 
	V = 1/kIn; 
	Pars[2] = nuIn; 
	Pars[3] = muIn; 
	Pars[4] = sigmaIn; 
	
	Pars[5] = ratioIn; 
	Pars[9] = mIn;
	Pars[6] = Pars[9]*deltaIn;

	//printf("k:%f nu:%f mu:%f sigmaIn:%f ratio:%f delta:%f m:%f \n",Pars[1],Pars[2],Pars[3],Pars[4],Pars[5],deltaIn,Pars[9]); 
	//getc(stdin);
	


	lambda = r = FxdPars[1];
	if(Evoln == 1){
		if(sFix == 1){
			s = 1;
		}else{
			s = FxdPars[2];
		}
		phi = FxdPars[3];
		
		if(gammaFix !=1){
			gamma = FxdPars[4];
		}else{
			gamma = 0.0;
		}
		
		if(sigmaFix!=1){
			sigmaLongTerm = FxdPars[5]; //TEMP
		}else{
			sigmaLongTerm = 0.0;							
		}
		if(bFix!=1){	
			b = FxdPars[6];
		}else{
			b = 1;
		}

		
		if(ModelVsData==1){	 //Old parameters are from AllMEG3B28O28YA.dat and G3A8O28YA.in, and have output in NiceCycles13July16.dat	
			
			r = 0.1921268;
			s = 1.218437;
			phi = 7.828333;
			gamma = 0.3204130;
			sigma = 0.0;
			b = 0.1132944;

		}
			
	}else{
		phi = FxdPars[2];
		if(gammaFix==1){
			gamma = 0;
		}else{
			gamma = FxdPars[3];
		}

		if(sigmaFix==1){
			sigmaLongTerm = 0.0;
		}else{
			sigmaLongTerm = FxdPars[4]; 
		}
		
		//printf("sigmaLongTerm:%f\n",sigmaLongTerm);
		//printf("FxdPars1:%f FxdPars2:%f FxdPars3:%f FxdPars4:%f FxdPars5:%f FxdPars6:%f\n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6]);
		getc(stdin);
		

		if(ModelVsData==1){	 //These parameters are from AllMNEG3A8O28YA.dat and G3A8O28YA.in, and have output in NiceCyclesNE13July16.dat	
			/*
        		a = 0.069629;
        		sigma =1.337532;
			lambda = 1.57;
			phi = 7.071491;
			b = 0.1259;
		
			Pars[1] = 0.7357414;
			Pars[2] = 0.1849719;
			Pars[3] = 0.1923509;
			Pars[4] = 0.0;
			Pars[5] = 0.01794659;	
			Pars[9] = 100;
			Pars[6] = 100*0.05647819;
			*/
			
			//5.376435 44.39841 0.2970993 1.480939 //AllD6CNE8O42ZC.dat best params based on l-hood
			lambda = 5.376435;
			phi = 44.39841;
			gamma = 0.29870993;
			sigma = 1.480939;

		}

	}

		//printf("FxdPars1:%f FxdPars2:%f FxdPars3:%f FxdPars4:%f FxdPars5:%f FxdPars6:%f\n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6]);
		//getc(stdin);


	if(NiceVerbose==1){
		//fp = fopen("NiceCyclesNE13July16.dat","w");
		//fp = fopen("DebugCyclesEvol30Dec19.dat","w");
		//fp = fopen("OneRlznD7CE8O42ZC.dat","w");
		//printf("r:%f s:%f phi:%f gamma:%f sigma:%f b:%f \n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6]);  getc(stdin);

		//fp = fopen("OneRlznN0CE8O42ZC.dat","w"); 
		fp = fopen("OneRlznN1CNE8O42ZC.dat","w");
	  	fclose(fp);
	}

	if(Evoln!=1){
		lamV = pow(lambda,V);
		Neqm = (lamV - 1)*lambda/(V*(1+phi-gamma)*(lambda-1));
		Neqm *= 1.01*(muIn/nuIn);
		Zeqm = phi*(lamV-1)/(V*(1+phi-gamma));
		Zeqm *= (muIn/nuIn);
		nuEqm = nuIn;
		
	}else{
		
	        double cx2 = 0.0;
	        double dx2 = 0.9999;
	        double tol2 = 1e-60;
	        iEqm = zbrentd3(nuEqmfx,s,b,V,r,cx2,dx2,tol2);
	        nuEqm = (1 - r*(1-iEqm))/(r*s*pow(1-iEqm,V+1));
	        double N1 = 1/pow(1-iEqm,V);
	        double N2 = 1 - pow(1-iEqm,V);
	        double N3 = (muIn/(nuEqm*V*iEqm));
	        double N4 = (1-gamma)/(1 - gamma + phi);
		Neqm = N1*N2*N3*N4;
		Zeqm = phi*Neqm*iEqm/(1-gamma);
	}
		
	/*
	double Ntemp = 0.1;
	double Ztemp =  0.1;
	double nutemp = 5.0;

	double Nhi = 1.5; double Nlo = -1;
	double dummyStart = gsl_ran_flat(r2,Nlo,Nhi); 
	Ntemp = pow(10.0,dummyStart);

	
	double Zhi = 1.5; double Zlo = -1;
	dummyStart = gsl_ran_flat(r2,Zlo,Zhi); 
	Ztemp = pow(10.0,dummyStart);

	double nuhi = 1.0; double nulo = 0.1;
	dummyStart = gsl_ran_flat(r2,nulo,nuhi); 
	nutemp = dummyStart;
	*/
	
	//printf("r:%f s:%f phi:%f gamma:%f b:%f \n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[6]);  getc(stdin);
	


	//N = oldInfN = oldN = 10.1*Neqm; Z = oldInfZ = oldZ = 10*Zeqm; //51.874;
	//InfN = N = oldInfN = oldN = Neqm; InfZ = Z = oldInfZ = oldZ = Zeqm; //51.874;
	InfN = N = oldInfN = oldN = 1.01*Neqm; 
	InfZ = Z = oldInfZ = oldZ = Zeqm; 
	if(Evoln==1){
		Infnu = nu = oldInfnu = oldnu = nuEqm;
	}else{
		Infnu = oldInfnu = oldnu = nu = nuIn;
	}	
	
	long nStor;
	nStor = 1;
	
	for(n=1;n<=maxn;n++){
	 
	  if(Dump!=1){

		Covars[1][1] = oldN;  //PASSING HOST AND PATHOGEN DENSITIES TO ODE SOLVER
		Covars[2][1] = oldZ;
		Pars[2] = oldnu;

		MaxDays[1] = 57.0; //70.0; //1000.0;  //1000 LOOKS WRONG, I THINK IT SHOULD BE 56? NUMBER OF DAYS IN AN EPIZOOTIC...
		
		for(i=1;i<2;i++)
		for(j=1;j<2;j++)
			for(k=0;k<1000;k++){ //this k is a numeric index, not the htg param, +10 is just to ensure that there are enough...

				double RndNm = gsl_ran_gaussian(r2,Pars[4]); 
				RandNums[i][j][k] = Pars[2]*exp(RndNm); //gauss;
				//printf("i:%d j:%d k:%d Pars2:%f Pars4:%f RandNums:%f \n",i,j,k,Pars[2],Pars[4],RandNums[i][j][k]); getc(stdin);
				
			}

		
		//MaxDays[1] = 100; //TEMP

					
		  // FI2 = Distd1(Pars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead); //OLD SHITTY RK45 THAT GREG WROTE

		
		FI = DistdRK45(Pars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);  //ODE SOLVER FOR SEIR MODEL
		//printf("N:%f Z:%f FI:%f\n",N,Z,FI); 
		//getc(stdin);

		
		FractInf = zbrentd2(fx,kParam,(oldInfnu/muIn)*oldInfN,(oldInfnu/muIn)*oldInfZ,cx,dx,tol);
		//printf("oldInfN:%f oldInfZ:%f FractInf:%f\n",oldInfN,oldInfZ,FractInf); 
		//getc(stdin);
			
		/*  //obsolete, predation model is gone
		if(aFix!=1){
			predation = 1 - (2*w*a*oldN/(a*a + oldN*oldN));
		}else{
			predation = 1;
		}
		*/
		predation = 1;

		//printf("oldN:%e FI:%e predation:%e predterm:%e\n",oldN,FI,predation,(2*w*a*oldN/(a*a + oldN*oldN)));

		if(Evoln<1){
			//printf("lambda:%f phi:%f gamma:%f sigmaLongTerm:%f \n",lambda,phi,gamma,sigmaLongTerm);getc(stdin);
			//printf("oldN:%f oldZ:%f\n",oldN,oldZ);
			N = lambda*oldN*(1-FI)*predation;
			Z = phi*oldN*FI + gamma*oldZ;
			nu = oldnu;
			//printf("oldN:%f oldZ:%f FI:%f \n",oldN,oldZ,FI);

			
			
		}else{

			//N = oldN*(1-FI)*(r + r*s*oldnu*pow((1-FI),V))*predation;
			//N = oldN*(1-FI)*(r + r*s*oldnu*pow((1-FI),V))*predation;
			//printf("r:%f s:%f phi:%f gamma:%f sigmaLongTerm:%f \n",r,s,phi,gamma,sigmaLongTerm);getc(stdin);

			N = oldN*(1-FI)*(r + r*s*oldnu*pow((1-FI),V));
			Z = phi*oldN*FI + gamma*oldZ;
			//printf("V:%f r:%f s:%F FI:%f oldN:%f N:%f phi:%f gamma:%f oldZ:%f Z:%f\n",V,r,s,FI,oldN,N,phi,gamma,oldZ,Z);
		
			//Paez et al 2017
			double Numnu = oldnu*pow((1-FI),(b*V)) + s*(b*V+1)*pow(oldnu,2)*pow((1-FI),(2*b*V));
			double Denomnu = 1+s*oldnu*pow((1-FI),(b*V));

			nu = Numnu/Denomnu;
			//printf("Numnu:%f Denomnu:%f nu:%f \n",Numnu,Denomnu,nu);
		}
		N += epsilon; 

		
		

		//if(sigmaFix!=1){ //STOCHASTICITY FOR NON-BURNOUT
			
		  //double sigmaLongTerm = 0.01;
		  gauss = gsl_ran_gaussian(r2,sigmaLongTerm); //MEAN 0, SD = SIGMA
		  N *= exp(gauss);
		//}
	
		//N += eps; //IMMIGRATION

		
		if(Evoln<1){
			InfN = lambda*oldInfN*(1-FractInf)*predation;
			InfZ = phi*oldInfN*FractInf + gamma*oldInfZ;
					
		}else{
			InfN = oldInfN*(1-FractInf)*(r + r*s*oldInfnu*pow((1-FractInf),V))*predation; //BURNOUT N
			InfZ = phi*oldInfN*FractInf + gamma*oldInfZ; //BURNOUT Z

		 //Paez et al. 2017
			double Numnu = oldInfnu*pow((1-FractInf),(b*V)) + s*(b*V+1)*pow(oldInfnu,2)*pow((1-FractInf),(2*b*V));
			double Denomnu = 1+s*oldInfnu*pow((1-FractInf),(b*V));
			Infnu = Numnu/Denomnu;
		}
		
		// if(sigmaFix!=1){ //STOCHASTICITY ONLY FOR BURNOUT?
		  //sigma = 0.0;
		  gauss = sigma*gsl_ran_gaussian(r2,sigmaLongTerm);
		  InfN *= exp(gauss);
		//}

		oldN = N; oldZ = Z; oldnu = nu;
		oldInfN = InfN; oldInfZ = InfZ; oldInfnu = Infnu; 
		//oldInfnu = nu;
		

		//THIS MODEL TENDS TO CRASH A LOT
		//IF IT CRASHES, IT NEEDS TO KNOW THAT IT CRASHED.
		//THESE HORRIBLE IF STATEMENTS CHECK FOR CRASHES
		if( (isnan(oldN)) || (isnan(oldZ)) ){
			 Dump = 1;
			 oldN = -1;
			 //printf("Dump isnan...\n");
			 //getc(stdin);
		}
		if( (isinf(oldN)) || (isinf(oldZ)) ){
			 Dump = 1;
			 oldN = -1;
			//printf("Dump isinf...\n");
			//getc(stdin);
		}
		if((oldN>1e6)||(oldN<0)){
			 Dump = 1;
			 oldN = -1;
			//printf("Dump N too big or too small...\n");
			//getc(stdin);
		}
	
	        if((oldZ>1e8)||(oldZ<0)){ 
			Dump = 1;
			oldN = -1;
			//printf("Dump N too big or too small...\n");
			//getc(stdin);
		}

		
	  } //if Dump

	  if(n>nMin){
	  	NStor[nStor] = log10(N);
		//printf("nMin:%ld n:%d nStor:%d N:%f NStor[nStor]:%f \n",nMin,n,nStor,N,NStor[nStor]);
		nStor++;
	  }
	 // getc(stdin);

	//printf("n:%d N:%e Z:%e nu:%f FI:%e \n",n,N,Z,nu,FI);
	//printf("n:%d N:%e Z:%e nu:%f FI:%e InfN:%e InfZ:%e Infnu:%e FractInf:%e \n",n,N,Z,nu,FI,InfN,InfZ,Infnu,FractInf);
	//printf("n:%d InfN:%e InfZ:%e FractInf:%e \n",n,InfN,InfZ,FractInf);
        //getc(stdin);




          if(NiceVerbose==1){
	    //printf("n:%f N:%f Z:%f nu:%f FI:%f \n",n,N,Z,nu,FI);
	    //getc(stdin);

	   	//fp = fopen("NiceCyclesNE13July16.dat","a");
		//fp = fopen("DebugCyclesEvol30Dec19.dat","a");  //TEMP
		//fp = fopen("OneRlznD7CE8O42ZC.dat","a");
		printf("%e %e %e %e %e %e %e\n",n,N,Z,nu,InfN,InfZ,Infnu);

		//fp = fopen("OneRlznN0CE8O42ZC.dat","a");
		fp = fopen("OneRlznN1CNE8O42ZC.dat","a");
	      	fprintf(fp,"%e %e %e %e %e %e %e\n",n,N,Z,nu,InfN,InfZ,Infnu);
		fclose(fp);
          }

	//printf("n:%d N:%e Z:%e nu:%f FI:%e InfN:%e InfZ:%e Infnu:%e FractInf:%e \n",n,N,Z,nu,FI,InfN,InfZ,Infnu,FractInf);


	
	} //generation

	//printf("n:%d N:%e Z:%e nu:%f FI:%e InfN:%e InfZ:%e Infnu:%e FractInf:%e \n",n,N,Z,nu,FI,InfN,InfZ,Infnu,FractInf);
	//getc(stdin);

	//printf("s:%f b:%f V:%f r:%f iEqm:%f nuEqm:%f Neqm:%f Zeqm:%f\n",s,b,V,r,iEqm,nuEqm,Neqm,Zeqm); 

	//printf("end of realization...\n"); 
	//getc(stdin); 

	free_f3tensor(ModelFI,1,MaxRlzns,1,25,0,25);
	free_f3tensor(ModelS,1,MaxRlzns,1,25,0,25);
	free_d3tensor(RandNums,1,MaxRlzns,1,25,0,1000);
	free_matrix(Covars,1,25,1,25);
	free_ivector(Split,1,100);
	free_ivector(MaxDays,1,100);

	if(NiceVerbose==1) exit(1); //TEMP
	

	if(Dump!=1){

		int NumStats;
	       if(ReturnTimeStat==1){
			NumStats = ReturnTime(nStor, NStor, Amp, Prd);
		} else{
			NumStats = AmpCalc(nStor, NStor, Amp, Prd);
		}

		
		//printf("NumStats:%d\n",NumStats);
		if(NumStats>0){
	  		double AvgAmp = 0; 
	  		double AvgPrd = 0;
	  		for(i=0;i<NumStats;i++){
		
		  		AvgAmp += Amp[i];
		  		if(i>0) AvgPrd += Prd[i]; //The first period is -1
		  		//printf("i:%d AvgAmp:%f Amp:%f AvgPrd:%f Prd:%f \n",i,AvgAmp,Amp[i],AvgPrd,Prd[i]); getc(stdin);

	  		}
	  		AvgAmp /= NumStats;
	  		AvgPrd /= (NumStats-1); //The first period is -1

	  		//printf("AvgAmp:%f AvgPrd:%f\n",AvgAmp,AvgPrd);
	
	  		double LH = LHoodFunc(AvgPrd,AvgAmp);
			//printf("AvgAmp:%f AvgPrd:%f LH:%e\n",AvgAmp,AvgPrd,LH); 
			//getc(stdin);

	  		return(LH);
		}else{ //LH calc no good

	  		return(-1e100);
		}//LH calc no good
	}else{ //Dumped
		return(-1e100);
	}
	
} //End of LHood

int GetDim(int DataSet, int Plot){
	


	
       return 1;
}





int main(void){

  FILE *fp, *fp2, *fp3;
int ParamChoice2, ParamChoice3;
 //float LHoodStor[100][5000];
 int LHCount,MaxLHCount;
 float BestParamStor[25];
 int Trial;
 int Stop;


 ParamStruct ParamVals;
 double RandNums[100];
 int i, ii;
 //double test;
 double nubar, sigma;
 int Rlzn, MaxRlzns = 1;
 double AvgLHood = 0;
 double AvgLH, OldAvgLH;
 int boot, bootGood, maxBoots = 1000; //25
 int minBootGood = 25;
 int Cancel;
 size_t calls = 1;
 double Gaussian,dummy;
 double LHoodPerPop, OldLHoodPerPop, NewOldLHoodPerPop;
 double temp, NewPost,OldPost,NewPropAdj,OldPropAdj,Criterion;

 int Pop;
 int ParamChoice;
 double Itn, MaxItn;
 
 double FxdPars[30],CurrentPars[30];

 int Profile = 0;

//ParamVals.FxdPars[0] = 5; //the Plot Number
//ParamVals.FxdPars[1] =  100; // pow(10,9.2); //large k //k
ParamVals.FxdPars[1] = 0.25; //2e5; //0.25; //200000.0; //1.19602; //500; //0.1; // pow(10,-0.7); small k
double kHtg = ParamVals.FxdPars[1];

ParamVals.FxdPars[2] = 1.0; //0.67732; //0.89811; //0.75323;

 nubar = ParamVals.FxdPars[2];
 double mu = 0.6;
 //mu = 0.6;

 double InitScale = 1.25, VarScale;



 ParamVals.FxdPars[3] = mu; //0.39268; //0.38678; //M:0.525; //mu

 ParamVals.FxdPars[4] = 0.0; //0.0; //0.3; //M:0.00398; //sigma
 //sigma = ParamVals.FxdPars[4];
 ParamVals.FxdPars[5] = 0.0; //0.01303; //0.71616; //1.17653; //2.75; //-2; //RATIO KEEP THIS AT -2 FOR OTVOS/SHEP AND FOR MOREAU  //in fact, this has been at -2.8 for O/S and M
 double ratio = 0.013;
 
ParamVals.FxdPars[6] = -2; //0.06417; //0.07794; //1.0/12.0; //(10.0/12.0); //M:15; //tau, speed of kill
 ParamVals.FxdPars[7] = 0.1; //0.01; //0.25; //h

 ParamVals.FxdPars[8] = 1; //corrln interval

 ParamVals.FxdPars[9] = 55; //15.0; //Number of E stages
 ParamVals.FxdPars[10] = 0.5; //Obsvn Error Var, actually delta
 ParamVals.FxdPars[11] = DataSet; //DataSet 1 is GM5 (Obsolete),8 is OtvosCon (Obsolte), 2 is GM5 beta-bnml (CURRENT), 12 is Otvos bnml, 13 is Otvos/Shep Con beta bnml (DON'T USE), 
//14 is Moreau con beta bnml, 15 is Otvos/Shep fixed NO 13!!

 size_t dim = GetDim(ParamVals.FxdPars[11],1)/ParamVals.FxdPars[8];
 int NoGo;


/*
 ParamVals.FxdPars[12] =  -3;//FI1
 ParamVals.FxdPars[13] =  -3; //FI2
 ParamVals.FxdPars[14] =  -3; //FI3
 ParamVals.FxdPars[15] =  -3; //FI4
 ParamVals.FxdPars[16] =  -3; //FI5
*/

 ParamVals.FxdPars[12] =  -4;//FI1
 ParamVals.FxdPars[13] =  -3.8; //FI2
 ParamVals.FxdPars[14] =  -3.3; //FI3
 ParamVals.FxdPars[15] =  -4; //FI4
 ParamVals.FxdPars[16] =  -3.9; //FI5
 ParamVals.FxdPars[17] =  -3; //FI5

 ParamVals.FxdPars[20] = 1.88480; //4.05330; //3.84827; //4.25; //gamma - THIS IS FOR THE PSEUDO LHOOD or the beta binomial inverse variance

for(i=0;i<=30;i++)
	FxdPars[i] = 0.0;
 for(i=0;i<=20;i++){
	FxdPars[i] = ParamVals.FxdPars[i];
	CurrentPars[i] = FxdPars[i];
 }

int NumPops[25],ParamCount[25];
//int DataSet = ParamVals.FxdPars[11];


NumPops[1] = NumPops[2] = 5; //Woods and Elk
NumPops[1] = NumPops[2] = 8; //Woods and Elk TEMP

NumPops[13] = 5; //Otvos/Shep control pseudo/bnml/beta bnml
NumPops[14] = 3; //Moreau control pseudo/bnml/beta bnml
NumPops[15] = 6; //Otvos-Shep fixed up 
NumPops[16] = 7; //this is Otvos All Con Bnml

/*
ParamCount[2] = 8; //Woods and Elk, beta binomial, obsvn fixed
ParamCount[13] = 11; //Otvos/Shep Control beta binomial, no ratio, no obsvn error
ParamCount[14] = 9; //Otvos/Shep Control beta binomial, no ratio, no obsvn error
ParamCount[15] = 12; //Otvos/Shep Control beta binomial, no ratio, no obsvn error, fixed up ICK!!
ParamCount[16] = 7; ////this is Otvos All Con Bnml
*/

//ParamCount[2] = 5; //lambda, r, phi, gamma, SEIR params


	if(Evoln<1){
		ParamCount[2] = 4;
	}else{
		//ParamCount[2] = 5; //lambda, r, phi, gamma, SEIR params
		ParamCount[2] = 6; //lambda, r, phi, gamma, SEIR params

	}


 int PC, PC2;
 int ParamNum[30]; //These are only for looping through the parameters

 for(PC=1;PC<=16;PC++)
 	ParamNum[PC] = PC; //GM bnml, no stage

//GM beta binomial
if(DataSet==2){
	ParamNum[7] = 9;
	ParamNum[8] = 20;
}

/*For Otvos/Shep Con beta binomial - NO ratio, no obsvn error */
if(DataSet==13){
 ParamNum[5] = 6;
 ParamNum[6] = 12;
 ParamNum[7] = 13;
 ParamNum[8] = 14;
 ParamNum[9] = 15;
 ParamNum[10] = 16;
 ParamNum[11] = 20;
}

/* Moreau Con beta binomial - NO ratio, no obsvn error */
if(DataSet==14){
 ParamNum[5] = 6;
 ParamNum[6] = 12;
 ParamNum[7] = 13;
 ParamNum[8] = 14;
 ParamNum[9] = 20;
}

/*For Otvos/Shep Con beta binomial - NO ratio, no obsvn error */
//if((ParamVals.FxdPars[11]==15)||(ParamVals.FxdPars[11]==16)){ //15 is the data set number
if(ParamVals.FxdPars[11]>=15){ //15 is the data set number
 ParamNum[5] = 6;
 ParamNum[6] = 9;
 ParamNum[7] = 20;
}

double res, err, oldres;
double AvgL[50], Var[50];
     
       double xl[300];
       double xu[300];
     
       //const gsl_rng_type *T; //don't seem to be used...
       //gsl_rng *r, *pr;

float lambda;

//T = gsl_rng_default;
//r = gsl_rng_alloc (T);  //These 2 really aren't used, screw'em

Verbose = 0;
Verbose2 = 0;
//calls = 5;//50; //this is now data-set specific, below

float BestLHood = -1e50;
float BestParam;
int BestCrashOut;
double Variance,VarBest;
int MaxTrial = 1000;


const gsl_rng_type *T2;
//gsl_rng *r2;
//gsl_rng_env_setup(); //you do this down below, AFTER seed, so don't do it here



   long seed;
       srand((unsigned) time(NULL));
      // seed = -rand();
    seed = time(NULL)*(int)getpid(); 
    //seed = -1;
    
   

       
//Monte function here G() //THIS IS OLD
       //gsl_monte_function G = { &LHood, dim, &ParamVals }; //MOVED
      
       for (i=0; i<300;i++) {xl[i]=0;xu[i]=1;}
       gsl_rng_env_setup ();

       gsl_rng_default_seed = seed;
       //pr = gsl_rng_alloc(gsl_rng_mt19937);

T2 = gsl_rng_default;
r2 = gsl_rng_alloc(T2);

if(Profile==1){
	MaxTrial = 1;
	MaxItn = 2;
}

int Accept[30], ParamItn[30];
int PropDist[30], PriorDist[30];

//Reading in SEIR parameters
	//double nuIn[5000],muIn[5000],ratioIn[5000],deltaIn[5000],mIn[5000];
	double *kIn, *nuIn, *muIn, *sigmaIn, *ratioIn, *deltaIn, *mIn; 
	

	const NumSets = 10000000;
	kIn = dvector(0,NumSets);
	nuIn = dvector(0,NumSets);
	muIn = dvector(0,NumSets);
	sigmaIn = dvector(0,NumSets);
	ratioIn = dvector(0,NumSets);
	deltaIn = dvector(0,NumSets);
	mIn = dvector(0,NumSets);

	double Error;
	
	//fp = fopen("Param4EB.in","r"); //GM, but check to make sure rows are ok!
	//fp = fopen("PrmOAOC5FThin1e2.in","r"); //GM, but check to make sure rows are ok!
	//fp = fopen("PrmGx5O5JThin1e2.in","r"); //GM, but check to make sure rows are ok!
       //fp = fopen("PrmGx5O8CThin1e2.in","r"); //GM, but check to make sure rows are ok!
 	
	//fp = fopen("PrmGx5O2FThin1e2.in","r"); //GM, but check to make sure rows are ok!

	//fp = fopen("G3A8O28YA.in","r"); //GM, but check to make sure rows are ok!

/* 17 December 2019 */

	fp = fopen("GA8O42ZC.in","r"); //GM //Priors on everything
	//fp = fopen("GA8O15ZC.in","r"); //Priors on nubar and k only
	//fp = fopen("GA8P0ZC.in","r");
	//printf("file opened...\n");


	i = 0;
	
	long SEIRParamNum = 0;
	double junkIn; //represents sigma
	while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&junkIn,&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error)!=EOF){
	//while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&junkIn,&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error)!=EOF){
		
		if(kHtg>1e5) kIn[SEIRParamNum] = 2e6;
				
	/* Old
	//while(fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum])!=EOF){ //mFix
        //while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&Error)!=EOF){
		//printf("SEIRParamNum:%d sigmaIn:%f\n",SEIRParamNum,sigmaIn[SEIRParamNum]); 
	//mIn[SEIRParamNum] = 100.0; //Ditched as of 17 Dec 2019

	//while(fscanf(fp,"%f %f %f %f %f %f %f\n",&k,&nubar,&mu,&ratio,&delta,&m,&Error)!=EOF){
	*/
	
		//sigmaIn[SEIRParamNum] = 0.0;
		
		//printf("%d %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],mIn[SEIRParamNum]);

		//getc(stdin);
		SEIRParamNum++;	
	}
	fclose(fp);
	SEIRParamNum--;	
	//printf("SEIRParamNum:%d\n",SEIRParamNum);
	//getc(stdin);
	//exit(1);

/****************** Priors and Proposals *********************/
double PropMean[30], PropVar[30], PriorMean[30], PriorVar[30], PriorUnif[30], InitFlag[30];
double MaxParam[30];
double PropMean2[30], PropVar2[30], MixFract[30], MixRand;
int ThinStop;

//printf("DataSet:%d\n",DataSet);

 if(DataSet==2){



  ThinStop = 10;
  //ThinStop = 1; //TEMP

  MaxItn = 1e9;
  calls = 100; //25; //50;
 
  for(i=1;i<30;i++){
    PriorMean[i] = -1; PriorVar[i] = 1e3;    
    PropVar[i] = 0.4; //0.4
    PriorUnif[i] = 1;
    PriorDist[i] = 1; //Uniform
    PropDist[i] = 3; //Log-normal
    Accept[i] = 0;
    ParamItn[i] = 0;
    PriorVar[i] = 10000;
    PropMean2[i] = -1; PropVar2[i] = 2.0;
    MixFract[i] = 0.01;
    InitFlag[i] = 1;
  }
  PriorUnif[3] = 0;

 if(Evoln==1){
//Here are the proposals for r
//Uninformative

   PropMean[1] = 0.2; //log(40.0); //log(25.0); //r
   PropVar[1] =  0.1; //0.25; //0.25; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PropDist[1] = 2; //normal-log

   PriorDist[1] = 1;
   PriorVar[1] = 1.0; //r has to be less than 1
   MaxParam[1] = 1.0;



//Informative r

/*
//SF  1:21 (0.2, 1.9), rF  0:22 (0.05, 0.51) //page 5, PaezEtAl17SI (Appendix B)
   PriorMean[1] = log(0.22); //log(40.0); //log(25.0); //r, see NatSel16, SI, page 52
   double VarTemp = log(0.22) + log(0.46/2);
   PriorVar[1] =  log(0.46/2.0); //0.25; //0.25; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PriorDist[1] = 2; //normal-log
   MaxParam[1] = 1.0;
*/


// Uninformative r AGAIN?  WHY? Ugh
/*
   PriorDist[1] = 1;
   PriorVar[1] = 1.0; //r has to be less than 1
   MaxParam[1] = 1.0;
*/
	
 
//Here are the proposals for s WHICH HAS BEEN SCALED AWAY AND THEREFORE NO LONGER EXISTS
//Uninformative	

  PropDist[2] = 2; //s - slope
  PropMean[2] = 0.5; //-3; //-2.5; //log(0.15); //-4; //MGH: 1, cuz MGF:2.5; //use 5 for uinformative...//3; //-1; //log(FxdPars[2]);  //note that we are here nuking the log-scale for nubar
  PropVar[2] = 1.0; //0.75; //0.5; //0.25; //0.1; //0.025; //0.05; //0.25; //0.1; //0.1; //0.5; //1.5; //0.5; //1.5; 
  
  PriorDist[2] = 1;
  PriorVar[2] = 10000; //a guess...
  MaxParam[2] = 10000;


/*
//SF  1:21 (0.2, 1.9), rF  0:22 (0.05, 0.51) //page 5, PaezEtAl17SI (Appendix B)
  PropDist[2] = 2; //s - slope
  PropMean[2] = log(1.21);  
  PropVar[2] = log(1.7/2.0);
  MaxParam[2] = 1e4;
  
  PriorDist[2] = 2; //s - slope
  PriorMean[2] = log(1.21);  
  PriorVar[2] = log(1.7/2.0);
*/

 

 

//phi
//Uninformative Prior

  PropDist[3] =  2; //phi
  PropMean[3] = log(10); //-0.89; //-1.25; //-4; //log(0.3); //log(0.41); //SOK prior's only (MGH) 2.0; //MGF: log(20); //uninf: log(100); //log(0.41);
  PropVar[3] = 0.5; //0.25; //1.0; //0.75; //0.5; //0.25; //0.05; //0.5; //1.25; //0.5; //1.25; 
  PriorVar[3] = 1e5; //phi bigger than 1000 is hard to imagine BUT WITH OVERDISPERSION YOU CAN RUN INTO BIG TROUBLE
  PriorDist[3] = 1;
  MaxParam[3] = 1e5; 

//Informative Prior
  //Priors
/*
  PriorDist[3] = 2;
  PriorMean[3] = 2.02;
  PriorVar[3] = 0.48; //phi bigger than 100 is hard to imagine
  //Proposals
  PropDist[3] =  2; //phi
  PropMean[3] = 2.02; //-0.89; //-1.25; //-4; //log(0.3); //log(0.41); //SOK prior's only (MGH) 2.0; //MGF: log(20); //uninf: log(100); //log(0.41);
  PropVar[3] = 0.48; //1.0; //0.75; //0.5; //0.25; //0.05; //0.5; //1.25; //0.5; //1.25; 
  MaxParam[3] = 1000.0;
*/




//4, which is a, which is the density at which predation is highest - THIS IS  b in DWYER ET AL. 2004! This is Paez et al. 2017
//NO, now it's gamma
  //Here are the proposals
  PropMean[4] = -1.5; //log(0.25); //log(0.3); //gamma 
  PropVar[4] = 0.5; //0.35; //0.2; //0.1; //P thru S; GMOFL: 1.0; 
  PropDist[4] = 1.0; //
  MaxParam[4] = 100.0; 
  
 
  //Uninformative gamma   
 PriorDist[4] = 1;
 PriorVar[4] = 1.0; //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
 MaxParam[4] = 1.0;
 // PriorUnif[4] = 1; //no longer in use?

//Informative on predation density
  //PriorDist[4] = 2;
  //PriorMean[4] = -2.732795114;
  //PriorVar[4] = 0.479460599;
   

     //5, which is sigma
 PriorDist[5] = 1;
 PriorVar[5] = 10.0; //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
 MaxParam[5] = 10.0;
 PriorUnif[5] = 1; //no longer in use?

     //6, which is heritability b
//Uninformative b

PriorDist[6] = 1;
PriorVar[6] = 1.0; //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
MaxParam[6] = 1.0;
//PriorUnif[6] = 1; //no longer in use?



//Informative b 
/*
//b = 0.1259;
//#2.5%      10%       25%        50%        60%        70%        75%        80%      97.5%
//#0.0013   0.0193    0.0639    0.1259    0.1542    0.1879    0.2078    0.2322   0.4771
  PriorDist[6] = 2;
  PriorMean[6] = log(0.1259);
//(log(0.1259) - log(0.4771))/2 = -0.6661191; 
  PriorVar[6] = 0.6661191;
  MaxParam[6] = 1.0;
*/

 

  } else { //No Evoln below

  
   //Uninformative - lambda

   //Here are the proposals
   PropMean[1] = 2.0; //1.5; //3.22; //log(40.0); //log(25.0); //lambda
   PropVar[1] =  0.65; //0.25; //0.25; //0.25; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PropDist[1] = 2; //normal-log
   PriorVar[1] = log(1e4); //max reasonable no. eggs/egg mass
   PriorDist[1] = 1;
   MaxParam[1] = log(1e4);

  

    //Informative - lambda - Look in ElkEtAl96DensityNumbers
   //Here are the proposals
/*
   PropMean[1] = -2.5; //log(3.0); //3.22; //log(40.0); //log(25.0); //lambda
   PropVar[1] =  1.0; //log(1.75); //0.25; //0.25; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PropDist[1] = 2; //normal-log
   PriorDist[1] = 2;
   PriorMean[1] = 1.87; //See ElkEtAl96DensityNumbers.xls
   PriorVar[1] = 0.53; //
   MaxParam[1] = log(1e4);
*/
			
 


//Informative Prior on phi
  //Priors
/*
  PriorDist[2] = 2;
  PriorMean[2] = 2.02;
  PriorVar[2] = 10.0; 
  MaxParam[2] = log(1e4);


  //Proposals
  PropDist[2] =  2; //phi
  PropMean[2] = 2.02; //-0.89; //-1.25; //-4; //log(0.3); //log(0.41); //SOK prior's only (MGH) 2.0; //MGF: log(20); //uninf: log(100); //log(0.41);
  PropVar[2] = 0.48; //1.0; //0.75; //0.5; //0.25; //0.05; //0.5; //1.25; //0.5; //1.25; 
*/
 //Uninformative prior on phi
   PropMean[2] = 2.0; //1.5; //3.22; //log(40.0); //log(25.0); //lambda
   PropVar[2] =  0.65; //0.25; //0.25; //0.25; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PropDist[2] = 2; //normal-log
   PriorDist[2] = 1;
   PriorVar[2] = 100; //max reasonable no. eggs/egg mass
   MaxParam[2] = 1e5;



//3, which is gamma
  //Here are the proposals
  PropMean[3] = -1.5; //log(0.25); //log(0.3); //gamma 
  PropVar[3] = 0.5; //0.35; //0.2; //0.1; //P thru S; GMOFL: 1.0; 
  PropDist[3] = 1.0; //
  MaxParam[3] = 1000.0; 
  
  
     //Uninformative
   PriorDist[3] = 1;
   PriorVar[3] = 1.0; //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
   MaxParam[3] = 1.0;
 // PriorUnif[4] = 1; //no longer in use?

  //Informative - old, for "a"
  /*
  PriorDist[3] = 2;
  PriorMean[3] = -2.732795114;
  PriorVar[3] = 0.479460599;
  */
   

     //4, which is sigma
 PriorDist[4] = 1;
 PriorVar[4] = 100.0; //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
 MaxParam[4] = 1000.0;
 PriorUnif[4] = 1; //no longer in use?
  




  } //No Evoln

 }//GM priors



 /* OS from here on */
 if(DataSet>=15){


  ThinStop = 1000;
  MaxItn = 1e9;
  calls = 25;
  for(i=1;i<30;i++){
    PriorMean[i] = -1; PriorVar[i] = 1e3;    
    PropVar[i] = 0.4; //0.4
    PriorUnif[i] = 1;
    PriorDist[i] = 1; //Uniform
    PropDist[i] = 3; //Log-normal
    Accept[i] = 0;
    ParamItn[i] = 0;
    PriorVar[i] = 1000;
    InitFlag[i] = 1;
  }
  //PriorUnif[3] = 0;
//PriorMean[3] = 0.41; PriorVar[3] = 1.0; //this is Emma's mu, more or less, on linear scale

  PropMean[1] = 5.0; //1.5; //lambda
  PropVar[1] = 0.25; //0.1; //G: 0.1; //C: 0.25 A: 1.0?  B:0.5; //0.25; //0.4;
  PropDist[1] = 2;
 
  //printf("PropMean1:%f\n",PropMean[1]);
  //PropVar[1] = 0.75; //0.75
  PriorVar[1] = 1000; //for uniform case
  

  PropDist[2] = 2; //r normal
  PropMean[2] = -2; //log(10); //C: log(2.0) A? 2.0; //log(FxdPars[2]);  
  PropVar[2] = 0.25; //0.1; //G: 0.1; //C:0.5 B: 1.0? A: 0.5;//0.4
  PriorVar[2] = 1000;
  PriorDist[2] = 1;
  //PropMean[2] = 0.5;
  // PropVar[2] = 0.5;
  
  //Informative
  /*
  PriorDist[2] = 2; //normal
  PriorMean[2] = log(0.5); //From Dwyer '90
  PriorVar[2] = 0.4; //From Dwyer '90 (sort of)
  */

 

  //PriorUnif[3] = 0;
  //Uninformative
  PriorVar[3] = 1e5; //Uniform Prior 
  PriorDist[3] = 1; //uniform
  //Informative
  /*
  PriorDist[3] = 2; //log-normal 
  PriorMean[3] = log(0.41);  //Emma's data 
  PriorVar[3] = 0.01; //Emma's data 
  */
  
  PropDist[3] = 2; //Normal
  PropMean[3] = -0.5; //0.2; //2.0; //B:log(10)? //A: 2.0; //log(0.41); //Normal
  PropVar[3] = 0.5; //0.25; //0.1; 
   
  PropMean[4] = 0.0; //gamma
  PropVar[4] = 0.25; //0.1;
  PropDist[4] = 2;
  PriorVar[4] = 24;

  PropMean[5] = log(FxdPars[5]);  //SEIR from her on down...
  PriorVar[5] = 1000;
 
  PropDist[6] = 2;
  PropMean[6] = log(0.02); //A?:-2;
  PropVar[6] =  0.1; //G: 0.1; //A?:0.76;

  PriorDist[6] = 1;
  PriorVar[6] = 10.0; 


  PropMean[9] = 1000; //10; //25; //25; //28; //15; //This is x in the neg binom, because this is the number of exposed stages, so it has to be an integer
  PropVar[9] = 0.995; //0.5; //0.5; //This is p in the neg binom
  PropDist[9] = 7; //neg binom
  PriorVar[9] = 125; //75;
  PriorDist[9] = 1;
  
  //Informative: (from Arietta)
  //PriorDist[9] = 2; //Normal on a log scale
  //PriorMean[9] = 3.305926445; //3.969886639; //3.305926455;
 // PriorVar[9] = 0.048050982; //0.074451582; //0.048050982;
 

  // PropMean[17] = 0.01;
  
  PropDist[20] = 2.0; //A? Nothing?
  PropMean[20] = log(1.0); //A: 1.0; //0.7 log(FxdPars[20]); 
  PropVar[20] = 0.1; //G: 0.1; //A?: 0.4;
  PriorVar[20] = 1000;

 }//OS priors and post's



	char testbuff[128];
	char param_3[100];   
	
//L - long term, G - gypsy moth, M4E is the source of the samples
//N - No evolution, G - gypsy moth M4E is the source
//D - DFTM
//M4E2 is a trial from June 2014
//M4EX means param set 4E, and X means updated LHFunc and data
//M4EY like M4EX, but with amps and periods calculated as times above the mean, instead of discretized second deriviatives
//LG2 has informative prior on phi
//LG3 has informative priors on both phi and lambda (probably only for no evolution)
//L2 has maxn = 200 instead of 100, ditto for N2; L3 is maxn = 1000; L4 is 500
//x used to be an update of the Woods and Elk data, but now let's say it's how peaks are calculated - as a bend in the data
//y is then peaks calculated using the return time

//MEGA did not have Skaller et al, MEGB DOES have Skaller et al.  Baseline is 25 calls, MEG100B8O23ZA is 100 calls 
//GC fits w and sigma
//B2 has an informative prior on heritability
//G4 has informative priors on everything
/******/
//17 Dec 2019
//E for Evolution, N for No Evolution, we're only doing gypsy moth, so never mind the G4
//LS for Line Search, M for MCMC
//8042YC is the params from fitting to the Woods and Elk Data 
//DLS is Dwyer lab models: Dwyer et al 2000 for no evoln, Elderd et al 2008 for evoln
//D2LS is uninformative priors
//D3 has sigmaFix != 1, and uses the best-fit sigma values
//D3B like D3, but with minBootgood = 50, maxBoots 100 in Line Search; minBootGood = maxBoots = 25 in MCMC; VarScale = 1.2
//D3C like D3B, but with minBootGood = maxBoots = 50 in MCMC; VarScale = 1.2
//D3DLS has minBootGood = 25, maxBoots 25 in Line Search, sigmaFix != 1; D4D like D3D, but sigmaFix = 1
//D4E like D4D, except maxn = 400 instead of 300
//
//On 30 January 2019, code was changed so that it calculated LH's over the entire posterior from the Epiz model fits, instead of booting.
//Previous D series is superseded.
//D is still Dwyer-lab-models, and NE is still No Evolution  
//D3 is Sigma NOT fixed; D4 is Sigma FIXED at 0.01; D4B is epsilon = 1e-5; D4C is like D4B, but uses SEIRParamNum instead of bootGood to calc AvgLH;
//D4D like D4C, but epsilon = 0.0; D5C like D4C, except MaxLSItn = 3; D6C like D5C, except MaxLSItn = 2 (no, it's also 3, whoops).
//D7C like D6C, except sigmaFix != 0; D8C like D6C, except b!=1; D9C like D8C, except also sigmaFix!=1

//Previous D series unintentionally used a bogus informative prior on gamma, not to mention on b.  Darn it.
//D2 has sigmaFix = 1, bFix = 1, uninformative prior on gamma; D3 like D2, but sigmaFix !=1; D4 has sigmaFix = 1, bFix != 1; 
//D5 like D4, but informative prior on b; D6 like D5, but sigmaFix == 0; D7 has informatives on r, s and b, but sigmaFix = 1; D8 like D7, but sigmaFix = 0
//Due to an error, D7 turned back into sigmaFix = 0 (I think, at least for Evoln = 1)
//testD7 does not re-do LH calc for old param set when sigmaFix == 1; D8 has no informative priors
//E as first letter is after I realized that s can be scaled away...UGH. F like E, except that epsilon = 0 (because gamma = 0 appeared to be the best fit)
//G has gamma fixed at gamma = 0.0, epsilon = 1e-5; H calculates lhood 2x, as in the stoch models, but even if stoch=0
//7 is actually r uninformative, b informative, NOT r informative; 8 is now all uninformative evoln model
//H8CE8O15ZC has no stoch, but gamma not zero; H9CE8O15ZC has no stoch, gamma = 0
//I has epsilon = 1e-6; J has epsilon 1e-4; H9B is a re-do, based on a possible error in NewDLPCA.R
//L has revised periods and amps, uses epsilon = 1e-5, b uninformative, epsilon = 1e-5 
//L8 has no stoch, gamma = 0; L7 has stoch, gamma not zero, NE; L6 is evoln WITH stoch, and with gamma = 0
//L9 like L8, but using all of the Line Search output, because L8 did not mix...
//L0 is L6, but keeping a bigger sample (fraction = 1.0) of the PCA output; L1 is L6, but with smaller sample (fraction = 0.25) of the PCA output
//M series is with an informative prior on phi; M0 has no stoch, gamma = 0
//N series uses return times for data, with the addition of Skaller, for which there are no return times, sadly
//N0 is evoln with no stoch, N1 is no evoln with stoch, and with gamma > 0; N2 is N0 but using the entire PCA output to make proposals
//N3 is N0 (or N2) but with informative prior on b, uses entire PCA output; N4 is N3 but with wider bounds on LS (N3 is basically dead...)
//N5 is evoln model with stochasticity, but with gamma = 0, but best fit turned out to be with very low stochasticity, so LS only 
//N6 is with SEIR parameters from GA8O15ZC, for which the priors were only on nubar and k, with no stoch; N7 like N6 but with stochasticity
//N8 is no evoln
//N2e4 is N2 with eps = 1e-4; N2e6 is with eps = 1e-6; N2e55 is with epsilon = 5e-5 (because 1e-4 sucks!); N2e56 is epsilon = 5e-6
//N2e7 is N2 with eps = 1e-7
//N8e4 should really have been N1e4; N1e6 is N1 with epsilon = 1e-6
//N1e5CNE8O42ZC is a repeat of N1CNE8O42ZC, but with a re-done line search with broader bounds
//N4e5 is a repeat of N4, but with a re-done line search; N4e56 has epsilon = 5e-6
//N2B is a repeat of N2, but with line search limits based on MCMC from N4e5
//N7 is N6 with stoch; by accident, N7CLSE8O42ZC is actually N7CLSE8O15ZC.  N7 already existed, so switched to N8
//N6CE8O15ZC has nans in the PrdAmpStor, so N6B is a re-do

	
	char *test[] =
	  {
	    "N6BCE8O15ZC", /* MCMC */
	    "AccGM3C", /* NOT IN USE - this used to record the acceptance rate */
	    "N6BCLSE8O15ZC", /* Line Search Output */
	    "N6BCE8O15ZC", /* Priors and Seeds for MCMC only */
	  };

    char testbuff2[128];
    char buffer0[128],buffer1[128],buffer2[128],buffer3[128];
    char bufferB[128];
    char fname;
    int pid;
    pid = getpid();
//printf("pid:%d\n",pid);

	char *strFileType = ".dat";
	strcpy(buffer0, test[0]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer0, bufferB);
	strcat(buffer0,strFileType); 
	test[0] =  buffer0;

	strcpy(buffer1, test[1]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer1, bufferB);
	strcat(buffer1,strFileType); 
	test[1] =  buffer1;

	strcpy(buffer2, test[2]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer2, bufferB);
	strcat(buffer2,strFileType); 
	test[2] =  buffer2;

	char *strFileType2 = ".par";
	strcpy(buffer3, test[3]);
     	sprintf(bufferB, "%d", pid);
     	strcat(buffer3, bufferB);
	strcat(buffer3,strFileType2); 
	test[3] =  buffer3;



		if(LineSearch==0){
			 			
			//Write proposal distributions etc to a file
			strcpy(testbuff,test[3]);
		  	fp2 = fopen(testbuff,"w");
			fprintf(fp2,"******************************************************************************\n");
			fprintf(fp2,"Seed:%d Boots:%d minBootGood:%d MAXN:%d ReturnTime:%d Evoln:%d\n",seed,maxBoots,minBootGood,MAXN,ReturnTimeStat,Evoln);
			fprintf(fp2,"DataSet:%d Seed:%d sigmaFix:%d bFix:%d InitScale:%f kHtg:%f\n\n\n",DataSet,seed,sigmaFix,bFix,InitScale,kHtg);
			fprintf(fp2,"PC   ParamNum PropDist   PropMean       PropVar         PropMean2            PropVar2	PriorDist   PriorMean       PriorVar    \n"); 
			for(PC=1;PC<=ParamCount[DataSet];PC++){
				ParamChoice = ParamNum[PC];
				fprintf(fp2,"%d     \t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\n",PC,ParamNum[PC],PropDist[ParamChoice],PropMean[ParamChoice],PropVar[ParamChoice],PropMean2[ParamChoice],PropVar2[ParamChoice],MixFract[ParamChoice],PriorDist[ParamChoice],PriorMean[ParamChoice],PriorVar[ParamChoice]);
			}
			fprintf(fp2,"*******************************************************************************\n");
			fclose(fp2);
		}


//*************************** LineSearch

if(LineSearch==1){

printf("Starting line search...sampSize:%d sigmaFix:%d\n",sampSize,sigmaFix);
	
   double Upper[30], Lower[30], Jump[30];
   int PC; //Parameter Choice
   double LStemp;
   int LSItn, MaxLSItn = 3;
   int NumParams;
   
   
 if(Evoln!=1){

	/* //k large
	Lower[1] = 0.001; Upper[1] = 0.022; Jump[1] = 0.005; //lambda
	//Lower[2] = 1.6; Upper[2] = 2.5; Jump[2] = 0.25; //phi
	Lower[2] = 0.25; Upper[2] = 7.0; Jump[2] = 1.0; //phi
	Lower[3] = -4; Upper[3] = 0; Jump[3] = 0.5; //gamma
	*/

	Lower[1] = 0.1; Upper[1] = 4.0; Jump[1] = 0.75; //lambda
	//Lower[1] = 0.5; Upper[1] = 2.5; Jump[1] = 0.5; //lambda

 	//Lower[2] = -1; Upper[2] = 2.5; Jump[2] = 0.25; //phi
	Lower[2] = 0; Upper[2] = 7; Jump[2] = 1; //phi

	//Lower[3] = -2.5; Upper[3] = 0.0; Jump[3] = 0.5; //gamma
	//Lower[3] = -2.5; Upper[3] = -1.0; Jump[3] = 0.5; //gamma
	Lower[3] = -2.5; Upper[3] = -0.5; Jump[3] = 0.5; //gamma

	//Lower[4] = 0.25; Upper[4] = 0.75; Jump[4] = 0.25; //sigmaLongTerm
	Lower[4] = 0.4; Upper[4] = 1.0; Jump[4] = 0.1; //sigmaLongTerm
	//Lower[4] = 0.35; Upper[4] = 0.65; Jump[4] = 0.1; //sigmaLongTerm

	
	NumParams = 3; //2;
	if(sigmaFix==1){
		 FxdPars[5] = 0.0;
	}else{
		NumParams++;
	}
		

 }else{ //with evoln

	
	////SF  1:21 (0.2, 1.9), rF  0:22 (0.05, 0.51) //page 5, PaezEtAl17SI (Appendix B)
	//Lower[1] = -0.5; Upper[1] = -0.1; Jump[1] = 0.1; //r 
	//Lower[1] = -2; Upper[1] = -0.1; Jump[1] = 0.2; //r informative 
	Lower[1] = -0.6; Upper[1] = -0.01; Jump[1] = 0.05; //r //uninformative

	//Lower[1] = log(0.22) - 0.15 ; Upper[1] = log(0.22) + 0.15; Jump[1] = 0.075; //r informative

 	Lower[2] = 1; Upper[2] = 8.0; Jump[2] = 1.0; //s uninformative [actually, s is dead]
	//Lower[2] = log(1.21) - 0.15; Upper[2] = log(1.21) + 0.15; Jump[2] = 0.075; //s

 	
	//Lower[3] = 2.0; Upper[3] = 6.0; Jump[3] = 1.0; //phi - uninformative
	//Lower[3] = 0; Upper[3] = 4.0; Jump[3] = 1.0; //phi - informative
	Lower[3] = 2.0; Upper[3] = 6.0; Jump[3] = 1.0; //phi - with b informative

	//Lower[4] = -2.5; Upper[4] = -0.5; Jump[4] = 0.5; //gamma
	//Lower[4] = -3.25; Upper[4] = -2.75; Jump[4] = 0.5; //gamma
	//Lower[4] = -4.5; Upper[4] = -1.5; Jump[4] = 0.5; //gamma
	//Lower[4] = -2.5; Upper[4] = -1.5; Jump[4] = 0.25; //gamma
	//Lower[4] = -2.9; Upper[4] = -2.5; Jump[4] = 0.1; //gamma
	Lower[4] = -5; Upper[4] = -1; Jump[4] = 2; //gamma 


	//Lower[5] = -3; Upper[5] = 2.0; Jump[5] = 1; //sigmaLongTerm
	Lower[5] = -9; Upper[5] = 1; Jump[5] = 1; //sigmaLongTerm

	//Lower[6] = -3.5; Upper[6] = -1.5; Jump[6] = 0.5; //b with informative prior
	//Lower[6] = -7; Upper[6] = -1; Jump[6] = 1.0; //b with uninformative prior
	Lower[6] = -7.0; Upper[6] = -3; Jump[6] = 1.0; //b with uninformative prior

	NumParams = 3; //4; //5 //3;
	/*
	if(sigmaFix==1){
		 SigmaLongTerm = 0.0; //Really? 
	}
	else{
		NumParams++;
	}
	*/
	if(sFix != 1){
		NumParams++;
	}else{
		FxdPars[2] = 1;
	}
	if(sigmaFix!=1){
		 NumParams++;
	}else{
		FxdPars[5] = 0.0;
	}
	
	if(bFix==1){
		 FxdPars[6] = 1.0;
	}else{
		 NumParams = 6;
	}


 }


 //printf("NumParams:%d\n",NumParams);// getc(stdin);

   maxBoots = 500; //used to be 25
   minBootGood = 50;
   long boot;



   bootGood = 0; 
   Cancel = 0; 
   int LSflag = 0;
   double OldPost = -1e10, NewPost;
   double paramtemp;



   Cancel = 1;
   do{
	//printf("starting do-until loop.\n");
   			for(PC=1;PC<=NumParams;PC++){
				
				NoGo = 0;
				//printf("PC:%d\n",PC);
				if(Evoln==1){
					if((sFix==1)&&(PC==2)) NoGo = 1;
					if((gammaFix==1)&&(PC==4)) NoGo = 1;
					if((sigmaFix==1)&&(PC==5)) NoGo = 1;
					if((bFix==1)&&(PC==6)) NoGo = 1;
					
				} else{
					if((gammaFix==1)&&(PC==3)) NoGo = 1;
					if((sigmaFix==1)&&(PC==4)) NoGo = 1;
				}
				
				//printf("PC:%d NoGo:%d \n",PC,NoGo);
			
				if(NoGo!=1){
				 
					double randtemp = gsl_rng_uniform(r2); 
					 
					LStemp = Lower[PC] + (Upper[PC] - Lower[PC])*randtemp;
					FxdPars[PC] = exp(LStemp); 
					//printf("PC:%d randtemp:%f Lower:%f Upper:%f FxdPars[PC]:%f \n",PC,randtemp,Lower[PC],Upper[PC],FxdPars[PC]); 
					//getc(stdin);
				}

   			} //PC loop


			AvgLH = 0.0;
			bootGood = 0;
			//FxdPars[1] = 0.1432; FxdPars[2] = 14.13; FxdPars[3] = 2.2689; FxdPars[4] = 0.12811; //D4DE8O42ZC good fit
			//FxdPars[1] = 0.1943335; FxdPars[2] = 1.214359; FxdPars[3] = 7.574822; FxdPars[4] = 0.3279860; FxdPars[6] = 0.1090024; //testAllD7CE8O42ZC
			//0.1943335 1.214359 7.574822 0.3279860 0.01303 0.1090024
			
			//maxBoots = 1;
   			//for(boot=0;boot<maxBoots;boot++){
			for(boot=0;boot<SEIRParamNum;boot++){
			//for(boot=0;boot<50;boot++){


				//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
				int Index = boot;   
												
				res = LHood(FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);
				
				//printf("Boot loop boot:%d Index:%d res:%e sigmaIn:%f\n",boot,Index,res,sigmaIn[Index]);
				//exit(1);

				 //If res is a really bad likelihood, it will effectively be adding zero to AvgLH
				if(res>-1e10){
		 			bootGood++;
		 			AvgLH += exp(res);
		 			//printf("Initial boot:%d bootGood:%d Index:%d res:%e AvgLH:%e\n",boot,bootGood,Index,res,log(AvgLH/(boot+1))); //getc(stdin);
					//if(res>-16) exit(1);
				}
				

   			} //for(boot

			
			if((isnan(AvgLH)==1)||(bootGood<1)){
				Cancel = 1; 
				//printf("Didn't work...bootGood:%d AvgLH:%e\n",bootGood,AvgLH); getc(stdin);
			}else{
				Cancel = 0;
				AvgLH /= SEIRParamNum;
			}
			
			//printf("bootGood:%d Cancel:%d AvgLH:%f\n",bootGood,Cancel,log(AvgLH)); 
			//getc(stdin);
   }while(Cancel==1);

   for(LSItn=1;LSItn<=MaxLSItn;LSItn++){
     for(PC=1;PC<=NumParams;PC++){

      paramtemp = FxdPars[PC]; 

      
      for(LStemp=Lower[PC];LStemp<=Upper[PC];LStemp+=Jump[PC]){
				
			FxdPars[PC] = exp(LStemp);
			NoGo = 0;
			if(Evoln==1){
					if((sFix==1)&&(PC==2)) NoGo = 1;
					if((gammaFix==1)&&(PC==4)) NoGo = 1;
					if((sigmaFix==1)&&(PC==5)) NoGo = 1;
					if((bFix==1)&&(PC==6)) NoGo = 1;
					
			} else{
					if((gammaFix==1)&&(PC==3)) NoGo = 1;
					if((sigmaFix==1)&&(PC==4)) NoGo = 1;
			}

		if(NoGo!=1){
			//printf("PC:%d LStemp:%f FxdPars:%f \n",PC,LStemp,FxdPars[PC]);

			AvgLH = 0.0;
			Cancel = 0; 
			bootGood = 0;
			//FxdPars[1] = 0.2051665; FxdPars[2] = 1.185351; FxdPars[3] = 20.69203; FxdPars[4] = 0.3115007; FxdPars[6] = 0.1114308;
			//0.2051655	1.185351	20.69203	0.3115007	0.1114308
			//printf("Fxd1:%f Fxd2:%f Fxd3:%f Fxd4:%f FxdPars[6]:%f\n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[6]);
			
   			//for(boot=0;boot<maxBoots;boot++){
			for(boot=0;boot<SEIRParamNum;boot++){

				//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
				int Index = boot;  
   				res = LHood(FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);
				//printf("Next boot:%d Index:%d res:%e \n",boot,Index,res);

				
				if(res>-1e10){
		 			bootGood++; 
		 			AvgLH += exp(res);
		 			//printf("Ready to exit boot:%d bootGood:%d Index:%d res:%e AvgLH:%e\n",boot,bootGood,Index,res,log(AvgLH/bootGood));
					//exit(1);
				}

				

   			} //for(boot


			if((isnan(fabs(AvgLH))==1)||(bootGood<1)){
				Cancel = 1;
			       //printf("Cancelled...\n"); 
			}else{
				AvgLH /= SEIRParamNum;
				
				AvgLH = log(AvgLH);
				
				NewPost = AvgLH;
				
			}
			

			//printf("Cancel:%d NewPost:%f\n",Cancel,NewPost);
	

			if(Cancel<1){


			  int PC2;
			  
			  
			  for(PC2=1;PC2<=NumParams;PC2++){
			    NoGo = 0;
			    if(Evoln==1){
					if((sFix==1)&&(PC2==2)) NoGo = 1;
					if((gammaFix==1)&&(PC2==4)) NoGo = 1;
					if((sigmaFix==1)&&(PC2==5)) NoGo = 1;
					if((bFix==1)&&(PC2==6)) NoGo = 1;
					//printf("PC2:%d NoGo:%d \n",PC2,NoGo);
					
				} else{
					if((gammaFix==1)&&(PC2==3)) NoGo = 1;
					if((sigmaFix==1)&&(PC2==4)) NoGo = 1;
			     }


			    if(NoGo!=1){
		 		
							   
		         	float tmp2;
		         	switch(PriorDist[PC2]){
					case 1 :
			  			if(log(FxdPars[PC2])>PriorVar[PC2]){
			       				Stop = 1;
			    				NewPost = -1e50;
			  			} else {
			    				tmp2 =  gsl_ran_flat_pdf((FxdPars[PC2]),0,PriorVar[PC2]);
						}
			  			break;
					case 2 :tmp2 = gsl_ran_gaussian_pdf(log(FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // a slightly bogus kind of normal, really sort log normal
				       		break;
					case 3 :tmp2 = gsl_ran_lognormal_pdf(FxdPars[PC2],PriorMean[PC2],PriorVar[PC2]);
				       		break;
	                		case 4 :tmp2 = gsl_ran_gaussian_pdf((FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // more like a real normal...
	   		               		break;
	                		default : printf("No prior distribution specified. Case:%d PriorDist:%d Bailing.\n",PC2,PriorDist[PC2]); getc(stdin); exit(1); break;
		      	 	} //End of switch

					

				//printf("PC2:%d FxdPars:%f tmp2:%e \n",PC2,FxdPars[PC2],tmp2);// getc(stdin);
		 
			 	if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
									
  					Cancel = 1;
					NewPost = -1e10;
			 	} else{
			    		NewPost += log(tmp2);
					//printf("PC2:%d FxdPars:%f tmp2:%e NewPost:%f\n",PC2,FxdPars[PC2],tmp2,NewPost); getc(stdin);
				
			 	}
			 	if((isnan(NewPost)!=0)||(isinf(NewPost)!=0)){
				
					Cancel = 1;
					NewPost = -1e10;
			 	}
				
						


		      	    } //NoGo
			  } //PC2

			  //printf("After Prior adjustment...NewPost:%f.  Hit return.\n",NewPost); getc(stdin);

			  

				if(LSflag<1){ //First time through
					OldPost = NewPost;
					OldAvgLH = AvgLH;
					FxdPars[PC] = exp(LStemp);
					LSflag = 1; 
				
				} else{
				
					if(NewPost>OldPost){

		
						//printf("PC:%d old:%f FxdPars:%f NewPost:%f OldPost:%f\n",PC,paramtemp,FxdPars[PC],NewPost,OldPost); 

						OldPost = NewPost;
						OldAvgLH = AvgLH;

						paramtemp = exp(LStemp);

						
					}//Post>BestPost
				
				} //LSFlag

				//printf("LSItn:%d PC:%d LStemp:%6.3f Fxd1:%f Fxd2:%f Fxd3:%f Fxd4:%f Fxd5:%f Fxd6:%f OldPost:%f NewPost:%f\n",LSItn,PC,LStemp,FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6],OldPost,NewPost);
				
			} else{ 
				//LSflag = 0;
				
			        //printf("Cancelled...\n");

			} //Cancel

			

       } //LStemp
       FxdPars[PC] = paramtemp;
     } // if NoGo!=1
     } //for PC
	

   } //for LSItn

   //for(PC=1;PC<=4;PC++) printf("PC:%d FxdPars:%f \n",PC,FxdPars[PC]);
  // getc(stdin);

   	if((isnan(OldPost)==0)&&(isinf(OldPost)==0)){
		strcpy(testbuff,test[2]);
		fp = fopen(testbuff,"a");
		fprintf(fp,"%e %e %e %e %e %e %e %e\n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6],(OldAvgLH),OldPost);	
		fclose(fp);
		
	}
	
   exit(1);

} //if(LineSearch)


//***************** End of Linesearch

//*************************** NiceVerbose Plotting Function, which calls lhood to get one realization plotted
//  AS OF 12 JULY 2016, THIS IS REALLY MORE OF A MODEL VS DATA CODE THAN A DIC CALCULATION.  DIC IS OBSOLETE, SO DIC CALCULATION IS KIND OF POINTLESS...


if(NiceVerbose==1){
	printf("NiceVerbose...\n");

	double Dbar;

	
	char testbuffA[128];
    	char bufferA[128];

	//char *test2 = "MeanN4G2Ax5O16XA";
	char *test2 = "DLModelVsDataGA8O28YA";

	char *test3;
	char *strFileType = ".med";
	strcpy(bufferA, test2);
	strcat(bufferA,strFileType); 
	test3 =  bufferA;

	strcpy(testbuffA,test3);
	//fp = fopen(testbuffA,"r");

	printf("just about to open DLModelVsData file\n");

	//fp = fopen("DLModelVsDataGA8O28YA.med","r");
	//fp = fopen("N0CEAddByHandGA8O42ZC.med","r");

	
	//printf("Evoln:%d file opened successfully...\n",Evoln);

	//r, s, phi, gamma, sigma, b
	//actually reading in - r, phi, b
	double mpsrf, Last;
	if(Evoln==1){
		printf("here we are...\n");
		printf("FxdPars:%f\n",FxdPars[1]);
		double junk;
	//&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&junkIn,&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error
		//Obsolete: &kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&Error
		//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&FxdPars[1],&FxdPars[2],&FxdPars[3],&FxdPars[4],&kIn[0],&nuIn[0],&muIn[0],&junk,&ratioIn[0],&deltaIn[0]);
		//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&FxdPars[1],&FxdPars[3],&FxdPars[6],&FxdPars[4],&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&junkIn,&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error);

		//mIn[0] = 100;
		//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&FxdPars[1],&FxdPars[2],&FxdPars[3],&FxdPars[4],&kIn[0],&nuIn[0],&muIn[0],&ratioIn[0],&deltaIn[0],&mIn[0],&Dbar,&mpsrf,&Last);
		//printf("%f\n",FxdPars[1]);
		
		//N0CEAddByHandGA8O42ZC.med
		FxdPars[1] = 0.5950393; //r
		FxdPars[3] = 133.8007;	//phi
		FxdPars[6] = 0.01766681; //b

		//L9CE8O42ZC.med
		/*
		FxdPars[1] = 0.6651729;
		FxdPars[3] = 21.57414;	
		FxdPars[6] = 0.04470024;
		*/
		//0.6651729	21.57414	0.04470024

		
	}else{
		//N1CNE8O42ZC
				 

		FxdPars[1] = 5.257213; //lambda
		FxdPars[2] = 29.42765;	//phi
		FxdPars[3] = 0.1904699; //gamma
		FxdPars[4] = 1.656735; //sigma	
	//5.257213	29.42765	0.1904699	1.656735
		
		//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&FxdPars[1],&FxdPars[3],&FxdPars[4],&kIn[0],&nuIn[0],&muIn[0],&ratioIn[0],&deltaIn[0],&mIn[0],&Dbar,&mpsrf,&Last);
	}
	//fclose(fp);

	//GA8O42ZC SEIR parameters
	kIn[0] = 0.7362992; nuIn[0] = 0.4609148; muIn[0] = 0.4103536; //0.9859989 stoch? 
	ratioIn[0] = 0.01794859; deltaIn[0] = 0.06424359; mIn[0] = 27.7875; //3.342443 error	
	res = LHood(FxdPars,kIn[0],nuIn[0],muIn[0],sigmaIn[0],ratioIn[0],deltaIn[0],mIn[0]); //calculates lhood and generates output

	exit(1);


} //if(0)

//*************************** End of NiceVerbose


///******************* WAIC calculation **********************////


if(WAIC==1){ //Turn on WAIC calculation

printf("starting WAIC...\n"); getc(stdin);

			if(PrdAmpStor==1){
				GlobalAvgAmp = 0.0;
				GlobalAvgPrd = 0.0; 
				//fp2 = fopen("PrdAmpStorH9CE8O15ZC.dat","w");
				//fp2 = fopen("PrdAmpStorNoSigD6CNE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN2CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN0CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN1CNE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN2CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZCtmp.dat","w");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZCtmp.dat","w");
				//fp2 = fopen("PrdAmpStorN6CE8O42ZCtmp.dat","w");
				//fp2 = fopen("PrdAmpStorN8CNE8O15ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN1CNE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN2CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN4e5CE8O42ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN6CE8O15ZC.dat","w");
				//fp2 = fopen("PrdAmpStorN6BCE8O15ZC.dat","w");
				fp2 = fopen("PrdAmpStorN1e5CNE8O42ZCnoSig.dat","a");


				fclose(fp2);
				
			}


	
//mFix = 0;
maxBoots = 500;
//printf("Doing WAIC... maxBoots:%d\n",maxBoots); getc(stdin);

	FILE *fp4;
	int samp, sampSize = 5;


	
/*
	int weekNum;
	for(Pop=0;Pop<30;Pop++){
		for(weekNum=0;weekNum<100;weekNum++){

			PointwiseLH[Pop][weekNum] = -1.0e10;
			PointwiseVar[Pop][weekNum] = -1.0e10;
			PointwiseMean[Pop][weekNum] = -1.0e10;

		}


	}
*/

	double *kIn, *nuIn, *muIn, *sigmaIn, *ratioIn, *deltaIn, *mIn, *Error;
	double *rIn, *sIn, *phiIn, *wIn, *sigIn, *bIn, *gammaIn;
	const NumSets = 10000000;
	kIn = dvector(0,NumSets);
	nuIn = dvector(0,NumSets);
	muIn = dvector(0,NumSets);
	sigmaIn  = dvector(0,NumSets);
	ratioIn = dvector(0,NumSets);
	deltaIn = dvector(0,NumSets);
	mIn = dvector(0,NumSets);
	Error = dvector(0,NumSets);

	

	rIn = dvector(0,NumSets);
	sIn = dvector(0,NumSets);
	phiIn = dvector(0,NumSets);
	gammaIn = dvector(0,NumSets);
	wIn = dvector(0,NumSets);
	
	sigIn = dvector(0,NumSets);
	bIn = dvector(0,NumSets);

	double dummy;


/*************** Change file name here ******************/

	
	//printf("just about to read in the SEIR parameters file...\n");
	//fp = fopen("G3A8O28YA.in","r");
	fp = fopen("GA8O42ZC.in","r");
	//fp = fopen("GA8P0ZC.in","r");
	//fp = fopen("GA8O15ZC.in","r");

	//printf("done opening...\n");

	i = 0;
	long MaxSEIRParamNum, SEIRParamNum = 0;
			
		if(kHtg>1e5){
			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf \n",&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				kIn[SEIRParamNum] = 2e5;

				//printf("Small k %d %f %f %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],mIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			
				SEIRParamNum++;	
			}
		}else{
			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf \n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				//printf("Small k %d %f %f %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],mIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			
				SEIRParamNum++;
			 }
		}	

	
	
		
	fclose(fp);

	MaxSEIRParamNum = SEIRParamNum;

	//exit(1);

	//fp = fopen("MNEG48O28YA.in","r");
	//fp = fopen("MNEG48O28YA.in","r");
	//fp = fopen("MNE8O42YC.in","r");
	//fp = fopen("ME8O42YC.in","r");
	//fp = fopen("D3NE8O42ZC.in","r");
	//AFTER 17 Dec 2019
	//fp = fopen("testD7CE8O42ZC.in","r");
	//fp = fopen("D6CNE8O42ZC.in","r");
	//fp = fopen("D6CNE8P0ZC.in","r");
	//fp = fopen("D7CE8O15ZC.in","r");
	//fp = fopen("D7CNE8O15ZC.in","r");
	//fp = fopen("G7CE8O42ZC.in","r");
	//fp = fopen("H8CE8O42ZC.in","r");
	//fp = fopen("H8CE8O42ZC.in","r");
	//fp = fopen("I9CE8O42ZC.in","r");
	//fp = fopen("I7CNE8O42ZC.in","r");
	//fp = fopen("I9BCE8O42ZC.in","r");
	//fp = fopen("J9CE8O42ZC.in","r");
	//fp = fopen("J7CNE8O42ZC.in","r");
	//fp = fopen("L6CNE8O42ZC.in","r");
	//fp = fopen("L9CE8O42ZC.in","r"); //200 samples
	//fp = fopen("LB9CE8O42ZC.in","r"); //321 samples
	//fp = fopen("N0CE8O42ZCtmp.in","r"); //221 samples
	//fp = fopen("N1CNE8O42ZCtmp.in","r"); //200 samples
	//fp = fopen("N0CE8O42ZC.in","r"); 
	//fp = fopen("N2CE8O42ZC.in","r"); 
	//fp = fopen("N1CNE8O42ZC.in","r");
	//fp = fopen("N4CE8O42ZCtmp.in","r"); 
	//fp = fopen("N4CE8O42ZCtmp.in","r"); 
	//fp = fopen("N6CE8O15ZCtmp.in","r");
	//fp = fopen("N8CNE8O15ZC.in","r");
	//fp = fopen("N4CE8O42ZC.in","r"); 
	//fp = fopen("N2e6CE8O42ZC.in","r");
	//fp = fopen("N2e7CE8O42ZCtmp.in","r"); 
	//fp = fopen("N2e56CE8O42ZCtmp.in","r"); 
	//fp = fopen("N2e55CE8O42ZCtmp.in","r"); 
	//fp = fopen("N2e6CE8O42ZC.in","r"); 
	//fp = fopen("N2e4CE8O42ZC.in","r");
	//fp = fopen("N2e55CE8O42ZC.in","r");
	//fp = fopen("N2e56CE8O42ZC.in","r"); 
	//fp = fopen("N8e4CNE8O42ZCtmp.in","r");
	//fp = fopen("N8e4CNE8O42ZCtmp.in","r");
	//fp = fopen("N1e5CNE8O42ZCtmp.in","r");
	//fp = fopen("N2e5CE8O42ZCtmp.in","r"); 
	//fp = fopen("N2e5CE8O42ZCtmpB.in","r"); 
	//fp = fopen("N2e5CE8O42ZCtmpC.in","r"); 
	fp = fopen("N1e5CNE8O42ZC.in","r"); // **********
	//fp = fopen("N2e5CE8O42ZC.in","r");
	//fp = fopen("N2e7CE8O42ZC.in","r"); 
	//fp = fopen("N1e6CNE8O42ZC.in","r");
	//fp = fopen("N8e4CNE8O42ZC.in","r"); 
	//fp = fopen("N4e5CE8O42ZCtmp.in","r");
	//fp = fopen("N1e3CNE8O42ZCtmp.in","r"); 
	//fp = fopen("N4CE8O42ZCtmp.in","r"); 
	//fp = fopen("N4CE8O42ZCtmpB.in","r");
	//fp = fopen("N4CE8O42ZCtmpC.in","r"); 
	//fp = fopen("N1e3CNE8O42ZC.in","r"); 
	//fp = fopen("N4CE8O42ZC.in","r"); 
	//fp = fopen("N4e5CE8O42ZCtmp.in","r"); 
	//fp = fopen("N4e5CE8O42ZC.in","r"); 
	//fp = fopen("N2BCE8O42ZCtmpA.in","r");
	//fp = fopen("N2BCE8O42ZC.in","r"); 
	//fp = fopen("N6CE8O15ZC.in","r"); 
	//fp = fopen("N6BCE8O15ZCtmpB.in","r"); 
	//fp = fopen("N6BCE8O15ZCtmpC.in","r"); 
	//fp = fopen("N2e4CE8O42ZC.in","r"); 
	//fp = fopen("N6BCE8O15ZC.in","r");
  	
	i = 0;
	long MaxDLParamNum, DLParamNum = 0;
	
	if(Evoln==1){
			//while(fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&rIn[DLParamNum],&sIn[DLParamNum],&phiIn[DLParamNum],&gammaIn[DLParamNum],&sigIn[DLParamNum],&bIn[DLParamNum])!=EOF){
			while(fscanf(fp,"%lf %lf %lf \n",&rIn[DLParamNum],&phiIn[DLParamNum],&bIn[DLParamNum])!=EOF){
	
				if(sigmaFix==1) sigIn[DLParamNum] = 0.0;
				if(sFix==1) sIn[DLParamNum] = 1.0;
				if(gammaFix==1) gammaIn[DLParamNum] = 0.0;
			
				//printf("%d %f %f %f %f %f %f\n",DLParamNum,rIn[DLParamNum],sIn[DLParamNum],phiIn[DLParamNum],gammaIn[DLParamNum],sigIn[DLParamNum],bIn[DLParamNum]);
				printf("With Evoln %d %f %f %f \n",DLParamNum,rIn[DLParamNum],phiIn[DLParamNum],bIn[DLParamNum]);
				//getc(stdin);
				
				//getc(stdin);
				DLParamNum++;	
			}
	} else{

			
			while(fscanf(fp,"%lf %lf %lf %lf\n",&rIn[DLParamNum],&phiIn[DLParamNum],&gammaIn[DLParamNum],&sigIn[DLParamNum])!=EOF){

				sIn[DLParamNum] = 0.0;
				if(sigmaFix==1) sigIn[DLParamNum] = 0.0;

				printf("%d %f %f %f %f \n",DLParamNum,rIn[DLParamNum],phiIn[DLParamNum],gammaIn[DLParamNum],sigIn[DLParamNum]);

			
				
				//if(DLParamNum==0){
					//printf("DLParamNum:%d r:%f phi:%f gamma:%f sigma:%f\n",DLParamNum,rIn[DLParamNum],phiIn[DLParamNum],gammaIn[DLParamNum],sigIn[DLParamNum]); 
					 //getc(stdin);
				//}
				
				DLParamNum++;	
			}
	} //Evoln

	
		
	fclose(fp);

	//printf("DLParamNum:%d \n",DLParamNum);
	//getc(stdin); 
	//exit(1);


	MaxDLParamNum = DLParamNum;


	for(i=0;i<100;i++){
		PopWiseMean[i] = 0.0;
		PopWiseM2Log[i] = 0.0;
		PopWiseMeanLog[i] = 0.0;
		PopWiseVarLog[i] = 0.0; //unnecessary?
		PopWiseNumSamp = 0;
	}

//DLParamNum = 0;


for(DLParamNum = 0;DLParamNum<MaxDLParamNum;DLParamNum++){


	FxdPars[1] = rIn[DLParamNum]; 
	
	if(Evoln==1){
		FxdPars[2] = sIn[DLParamNum]; 
		FxdPars[3] = phiIn[DLParamNum];
		FxdPars[4] = gammaIn[DLParamNum];
		FxdPars[5] = sigIn[DLParamNum];
		FxdPars[6] = bIn[DLParamNum];
	}else{
		FxdPars[2] = phiIn[DLParamNum];
		FxdPars[3] = gammaIn[DLParamNum];
		FxdPars[4] = sigIn[DLParamNum];

	}

	//printf("DLParamNum:%d sig:%f FxdPars4:%f \n",DLParamNum,sigmaIn[DLParamNum],FxdPars[4]);
	

//Calculate the likelihood	
	AvgLH = 0.0;
	OldAvgLH = 0.0;
	Cancel = 0; 
	bootGood = 0;
	maxBoots = 50;

		
			if(PrdAmpStor==1){
				GlobalAvgAmp = 0.0;
				GlobalAvgPrd = 0.0; 
			}
			//for(boot=0;boot<maxBoots;boot++){
			for(boot=0;boot<SEIRParamNum;boot++){
			
				//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
				int Index = boot; 

				//printf("boot:%d\n",boot);

				res = LHood(FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);
				//printf("Inside boot loop Fxd1:%f Fxd2:%f Fxd3:%f Fxd4:%f Fxd6:%f res:%f\n",FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[6],res);
				//getc(stdin);

				PopWiseNumSamp++; //do it here, because otherwise you will increase it even though the parameter set gave only nan
								
				if(res>-1e10){
		 			 
					//if(DLParamNum==1) printf("res:%e\n",res);
		 			OldAvgLH += exp(res);
					//PopWiseNumSamp++; //NO LONGER TRUE: do it here, because otherwise you will increase it even though the parameter set gave only nan
					// getc(stdin);
					bootGood++;
					//printf("DLParamNum:%d boot:%d bootGood:%d res:%f OldAvgLH:%e PopWiseNumSamp:%d AvgAmp:%f AvgPrd:%f\n",DLParamNum,boot,bootGood,res,log(OldAvgLH/(boot+1)),PopWiseNumSamp,GlobalAvgAmp/bootGood,GlobalAvgPrd/bootGood);
							 			
				}

			
   			} //for(boot
			if(PrdAmpStor==1){
				//fp2 = fopen("PrdAmpStorH9CE8O15ZC.dat","a");
				//fp2 = fopen("PrdAmpStorNoSigD6CNE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN0CE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN1CNE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN2CE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZCtmp.dat","a");
				//fp2 = fopen("PrdAmpStorN6CE8O15ZCtmp.dat","a");
				//fp2 = fopen("PrdAmpStorN8CNE8O15ZC.dat","a");
				///fp2 = fopen("PrdAmpStorN4CE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN1CNE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN2CE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN4CE8O42ZCtmpC.dat","a");
				//fp2 = fopen("PrdAmpStorN4e5CE8O42ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN6CE8O15ZC.dat","a");
				//fp2 = fopen("PrdAmpStorN6BCE8O15ZC.dat","a");
				fp2 = fopen("PrdAmpStorN1e5CNE8O42ZCnoSig.dat","a");

				
				fprintf(fp2,"%f %f \n",GlobalAvgAmp/bootGood,GlobalAvgPrd/bootGood);
				printf("printing Prd and Amp...Prd:%f Amp:%f\n",GlobalAvgAmp/bootGood,GlobalAvgPrd/bootGood);
				fclose(fp2);
			}


			if(isnan(AvgLH)==1){
				//if(DLParamNum==1) printf("Ditch...\n");
				Cancel = 1; 
			}else{


				//if(DLParamNum==1) printf("OldAvgLH:%e\n",OldAvgLH);
				
				//OldAvgLH /= bootGood;  //used to do only a sub-sample
				//printf("SEIRParamNum:%d\n",SEIRParamNum); getc(stdin);
				OldAvgLH /= SEIRParamNum;  //currently doing the entire sample


				//if(DLParamNum==1) printf("divide bootGood OldAvgLH:%e\n",OldAvgLH);

				OldAvgLH = log(OldAvgLH);

				//if(DLParamNum==1) printf("log OldAvgLH:%e\n",OldAvgLH);
			
				OldPost = (OldAvgLH);
				OldLHoodPerPop = OldAvgLH;

				

				
			}
			//printf("DLParamNum:%d FXd1:%f Fxd2:%f FXd3:%f Fxd4:%f Fxd5:%f Fxd6:%f bootGood:%d OldAvgLH:%f\n",DLParamNum,FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6],bootGood,OldAvgLH);
			//getc(stdin);


 	} //DLParamNum

	fp = fopen("DLWAIC2020.dat","a");

	double pWAIC2B = 0.0;
	double TotPopWiseMeanLog = 0.0;
	printf("PopWiseNumSamp:%d\n",PopWiseNumSamp);
	for(samp=0;samp<sampSize;samp++){
				//printf("samp:%d PopWiseMean:%e PopWiseMeanLog:%e PopWiseM2Log:%e\n",samp,PopWiseMean,PopWiseMeanLog[samp],PopWiseM2Log[samp]); 

				printf("samp:%d PopWiseMeanLog:%f PopWiseM2Log:%e\n",samp,PopWiseMeanLog[samp],PopWiseM2Log[samp]);//  getc(stdin);
				pWAIC2B += (PopWiseM2Log[samp]/PopWiseNumSamp) - (PopWiseMeanLog[samp]/PopWiseNumSamp)*(PopWiseMeanLog[samp]/PopWiseNumSamp);
				TotPopWiseMeanLog += log(PopWiseMean[Pop]/PopWiseNumSamp);
	}
	double AltWAIC = -2*(TotPopWiseMeanLog) + 2*(pWAIC2B);

	printf("TotPopWiseMeanLog:%f pWAIC2B:%f AltWAIC:%f \n",TotPopWiseMeanLog,pWAIC2B,AltWAIC);

	fprintf(fp,"PA_N1e5CNE8O42ZCnoSig TotPopWiseMeanLog:%f pWAIC2B:%f AltWAIC:%f Sample Size:%d\n",TotPopWiseMeanLog,pWAIC2B,AltWAIC,MaxDLParamNum);
	fclose(fp);

	//printf("hit return to continue...\n");
	//getc(stdin);
	
	free_dvector(rIn,0,NumSets);
	free_dvector(sIn,0,NumSets);
	free_dvector(phiIn,0,NumSets);
	free_dvector(gammaIn,0,NumSets);
	free_dvector(wIn,0,NumSets);
	free_dvector(sigIn,0,NumSets);
	free_dvector(bIn,0,NumSets);
	

	exit(1);

}

///******************* End of WAIC Calculation ***********************////




/*************** MCMC *****************/

///////////////////////////////   Reading in PCA Results ///////////////////////////


minBootGood = 50; maxBoots = 50; //No longer relevant, I believe
printf("Starting MCMC, no need to hit return. minBootGood:%d maxBoots:%d \n",minBootGood,maxBoots);

        double Coefficients[30][30];
	double Scale[30],Center[30],SD[30];
	int j;
	int ParCnt2 = ParamCount[DataSet];

	if(Evoln==1){
		if(sFix==1) ParCnt2--;
		if(gammaFix==1) ParCnt2--;
		if(sigmaFix==1) ParCnt2--;
		if(bFix==1) ParCnt2--;
	} else{
		if(gammaFix==1) ParCnt2--;
		if(sigmaFix==1) ParCnt2--;
	}

	//fp = fopen("D7CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("D6CLSNE8P0ZCRotations.txt", "r");
	//fp = fopen("D7CLSE8O15ZCRotations.txt", "r");
	//fp = fopen("D7CLSNE8O15ZCRotations.txt", "r");
	//fp = fopen("D8CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("E8CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("G7CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("G8CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("H8CLSE8O15ZCRotations.txt", "r");
	//fp = fopen("H9CLSE8O15ZCRotations.txt", "r");
	//fp = fopen("I9CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("I7CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("J9CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("J7CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("L8CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("L7CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("L8CLSE8O42ZCRotations2.txt", "r");
	//fp = fopen("L6CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("L0CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("L1CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("M0CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N0CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N1CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("N2CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N3CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N6CLSE8O15ZCRotations.txt", "r");
	//fp = fopen("N7CLSE8O15ZCRotations.txt", "r");
	//fp = fopen("N8CLSNE8O15ZCRotations.txt", "r");
	//fp = fopen("N2e4CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e6CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e6CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e56CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e7CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N8e4CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("N1e6CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("N1e5CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N4e6CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N4e56CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N1e3CLSNE8O42ZCRotations.txt", "r");
	//fp = fopen("N4CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N2BCLSE8O42ZCRotations.txt", "r");
	//fp = fopen("N7CLSE8O15ZCRotations.txt", "r");
	fp = fopen("N6BCLSE8O15ZCRotations.txt", "r");



	//printf("Rotations...\n");
	for(i = 0;i<ParCnt2;i++){
	  for(j = 0;j<ParCnt2;j++){
	    	fscanf(fp, "%lf\t", &Coefficients[i][j]);
		//printf("%lf \n", Coefficients[i][j]); 
	  }
	  fscanf(fp,"\n");
        //printf("\n");
	}
	fclose(fp);
	//getc(stdin); 

	
	//fp = fopen("D7CLSE8O42ZCScale.txt", "r");
	//fp = fopen("D6CLSNE8P0ZZCScale.txt", "r"); //spelling problem here...
	//fp = fopen("D7CLSE8O15ZCScale.txt", "r");
	//fp = fopen("D7CLSNE8O15ZCScale.txt", "r");
	//fp = fopen("D8CLSE8O42ZCScale.txt", "r");
	//fp = fopen("E8CLSE8O42ZCScale.txt", "r");
	//fp = fopen("G7CLSE8O42ZCScale.txt", "r");
	//fp = fopen("G8CLSE8O42ZCScale.txt", "r");
	//fp = fopen("H8CLSE8O15ZCScale.txt", "r");
	//fp = fopen("H9CLSE8O42ZCScale.txt", "r");
	//fp = fopen("I9CLSE8O42ZCScale.txt", "r");
	//fp = fopen("I7CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("L8CLSE8O42ZCScale.txt", "r");
	//fp = fopen("L7CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("L8CLSE8O42ZCScale2.txt", "r");
	//fp = fopen("L6CLSE8O42ZCScale.txt", "r");
	//fp = fopen("L0CLSE8O42ZCScale.txt", "r");
	//fp = fopen("L1CLSE8O42ZCScale.txt", "r");
	//fp = fopen("M0CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N0CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N1CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("N2CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N3CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N6CLSE8O15ZCScale.txt", "r");
	//fp = fopen("N7CLSE8O15ZCScale.txt", "r");
	//fp = fopen("N8CLSNE8O15ZCScale.txt", "r");
	//fp = fopen("N2e4CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N2e6CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N2e56CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N2e7CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N8e4CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("N1e6CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("N1e5CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N4e6CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N4e56CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N1e3CLSNE8O42ZCScale.txt", "r");
	//fp = fopen("N4CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCScale.txt", "r");
	//fp = fopen("N2BCLSE8O42ZCScale.txt", "r");
	//fp = fopen("N7CLSE8O15ZCScale.txt", "r");
	fp = fopen("N6BCLSE8O15ZCScale.txt", "r");

	//printf("Scales...\n");
	for (i=0; i<ParCnt2; i++){
	  fscanf(fp, "%lf ", &Scale[i]);
	  //printf("i:%d Scale:%e\n", i,Scale[i]); //TEMP
	 }
	fclose(fp);
	//getc(stdin); 

	//fp = fopen("D7CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("D6CLSNE8P0ZCCenter.txt", "r");
	//fp = fopen("D7CLSE8O15ZCCenter.txt", "r");
	//fp = fopen("D8CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("E8CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("G7CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("G8CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("H8CLSE8O15ZCCenter.txt", "r");
	//fp = fopen("H9CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("I9CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("I7CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("J9CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("L8CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("L7CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("L8CLSE8O42ZCCenter2.txt", "r");
	//fp = fopen("L6CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("L0CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("L1CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("M0CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N0CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N1CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("N2CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N3CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N6CLSE8O15ZCCenter.txt", "r");
	//fp = fopen("N8CLSNE8O15ZCCenter.txt", "r");
	//fp = fopen("N2e4CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N2e6CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N2e56CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N2e7CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N8e4CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("N1e6CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("N1e5CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N4e6CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N4e56CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N1e3CLSNE8O42ZCCenter.txt", "r");
	//fp = fopen("N4CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N2BCLSE8O42ZCCenter.txt", "r");
	//fp = fopen("N7CLSE8O15ZCCenter.txt", "r");
	fp = fopen("N6BCLSE8O15ZCCenter.txt", "r");


	//printf("Center...\n");
	for (i=0; i<ParCnt2; i++){
	   fscanf(fp,"%lf ", &Center[i]); //TEMP
	  // printf("i:%d Center:%lf\n", i,Center[i]);
	}
	fclose(fp);
	//getc(stdin); 

	//fp = fopen("D7CLSE8O42ZCsd.txt", "r");
	//fp = fopen("D6CLSNE8P0ZCsd.txt", "r");
	//fp = fopen("D7CLSE8O15ZCsd.txt", "r");
	//fp = fopen("D8CLSE8O42ZCsd.txt", "r");
	//fp = fopen("E8CLSE8O42ZCsd.txt", "r");
	//fp = fopen("G7CLSE8O42ZCsd.txt", "r");
	//fp = fopen("G8CLSE8O42ZCsd.txt", "r");
	//fp = fopen("H8CLSE8O15ZCsd.txt", "r");
	//fp = fopen("H9CLSE8O42ZCsd.txt", "r");
	//fp = fopen("I9CLSE8O42ZCsd.txt", "r");
	//fp = fopen("I7CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("J9CLSE8O42ZCsd.txt", "r");
	//fp = fopen("L8CLSE8O42ZCsd.txt", "r");
	//fp = fopen("L7CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("L8CLSE8O42ZCsd2.txt", "r");
	//fp = fopen("L6CLSE8O42ZCsd.txt", "r");
	//fp = fopen("L0CLSE8O42ZCsd.txt", "r");
	//fp = fopen("L1CLSE8O42ZCsd.txt", "r");
	//fp = fopen("M0CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N0CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N1CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("N2CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N3CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N6CLSE8O15ZCsd.txt", "r");
	//fp = fopen("N8CLSNE8O15ZCsd.txt", "r");
	//fp = fopen("N2e4CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N2e6CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N2e56CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N2e56CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N8e4CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("N1e6CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("N1e5CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("N2e5CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N4e6CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N4e56CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N1e3CLSNE8O42ZCsd.txt", "r");
	//fp = fopen("N4CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N4e5CLSE8O42ZCsd.txt", "r");
	//fp = fopen("N2BCLSE8O42ZCsd.txt", "r");
	//fp = fopen("N7CLSE8O15ZCsd.txt", "r");
	fp = fopen("N6BCLSE8O15ZCsd.txt", "r");



	//printf("SD...\n");
	 for (i=0; i<ParCnt2;i++){
	   fscanf(fp, "%lf\n", &SD[i]);
	   //printf("i:%d SD:%lf\n",i,SD[i]);
	}
	fclose(fp);
	//printf("Done...\n");
	//getc(stdin); 


	//printf("All Done with PCA, hit return...\n"); getc(stdin); 
	//printf("All Done with PCA\n"); 

	//exit(1);


	for(i=1;i<=4;i++) FxdPars[i] = CurrentPars[i] = 0.0;



	
///////////////////////   Done reading in PCA Results /////////////////////////////////


///// Dot product function /////////

	
	
float DotProduct (int Length, double *Holder, double *PCA)
{

  double answer = 0;
  
  int i;
  for(i=0;i<Length;i++){
    answer += PCA[i]*Holder[i];
   //printf("i:%d Holder:%f PCA:%f answer:%f \n",i,Holder[i],PCA[i],answer);
    
  }
  return(answer);

}


///////////////////////  Initial Parameter Set ///////////////////

//Draw a complete set of parameters on the PCA scale
 double PCA[30], OldPCA[30];
 VarScale = 1.2;
 

 int PCAParCnt = 0;
  for(i=0;i<ParCnt2;i++){
		NoGo = 0;
		
		
		if(Evoln==1){
			if((sFix==1)&&(i==2)) NoGo = 1;
			if((gammaFix==1)&&(i==4)) NoGo = 1;
			if((sigmaFix==1)&&(i==5)) NoGo = 1;
			if((bFix==1)&&(i==6)) NoGo = 1;
		} else{
			if((gammaFix==1)&&(i==3)) NoGo = 1;
			if((sigmaFix==1)&&(i==4)) NoGo = 1;
		}
		

		if(NoGo==0){ OldPCA[PCAParCnt] = PCA[PCAParCnt] = gsl_ran_gaussian (r2, VarScale*SD[PCAParCnt]);}
		
		if(NoGo!=1){
			 //PCAParCnt allows you to skip past parameters that are not being varied
			 //printf("PCAParCnt:%d PCA:%lf SD:%f\n",PCAParCnt,PCA[PCAParCnt],SD[PCAParCnt]);
			 PCAParCnt++;
		}
  }


ParCnt2 = PCAParCnt;


//Now convert the first set to the original, linear scale

int DontLHood = 0;
//Inflate = InitScale;


double Holder[30];



LHCount = 1;
OldPost = 0;
int Cancel2 = 0;
 
  
//Get parameters onto original scale
	PCAParCnt = 0; OldPropAdj = 0.0;
 	for(PC=1;PC<=ParamCount[DataSet];PC++){
		  NoGo = 0;
		  ParamChoice = ParamNum[PC];


		  if(Evoln==1){
			if((sFix==1)&&(ParamChoice==2)) NoGo = 1;
			if((gammaFix==1)&&(ParamChoice==4)) NoGo = 1;
			if((sigmaFix==1)&&(ParamChoice==5)) NoGo = 1;
			if((bFix==1)&&(ParamChoice==6)) NoGo = 1;
		  } else{
			if((gammaFix==1)&&(ParamChoice==3)) NoGo = 1;
			if((sigmaFix==1)&&(ParamChoice==4)) NoGo = 1;
		  }

		  //printf("NoGo:%d PC:%d ParamNum:%d ParamChoice:%d ParamCount:%d \n",NoGo,PC,ParamNum[PC],ParamChoice,ParamCount[DataSet]);
		  //getc(stdin);
		  		  
		  if(NoGo!=1){

				
				for (j=0;j<ParCnt2; j++) Holder[j] = Coefficients[PCAParCnt][j]; //Holder just holds the jth column of Coefficients ALL AT ONCE 
		
    
				   double dummyDot = DotProduct(ParCnt2,Holder,PCA);				   
				   temp = exp(dummyDot*VarScale*Scale[PCAParCnt]+Center[PCAParCnt]);

				  
				   if(temp<0) temp = 1e-10;

				    ParamVals.FxdPars[PC] = temp;
				    FxdPars[PC] = temp;
			            CurrentPars[PC] = temp; //Only the first time...I think...
				  //printf("ParamChoice:%d dummyDot:%f VarScale:%f Scale:%f Center:%f temp:%f\n",ParamChoice,dummyDot,VarScale,Scale[PCAParCnt],Center[PCAParCnt],temp);
				  // getc(stdin);
				
				    PCAParCnt++;

		  } //NoGo 
	} //PC


//Calculate the likelihood	
	AvgLH = 0.0;
	Cancel = 0; 
	
	//maxBoots = 25;


   			

			bootGood = 0; boot = 0;
			//while( (bootGood<minBootGood) && (boot<maxBoots) ) {
			//while(boot<maxBoots) {
			for(boot=0;boot<SEIRParamNum;boot++){
 				//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
				int Index = boot; 

				//printf("lambda:%f phi:%f gamma:%f \n",FxdPars[1],FxdPars[2],FxdPars[3]); 
				//printf("k:%f nu:%f mu:%f sigma:%f ratio:%f delta:%f m:%F \n",kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]); 

   				res = LHood(FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);

				//printf("boot:%d res:%f\n",boot,res); 
				if(res>-1e10){
		 			bootGood++; 
		 			OldAvgLH += exp(res);
					//printf("Initial LH calc:%d LH:%f\n",boot,res); //getc(stdin);  getc(stdin);
		 			
				}

				//boot++;
		

				

   			} //for(boot


			if(isnan(AvgLH)==1){
				Cancel = 1; 
			}else{

				
				//OldAvgLH /= bootGood;
				OldAvgLH /= SEIRParamNum;
				OldAvgLH = log(OldAvgLH);
				OldPost = (OldAvgLH);
				OldLHoodPerPop = OldAvgLH;

				//printf("OldPost:%f\n",OldPost);
			}

	//printf("MCMC Initial LH calc finished...boot:%d bootGood:%d OldAvgLH:%f OldPost:%f\n",boot,bootGood,OldAvgLH,OldPost,OldPost); 
	//getc(stdin);
  
	//exit(1);




//////////// Adjusting for prior /////////////
		for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		  
		      
		  NoGo = 0;
		  if(Evoln==1){
			if((sFix==1)&&(PC2==2)) NoGo = 1;
			if((gammaFix==1)&&(PC2==4)) NoGo = 1;
			if((sigmaFix==1)&&(PC2==5)) NoGo = 1;
			if((bFix==1)&&(PC2==6)) NoGo = 1;
		  } else{
			if((gammaFix==1)&&(PC2==3)) NoGo = 1;
			if((sigmaFix==1)&&(PC2==4)) NoGo = 1;
		  }

		
		  
		  if(NoGo!=1){
		 		
	               float tmp2;
		        switch(PriorDist[PC2]){
			case 1 :
			  if(FxdPars[PC2]>MaxParam[PC2]){
			    
			    Stop = 1;
			    NewPost = 0;
			    //printf("PC2:%d FxdPars:%f MaxParam[PC2]:%f\n",PC2,FxdPars[PC2],MaxParam[PC2]); getc(stdin);
			    

			  } else {
			    //printf("PC2:%d PriorVar:%f FxdPars:%f\n",PC2,PriorVar[PC2],FxdPars[PC2]); getc(stdin);
			    tmp2 =  (gsl_ran_flat_pdf(FxdPars[PC2],0,PriorVar[PC2]));
				


			  }
			  break;
			case 2 :tmp2 = gsl_ran_gaussian_pdf(log(FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // a slightly bogus kind of normal, really sort log normal

			
				break;

			case 3 :  tmp2 = (gsl_ran_lognormal_pdf(FxdPars[PC2],PriorMean[PC2],PriorVar[PC2]));

					 
					 break;
	

			 case 4 :tmp2 = gsl_ran_gaussian_pdf((FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // more like a real normal...

				
				 break;
	                 default : printf("No prior distribution specified.  Bailing.\n"); getc(stdin); exit(1); break;
		      	} 
			if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
				Cancel = 1;
				
				tmp2 = 1e-300;
			} else{
				OldPost += log(tmp2);

				//printf("FxdPars:%f PC2:%d OldPost:%f tmp2:%f\n",FxdPars[PC2],PC2,OldPost,tmp2); getc(stdin);
				
			}
			if((isnan(OldPost)!=0)||(isinf(OldPost)!=0)){
				Cancel = 1;
								OldPost = -1000;
								
			}
			



		      } //ParamChoice2!=4, we're not doing sigma
	
		    } //PC2


		    //printf("OldPost:%f\n",OldPost); getc(stdin);


//////////////// End of Initial Parameter Set



Itn = 1;
Stop = 0;
int Thin = 0;

 //printf("Just before while statement...\n");

//Metropolis-Hastings steps
 while(Itn<MaxItn){

	//printf("Itn:%f\n",Itn); getc(stdin);
		
		

	
	int PCAParCnt = 0;
	
  	for(i=0;i<ParCnt2;i++){
			
		
		

		NoGo = 0;
		if(Evoln==1){
			if((sFix==1)&&(i==2)) NoGo = 1;
			if((gammaFix==1)&&(i==4)) NoGo = 1;
			if((sigmaFix==1)&&(i==5)) NoGo = 1;
			if((bFix==1)&&(i==6)) NoGo = 1;
		} else{
			if((gammaFix==1)&&(i==3)) NoGo = 1;
			if((sigmaFix==1)&&(i==4)) NoGo = 1;
		}
		if(NoGo!=1){

//Draw one parameter
			 
			PCA[PCAParCnt] = gsl_ran_gaussian (r2, SD[PCAParCnt]); 


			
	///////////// Given the new parameters, make a complete set
		

		int PCAParCnt2 = 0;
		
 		for(PC=1;PC<=ParamCount[DataSet];PC++){
			
		  NoGo = 0;
		  if(Evoln==1){
			if((sFix==1)&&(PC==2)) NoGo = 1;
			if((gammaFix==1)&&(PC==4)) NoGo = 1;
			if((sigmaFix==1)&&(PC==5)) NoGo = 1;
			if((bFix==1)&&(PC==6)) NoGo = 1;
		  } else{
			if((gammaFix==1)&&(PC==3)) NoGo = 1;
			if((sigmaFix==1)&&(PC==4)) NoGo = 1;
		  }
		
		//printf("PC:%d MaxParam:%f\n",PC,MaxParam[PC]);
		  
		  if(NoGo!=1){

			for (j=0;j<ParCnt2; j++) Holder[j] = Coefficients[PCAParCnt2][j]; //Holder just holds the jth column of Coefficients ALL AT ONCE
    

			double temp2;

			double dummyDot = DotProduct(ParCnt2,Holder,PCA);

				   
			temp = exp(dummyDot*Scale[PCAParCnt2]+Center[PCAParCnt2]);
			if(temp<0) temp = 1e-10;

			ParamVals.FxdPars[PC] = temp;
			FxdPars[PC] = temp;

				    
			if(FxdPars[PC]>MaxParam[PC]){
			    			
						//printf("PC:%d FxdPars:%f MaxParam:%f \n",PC,FxdPars[PC],MaxParam[PC]); 
			    			Stop = 1;
			    			NewPost = 0;
			  }
				    

			//printf("Stop:%d j:%d ParCnt2:%d ParamCount:%d PC:%d temp:%f MaxParam:%f\n",Stop,j,ParCnt2,ParamCount[DataSet],PC,temp,MaxParam[PC]); 
			//getc(stdin); 
			PCAParCnt2++;
				  	  					  
				    
		  } //NoGo 
		} //PC
		

		//printf("done with PC loop...\n");// getc(stdin);

	  ////////////// Done converting complete set to original axes ////////////



	 ////////////// Calculate likelihood ///////////////

		//printf("Stop:%d...hit return\n",Stop); 
		//getc(stdin);

		if(Stop<1){
			CrashOut = 0;
			NewOldLHoodPerPop = 0.0; 
			LHoodPerPop = 0.0;


       		        AvgLHood = 0.0;
		

					Cancel = 0;

					
 					
      					//First time: do it for proposed parameter values

					AvgLH = 0.0;
					bootGood = 0; 
					//boot = 0;
					//while((bootGood<minBootGood) && (boot<maxBoots)) {
					for(boot=0;boot<SEIRParamNum;boot++){

						//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
						int Index = boot;  
   						res = LHood(ParamVals.FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);
						//printf("New LHood boot:%d Index:%d res:%e\n",boot,Index,res);

						if(res>-1e10){
		 					bootGood++; 
		 					AvgLH += exp(res);
		 					//printf("New LHood boot:%d Index:%d res:%e\n",boot,Index,res);
						}

						//boot++;

   					} //for(boot

					//printf("boot:%d bootGood:%d not yet averaged OldAvgLH:%f\n",boot,bootGood,log(AvgLH)); 
					//printf("hit return...\n");getc(stdin);

					if((isnan(AvgLH)==1)||(bootGood<1)){
						Cancel = 1; 
					}else{
						Cancel = 0;
					}
					//printf("Cancel:%d \n",Cancel); getc(stdin);

					if(Cancel!=1){

						
					
						//AvgLH /= bootGood;
						AvgLH /= SEIRParamNum;
						AvgLH = log(AvgLH);
		   				NewPost = AvgLH;
						LHoodPerPop = AvgLH;
						//printf("NewPost:%f\n",NewPost);

						
					}else{
						NewPost = -1000;
					}
					
					
					//Now set ParamVals.FxdPars back to the OLD set of parameter values...
					for(PC=1;PC<=ParamCount[DataSet];PC++){
						
						ParamVals.FxdPars[PC] = CurrentPars[PC]; //Go back and re-calculate the likelihood again for the current parameters
						

					}
					//getc(stdin);

					
					
				       //and re-calculate the OLD likelihood, to avoid winner's curse
				       //but only for a stochastic model...
				      //if(sigmaFix!=1){
						OldAvgLH = 0.0;
						bootGood = 0; boot = 0;
						for(boot=0;boot<SEIRParamNum;boot++){
					
						//int Index = gsl_ran_flat(r2,0.0,SEIRParamNum); 
							int Index = boot; 
   							res = LHood(ParamVals.FxdPars,kIn[Index],nuIn[Index],muIn[Index],sigmaIn[Index],ratioIn[Index],deltaIn[Index],mIn[Index]);

							if(res>-1e10){
		 						bootGood++; 
		 						OldAvgLH += exp(res);
		 						//printf("Old LHood boot:%d Index:%d res:%e\n",boot,Index,res);
							}

						} //for(boot	
						//OldAvgLH /= bootGood;
						OldAvgLH /= SEIRParamNum;
						OldAvgLH = log(OldAvgLH);
						
				      // }
				
					//printf("bootGood:%d OldAvgLH:%f OldPost:%f\n",bootGood,(OldAvgLH),OldPost);
					//getc(stdin);
					//printf("OldAvgLH:%f\n",OldAvgLH);

								
					
					//Now set ParamVals.FxdPars back to the proposed set of parameters				
				      	for(PC=1;PC<=ParamCount[DataSet];PC++){
						
						ParamVals.FxdPars[PC] = FxdPars[PC]; //CurrentPars holds the old parameters, so go back to using ParamVals.FxdPars to  hold the proposed parameters

					}

					if((isnan(AvgLH)!=1)&&(bootGood>0)){
												
						NewOldLHoodPerPop = OldAvgLH;
						//printf("Before OldPost OldPost:%f OldAvgLH:%f OldLHoodPerPop:%f\n",OldPost,OldAvgLH,OldLHoodPerPop);
		    				OldPost += NewOldLHoodPerPop - OldLHoodPerPop;
						//printf("OldPost:%f NewOld:%f Old:%f \n",OldPost,NewOldLHoodPerPop,OldLHoodPerPop); 
						OldLHoodPerPop = NewOldLHoodPerPop;
					
					} else{
						NewOldLHoodPerPop = OldLHoodPerPop;
					}				


					
					

	
			
	
		/////////////// Doing the calculation to adjust (later) for the proposal ////////////////////

		NewPropAdj = log(gsl_ran_gaussian_pdf(PCA[PCAParCnt],SD[PCAParCnt])); 
		if(Itn>1)  VarScale = 1.0;  //WRONG?? SHOULD PROBABLY BE IF(ACCEPTCOUNT>0)
		OldPropAdj = log(gsl_ran_gaussian_pdf(OldPCA[PCAParCnt],VarScale*SD[PCAParCnt]));

		

		/////////////// Done doing the calculation to adjust (later) for the proposal



 
		//////////// Adjusting for prior /////////////
		for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		      
		  NoGo = 0;
		  if(Evoln==1){
			if((sFix==1)&&(PC2==2)) NoGo = 1;
			if((gammaFix==1)&&(PC2==4)) NoGo = 1;
			if((sigmaFix==1)&&(PC2==5)) NoGo = 1;
			if((bFix==1)&&(PC2==6)) NoGo = 1;
		 } else{
			if((gammaFix==1)&&(PC2==3)) NoGo = 1;
			if((sigmaFix==1)&&(PC2==4)) NoGo = 1;
		 }

		//printf("inside PC2 loop...\n");
		
		  
		  if(NoGo!=1){
		 		
	               float tmp2;
		        switch(PriorDist[PC2]){
			case 1 :
			  if(FxdPars[PC2]>MaxParam[PC2]){
			    
			    Stop = 1;
			    NewPost = -1000;
			    
			  } else {
			    tmp2 =  (gsl_ran_flat_pdf(FxdPars[PC2],0,PriorVar[PC2]));
			  }
			  break;
			case 2 :tmp2 = gsl_ran_gaussian_pdf(log(FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // a slightly bogus kind of normal, really sort log normal
		
				break;

			case 3 :  tmp2 = (gsl_ran_lognormal_pdf(FxdPars[PC2],PriorMean[PC2],PriorVar[PC2]));
				 
					 break;
	
			 case 4 :tmp2 = gsl_ran_gaussian_pdf((FxdPars[PC2])-PriorMean[PC2],PriorVar[PC2]); // more like a real normal...

				
				 break;
	                 default : printf("No prior distribution specified.  Bailing.\n"); getc(stdin); exit(1); break;
		      	} //End of switch

			if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
				Cancel = 1;
				
				tmp2 = 1e-300;
			} else{
			
				//printf("Here we are PC2:%d PriorVar:%f FxdPars:%f log(tmp2):%f NewPost:%f\n",PC2,PriorVar[PC2],FxdPars[PC2],log(tmp2),NewPost);
			    	NewPost += log(tmp2);
				//printf("log(tmp2):%f NewPost:%f\n",log(tmp2),NewPost); getc(stdin);
	

				
			}


			if((isnan(NewPost)!=0)||(isinf(NewPost)!=0)){
				Cancel = 1;
				NewPost = -1000;
						
			}
			
		      }  //NoGo	
	        } //PC2
	
		
		//printf("Itn:%f Fxd1:%f Fxd2:%f Fxd3:%f Fxd4:%f Fxd5:%f Curr1:%f Curr2:%f Curr3:%f AvgLH:%f OldAvgLH:%f NewPost:%f OldPost:%f\n",Itn,FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],CurrentPars[1],CurrentPars[2],CurrentPars[3],AvgLH,OldAvgLH,NewPost-NewPropAdj,OldPost-OldPropAdj);

		if(Itn>0){ 

		 
		  if(Cancel!=1){//Cancel if prior gives bad calc or lhood gives bad calc
		    

			if(OldPost>-900){
				Criterion = exp(NewPost - NewPropAdj -(OldPost-OldPropAdj));
			}else{
				Criterion = 100.0;
			}
			


			if((isnan(Criterion)!=0)||(isinf(Criterion)!=0)){
				Criterion = -1;
			}
			
			temp = gsl_rng_uniform(r2);
			
			if((Criterion>1)||(temp<Criterion)){  //accept


			
			  if(Criterion<1){
			  	//printf("Itn:%f Accept the hard way...\n",Itn);
			  }else{
				//printf("Itn:%f Accept the easy way...\n",Itn);
			  }
			
			  			
			  //Acceptance means changing the current parameters to the proposed parameters
			  for(PC2=1;PC2<=ParamCount[DataSet];PC2++){  //Isn't this whole loop pointless?  No, this is where FxdPars gets reset to CurrentPars.  Soon, ONE element of FxdPars will change, BUT ONLY ONE 
		  		 
		                CurrentPars[PC2] = FxdPars[PC2];
			     
			                
		  	  }



			  OldPCA[PCAParCnt] = PCA[PCAParCnt];
				
	
			  OldPost = NewPost;
			  OldLHoodPerPop = LHoodPerPop;
			  
			 

			   
			} else { //Criterion < 1 and temp > Criterion
				OldLHoodPerPop = NewOldLHoodPerPop;
				//printf("Itn:%f Do not accept...\n",Itn);
				
			}
		 } //if Cancel

		      
		 //Reset FxdPars to the current parameters, which may or may not have been updated	         	
		  for(PC2=1;PC2<=ParamCount[DataSet];PC2++){  //Isn't this whole loop pointless?  No, this is where FxdPars gets reset to CurrentPars.  Soon, ONE element of FxdPars will change, BUT ONLY ONE 
		  		 
		               			     
			        FxdPars[PC2] = CurrentPars[PC2];
				ParamVals.FxdPars[PC2] = CurrentPars[PC2];
				
				
			         
		  }


   
		  if(Thin==ThinStop){

			int prnt = 0;
			 strcpy(testbuff,test[0]);
		  	 fp2 = fopen(testbuff,"a");
		        
		        fprintf(fp2,"%f  ",Itn);
			 if(prnt==1) printf("%f  ",Itn);

			 for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		                 		     
			         
				 NoGo = 0;
                            				 
				 fprintf(fp2,"%e  ",FxdPars[PC2]);
				 if(prnt==1) printf("PC2:%d %e  ",PC2,FxdPars[PC2]);

			 }
			 
			 fprintf(fp2,"%e %e\n",OldLHoodPerPop,OldPost);
			 if(prnt==1) printf("%e %e\n",OldLHoodPerPop,OldPost-OldPropAdj);

			 if(OldLHoodPerPop<(OldPost-OldPropAdj)){ printf("Problem...\n"); getc(stdin);}
			

			 fclose(fp2); 
			 
			 Thin = 0;
		 }//Thin == ThinStop
			  
			
	} //If(Itn>0)
			
					          	  
	


		 Itn++;
		 Thin++;
		 //ParamItn[ParamChoice]++;
		//printf("Itn:%d ParamItn:%d\n",Itn,ParamItn[ParamChoice]);
		
		if(Verbose2==1){
			strcpy(testbuff,test[2]); //Stripped file no longer needed
		  	 fp3 = fopen(testbuff,"a");
		        fprintf(fp3,"At end of PC loop %d  ",Itn);
			 for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		            ParamChoice3 = ParamNum[PC2];
			     fprintf(fp3,"%f  ",FxdPars[ParamChoice3]); 
				 
			 }
			 fprintf(fp3,"\n");
			 fclose(fp3);
		}
		

		 

	
	         } //Stop 

			Stop = 0;
			PCAParCnt++;
		} //PCA loop?
 	}  //i??		




 



} //Itn










	



//gsl_rng_free (r);

free_dvector(kIn,1,NumSets); //check
free_dvector(nuIn,0,NumSets);
free_dvector(muIn,0,NumSets);
free_dvector(sigmaIn,0,NumSets);
free_dvector(ratioIn,0,NumSets);
free_dvector(deltaIn,0,NumSets);
free_dvector(mIn,0,NumSets);



/*
	rIn = dvector(0,NumSets);
	sIn = dvector(0,NumSets);
	phiIn = dvector(0,NumSets);
*/


 return 0;

}


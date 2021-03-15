int Fuck2;
int Verbose;
int Verbose2; 
int CrashOut=0; 
int ExtRlzns = 0;
int Count = 0;
int DataSet = 3; //16;  //3; //17; //2; //16;  //2 is woods & elk 5 plots, 3 is woods and elk 8 plots, 16 is Otvos et al, 17 is Karl Polivka, 18 is both Otvos and Shepherd, both control and treatment, 19 is O&S with tmt only, 13 is O&S with con only
 
int DIC = 0; 
int WAIC = 1; //Now go and look for Prm, the input file (really .in)
int LineSearch = 0;
int Figures = 0; //SEIR figures, along with Fract Inf.  The SEIR outpart part is in Distd1.c, go change the file name there.  The remainder is in main text and FigureGenerate.  IC's are in FigureGenerate.  
//Have to change 2 output file names for Figures
int ModelVsData = 0; //Just Fract Inf with multiple plots and stochasticty.  Found in LHood only.  FOR THIS TO WORK YOU SHOULD HAVE DIC = 1, ALSO CHECK OUTPUT FILE NAME!

double PointwiseLH[30][100], PointwiseVar[30][100], PointwiseMean[30][100];

float FractI2[20];

typedef struct {

	double FxdPars[30];
	

} ParamStruct;

#include <time.h>
#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include <omp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>	// GSL_SUCCESS ...
#include <gsl/gsl_odeiv.h>	// ODE solver
//#include <gsl/gsl_odeiv2.h>	// ODE solver - THE NEW ODE SOLVER, WHICH MIDWAY DOESN'T SUPPORT!  THE HORROR! THE HORROR!

#include <gsl/gsl_types.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>


#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//#define PI 3.141592654

//#include "round.c"

//#include "MultiSG9.auto.h"
#include "math.h"
#include "stdio.h"
#include "DataFuncs.c"

#include "fast_odesGD.h"

#include "Nrutil.c"
#include "LHoodFunc.c"

//#include "ran2.c"
long seed;


//#include "Distd1.c"
#include "DistdRK45.c"

gsl_rng *r2;






double FigureGenerate(double *RandNumsPass, int dim, double *FxdPars){


// double *AdjPars, double *lhoodVec, int DataSet){

	float *S0, *FractI;
	float ***Data, **Covars;
	float ***ModelFI,***ModelS;
	float ***SSE;
	double *lhoodTemp1, *lhoodTemp2;
	int *Split;
	//double ***RandNums;
	int *MaxDays;
	long seed = 18; 
	//double PriorMedian = 1, PriorSigma = 0.6, PostSigma = 0.6;
	//double **PostMedian, *MixRatio;
	
	//int N;
	int nuNum;
	double Output[100][20];

	double m_Best_LHood = 1e10; 

	double ax,bx,cx,xmin,gold;
	int nmin=0;
	int MaxWk = 9;
	int Rlzn, Wk;
	int MaxPlots = 5;
	//double  *FxdPars;
	double tau;
	double TOL = 1e-6;
	double junk;
	int i,j;
	double answer;
	int test;
	int DataInterval,DataType;
	int SData, FractDead;
	int WhichData = 1;
        double AdjPars[30];
        //int DataSet;
       //ParamStruct* Params;
	//Params = (ParamStruct*) ParamStuff; 
        int MaxRlzns = 1;
	double lhoodVec[1];


	int Plot = FxdPars[0];
	FILE *fp, *fp2; 

        double ***RandNums;
	RandNums = d3tensor(1,MaxRlzns,1,25,0,50);


        //fp = fopen("Junk2.dat","wt");
	//Inverse Gaussian cdf: double gsl_cdf_gaussian_Pinv (double P, double sigma)
	  /*@@@*/
	  double mean=FxdPars[2];
	  double sd=sqrt(FxdPars[4]);

	//sd = 0.0;

	
	//fp = fopen("RandNumsJunk.dat","w");
	for(i=0;i<dim;i++){
		
            //RandNums[1][Plot][i] = exp(mean+gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd));  //Greg's version 
	    if(sd>0){
		//The above line assumes that random numbers come from Monte miser, but that's no longer true, so have to gen them yrslf
	     	//RandNums[1][1][i] = exp(mean+gsl_ran_gaussian(r2,sd));  //WRONG - in this version, mean has been exponentiated twice.  Oops. 
		RandNums[1][1][i] = mean*exp(gsl_ran_gaussian(r2,sd));
	    }else{
		RandNums[1][1][i] = exp(mean); 
	    } 
		
	      //fprintf(fp,"%f\n",RandNumsPass[i]);
	     //printf("%e\n",i,RandNums[1][1][i]);

	}

	//printf("just after RandNums...\n");
	//fclose(fp);
	//exit(1);


	//FxdPars = dvector(1,100);
	S0 = vector(1,100);
      	FractI = vector(1,100);
	Covars = matrix(1,25,1,25);
	Data = f3tensor(0,10,0,20,0,50);
	ModelFI = f3tensor(1,MaxRlzns,1,25,0,25); //indices:Rlzns,plots,weeks
        ModelS = f3tensor(1,MaxRlzns,1,25,0,25); //indices:Rlzns,plots,weeks
	SSE = f3tensor(1,MaxRlzns,1,25,0,25);
	lhoodTemp1 = dvector(1,MaxRlzns);
	lhoodTemp2 = dvector(1,MaxRlzns);
	Split = ivector(1,100);
	MaxDays = ivector(1,100);

	int RandNumCount = 5;

	for(DataType=0;DataType<=10;DataType++)
		for(Plot=0;Plot<=20;Plot++)
			for(i=0;i<=50;i++)
				Data[DataType][Plot][i] = -2;

	//DataSet = (int) Params->FxdPars[11];

  
  	
  MaxPlots = 1;
  MaxDays[1] = 42; //42
  S0[1] = 250; //hi is 250, low is 70
  FractI[1] = 0.5; //0.5; //0.1, 0.5 goes with 140, 0.1 goes with 70
  Split[1] = 0;
  FxdPars[0] = 1;
  FxdPars[4] = 0.0; //sigma
  //FxdPars[7] = 1e-3;
	

//printf("Inside LHood MaxPlots:%d\n",MaxPlots);
		for(i=1;i<=MaxPlots;i++){
			Covars[1][i] = S0[i];
			Covars[2][i] = FractI[i];
			Covars[3][i] = Split[i];
		}

        for(i=0;i<30;i++){
		AdjPars[i] = FxdPars[i];
	}

	

//AdjPars[7] = 0.0005; //0.01; //0.25; //h  //I tried different values of h, but the effect was trivial

//test = Distd1(AdjPars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);

test = DistdRK45(AdjPars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);


//fp = fopen("ModelFIOut9Aug14lokLoDens.dat","a"); //Det version
fp = fopen("StochFIOut9Aug14lokHiDens.dat","a"); //Stoch version

for(i=1;i<=5;i++){
	 printf("i:%d ModelFI[i]:%f\n",i,ModelFI[1][1][i]);
	 fprintf(fp,"%d %f\n",i,ModelFI[1][1][i]);

}
fclose(fp);

//test = DDE_SSE9(AdjPars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);


	CrashOut += test;
        //fp = fopen("CrashOut.dat","a");
	//fprintf(fp,"%d\t%d\t%d\n",ExtRlzns,CrashOut,test);
	//fclose(fp);
	ExtRlzns++;




	ax = 1e-20;
	bx = 200;
	cx = 0.15;

	WhichData = 1;
        MaxPlots = AdjPars[0];
	

	/*
	if((DataSet!=2)&&(DataSet<12)){
		test = SSECalc(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,SSE,WhichData);
        	xmin = AdjPars[10];
		answer = LHoodFunc2(MaxRlzns,MaxPlots,MaxWk,SSE,xmin,lhoodTemp1);
	} else {
//for straight-up binomial, use just the next line
		//answer = LhoodFuncBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1); //straight up binomial 
//for the pseudo-likelihood, just the next TWO lines
		//test = SSECalcPseudo(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,SSE); //Pseudo-likelihood binomial 
		//answer =LHoodFuncPseudo(MaxRlzns,MaxPlots,MaxWk,Data,SSE,AdjPars[10],AdjPars[20],lhoodTemp1); //Pseudo-lhood 
//for beta-binomial, use the next line
  		answer = LhoodFuncBetaBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1,AdjPars[20],AdjPars[10]); //straight up binomial, 
//20 is gamma, 10 is obsvn error, which in this case is "delta"
	}
	*/


	answer = LhoodFuncBetaBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1,AdjPars[20],AdjPars[10]);


	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
		lhoodVec[Rlzn] = lhoodTemp1[Rlzn];
		//if(SData==1) lhoodVec[Rlzn] += lhoodTemp2[Rlzn];
	}		



	free_vector(S0,1,100); //check
	free_vector(FractI,1,100); //check
	free_f3tensor(Data,0,10,0,20,0,50);
	free_matrix(Covars,1,25,1,25); //check
	free_ivector(Split,1,100); //check
	free_ivector(MaxDays,1,100); //check
	free_f3tensor(ModelFI,1,MaxRlzns,1,25,0,25); //check
	free_f3tensor(ModelS,1,MaxRlzns,1,25,0,25); //check
	free_f3tensor(SSE,1,MaxRlzns,1,25,0,25); //check
	free_dvector(lhoodTemp1,1,MaxRlzns); //check
	free_dvector(lhoodTemp2,1,MaxRlzns); //check
	free_d3tensor(RandNums,1,MaxRlzns,1,25,0,50); //check
	

	//if(lhoodVec[1]>-250)
		return exp(lhoodVec[1]);
	//else
		//return exp(-250);
}






double LHoodRK45(double *RandNumsPass, size_t dim, void *ParamStuff){


	

// double *AdjPars, double *lhoodVec, int DataSet){

	float *S0, *FractI;
	float ***Data, **Covars;
	float ***ModelFI,***ModelS;
	float ***SSE;
	double *lhoodTemp1, *lhoodTemp2;
	int *Split;
	//double ***RandNums;
	int *MaxDays;
	long seed = 18; 
	//double PriorMedian = 1, PriorSigma = 0.6, PostSigma = 0.6;
	//double **PostMedian, *MixRatio;
	
	//int N;
	int nuNum;
	double Output[100][20];

	double m_Best_LHood = 1e10; 

	double ax,bx,cx,xmin,gold;
	int nmin=0;
	int MaxWk = 9;
	int Rlzn, Wk;
	int MaxPlots = 5;
	double  *FxdPars;
	double tau;
	double TOL = 1e-6;
	double junk;
	int i,j;
	double answer;
	int test;
	int DataInterval,DataType;
	int SData, FractDead;
	int WhichData = 1;
        double AdjPars[30];
        //int DataSet;
                ParamStruct* Params;
	Params = (ParamStruct*) ParamStuff; 
        int MaxRlzns = 1;
	double lhoodVec[1];


	int Plot = Params->FxdPars[0];
	
	FILE *fp, *fp2; 

        double ***RandNums;
	RandNums = d3tensor(1,MaxRlzns,1,25,0,50);


        //fp = fopen("Junk2.dat","wt");
	//Inverse Gaussian cdf: double gsl_cdf_gaussian_Pinv (double P, double sigma)
	  /*@@@*/
	  double mean=Params->FxdPars[2];
	  double sd=sqrt(Params->FxdPars[4]);

	

	//fp = fopen("RandNumsJunk.dat","w");
	for(i=0;i<dim;i++){
             //RandNums[1][Plot][i] = exp(mean+gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd)); //WRONG - in this version, mean has been exponentiated twice.  Oops.
	     RandNums[1][Plot][i] = mean*exp(gsl_cdf_gaussian_Pinv(RandNumsPass[i],sd)); 
	      //fprintf(fp,"%f\n",RandNumsPass[i]);
	}
	//fclose(fp);
	//exit(1);


	FxdPars = dvector(1,100);
	S0 = vector(1,100);
      	FractI = vector(1,100);
	Covars = matrix(1,25,1,25);
	Data = f3tensor(0,10,0,20,0,50);
	ModelFI = f3tensor(1,MaxRlzns,1,25,0,25); //indices:Rlzns,plots,weeks
        ModelS = f3tensor(1,MaxRlzns,1,25,0,25); //indices:Rlzns,plots,weeks
	SSE = f3tensor(1,MaxRlzns,1,25,0,25);
	lhoodTemp1 = dvector(1,MaxRlzns);
	lhoodTemp2 = dvector(1,MaxRlzns);
	Split = ivector(1,100);
	MaxDays = ivector(1,100);

	int RandNumCount = 5;

	for(DataType=0;DataType<=10;DataType++)
		for(Plot=0;Plot<=20;Plot++)
			for(i=0;i<=50;i++)
				Data[DataType][Plot][i] = -2;

	//DataSet = (int) Params->FxdPars[11];


	switch(DataSet){
			case 1 : MaxPlots = swplots5(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 2 : MaxPlots = swplots5Bnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead);  
				 break;
			case 3 : MaxPlots = swplots8Bnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead);  
				 break;

			//case 3:  MaxPlots = Mandarte(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 4:  MaxPlots = ShepherdControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 5:  MaxPlots = ShepherdAll(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 6:  MaxPlots = ShepherdSpray(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 7:  MaxPlots = ShepherdAerial(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 8:  MaxPlots = OtvosControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 9:  MaxPlots = OtvosAll(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 10:  MaxPlots = OtvosSpray(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 11:  MaxPlots = MoreauControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 12:  MaxPlots = OtvosConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 13:  MaxPlots = OtvosShepConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 14:  MaxPlots = MoreauConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 15:  MaxPlots = OtvosShepConBnml2(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 16:  MaxPlots = OtvosAllBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 17:  MaxPlots = KPBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 18:  MaxPlots = OtvosShepAllBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 19:  MaxPlots = OtvosShepTmtBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;

	}


		
		for(i=1;i<=MaxPlots;i++){
			Covars[1][i] = S0[i];
			Covars[2][i] = FractI[i];
			Covars[3][i] = Split[i];
			
		}
	
		
        for(i=0;i<30;i++){
		AdjPars[i] = Params->FxdPars[i];
		//printf("i:%d Params->FxdPars[i]:%f\n",i,Params->FxdPars[i]);
	}
	// getc(stdin);

Plot = Params->FxdPars[0];



test = DistdRK45(AdjPars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);

int TotWeeks = ((int) MaxDays[Plot]/7) - 1;
if(ModelVsData==1){
	//fp = fopen("r2ModelFIOutGA8O42ZC.dat","a");
	fp = fopen("r2ModelFIOutGA8P0ZC.dat","a");
	//fp = fopen("r2ModelFIOutGA8O15ZC.dat","a");


	fprintf(fp,"%d %d %f %f %f\n",Plot,0,0.0,0.0,S0[Plot]);
	fprintf(fp,"%d %d %f %f %f\n",Plot,1,FractI[Plot],0.0,S0[Plot]);
	fprintf(fp,"%d %d %f %f %f\n",Plot,2,FractI2[Plot],0.0,S0[Plot]);
	
	if(Plot==8) TotWeeks++;
	for(i=1;i<=TotWeeks;i++){
		double dummy = ((float) Data[1][Plot][i])/((float) Data[2][Plot][i]);

	 	
		//if(Plot<8){
	 		fprintf(fp,"%d %d %f %f %f\n",Plot,i+2,dummy,ModelFI[1][Plot][i],S0[Plot]);
		//}else{
			//printf("MvD Plot:%d i:%d ModelFI[i]:%f\n",Plot,i+2,ModelFI[1][Plot][i]);
			//fprintf(fp,"%d %d %f %f %f\n",Plot,i+1,dummy,ModelFI[1][Plot][i],S0[Plot]);
		//}
			

	}
	fclose(fp);
	
}




//test = DDE_SSE9(AdjPars,Covars,ModelFI,ModelS,RandNums,MaxDays,MaxPlots,DataInterval,FractDead);


	CrashOut += test;
        //fp = fopen("CrashOut.dat","a");
	//fprintf(fp,"%d\t%d\t%d\n",ExtRlzns,CrashOut,test);
	//fclose(fp);
	ExtRlzns++;




	ax = 1e-20;
	bx = 200;
	cx = 0.15;

	WhichData = 1;
        MaxPlots = AdjPars[0];
	
/*
	if((DataSet!=2)&&(DataSet<12)){
		test = SSECalc(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,SSE,WhichData);
        	xmin = AdjPars[10];
		answer = LHoodFunc2(MaxRlzns,MaxPlots,MaxWk,SSE,xmin,lhoodTemp1);
	} else {
//for straight-up binomial, use just the next line
		//answer = LhoodFuncBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1); //straight up binomial 
//for the pseudo-likelihood, just the next TWO lines
		//test = SSECalcPseudo(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,SSE); //Pseudo-likelihood binomial 
		//answer =LHoodFuncPseudo(MaxRlzns,MaxPlots,MaxWk,Data,SSE,AdjPars[10],AdjPars[20],lhoodTemp1); //Pseudo-lhood 
//for beta-binomial, use the next line
  		answer = LhoodFuncBetaBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1,AdjPars[20],AdjPars[10]); //straight up binomial, 
//20 is gamma, 10 is obsvn error, which in this case is "delta"
	}
*/

	answer = LhoodFuncBetaBnml(MaxRlzns,MaxPlots,MaxWk,Data,ModelFI,lhoodTemp1,AdjPars[20],AdjPars[10]); 



	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
		lhoodVec[Rlzn] = lhoodTemp1[Rlzn];
		//if(SData==1) lhoodVec[Rlzn] += lhoodTemp2[Rlzn];
	}		



		
	free_dvector(FxdPars,1,100); //check
	free_vector(S0,1,100); //check
	free_vector(FractI,1,100); //check
	free_f3tensor(Data,0,10,0,20,0,50);
	free_matrix(Covars,1,25,1,25); //check
	free_ivector(Split,1,100); //check
	free_ivector(MaxDays,1,100); //check
	free_f3tensor(ModelFI,1,MaxRlzns,1,25,0,25); //check
	free_f3tensor(ModelS,1,MaxRlzns,1,25,0,25); //check
	free_f3tensor(SSE,1,MaxRlzns,1,25,0,25); //check
	free_dvector(lhoodTemp1,1,MaxRlzns); //check
	free_dvector(lhoodTemp2,1,MaxRlzns); //check
	free_d3tensor(RandNums,1,MaxRlzns,1,25,0,50); //check
	

	//if(lhoodVec[1]>-250)
		return exp(lhoodVec[1]);
	//else
		//return exp(-250);
} //LHoodRK45



int GetDim(int DataSet, int Plot){
	
        float S0[30], FractI[30];
	float ***Data;
	int  Split[30];
	int  MaxDays[30];	
	int MaxWk, DataInterval;
	int SData, FractDead;
	int MaxPlots;

	Data = f3tensor(0,10,0,20,0,50);

	
	switch(DataSet){
			case 1 : MaxPlots = swplots5(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 2 : MaxPlots = swplots5Bnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 3 : MaxPlots = swplots8Bnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;

			//case 3:  MaxPlots = Mandarte(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 4:  MaxPlots = ShepherdControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 5:  MaxPlots = ShepherdAll(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 6:  MaxPlots = ShepherdSpray(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 7:  MaxPlots = ShepherdAerial(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 8:  MaxPlots = OtvosControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 9:  MaxPlots = OtvosAll(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 10:  MaxPlots = OtvosSpray(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 11:  MaxPlots = MoreauControl(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); break;
			case 12:  MaxPlots = OtvosConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 13:  MaxPlots = OtvosShepConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
			
break;
			case 14:  MaxPlots = MoreauConBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
			
break;
			case 15:  MaxPlots = OtvosShepConBnml2(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
			
break;
			case 16:  MaxPlots = OtvosAllBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 17:  MaxPlots = KPBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 18:  MaxPlots = OtvosShepAllBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;
			case 19:  MaxPlots = OtvosShepTmtBnml(S0,FractI,Data,Split,MaxDays,&MaxWk,&DataInterval,&SData,&FractDead); 
break;


	}

	free_f3tensor(Data,0,10,0,20,0,50);
       return MaxDays[Plot];
}




int main(void){


//int main(int argc, char **argv){


int iii;
	
  for(iii=0;iii<20;iii++)
	FractI2[iii] = 0.0;

 if((DataSet==2)||(DataSet==3)){
  FractI2[1] = 0.286;  //Fraction hatching infected, GM data, plot 1, 1983, second week
  //FractI2[2] = 0.0885;
  FractI2[4] = 0.134;  //Fraction hatching infected, GM data, plot 10, 1984, second week
 }

  FILE *fp, *fp2, *fp3;
int ParamChoice2, ParamChoice3;
 //float LHoodStor[100][5000];
 int LHCount,MaxLHCount;
 float BestParamStor[25];
 int Trial;
 int Stop;


 ParamStruct ParamVals;
 double RandNums[100];
 int i;
 //double test;
 double nubar, sigma;
 int Rlzn, MaxRlzns = 1;
 double AvgLHood = 0;
 size_t calls = 1;
 double Gaussian,dummy;
 double LHoodPerPop, OldLHoodPerPop, NewOldLHoodPerPop;
 double temp, NewPost,OldPost,NewPropAdj,OldPropAdj,Criterion;
 int Cancel, Cancel2;

 double LHoodPerPop2;
 int Pop;
 int ParamChoice;
 double Itn, MaxItn;
 
 double FxdPars[30],CurrentPars[30];

 int Profile = 0;

//ParamVals.FxdPars[0] = 5; //the Plot Number
ParamVals.FxdPars[1] = 0.25; //2.0e6; //0.25; //2.0e6;  //k
double kHtg = ParamVals.FxdPars[1];
 //ParamVals.FxdPars[2] = exp(-0.25); //M: -1.2; nubar - see note below
 ParamVals.FxdPars[2] = 0.67732; //0.89811; //0.75323;
//REMEMBER: the Monte Carlo integration code is drawing values of what is really log(\bar{\nu}).  
//As a result, the parameter written as "nubar" must range from negative infinity to positive infinity. 
// This is not the same situation as k, or as
//the FI's, because THOSE parameters are positive as far as the MC integration code and DDE7 are concerned, they're only negative for 
//the purposes of searching a grid.
 nubar = ParamVals.FxdPars[2];
 double mu = 0.6;
 //mu = 0.6;
int ratioFix = 0;
double ratioFixVal = 1.0; //exp(-2.78);
int mFix = 0;
double mFixVal = 100; //exp(3.226452); //100
int deltaFix = 0;
double deltaFixVal = exp(-2.665037);

  //Boostrap of KP Cool				
  //Avg log 1/delta: -2.665037 SD log 1/delta: 1.781641e-05 (actually it must be delta itself);				
  //Avg log m: 3.226452 SD log m: 0.0001195957 	

 double InitScale = 1.25, Inflate;

 ParamVals.FxdPars[3] = mu; //0.39268; //0.38678; //M:0.525; //mu

 ParamVals.FxdPars[4] = 0.25; //0.25; //0.0; //0.3; //M:0.00398; //sigma
 sigma = ParamVals.FxdPars[4];

 ParamVals.FxdPars[5] = 0.0; //0.01303; //0.0; //0.01303; //0.71616; //1.17653; //2.75; //-2; //RATIO KEEP THIS AT -2 FOR OTVOS/SHEP AND FOR MOREAU  //in fact, this has been at -2.8 for O/S and M
 double ratio = 0.013; //-1; //0.013;
 
ParamVals.FxdPars[6] = -2; //0.06417; //0.07794; //1.0/12.0; //(10.0/12.0); //M:15; //tau, speed of kill
 ParamVals.FxdPars[7] = 0.1; //0.1; //0.01; //0.25; //h

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


//NumPops[1] = NumPops[2] = 5; //Woods and Elk
//NumPops[1] = NumPops[2] = 8; //Woods and Elk TEMP
NumPops[1] = NumPops[2] = 5; //Woods and Elk TEMP
NumPops[3] = 8;

NumPops[13] = 5; //Otvos/Shep control pseudo/bnml/beta bnml
NumPops[14] = 3; //Moreau control pseudo/bnml/beta bnml
NumPops[15] = 6; //Otvos-Shep fixed up 
NumPops[16] = 7; //this is Otvos All Con Bnml
NumPops[17] = 3;
NumPops[18] = 15;
NumPops[19] = 10;

ParamCount[3] = ParamCount[2] = 8; //Woods and Elk, beta binomial, obsvn fixed
//ParamCount[13] = 11; //Otvos/Shep Control beta binomial, no ratio, no obsvn error OLD
ParamCount[13] = 8; //Otvos/Shep Control beta binomial
ParamCount[14] = 9; //Otvos/Shep Control beta binomial, no ratio, no obsvn error
ParamCount[15] = 12; //Otvos/Shep Control beta binomial, no ratio, no obsvn error, fixed up ICK!!
//ParamCount[16] = 7; ////this is Otvos All Con Bnml - This was when there was no ratio 
ParamCount[16] = 8;
ParamCount[17] = 8;
ParamCount[18] = 8;
ParamCount[19] = 8;



 int PC, PC2;
 int ParamNum[30]; //These are only for looping through the parameters

 for(PC=1;PC<=16;PC++)
 	ParamNum[PC] = PC; //GM bnml, no stage

//GM beta binomial
if((DataSet==2)||(DataSet==3)){
	ParamNum[7] = 9;
	ParamNum[8] = 20;
}

/*For Otvos/Shep Con beta binomial - NO ratio, no obsvn error */
/*
if(DataSet==13){
 ParamNum[5] = 6;
 ParamNum[6] = 12;
 ParamNum[7] = 13;
 ParamNum[8] = 14;
 ParamNum[9] = 15;
 ParamNum[10] = 16;
 ParamNum[11] = 20;
}
*/

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
/* It used to be that ratio was not used in OS, but that's no longer true
if(ParamVals.FxdPars[11]>=15){ //15 is the data set number
 ParamNum[5] = 6;
 ParamNum[6] = 9;
 ParamNum[7] = 20;
}
*/

if(ParamVals.FxdPars[11]==13){ //This uses ratio, 9 is OtvosAll, 13 is OtvosShepConAll
 ParamNum[7] = 9;
 ParamNum[8] = 20;
}


if(ParamVals.FxdPars[11]>=15){ //This uses ratio
 ParamNum[7] = 9;
 ParamNum[8] = 20;
}

double res, res2, err, oldres;
double AvgL[50], Var[50];
     
       double xl[300];
       double xu[300];
     
       //const gsl_rng_type *T; //don't seem to be used...
       //gsl_rng *r, *pr;

float lambda;

//T = gsl_rng_default;
//r = gsl_rng_alloc (T);  //These 2 really aren't used, screw'em


/* Start here for MCMC*/
Verbose = 0;
Verbose2 = 0;
//calls = 100;//50; //this is now data-set specific, below

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

/****************** Priors and Proposals *********************/
double PropMean[30], PropVar[30], PriorMean[30], PriorVar[30], PriorUnif[30], InitFlag[30];
double PropMean2[30], PropVar2[30], MixFract[30], MixRand;
int ThinStop;

  double MaxParam[30];
  MaxParam[1] = 1e3; //Narrow is 1e3, else 1e6; //k
  MaxParam[2] = 1e2; //Narrow is 1e2, else 1e5 //50.0; //nubar
  MaxParam[3] = 1e2; //Narrow is 1e2, else 1e5 //100.0; //10.0; //mu
  MaxParam[4] = 1e2; //Narrow is 1e2, else 1e4 //10.0; //sigma
  MaxParam[5] = 1e2; //Narrow is 1e2, else 1e3 //10.0; //ratio
  MaxParam[6] = 1; //Narrow is 1, else 1e2 //5.0; //0.5; // 1/SOK
  MaxParam[9] = 1e2;  //Narrow is 1e2, else 1e4 //100.0; //m
  MaxParam[20] = 1e3;


 if((DataSet==2)||(DataSet==3)){


  MaxItn = 1e9;
  calls = 100; //250; //25; //50; 
  ThinStop = 100;
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


//So far, k's prior is always uninformative
//No longer true as of 21 August 2014

// k
//Here are the proposals

   PropMean[1] = 2; //2; //7 for uninformative; //1; //k on a log scale
   PropVar[1] =  1.5; //2.0; // 1.5; 2.5; //MN2A: 5.0 OFN: 1.5; //OFM: 1;
   PropDist[1] = 2; //normal-log

   
   //k Uninformative
/*
   PriorDist[1] = 1;
   PriorVar[1] = (MaxParam[1]);
*/



//k Informative

   PropDist[1] = 2; //normal-log
   PriorDist[1] = 2;
   PriorMean[1] = -0.320; //0.263; //-0.320; //0.263; //-0.831;// //-0.821784; //Data are in ROWO03ANoConB.csv, code is BootROWO.R
   PriorVar[1] =  0.384; //2.01; //0.384; //2.01; //0.2536864;

   
   PropMean[1] = PriorMean[1]; //7 for uninformative; //1; //k on a log scale
   PropVar[1] =  PriorVar[1];  //2.5; //MN2A: 5.0 OFN: 1.5; 

//log k: -0.8290467 SD k: 0.262467  ROWO03 1000 boots - Call this R
//mean log k: -0.8314351 SD k: 0.2579516  ROWO03 3000 boots
//log k: 0.262740349	2.010233473  Bret - Z? and Q
//log k: -0.319499557 sd 0.384211223 Bret, with 2 outliers removed - Y : this is from ReCalc.Nubar.OutliersGone.R  (see Manuscripts/EnvtlStochEpiMs/BretData)


			
 //nubar
//Uninformative
/*
  PropDist[2] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[2] = -5; //-4; //-3; //-2.5; //-4; //MGH: 1, cuz MGF:2.5; //use 5 for uinformative...//3; //-1; //log(FxdPars[2]);  //note that we are here nuking the log-scale for nubar
  PropVar[2] = 2.0; //1.5; //1.25;  //1.0; //1.5; //1.25; //0.25; //0.4; //1.0; //0.5; //1.5; //0.5; //1.5;
  PriorDist[2] = 1; 
  PriorVar[2] = (MaxParam[2]); //3.0; //GMOFP: 1e1;
*/


 //Informative 

  PriorDist[2] = 2; 
  PriorMean[2] = -0.775; //-1.687775; //Bret Classic //-0.776; //Bret htg
  PriorVar[2] = 0.538; //0.3116; //Bret Classic // 0.538; //Bret htg
  PropDist[2] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[2] =  PriorMean[2]; //0.22;  This is from ROWO3C //-4; //MGH: 1, cuz MGF:2.5; //use 5 for uinformative...//3; //-1; //log(FxdPars[2]);  //note that we are here nuking the log-scale for nubar
  PropVar[2] = PriorVar[2]; //0.38; This is from ROWO3C //1.5; //1.25; //0.25; //0.4; //1.0; //0.5; //1.5; //0.5; //1.5;

//-1.687775 SD: 0.3115552  //ONLY WITH k big!!  These are from ReCalc.NubarClassic.OutliersGone.R (see Manuscripts/EnvtlStochEpiMs/BretData)

//mean log nu: 0.2334692 SD nu: 0.3958792  ROWO
//mean log nu: 0.236289 SD nu: 0.3868806 ROO03 3000 boots
//mean: -0.850658 SD: 0.574246   Bret
//mean: -0.7757777 SD: 0.5376314 Bret with 2 outliers removed
 
  
//nubar mean with C = 0:  -1.686655
//> sd(log(nu.barNok))/sqrt(length(nu.barNok)-1)
//[1] 0.06259837


/* Next set is for Mix Fract only */
  PropMean2[2] = -4.5; 
  PropVar2[2] = 0.2; 
  MixFract[2] = 0.25;
 

//mu
//Uninformative Prior

  PriorVar[3] = 10000; //100; //Uniform Prior 
  PriorDist[3] = 1;
  PropDist[3] =  2; //5; //Double normal is 5
  PropMean[3] = -1; //0.0; //-1; //-3.0; //0.89; //-1.25; //-0.89; //-0.89; //-1.25; //-4; //log(0.3); //log(0.41); //SOK prior's only (MGH) 2.0; //MGF: log(20); //uninf: log(100); //log(0.41);
  PropVar[3] =  (MaxParam[3]); //1.25; //1.5; //2.0 //1.5; //1.25; //0.75; //0.5; //0.5; //1.25; //0.5; //1.25; 


//Informative Emma Prior
/*
  PriorDist[3] = 2;
  PriorMean[3] = log(0.41);
  //PriorVar[3] = 0.52; //95% CI on Emma's data is 0.24 to 0.51, a range of 0.27.  sqrt(0.27) is 0.52, roughly THIS IS BASED ON THE NON-LOG SCALE (it's wrong)
  PriorVar[3] = 0.235; //this is based on the log scale
  PropDist[3] =  2; //5; //Double normal is 5
  PropMean[3] = log(0.41); //-3; //0.89; //-1.25; //-0.89; //-0.89; //-1.25; //-4; //log(0.3); //log(0.41); //SOK prior's only (MGH) 2.0; //MGF: log(20); //uninf: log(100); //log(0.41);
  PropVar[3] = 0.235; //2.0; //1.5; //1.25; //0.75; //0.5; //0.5; //1.25; //0.5; //1.25; 
*/

  //Uninformative sigma - log-normal proposal
  PropDist[4] = 2; // //sort of Gaussian/log, actually //4 is exponential
  PropMean[4] = -1; //0.0; //-1; //-3.0; //R: -5 GMOFQ: -2.5 GMOFL: log(0.25); //GMOFI:5e-3; //2.5; //this is for sigma s
  PropVar[4] = 0.25; //0.2; //0.5; //0.25; //0.5; //P thru S; GMOFL: 1.0; 
    
  //Double-normal or half-normal
  //PropDist[4] = 5; //Double-normal, or half-normal, they're the same really
  //PropVar[4] = 0.25; //0.75; //0.5;
     
  PriorDist[4] = 1;
  PriorVar[4] = (MaxParam[4]); //0.1; //GMOFP: 0.5 //GMOFO and earlier: 1000.0;
  PriorUnif[4] = 1; //no longer in use?

//Ratio

//Next is ratio

//Uninformative prior - ratio

  PropMean[5] = -4; //-3; //-4; //-3; //-4; //log(FxdPars[5]);  //ratio?
  PropVar[5] = 0.05; //0.25; //0.05; //0.025; //0.05; //0.05; //0.025; //0.05; //0.025; //0.05; //0.2; //0.05; //0.5; // GMOFN: 0.2 GMOFI: 0.2 //0.005;
  PropDist[5] = 2;
  PriorDist[5] = 1; //Uniform
  PriorVar[5] = (MaxParam[5]);



//Informative prior - ratio
/*
  PropMean[5] = -4; //log(FxdPars[5]);  //ratio?
  PropVar[5] = 0.1; //0.04; //0.025; //0.05; //0.025; //0.05; //0.2; //0.05; //0.5; // GMOFN: 0.2 GMOFI: 0.2 //0.005;
  PropDist[5] = 2;
  PriorDist[5] = 2;
  PriorMean[5] = -4.022;  //-2.295416836; //0.0102; //(1.0/10.0); //1.0/12.0; 
  PriorVar[5] = 0.0398; //there is a piece of R code that generated this, in OBCounts directory
*/


  //Uninformative prior - 1/SOK

  PriorDist[6] = 1;
  PriorVar[6] = 10.0;
  // Single proposal
  PropDist[6] = 2;
  PropMean[6] = -3.25; //-3.315; //log(0.1); //-2.75; 
  PropVar[6] = (MaxParam[6]); //0.25; //0.5; //0.1; //0.25; //0.1; //0.05; //0.012; //0.025; //0.05; //0.015; //0.25; //0.015; 


  //Informative - 1/SOK - from Arietta
/*
  PriorDist[6] = 2;
  PriorMean[6] = -2.749462039;  //-2.295416836; //0.0102; //(1.0/10.0); //1.0/12.0;
  PriorVar[6] = 0.012081358; //0.02 is a made-up value, to make the chains mix.  The real value is apparently: 0.012081358; //See: mEstimatesFromAriettasData
  PropDist[6] = 2;
  PropMean[6] = -2.75; //log(0.1); //-2.75; 
  PropVar[6] = 0.05; //0.012; //0.025; //0.05; //0.015; //0.25; //0.015;
*/


//9 is for m
   /* Neg binom */
  /*
  PropMean[9] =  5500; //50; //50; //15; //25; //15;  //5 if uninf prior //28; //15; //This is x in the neg binom, because this is the number of exposed stages, so it has to be an integer
  PropVar[9] = 0.995; //0.45; //0.75; //0.5; //0.25; //0.25; //0.5; //This is p in the neg binom
  PropDist[9] = 7; //neg binom
  */
  
  //Uninformative - 1/SOK

  PriorDist[9] = 1;
  PriorVar[9] = MaxParam[9];
  PropDist[9] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[9] = 3.0; //3.5; //3.0;  //2.0; //0.5; //1.8;//1.5 //2.0; //3.3; 
  PropVar[9] = 0.75; //0.4; //0.1; //0.05; //0.5; // 0.05; 

    
 //Informative from Arietta for m, number of classes
/*
  PriorDist[9] = 2; //Normal on a log scale
  PriorMean[9] = 3.305926445; //3.969886639; //3.305926455;
  PriorVar[9] = 0.048050982; //0.074451582; //0.048050982;
  PropDist[9] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[9] = 3.3;//1.5 //2.0; //3.3; 
  PropVar[9] = (MaxParam[9]); //0.1; //0.05; //0.5; // 0.05;
*/
  
  PropMean[20] = 1.0; //log(FxdPars[20]); 

  PropVar[20] = 0.4; //0.4; //0.4;
  PriorDist[20] = 1;
  PriorVar[20] = (MaxParam[20]);

 }//GM priors



 /*******************  OS from here on **********************/
 if((DataSet==13)||(DataSet>=15)){


  //ThinStop = 1000;

  if(sigma<=0){
  	ThinStop = 1000;
  }else{
	ThinStop = 100;
  }

  MaxItn = 1e9;
  calls = 100; //10; //25;
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

  
  
  PropMean[1] = 1.0; //-5.0; //1.0; //-1.0; //-3.0; //1.0; //-4; //1.0; //1.5; //A, I think: 2.0; //0.5; //0.5; //k
  PropVar[1] = 0.5; //1.0; //1.5; //0.25; //1.5; //0.1; //G: 0.1; //C: 0.25 A: 1.0?  B:0.5; //0.25; //0.4;
  PropDist[1] = 2;
  
  
  /*
  PropDist[1] = 8; //Mixed Distribution
  PropMean[1] = 0.5; //1.0
  PropVar[1] = 0.5; //1.5
  PropMean2[1] = 3.0;
  PropVar2[1] = 1.0;
  MixFract[1] = 0.15; //Fraction of time that you choose second distribution
  */
  
// Uninformative
/*
  PriorVar[1] = 1000; //1000; //for uniform case
  PriorDist[1] = 1;
*/


  //Informative - Joe's experiment - New Mexico strain
 // PriorMean[1] = -3.155;
 // PriorVar[1] = 0.435*0.435;
 // PriorDist[1] = 2;

   //Testing intermediate priors on k
  PriorMean[1] = log(10);
  PriorVar[1] = 4.0;
  PriorDist[1] = 2;

/*
k:
mean: -3.155
sd: 0.435
*/

  /*
  PropDist[2] = 8; //Mixed Distribution
  PropMean[2] = 0.0; //2.0; //1.0
  PropVar[2] = 0.75; //0.75; //1.5
  PropMean2[2] = -5;
  PropVar2[2] = 0.25;
  MixFract[2] = 0.75;
  */


  //Uninformative
 /*
  PriorVar[2] = 1000; //1000;
  PriorDist[2] = 1;
 */
 
  

//From Joe's experiment
//Overall
/*
PriorDist[2] = 2; 
PriorMean[2] = -4.87; 
PriorVar[2] = 0.776;
*/

//Informative - Joe's experiment, New Mexico isolate
//mean: -2.410
//sd: 0.996


PriorDist[2] = 2; 
PriorMean[2] = -2.410; 
PriorVar[2] = 0.996*0.996;




   

  //Informative 5th's, from Dwyer '90
/*
  PriorDist[2] = 2;
  PriorMean[2] = 0.999876151;
  PriorVar[2] = 0.627642699;
*/


  
  //Informative
  /*
  PriorDist[2] = 2; //normal
  PriorMean[2] = log(0.5); //From Dwyer '90
  PriorVar[2] = 0.4; //From Dwyer '90 (sort of)
  */

  //PriorUnif[3] = 0;
  //Uninformative
// Decay
/*
  PriorVar[3] = 100; //1000; //Uniform Prior 
  PriorDist[3] = 1; //uniform
  PropDist[3] = 2; //log Normal
  PropMean[3] = -3; //-2.5; //-1.0; 
  PropVar[3] =  0.5; //1.0; //0.1; //0.25; //0.025; // 0.1; 
*/


  //Informative

  PriorDist[3] = 2; //log-normal 
  PriorMean[3] = -1.28; //log(0.41);  //Karl's data 
  PriorVar[3] = 0.5; //0.01; //Karl's data
  PropDist[3] = 2; //log Normal
  PropMean[3] = -1.28; //-1.0; 
  PropVar[3] =  0.5; //0.1; //0.25; //0.025; // 0.1; 


//Uninformative sigma
   
  PropDist[4] = 2; //log Normal
  PropMean[4] = -2; //-1; //sigma
  PropVar[4] = 1.0; //0.75; //0.4;
  
  /*
  PropDist[4] = 8; //Mixed distribution
  PropMean[4] = -3.0;  
  PropVar[4] = 0.25;
  PropMean2[4] = 1.0;
  PropVar2[4] = 0.75;
  MixFract[4] = 0.85;
  */

  PriorDist[4] = 1;
  PriorVar[4] = 50; //24;


//Ratio
  //Proposals
  PropMean[5] = -6.5; //-4 //log(FxdPars[5]);  //ratio?
  PropVar[5] = 0.5; //0.25; //0.5; //0.01; //0.5; //0.025; 
  PropDist[5] = 2;


  /*
  PropDist[5] = 8; //Mixed distribution
  PropMean[5] = -4;  
  PropVar[5] = 0.75;
  PropMean2[5] = -6;
  PropVar2[5] = 0.25;
  MixFract[5] = 0.75;
  */

//Next is ratio

//Uninformative prior
/*
 PriorDist[5] = 1; //Uniform
 PriorVar[5] = 10;
*/


//Informative - Joe's version
PriorMean[5] = -3.387; 
PriorVar[5] = 0.867;
PriorDist[5] = 2;



 //Informative, see RatioPriorDFTM.R in GregStuff/R Code/OBCounts2009And2015
/*
  PriorDist[5] = 2;
  PriorMean[5] = -2.78; 
  PriorVar[5] = 0.869;
*/

  
  

//Proposal on 1/SOK

  //Uninformative
  
  PriorDist[6] = 1;
  PriorVar[6] = 10.0;
  PropDist[6] = 2;
  PropMean[6] = -3; //log(0.1); //-2.75; 
  //for informative priors, PropVar of 0.25 gives inf's and nans, whereas 0.01 gives no progress (?)
  PropVar[6] = 0.5; //0.23; //0.05; //0.012; //0.025; //0.05; //0.015; //0.25; //0.015; 
  

//Informative
 //Informative from DFTMSOKSummary.xls in R Code/SOKGamma
/*
  PriorDist[6] = 2;
  PriorMean[6] = -2.67;  //-2.295416836; //0.0102; //(1.0/10.0); //1.0/12.0; 
  PriorVar[6] = 1.78e-5; 
  
  PropDist[6] = 2;
  PropMean[6] = -2.67; //-2.43; //log(0.1); //-2.75; 
  //for informative priors, PropVar of 0.25 gives inf's and nans, whereas 0.01 gives no progress (?)
  PropVar[6] = 1.78e-5; //0.23; //0.05; //0.012; //0.025; //0.05; //0.015; //0.25; //0.015; 
*/
  //Boostrap of KP Cool				
  //Avg log 1/delta: -2.665037 SD log 1/delta: 1.781641e-05 				
  //Avg log m: 3.226452 SD log m: 0.0001195957 	
			





//Here is the proposal on m

/*
  PropMean[9] = 1000; //10; //25; //25; //28; //15; //This is x in the neg binom, because this is the number of exposed stages, so it has to be an integer
  PropVar[9] = 0.995; //0.5; //0.5; //This is p in the neg binom
  PropDist[9] = 7; //neg binom
*/


/*
  PriorVar[9] = 125; //75;
  PriorDist[9] = 1;
*/

  //Uninformative
  /*
  PriorDist[9] = 1;
  PriorVar[9] = 1000;
  PropDist[9] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[9] = 3.0; //1.5; //3.5; //1.5; //3.3; //1.5; //3.3; //2.0; //3.3; 
  PropVar[9] = 0.75; //0.5; //0.25; //0.1; //0.4; //0.1; //0.05; //0.5; // 0.05; 
  */
  

 //Informative from DFTMSOKSummary.xls in R Code/SOKGamma
  
  PriorDist[9] = 2; //Normal on a log scale
  PriorMean[9] = 3.23; //3.5; //3.305926445; //3.969886639; //3.305926455;
  PriorVar[9] = 1.20e-4; //0.243; //0.074451582; //0.048050982;
  PropDist[9] = 2; //8:Mixture of normal logs   2; //normal-log
  PropMean[9] = 3.5; //1.5; //3.3; //1.5; //3.3; //2.0; //3.3; 
  PropVar[9] = 0.25; //0.1; //0.4; //0.1; //0.05; //0.5; // 0.05;
  
  //Boostrap of KP Cool				
  //Avg log 1/delta: -2.665037 SD log 1/delta: 1.781641e-05 				
  //Avg log m: 3.226452 SD log m: 0.0001195957 	 
  

 

  // PropMean[17] = 0.01;
  

  PropMean[20] = -1; //1; //-2; //log(FxdPars[20]); 
  PropVar[20] = 1.5; //0.4; //0.15; //0.4; //0.4; //0.4;
  //PropDist[20] = 2; 

  PriorDist[20] = 1;
  PriorVar[20] = 10.0;

 }//OS priors and post's


//12 June 2013
//
//GMN2A is all uniformative priors, big k,sigma = 0 (it got messed up in the end); 2B is like 2A, only lower means on prop nubar and mu; 2C has even lower man on prop nubar
//GMN2D like 2C, only lower mean on prop nubar, higher prop variances all around
//GMM2A is all uninformative, but k small, sigma = 0; 
//GMN3A has informative prior on SOK, from Arietta's data; GMM3A like GMN3A but k is small; 3B has lower proposal mean on m, 15 instead of 28; 3C has higher prop var on m, otherwise like 3B; 3D is 3C
//GMM4A has Arietta prior on m instead of SOK; GMM4B like 4A, but adjusted prior on m
//GMM5A has Arietta priors on both SOK and m; 6 has prior on mu, 7 has priors on SOK, m and mu
//P means sigma > 0, k big
//Most OS are actually ONLY Otvos
//SO in contrast is Otvos AND Shepherd, including both controls and treatments; ST is just treatments, SC is just controls
//S2O is with Shepherd controls using max density for whole season as initial density, instead of the average of the first few weeks.
//OAO is Otvos All, stochastic, small k
//ROAO is Otvos All, stochastic small k, NO RATIO


//Gx has adjusted data, 5 is the same 5 plots as in Dwyer et al. '97
//Gy is like Gx, except using lab EM data for initial fraction infected in plot 1
//G25y uses 25 calls instead of 50
//Gx5O8 has informative prior on nu, Gx5O9 has informative prior on k, Gx5O10 has informative priors on both.  Priors are from data in ROWO03ANoConB2.csv
//Gx5O11 is like 10 but with informative prior on mu as well, while 12 also has an informative prior on ratio
//Gx5O5JB is Gx5O5J but on long queue


//Gx5O14XA uses nubar values from Elderd et al. 08 as priors, and so on for other X's
//G2 is 250 calls, G3 is 500 calls, G4 is 100 calls
//G3 is 400 calls starting 16 Dec 2015
//y has a bug fixed: if(ratio<=5) NoGo for prior calc (ugh)
//Y has a fixed delay of m = 250, Z is m = 100
//No, Y is Bret's bugs in bags data, less two outliers
//R after the number means priors from RoWo03A, Q is from Bret, Z is from Bret with fixed m = 100, Y is from Bret with fixed m = 100, and two outliers removed from k data
//Gz8 and GzLS8 have fast_odesGD.h with a fix: no more setting initial pathogen density to its initial value if it exceeds the initial value of host density
//31 is nubar, mu, m and delta; //32 has m and delta fixed at KPCool values, nubar mu informative; 33 has m and delta fixed at KPCool values, nubar, mu and ratio informative
//O + Z is actually Karl's SOK, fixed
//OR is Otvos all, but with the initial host density calculated for the first week for which there is a fraction infected value, and ratio = 1 for controls
//GA8 has a fix.  nubar was getting exponentiated 2X, which in this version is fixed.  Search on Pinv from top of file to find the key lines in FigureGenerate and LHoodRK45.  50 calls for LS, 100 for MCMC

//1 is all uninformative
//2 is uninformative priors, 3 has only m?  4 has only SOK?  5 has both m and SOK, 7 has m, mu and SOK, 10 has m, SOK, nubar and C; 9 like 10 but informative on C not nubar; 8 like 9 but informative on nubar not C
//15 has informative priors ONLY on k and nubar, 14 is informative on nu not k, 13 is informative on k not nu, 18 is informative on k, nu, mu and ratio, 17 is informative k, nu, mu; 16 is nu and mu
//11 is mu only
//GA8O23ZA accidentally had a prior only on ratio, let's say 34 is ratio only
//19 is k and mu, 20 is nu, mu and delta, 22 is m, delta and mu, 21 is delta and mu
//22XSF has informative on mu, 1/SOK (delta) and m; 21 has delta, mu and ratio; 23 has informative on mu and ratio; 24 has informative on mu and k; 25 has informative m, k and mu; 30 is everything informative but nuisances params
//29 is k, nubar, mu and ratio informative
//26 has informative on mu, 1/SOK, m and ratio; 27 is mu, nubar, k; 28 is k, nubar, ratio

//35 is informative on ratio, mu, m, and delta (not k, nubar, sigma, gamma); 36 is informative nubar, mu and ratio; 37 is informative k, nubar, mu and ratio; 38 is nu and ratio; 39 is k, mu, ratio and delta; 
//40 is k and ratio (and m,sorta), 41 is k, mu, ratio, m fixed
//R means ratio is fixed at 1; RZ has ratio at 1, m = 100; RY has ratio at 1, m = 10;
//J means priors from Joe; T1 is testing a prior on k with mean 10, var 4
//Y is Bret transmission data with 2 outliers removed.   YA is VarScale = 1.1, YB is VarScale = 1.5
//YC is 2 Dec 2019, with VarScale = 1.2
//GA8O29YC is all informative except for ratio, 1/SOK, and m is fixed
//GA8O1YC is all uninformative, again m fixed
//GA8O42YC is all informative, and not even m fixed
//GA8O0YC is all uninformative, and not even m fixed
//GA8P0YC is high k, all uninformative, m not fixed
//GA8O5YC is low k, informative on SOK and m only, m not fixed
//The Y series seems to have been with incorrect priors on nubar, Z therefore has the proper Bret priors
//GA8O15YC is informative on k and nubar only, m not fixed
//GA8O42ZC is all informative, and not even m fixed
//GA8O15ZD is GA8O15ZC but with narrower uniform priors.


	char testbuff[128];
	char param_3[100];   
	
	char *test[] =
	  {
	    "GA8O15ZD", //"GA8O42YC", //"GA8P16YB",  //"GA8P11RYA", //"OzO32ZA", //"O4P11A", //"G4z8O29QC", //"Gy8P22XK", //"OAO1C", Gx5O13XA //"Gx5P16A", //"OAPC5J", //"KPM7B", //"OSM7I", //"G8O9H", //"GMOFZ", //"GMOFP", //"OSB", //"GMMFJ", //"GMMFI", //"GMWB",
	    "GA8O15ZD", //"GLSA8O42YC", //"GLSA8P16YB", //"GLSA8P11RYA", //"OLSP11B",//"GLS8O30RA", // Line search
	    "Null", //"GMPWB", //strips off L-hood scores, handy for looking at Bayesian results - OBSOLETE
	    "GA8O15ZD", //"GA8O42YC", //seeds and priors
//"GA8P16YB", //"GA8P11RYA", //"OzO32ZA", //"OLSP11B", //"GLS8O30RASeed", //"Gy8P22XKSeed", //"Gy8O14XASeed", //"OAO1CSeed", //"Gx5P16ASeed", //"OAPC5JSeed", //"KPM7BSeed", //"OSM7ISeed", //"G8O9HSeed"
	  };

    char testbuff2[128];
    char buffer0[128],buffer1[128],buffer2[128],buffer3[128];
    char bufferB[128];
    char fname;
    int pid;
    pid = getpid();


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


		 			
			//Write proposal distributions etc to a file
			if(LineSearch!=1){
				strcpy(testbuff,test[3]);
		  		fp2 = fopen(testbuff,"w");
				fprintf(fp2,"******************************************************************************\n");
				fprintf(fp2,"Seed:%d Calls:%d mFix:%d\n",seed,calls,mFix);
				fprintf(fp2,"DataSet:%d Seed:%d sigma:%f calls:%d InitScale:%f kHtg:%f\n\n\n",DataSet,seed,sigma,calls,InitScale,kHtg);
				fprintf(fp2,"PC   ParamNum PropDist   PropMean       PropVar         PropMean2            PropVar2	PriorDist   PriorMean       PriorVar    \n"); 
				for(PC=1;PC<=ParamCount[DataSet];PC++){
					ParamChoice = ParamNum[PC];
					fprintf(fp2,"%d     \t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\n",PC,ParamNum[PC],PropDist[ParamChoice],PropMean[ParamChoice],PropVar[ParamChoice],PropMean2[ParamChoice],PropVar2[ParamChoice],MixFract[ParamChoice],PriorDist[ParamChoice],PriorMean[ParamChoice],PriorVar[ParamChoice]);
				}
			
			

				fprintf(fp2,"*******************************************************************************\n");
				fclose(fp2);
			}


			

 



double Dbar;
//***************************Just for figure generation

//printf("Just before if(Figures==1)...\n");
if(Figures==1){

ParamVals.FxdPars[0] = Pop; //the Plot Number

FxdPars[7] = 0.0005; //It's easier to have this line in the actual function (no it turns out not...) - this is the time step for ODE solver
fp = fopen("FigureOutSEIR.dat","a");
fclose(fp);

//Next 2 lines are obsolete
//fp = fopen("MeanGx5O5J.dat","r");
//fp = fopen("MeanGx5P5H.dat","r"); 
fp = fopen("FigGenParams.dat","r"); //Read in the parameters here - really this one is MeanGx5O5J.dat



fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20],&Dbar);
printf("k:%f nu:%f mu:%f sigma:%f ratio:%f m:%f delta:%f DBar:%f \n",ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],Dbar);
//printf("Hit enter to proceed...\n");
//getc(stdin);

for(i=1;i<=20;i++){
	 
	 FxdPars[i] = ParamVals.FxdPars[i];
	 
}

printf("Dstar:%e\n",Dbar);
printf("ratio?:%f\n",ParamVals.FxdPars[5]);
fclose(fp);


       				dim = 100; calls = 1;
					
					res = FigureGenerate(RandNums,dim,FxdPars);

printf("And now, exit...\n");
exit(1);

}

//**************** Figure generation end

///******************* Just for WAIC calculation **********************////



if(WAIC==1){ //Turn on WAIC calculation

//mFix = 0;
printf("Doing WAIC: mFix = %d \n",mFix);// getc(stdin);

FILE *fp4;


	

	int weekNum;
	for(Pop=0;Pop<30;Pop++){
		for(weekNum=0;weekNum<100;weekNum++){

			PointwiseLH[Pop][weekNum] = -1.0e10;
			PointwiseVar[Pop][weekNum] = -1.0e10;
			PointwiseMean[Pop][weekNum] = -1.0e10;

		}


	}

	double *kIn, *nuIn, *muIn, *sigmaIn, *ratioIn, *deltaIn, *mIn, *Error;
	const NumSets = 10000000;
	kIn = dvector(0,NumSets);
	nuIn = dvector(0,NumSets);
	muIn = dvector(0,NumSets);
	sigmaIn  = dvector(0,NumSets);
	ratioIn = dvector(0,NumSets);
	deltaIn = dvector(0,NumSets);
	mIn = dvector(0,NumSets);
	Error = dvector(0,NumSets);

/*************** Change file name here ******************/
	//fp = fopen("GA8O41YA.in","r");

	//fp = fopen("GA8P23ZA.in","r");

	//fp = fopen("GA8O1YC.in","r");
	//fp = fopen("GA8O42YC.in","r");
	//fp = fopen("GA8O5YC.in","r");
	//fp = fopen("GA8P0YC.in","r");
	//fp = fopen("GA8O15ZC.in","r");
	//fp = fopen("GA8O42ZC.in","r");
	fp = fopen("GA8O15ZC4.in","r");
	//fp = fopen("GA8O15ZD4.in","r");

	i = 0;
	long MaxSEIRParamNum, SEIRParamNum = 0;

	
	if(mFix==0){
	
		if(kHtg<2e5){
			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("Small k %d %f %f %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],mIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			

				SEIRParamNum++;	
			}
		}else{

			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&mIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("Big k %d %f %f %f %f %f %f %f\n",SEIRParamNum,nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],mIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			

				SEIRParamNum++;	
			}



		}

	} else{

		printf("Just before while statement...\n");


		if(kHtg<2e5){
			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf \n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("You are here Small k %d %f %f %f %f %f %f  %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			

				SEIRParamNum++;	
			}
		}else{

			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf \n",&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("Big k %d %f %f %f %f %f %f\n",SEIRParamNum,nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);
			

				SEIRParamNum++;	
			}



		}

		/*

		if(deltaFix==0){
			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&deltaIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("%d %f %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);

				SEIRParamNum++;	
			}
		}else{

			while(fscanf(fp,"%lf %lf %lf %lf %lf %lf\n",&kIn[SEIRParamNum],&nuIn[SEIRParamNum],&muIn[SEIRParamNum],&sigmaIn[SEIRParamNum],&ratioIn[SEIRParamNum],&Error[SEIRParamNum])!=EOF){
			
				printf("%d %f %f %f %f %f %f %f\n",SEIRParamNum,kIn[SEIRParamNum],nuIn[SEIRParamNum],muIn[SEIRParamNum],sigmaIn[SEIRParamNum],ratioIn[SEIRParamNum],deltaIn[SEIRParamNum],Error[SEIRParamNum]); 
				//getc(stdin);

				SEIRParamNum++;	
			}


		}
		*/
	}
	

	//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[20]);;


		
	fclose(fp);

	MaxSEIRParamNum = SEIRParamNum;
	printf("SEIRParamNum:%d MaxSEIRParamNum:%d \n",SEIRParamNum,MaxSEIRParamNum); getc(stdin);

	//exit(1);



for(i=1;i<=20;i++){
	 
	 FxdPars[i] = ParamVals.FxdPars[i];
	 
}


		calls = 5000;


 
		double PopWiseMean[100], PopWiseM2Log[100], PopWiseVarLog[100], PopWiseMeanLog[100];
		for(i=0;i<100;i++){
			PopWiseMean[i] = 0.0;
			PopWiseM2Log[i] = 0.0;
			PopWiseMeanLog[i] = 0.0;
			PopWiseVarLog[i] = 0.0; //unnecessary?
		}
		

		for(SEIRParamNum=0;SEIRParamNum<MaxSEIRParamNum;SEIRParamNum++){

			if(kHtg<2e5){
				ParamVals.FxdPars[1] = kIn[SEIRParamNum];
			}else{
				ParamVals.FxdPars[1] = 2e5;
			}

			

				
			ParamVals.FxdPars[2] = nuIn[SEIRParamNum];
			ParamVals.FxdPars[3] = muIn[SEIRParamNum];
			ParamVals.FxdPars[4] = sigmaIn[SEIRParamNum];
			//ParamVals.FxdPars[5] = ratioIn[SEIRParamNum];
			//ParamVals.FxdPars[6] = deltaIn[SEIRParamNum];
			if(ratioFix==1){
				ParamVals.FxdPars[5] = ratioFixVal;
			}else{
				ParamVals.FxdPars[5] = ratioIn[SEIRParamNum];
			}
			if(deltaFix==1){
				ParamVals.FxdPars[6] = deltaFixVal;
			}else{
				ParamVals.FxdPars[6] = deltaIn[SEIRParamNum];

			}
			if(mFix==1){
				ParamVals.FxdPars[9] = mFixVal;
			}else{
				ParamVals.FxdPars[9] = mIn[SEIRParamNum];
			}
			
			ParamVals.FxdPars[20] = Error[SEIRParamNum];

			//printf("SEIRParamNum:%d k:%f nubar:%f mu:%f sigma:%e ratio:%e delta:%f m:%f error:%f \n",SEIRParamNum,ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20]);
			//getc(stdin);

			LHoodPerPop = 0.0;
			
       			AvgLHood = 0.0;
			



			for(Pop=1;Pop<=NumPops[DataSet];Pop++){
			//for(Pop=2;Pop<=NumPops[DataSet];Pop++){ //TEMP
 			//for(Pop=1;Pop<=1;Pop++){ //TEMP

				
				if(ParamVals.FxdPars[4]>0){ //TEMP
					ParamVals.FxdPars[0] = Pop; //the Plot Number

       				dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];
					//Monte function G()
      					
					//gsl_monte_function G = { &LHood, dim, &ParamVals }; //Old

					
					gsl_monte_function G = { &LHoodRK45, dim, &ParamVals };

					
			
					gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
         				gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &res, &err);

					//printf("res:%e\n",res);
					//getc(stdin);
					
				       gsl_monte_miser_free(s);
	     					

				} else { //TEMP

					FxdPars[0] = Pop;

					dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];

					//res = LHoodSigZero(dim,FxdPars);
					res = 1e50;
					
		
					AvgL[Pop] = res;
					Var[Pop] = 0.0;

				} //TEMP
			
				LHoodPerPop += log(res);
				PopWiseMean[Pop] += res;
				PopWiseMeanLog[Pop] += log(res);
				PopWiseM2Log[Pop] += log(res)*log(res);
				

				//printf("Pop:%d PopWiseMean:%f LHoodPerPop:%f\n",Pop,PopWiseMean[Pop],LHoodPerPop);

				
				//fprintf(fp,"%d\t%f\n",Pop,LHoodPerPop);
				//fclose(fp);
				//if(LHoodPerPop<-300) LHoodPerPop = -200;

     
         	//display_results ("miser", res, err);





     
       				//gsl_rng_free (r);
			} //Pop loop

			

			


		




	} //SEIRParamNum


	

	double pWAIC2B = 0.0;
	double TotPopWiseMeanLog = 0.0;
	for(Pop=1;Pop<=NumPops[DataSet];Pop++){
				pWAIC2B += (PopWiseM2Log[Pop]/SEIRParamNum) - (PopWiseMeanLog[Pop]/SEIRParamNum)*(PopWiseMeanLog[Pop]/SEIRParamNum);
				TotPopWiseMeanLog += log(PopWiseMean[Pop]/SEIRParamNum);
				printf("PopWiseM2Log - PopWiseMeanLogSqrd:%f \n",PopWiseM2Log[Pop]/SEIRParamNum - ((PopWiseMeanLog[Pop]/SEIRParamNum)*(PopWiseMeanLog[Pop]/SEIRParamNum)));
	}
	double AltWAIC = -2*TotPopWiseMeanLog + 2*pWAIC2B;

	

	

	
	double lppd = 0;
	double pWAIC2 = 0;
	
	for(Pop=0;Pop<8;Pop++){
			for(weekNum=0;weekNum<7;weekNum++){

				//PointwiseLH[Pop][weekNum] = 0.0;
				if(PointwiseLH[Pop][weekNum]>-1e10) {

					
					double LogP = (PointwiseLH[Pop][weekNum]/(calls*MaxSEIRParamNum));
					lppd += LogP;
					double dummyM = ((PointwiseMean[Pop][weekNum]/((calls*MaxSEIRParamNum-1)))*(PointwiseMean[Pop][weekNum]/((calls*MaxSEIRParamNum-1))));  //First moment, squared
					double dummyV = (PointwiseVar[Pop][weekNum]/((calls*MaxSEIRParamNum)-1)); //2nd moment
					pWAIC2 +=  (dummyV - dummyM);
					//printf("Pop:%d weekNum:%d PointwiseLH:%f PointwiseMean:%f dummyM:%f dummyV:%f\n",Pop,weekNum,PointwiseLH[Pop][weekNum],PointwiseMean[Pop][weekNum]/(calls*MaxSEIRParamNum),dummyM,dummyV);
				}
			}


	} //for Pop

	printf("lppd:%f pWAIC2:%f WAIC:%f\n",lppd,pWAIC2,-2.0*lppd+2*pWAIC2);
	

	fp = fopen("WAIC2019.dat","a");
	printf("TotPopWiseMeanLog:%f pWAIC2B:%f AltWAIC:%f  \n",TotPopWiseMeanLog,pWAIC2B,AltWAIC);
	fprintf(fp,"GA8O15ZC4 calls:%d MaxSEIRParamNum:%d TotPopWiseMeanLog:%f pWAIC2B:%f AltWAIC:%f lppd:%f pWAIC2:%f WAIC:%f \n",calls,MaxSEIRParamNum,TotPopWiseMeanLog,pWAIC2B,AltWAIC,lppd,pWAIC2,WAIC);
	fclose(fp);
	exit(1);
	//printf("hit return to continue...\n");
	//getc(stdin);


		//int weekNum;
	

	double pD = -2.*Dbar - (-2*LHoodPerPop);
	double DIC = -2.0*Dbar + pD;
	printf("LHoodPerPop:%f  Dbar:%f pD:%f DIC:%f \n",LHoodPerPop,Dbar,pD,DIC);


	//fprintf(fp,"%f %f %f %f \n",LHoodPerPop,Dbar,pD,DIC);

       //	fclose(fp);

	
	free_dvector(kIn,0,NumSets); //check
	free_dvector(nuIn,0,NumSets);
	free_dvector(muIn,0,NumSets);
	free_dvector(sigmaIn,0,NumSets);
	free_dvector(ratioIn,0,NumSets);
	free_dvector(deltaIn,0,NumSets);
	free_dvector(mIn,0,NumSets);
	free_dvector(Error,0,NumSets);


	//getc(stdin);
	exit(1);
} //if(0)




///******************* End of WAIC Calculation ***********************////



//***************************Just for DIC calculation



if(DIC==1){ //Turn on DIC calculation

printf("Doing DIC:\n");

FILE *fp4;


	char testbuffA[128];
    	char bufferA[128];




	//char *test2 = "MeanGx5P2G";

	char *test2 = "BestParamsGA8P23ZA";  //NO!  THIS IS OLD!


	char *test3;
	char *strFileType = ".dat";
	strcpy(bufferA, test2);
	strcat(bufferA,strFileType); 
	test3 =  bufferA;
	strcpy(testbuffA,test3);
//fp = fopen(testbuffA,"r");

//with sigma
double mpsrf, Last;
//Next line is actual DIC
//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20],&Dbar,&mpsrf,&Last);;

double BestLHood, BestPost;
//Nex line is Models vs Data

printf("just about to read in parameters...\n");

//fp = fopen("ModelVsDataGA8P23ZA.med","r");
//fp = fopen("ModelVsDataGA8O29ZA.med","r"); //med = median

//fp = fopen("ModelVsDataGA8O42ZC.mean","r");
fp = fopen("ModelVsDataGA8P0ZC.mean","r");
//fp = fopen("ModelVsDataGA8O15ZC.mean","r");




printf("opened file successfully...\n");
int LineNum;

//fscanf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&LineNum,&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20],&BestLHood,&BestPost);

if(kHtg<1e6){
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20]);
}else{
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[4],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20]);
        ParamVals.FxdPars[1] = 2e6;
}

printf("k?:%f nu:%f mu:%f sigma:%f ratio:%f delta:%f m:%f \n",ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9]);
       getc(stdin);

//without sigma
//ParamVals.FxdPars[4] = 0.0; 
//fscanf(fp,"%lf %lf %lf %lf %lf %lf %lf %lf",&ParamVals.FxdPars[1],&ParamVals.FxdPars[2],&ParamVals.FxdPars[3],&ParamVals.FxdPars[5],&ParamVals.FxdPars[6],&ParamVals.FxdPars[9],&ParamVals.FxdPars[20],&Dbar,&mpsrf,&Last);


for(i=1;i<=20;i++){
	 
	 FxdPars[i] = ParamVals.FxdPars[i];
	 
}

printf("Dstar:%e\n",Dbar);
printf("ratio?:%f\n",ParamVals.FxdPars[5]);
fclose(fp);

	char *strFileType2 = ".out";
	strcpy(bufferA, test2);
	strcat(bufferA,strFileType2); 
	test3 =  bufferA;
	
	strcpy(testbuffA,test3);
	fp = fopen(testbuffA,"a");

	printf("k:%f nubar:%f mu:%f sigma:%e ratio:%e 1/SOK:%f MaxStages:%f error:%f Dbar:%f\n",(ParamVals.FxdPars[1]),ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],Dbar);
	fprintf(fp,"%f %f %f %f %f %f %f %f  ",(ParamVals.FxdPars[1]),ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20]);




		//printf("lambda:%f phi:%f gamma:%f k:%f nu:%f mu:%f ratio:%f delta:%f m:%f Dbar:%f\n",FxdPars[1],FxdPars[3],FxdPars[4],kIn[0],nuIn[0],muIn[0],ratioIn[0],deltaIn[0],mIn[0],Dbar);
		//fprintf(fp,"%f %f %f %f %f %f %f %f %f  ",ParamVals.FxdPars[1],ParamVals.FxdPars[3],ParamVals.FxdPars[4],kIn[0],nuIn[0],muIn[0],ratioIn[0],deltaIn[0],mIn[0]);





calls = 500;
printf("DataSet:%d NumPops:%d calls:%d \n",DataSet,NumPops[DataSet],calls); getc(stdin);

			LHoodPerPop = 0.0;
			
       			AvgLHood = 0.0;
			



			for(Pop=1;Pop<=NumPops[DataSet];Pop++){
			//for(Pop=2;Pop<=NumPops[DataSet];Pop++){ //TEMP
 			//for(Pop=1;Pop<=1;Pop++){ //TEMP


				if(ParamVals.FxdPars[4]>0){ //TEMP
					ParamVals.FxdPars[0] = Pop; //the Plot Number

       				dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];
					//Monte function G()
      					
					//gsl_monte_function G = { &LHood, dim, &ParamVals }; //Old

					
					gsl_monte_function G = { &LHoodRK45, dim, &ParamVals };

			
					gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
         				gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &res, &err);

					//printf("res:%e\n",res);
					//getc(stdin);
					
				       gsl_monte_miser_free(s);
	     					

				} else { //TEMP

					FxdPars[0] = Pop;

					dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];

					//res = LHoodSigZero(dim,FxdPars);
					res = 1e50;
					
		
					AvgL[Pop] = res;
					Var[Pop] = 0.0;

				} //TEMP
			
				LHoodPerPop += log(res);
				

				printf("Pop:%d log(res):%f LHoodPerPop:%f\n",Pop,log(res),LHoodPerPop);

				
				//fprintf(fp,"%d\t%f\n",Pop,LHoodPerPop);
				//fclose(fp);
				//if(LHoodPerPop<-300) LHoodPerPop = -200;

     
         	//display_results ("miser", res, err);





     
       				//gsl_rng_free (r);
			} //Pop loop

	double pD = -2.*Dbar - (-2*LHoodPerPop);
	double DIC = -2.0*Dbar + pD;
	printf("LHoodPerPop:%f  Dbar:%f pD:%f DIC:%f \n",LHoodPerPop,Dbar,pD,DIC);
	fprintf(fp,"%f %f %f %f \n",LHoodPerPop,Dbar,pD,DIC);

       //fprintf(fp,"%f %f %f %f %f\n",LHoodPerPop,Dbar,-2*LHoodPerPop  - 4*Dbar,mpsrf,Last);
	fclose(fp);


	//getc(stdin);
	exit(1);
} //if(0)

//******************* the above is just for DIC calculation


//************************************ Line Search **********************************//

if(LineSearch==1){


double Upper[30], Lower[30], Jump[30];
double paramtemp;
calls = 50; //250;





///////////// Gypsy Moth 

if(kHtg<2e5){
 Lower[1] = -2; Upper[1] = 2; Jump[1] = 0.5; //k 
}else{
	Lower[1] =  log(2.0e5); Upper[1] = log(2.0e5); Jump[1] = log(2.0e5); //k
}
 
 Lower[2] = -6; Upper[2] = 6.0; Jump[2] = 1.0; //nubar with small k, broad is -6 to 6 by 2
 Lower[3] = -4; Upper[3] = 3.5; Jump[3] = 1.5; //mu //Formerly -4,2,1
 Lower[4] = -2.0; Upper[4] = 3.25; Jump[4] = 0.5; //sigma
if(ratioFix!=1){
 Lower[5] = -7; Upper[5] = -1; Jump[5] = 1.0; //ratio
}else{
 Lower[5] = log(ratioFixVal); Upper[5] = 10.0; Jump[5] = 1e5;
}
 Lower[6] = -5.5; Upper[6] = -1.5; Jump[6] = 0.25; //delta

if(mFix!=1){

 	Lower[9] = 1.8; Upper[9] = 2.2; Jump[9] = 0.05; //m
}else{
 
	Lower[9] = log(mFixVal); Upper[9] = log(mFixVal); Jump[9] = log(mFixVal); //m
	//printf("Lower:%f Upper:%f Jump:%f\n",Lower[9],Upper[9],Jump[9]);
	//getc(stdin);
 
}

Lower[20] = 1.0; Upper[20] = 1.6; Jump[20] = 0.1; //gamma



////////// Tussock moth
/*

if(kHtg<2e5){
 Lower[1] = -2; Upper[1] = 5; Jump[1] = 1.0; //k
}else{
	Lower[1] =  log(2.0e5); Upper[1] = log(2.0e5); Jump[1] = log(2.0e5); //k
}

 Lower[2] = -9; Upper[2] = -2; Jump[2] = 1.0; //nubar with small k
 Lower[3] = -3; Upper[3] = 6; Jump[3] = 1.0; //mu
 Lower[4] = -3; Upper[4] = 3; Jump[4] = 0.5; //sigma
 //Lower[5] = -7.1; Upper[5] = -6.9; Jump[5] = 0.025; //ratio

if(ratioFix!=1){

 	Lower[5] = -6; Upper[5] = -1; Jump[5] = 1; //delta
}else{
 
	Lower[5] = log(deltaFixVal); Upper[5] = log(deltaFixVal); Jump[5] = 1e10; //delta
	//printf("Lower:%f Upper:%f Jump:%f\n",Lower[9],Upper[9],Jump[9]);
	//getc(stdin);
 
}

if(deltaFix!=1){

 	Lower[6] = -5; Upper[6] = -2; Jump[6] = 0.5; //delta
}else{
 
	Lower[6] = log(deltaFixVal); Upper[6] = log(deltaFixVal); Jump[6] = 1e10; //delta
	//printf("Lower:%f Upper:%f Jump:%f\n",Lower[9],Upper[9],Jump[9]);
	//getc(stdin);
 
}

if(mFix!=1){

 	Lower[9] = 0.25; Upper[9] = 5.0; Jump[9] = 1.0; //m
}else{
 
	Lower[9] = log(mFixVal); Upper[9] = log(mFixVal); Jump[9] = log(mFixVal); //m
	//printf("Lower:%f Upper:%f Jump:%f\n",Lower[9],Upper[9],Jump[9]);
	//getc(stdin);
 
}



Lower[20] = -1; Upper[20] = 2.0; Jump[20] = 0.2; //gamma

*/

/*
//EXIT Here
printf("here I am!\n");
fp = fopen("junk.dat","w");
fprintf(fp,"Huh?\n");
fclose(fp);
*/



//Informative priors require narrow ranges, so set them here
double Deflate[100];
Deflate[1] = 0.05;
Deflate[2] = 0.001;
Deflate[3] = 0.005;
Deflate[4] = 0.5; //Meaningless
Deflate[5] = 0.05;
Deflate[6] = 0.5;
Deflate[9] = 0.5; 
Deflate[20] = 0.25; //Meaningless
for(PC=1;PC<=ParamCount[DataSet];PC++){
	ParamChoice = ParamNum[PC];
	if(PriorDist[ParamChoice]>1){
		Upper[ParamChoice] = PriorMean[ParamChoice] + Deflate[ParamChoice]*PriorVar[ParamChoice];
		Lower[ParamChoice] = PriorMean[ParamChoice] -  Deflate[ParamChoice]*PriorVar[ParamChoice];
		Jump[ParamChoice] = (Upper[ParamChoice]-Lower[ParamChoice])/15.0;

		//printf("ParamChoice:%d Prior:%f Upper:%f Lower:%f Jump:%f \n",ParamChoice,PriorMean[ParamChoice],Upper[ParamChoice],Lower[ParamChoice],Jump[ParamChoice]);
		
	}
}

		
		
				strcpy(testbuff,test[3]);
		  		fp2 = fopen(testbuff,"w");
				fprintf(fp2,"******************************************************************************\n");
				fprintf(fp2,"Seed:%d Calls:%d \n",seed,calls);
				fprintf(fp2,"DataSet:%d Seed:%d sigma:%f calls:%d InitScale:%f kHtg:%f\n\n\n",DataSet,seed,sigma,calls,InitScale,kHtg);
				fprintf(fp2,"PC   ParamNum PropDist   PropMean       PropVar         PropMean2            PropVar2	PriorDist   PriorMean       PriorVar    \n"); 
				for(PC=1;PC<=ParamCount[DataSet];PC++){
					ParamChoice = ParamNum[PC];
					fprintf(fp2,"%d     \t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%f\t%f\n",PC,ParamNum[PC],PropDist[ParamChoice],PropMean[ParamChoice],PropVar[ParamChoice],PropMean2[ParamChoice],PropVar2[ParamChoice],MixFract[ParamChoice],PriorDist[ParamChoice],PriorMean[ParamChoice],PriorVar[ParamChoice]);
				}

				
			  
				fprintf(fp2,"PC   ParamNum Upper		 Lower     	Jump		Deflate\n"); 
				for(PC=1;PC<=ParamCount[DataSet];PC++){
					ParamChoice = ParamNum[PC];
					fprintf(fp2,"%d     \t%d\t%f\t%f\t%f\t%f\n",PC,ParamNum[PC],Upper[ParamChoice],Lower[ParamChoice],Jump[ParamChoice],Deflate[ParamChoice]);
				}


			
			

				fprintf(fp2,"*******************************************************************************\n");
				fclose(fp2);
				
	
//printf("After PriorMean[1]:%f\n",PriorMean[1]); getc(stdin); exit(1);		
	 

	


      //***** Initial parameters
    
double LStemp;
int LSItn, MaxLSItn = 10;
int LSflag = 0;
//calls = 100;

	for(PC=1;PC<=ParamCount[DataSet];PC++){
				
				NoGo = 0;
				ParamChoice = ParamNum[PC];

				if((ParamChoice==1)&&(kHtg>1000)){ NoGo = 1; ParamVals.FxdPars[ParamChoice] = kHtg;} 

				if((ParamChoice==5)&&(ratioFix==1)){ NoGo = 1; ParamVals.FxdPars[ParamChoice] = ratioFixVal;}

				if((ParamChoice==6)&&(deltaFix==1)){ NoGo = 1; ParamVals.FxdPars[ParamChoice] = deltaFixVal;}
                      
		      		if((ParamChoice==9)&&(mFix==1)){ NoGo = 1; ParamVals.FxdPars[ParamChoice] = mFixVal;}

				

				

			
				if(NoGo!=1){
				 
					double randtemp = gsl_rng_uniform(r2); 
					//double dummyrand =  (Upper[ParamChoice] - Lower[ParamChoice])*gsl_rng_uniform(r2); 
					LStemp = Lower[ParamChoice] + (Upper[ParamChoice] - Lower[ParamChoice])*randtemp;
					ParamVals.FxdPars[ParamChoice] = exp(LStemp); 
					//printf("ParamChoice:%d randtemp:%f Lower:%f Upper:%f ParamVals:%f \n",ParamChoice,randtemp,Lower[ParamChoice],Upper[ParamChoice],ParamVals.FxdPars[ParamChoice]); getc(stdin);
				}

				//printf("ParamChoice:%d ParamVals:%f \n",ParamChoice,ParamVals.FxdPars[ParamChoice]); //getc(stdin);

				
         
	} //for PC



//200000.000000 1.223367 131.999922 1.028017 0.036322 0.037458 4.032222 1.382479 -1.723326e+02 -1.000000e+50

//2.436363 152.622015 3.097691 0.170260 0.028121 1.604077 2.493068

/*	
ParamVals.FxdPars[1] = 200000.000000;
ParamVals.FxdPars[2] = 2.436363;
ParamVals.FxdPars[3] = 152.622015;
ParamVals.FxdPars[4] = 3.097691;
ParamVals.FxdPars[5] = 0.170260;
ParamVals.FxdPars[6] = 0.028121;
ParamVals.FxdPars[7] = 1.604077;
ParamVals.FxdPars[8] = 2.493068;
// -1.000000e+50 -1.000000e+50
*/

oldres = 1;
for(LSItn=1;LSItn<=MaxLSItn;LSItn++){
	for(PC=1;PC<=ParamCount[DataSet];PC++){
		
		ParamChoice = ParamNum[PC];

		Cancel = 0; Cancel2 = 0;
		for(LStemp=Lower[ParamChoice];LStemp<=Upper[ParamChoice];LStemp+=Jump[ParamChoice]){


	//printf("LSItn:%d PC:%d LStemp:%f k:%f nubar:%f mu:%f sigma:%f ratio:%f delta:%f m:%f gamma:%f LHoodPerPop:%f OldLHoodPerPop:%f NewPost:%f OldPost:%f\n",LSItn,PC,LStemp,ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],LHoodPerPop,OldLHoodPerPop,NewPost,OldPost);	
	//getc(stdin);	

			



        		
			 //Go back and re-calculate the likelihood again for the old parameter			
			LHoodPerPop = 0.0;
			NewOldLHoodPerPop = 0.0;			
			for(Pop=1;Pop<=NumPops[DataSet];Pop++){

									
					ParamVals.FxdPars[0] = Pop; //the Plot Number
			
					
       					dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];


					int ii;
					
					
					//Monte function G()
      					//gsl_monte_function G = { &LHood, dim, &ParamVals };
					gsl_monte_function G = { &LHoodRK45, dim, &ParamVals };

					paramtemp = ParamVals.FxdPars[ParamChoice];
					
					ParamVals.FxdPars[ParamChoice] = exp(LStemp);

					gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);
	


         				gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &res, &err);
					
					ParamVals.FxdPars[ParamChoice] = paramtemp;
					if(LSflag>0) gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &oldres, &err); //note oldres!
				       	
					gsl_monte_miser_free(s);
					
			             			

				//if(LSItn>1)
				 
				//getc(stdin); 

				
				
				//LHoodPerPop += log(res);
				//NewOldLHoodPerPop += log(oldres);

				//printf("res:%e oldres:%e LHoodPerPop:%e OldLHoodPerPop:%e NewOldLHoodPerPop:%e\n",res,oldres,LHoodPerPop,OldLHoodPerPop,NewOldLHoodPerPop); getc(stdin);

				if((isnan(res)!=0)||(isinf(res)!=0)) Cancel = 1;
                                if(Cancel!=1){
					LHoodPerPop += log(res);
				}else{
					LHoodPerPop = -1000;
				}


				//printf("LSflag:%d oldres:%e isnan:%d isinf:%d Cancel2:%d NewOldLHoodPerPop:%e \n",LSflag,log(oldres),isnan(oldres),isinf(oldres),Cancel2,NewOldLHoodPerPop);

				
				
				if(((isnan(oldres)!=0)||(isinf(oldres)!=0))&&(LSflag>0)) Cancel2 = 1;
                                if(Cancel2!=1){
					if(LSflag>0){ NewOldLHoodPerPop += log(oldres);
					}else{ NewOldLHoodPerPop = 0.0;
					}
				}else{
					NewOldLHoodPerPop = -1000;
				}
				
			       //printf("Pop:%d res:%e LHoodPerPop:%e\n",Pop,res,LHoodPerPop); getc(stdin);


				 
       				//gsl_rng_free (r);
			} //Pop
		
			
					

			NewPost = LHoodPerPop;
			if(LSflag>0){  //NOT first time through
				OldPost = OldPost - OldLHoodPerPop + NewOldLHoodPerPop;
				OldLHoodPerPop = NewOldLHoodPerPop;
			}else{
				OldLHoodPerPop = NewOldLHoodPerPop;
			}

	

		//Cancel = 0;
		
	if(Cancel<1)
		for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 

		      //printf("PC2:%d NewPost:%f\n",PC2,NewPost);
		      ParamChoice2 = ParamNum[PC2];
		      NoGo = 0;
		
		      if((PC2==1)&&(kHtg>1000)) NoGo = 1;
		      if((PC2==3)&&(mu<=0.0)) NoGo = 1;
		      if((PC2==4)&&(sigma==0.0)) NoGo = 1;
                      if((PC2==5)&&(ratio<0.0)) NoGo = 1;  //This is where it used to say <= 5
		      if((ParamChoice2==5)&&(ratioFix==1)) NoGo = 1;
		      if((ParamChoice2==6)&&(deltaFix==1)) NoGo = 1;
		      if((ParamChoice2==9)&&(mFix==1)) NoGo = 1;
		      

		
		    
			
	          		
		      if(NoGo!=1){
		 		
			
			//printf("ParamChoice2:%d ParamVals.FxdPars[ParamChoice2]:%f PriorVar[ParamChoice2]:%f \n",ParamChoice2,ParamVals.FxdPars[ParamChoice2],PriorVar[ParamChoice2]); getc(stdin);

		        float tmp2;
		        switch(PriorDist[ParamChoice2]){
				case 1 :
			  		if(log(ParamVals.FxdPars[ParamChoice2])>PriorVar[ParamChoice2]){
			       			Stop = 1;
						//printf("Stop:%d ParamChoice2:%d PriorDist:%d ParamVals:%f PriorVar:%f Upper:%f\n",Stop,ParamChoice2,PriorDist[ParamChoice2],ParamVals.FxdPars[ParamChoice2],PriorVar[ParamChoice2],Upper[ParamChoice2]); getc(stdin);
			    			NewPost = -1e50;
			  		} else {
			    			tmp2 =  gsl_ran_flat_pdf((ParamVals.FxdPars[ParamChoice2]),0,PriorVar[ParamChoice2]);
						//printf("ParamChoice2:%d tmp2:%e (ParamVals):%f PriorVar:%f\n",ParamChoice2,tmp2,(ParamVals.FxdPars[ParamChoice2]),PriorVar[ParamChoice2]); getc(stdin);
					}
			  		break;
				case 2 :tmp2 = gsl_ran_gaussian_pdf(log(ParamVals.FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // a slightly bogus kind of normal, really sort log normal
				       break;
				case 3 :tmp2 = gsl_ran_lognormal_pdf(ParamVals.FxdPars[ParamChoice2],PriorMean[ParamChoice2],PriorVar[ParamChoice2]);
				       break;
	                	case 4 :tmp2 = gsl_ran_gaussian_pdf((ParamVals.FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // more like a real normal...
	   		               break;
	                	default : printf("No prior distribution specified. Case:%d PriorDist:%d Bailing.\n",ParamChoice2,PriorDist[ParamChoice2]); getc(stdin); exit(1); break;
		      	} //End of switch

	
			//printf("ParamChoice2:%d PriorDist:%d Param:%e PriorMean:%e log(tmp2):%e NewPost:%e\n",ParamChoice2,PriorDist[ParamChoice2],log(ParamVals.FxdPars[ParamChoice2]),PriorMean[ParamChoice2],log(tmp2),NewPost); getc(stdin);
			if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
				
  				Cancel = 1;
			} else{
			    	NewPost += log(tmp2);
				
			}
			if((isnan(NewPost)!=0)||(isinf(NewPost)!=0)){
				
				Cancel = 1;
			}
			
		


		      } //NoGo
			//printf("end of PC2 loop...\n"); getc(stdin);
	
		    } //PC2
		  		

			
			if(Cancel<1){
				if(LSflag<1){ //First time through
					OldPost = NewPost;
					OldLHoodPerPop = LHoodPerPop;
					ParamVals.FxdPars[ParamChoice] = exp(LStemp);
					LSflag = 1; 
				
				} else{
				
					if(NewPost>OldPost){

		
		//printf("LSItn:%d k:%f nubar:%f mu:%f sigma:%f ratio:%f delta:%f m:%f gamma:%f LHoodPerPop:%f OldLHoodPerPop:%f NewPost:%f OldPost:%f\n",LSItn,ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],LHoodPerPop,OldLHoodPerPop,NewPost,OldPost);	
						//getc(stdin);
	


						OldPost = NewPost;
						OldLHoodPerPop = LHoodPerPop;
						
						ParamVals.FxdPars[ParamChoice] = exp(LStemp);

						
					}else{
						ParamVals.FxdPars[ParamChoice] = paramtemp;
					} //Post>BestPost
				
				} //LSFlag
			} else{ 
				LSflag = 0;
				

			} //Cancel

			
						

		
		} //LStemp
		
	}//PC
	//printf("LSItn:%d k:%f nubar:%f mu:%f sigma:%f ratio:%f delta:%f m:%f gamma:%f OldLHoodPerPop:%f OldPost:%f \n",LSItn,ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],OldLHoodPerPop,OldPost);	
	 
        


}//LSItn

printf("just before file print...\n");
printf("OldPost:%e\n",OldPost);


if((isnan(OldPost)==0)&&(isinf(OldPost)==0)){
		strcpy(testbuff,test[1]);
		//if(OldLHoodPerPop<-1e10) OldPost = -1e50;
 		fp = fopen(testbuff,"a");
		fprintf(fp,"%e %e %e %e %e %e %e %e %e %e \n",ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],OldLHoodPerPop,OldPost);	
		fclose(fp);
		//printf("%e %e %e %e %e %e %e %e %e %e \n",ParamVals.FxdPars[1],ParamVals.FxdPars[2],ParamVals.FxdPars[3],ParamVals.FxdPars[4],ParamVals.FxdPars[5],ParamVals.FxdPars[6],ParamVals.FxdPars[9],ParamVals.FxdPars[20],OldLHoodPerPop,OldPost);	
	
	}
	



		
	

	//getc(stdin);
	exit(1);
} 

//***************************** End of Line Search ***********************************//



	

//*************************** Beginning of MCMC *************************************//

 ///////////////////////////////   Reading in PCA Results ///////////////////////////

        double Coefficients[30][30];
	double Scale[30],Center[30],SD[30];
	int j;
	int ParCnt2 = ParamCount[DataSet];
	//printf("Initial ParCnt2:%d\n",ParCnt2);
	if(kHtg>1e5) ParCnt2--;
	if(ratioFix==1) ParCnt2--;
	if(deltaFix==1) ParCnt2--;
	if(mFix==1) ParCnt2--;
	

	//printf("ParCnt2:%d\n",ParCnt2); getc(stdin);


	//fp = fopen("GLSA8O29YARotations.txt", "r");
	//fp = fopen("GLSA8O1YCRotations.txt", "r");
	//fp = fopen("GLSA8O42YCRotations.txt", "r");
	//fp = fopen("GLSA8P0YCRotations.txt", "r");
	//fp = fopen("GLSA8O5YCRotations.txt", "r");
	//fp = fopen("GLSA8O42ZCRotations.txt", "r");
	//fp = fopen("GLSA8O15ZCRotations.txt", "r");
	fp = fopen("GLSA8O15ZDRotations.txt", "r");


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

	
	//fp = fopen("GLSA8O29YAScale.txt", "r");
	//fp = fopen("GLSA8O29YAScale.txt", "r"); //Should have been 801YC??
	//fp = fopen("GLSA8O42YCScale.txt", "r");
 	//fp = fopen("GLSA8P0YCScale.txt", "r");
	//fp = fopen("GLSA8O5YCScale.txt", "r");
	//fp = fopen("GLSA8O42ZCScale.txt", "r");
	//fp = fopen("GLSA8O15ZCScale.txt", "r");
	fp = fopen("GLSA8O15ZDScale.txt", "r");


	//printf("Scales...\n");
	for (i=0; i<ParCnt2; i++){
	  fscanf(fp, "%lf ", &Scale[i]);
	 // printf("i:%d Scale:%e\n", i,Scale[i]); //TEMP
	 }
	fclose(fp);
	//getc(stdin); 


	//fp = fopen("GLSA8O29YACenter.txt", "r");
	//fp = fopen("GLSA8O1YCCenter.txt", "r");
	//fp = fopen("GLSA8O42YCCenter.txt", "r");
	//fp = fopen("GLSA8P0YCCenter.txt", "r");
	//fp = fopen("GLSA8O5YCCenter.txt", "r");
	//fp = fopen("GLSA8O42ZCCenter.txt", "r");
	//fp = fopen("GLSA8O15ZCCenter.txt", "r");
	fp = fopen("GLSA8O15ZDCenter.txt", "r");



	for (i=0; i<ParCnt2; i++){
	   fscanf(fp,"%lf ", &Center[i]); //TEMP
	  //printf("i:%d Center:%lf\n", i,Center[i]);
	}
	fclose(fp);
	//getc(stdin); 


	//fp = fopen("GLSA8O29YAsd.txt", "r");
	//fp = fopen("GLSA8O1YCsd.txt", "r");
	//fp = fopen("GLSA8O42YCsd.txt", "r");
	//fp = fopen("GLSA8O5YCsd.txt", "r");
	//fp = fopen("GLSA8O42ZCsd.txt", "r");
	//fp = fopen("GLSA8O15ZCsd.txt", "r");
	fp = fopen("GLSA8O15ZDsd.txt", "r");


	 for (i=0; i<ParCnt2;i++){
	   fscanf(fp, "%lf\n", &SD[i]);
	 // printf("i:%d SD:%lf\n",i,SD[i]);
	}
	fclose(fp);
	//printf("Done...\n");
	//getc(stdin); 


	//printf("All Done with PCA, hit return...\n"); // getc(stdin); exit(1);

	if(kHtg>1e5) ParCnt2++;
	if(mFix==1) ParCnt2++;
	if(deltaFix==1) ParCnt2++;
	if(ratioFix==1) ParCnt2++;




	
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





//////////////// Initial parameter set /////////////////////


//Draw a complete set of parameters on the PCA scale
 double PCA[30], OldPCA[30];
 double VarScale = 1.2; //1.5; //1.1;
 

 int PCAParCnt = 0;
  for(i=0;i<ParCnt2;i++){
		NoGo = 0;
		
		if((i==1)&&(kHtg>=1e5)){NoGo = 1;}
		if((i==5)&&(ratioFix==1)){NoGo = 1;}
		if((i==6)&&(deltaFix==1)){NoGo = 1;} 
		if((i==7)&&(mFix==1)){NoGo = 1;} //i=7 corresponds to parameter 9, the number of classes
		if(NoGo==0){ OldPCA[PCAParCnt] = PCA[PCAParCnt] = gsl_ran_gaussian (r2, VarScale*SD[PCAParCnt]);}
		
		if(NoGo!=1){
			 //PCAParCnt allows you to skip past parameters that are not being varied
			 //printf("PCAParCnt:%d PCA:%lf SD:%f\n",PCAParCnt,PCA[PCAParCnt],SD[PCAParCnt]);
			 PCAParCnt++;
		}
  }


ParCnt2 = PCAParCnt;
//printf("PCAParCnt:%d\n",PCAParCnt);
//getc(stdin);



/*
First, draw proposals, put them in CurrentPars, FxdPars and ParamVals.FxdPars, and do it with a bigger variance.  CurrentPars is the current parameter set, FxdPars is the proposed parameter set, ParamVals.FxdPars is the 
proposed parameter set that 
gets passed to the monte function wrapper G.  Then in the first iteration, DO NOT propose any new values, just calculate the likelihood, the priors, and the posterior.  In some sense, the effect is that the first parameter
is not really proposed, but I believe that the code still calculates the adjustment to the posterior that is based on the proposal distribution.  Anyway, the effect of this first iteration is that it produces a posterior.
Then at the next iteration, you actually propose a new parameter value, NOT with a bigger variance, and the whole thing repeats.  THIS time, though, CurrentPars is not updated unless the new parameter set is accepted.      


*/



//Now convert the first set to the original, linear scale

LHCount = 1;
OldPost = 0;
Cancel2 = 0;
int DontLHood = 0;
//Inflate = InitScale;


double Holder[30];
 
 
 

	PCAParCnt = 0; OldPropAdj = 0.0;
 	for(PC=1;PC<=ParamCount[DataSet];PC++){
		  NoGo = 0;
		  ParamChoice = ParamNum[PC];

		  //if((ParamChoice==4)&&(sigma<=0)) {NoGo = 1; FxdPars[4] = CurrentPars[4] = ParamVals.FxdPars[4] = 0.0; }

		if((ParamChoice==1)&&(kHtg>=1e5)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = kHtg; }
		if((ParamChoice==5)&&(ratioFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = ratioFixVal; }
		if((ParamChoice==6)&&(deltaFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = deltaFixVal; }
		if((ParamChoice==9)&&(mFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = mFixVal; }


		  
		  //if((DataSet==15)&&(PC==1)) NoGo = 1;
		  
		  if(NoGo!=1){

				ParamChoice = ParamNum[PC];

				//printf("ParCnt2:%d\n",ParCnt2); getc(stdin);

				for (j=0;j<ParCnt2; j++) Holder[j] = Coefficients[PCAParCnt][j]; //Holder just holds the jth column of Coefficients ALL AT ONCE
    

				   double temp2;

		
				   double dummyDot = DotProduct(ParCnt2,Holder,PCA);
				   //printf("ParamChoice:%d\n",ParamChoice); getc(stdin);
				   temp = exp(dummyDot*VarScale*Scale[PCAParCnt]+Center[PCAParCnt]);

				  
				   if(temp<0) temp = 1e-10;


				    ParamVals.FxdPars[ParamChoice] = temp;
				    FxdPars[ParamChoice] = temp;
			            CurrentPars[ParamChoice] = temp; //Only the first time...I think...
				
				    //printf("ParamChoice:%d dummyDot:%f VarScale:%f Scale:%f Center:%f temp:%f\n",ParamChoice,dummyDot,VarScale,Scale[PCAParCnt],Center[PCAParCnt],temp);
					  					  
				    PCAParCnt++;

		  } //NoGo 
	} //PC
	//getc(stdin);

	 //200000.000000  0.009271  0.540657  0.843354  0.003188  0.064228  26.343295  3.248073  -156.427184 -77.441563

	/*
	ParamVals.FxdPars[1] = FxdPars[1] = CurrentPars[1] = 5.488e-1;
	ParamVals.FxdPars[2] = FxdPars[2] = CurrentPars[2] = 2.662;
	ParamVals.FxdPars[3] = FxdPars[3] = CurrentPars[3] = 3.339e-1;
	ParamVals.FxdPars[4] = FxdPars[4] = CurrentPars[4] = 2.7183;
	ParamVals.FxdPars[5] = FxdPars[5] = CurrentPars[5] = 2.26e-6;
	ParamVals.FxdPars[6] = FxdPars[6] = CurrentPars[6] = 6.9597e-2;
	ParamVals.FxdPars[9] = FxdPars[9] = CurrentPars[9] = 25.19;
	ParamVals.FxdPars[20] = FxdPars[20] = CurrentPars[20] = 1.8822;
	*/
	


	CrashOut = 0;
	NewOldLHoodPerPop = 0.0; 
	LHoodPerPop = 0.0;
        AvgLHood = 0.0;


			
	for(Pop=1;Pop<=NumPops[DataSet];Pop++){

				

					Cancel = 0;
					ParamVals.FxdPars[0] = Pop; //the Plot Number
								
					
       					dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];
					
					//Monte function G()
      					//First time: do it for proposed parameter values
					gsl_monte_function G = { &LHoodRK45, dim, &ParamVals };
					gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);


         				gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &oldres, &err);
				      
				 
				       gsl_monte_miser_free(s);
					
				
					
				if(oldres<=0){
					 Cancel = 1;
					 //printf("Cancel:%d res:%f\n",Cancel,res);
				}
				if((isnan(res)!=0)||(isinf(res)!=0)){
					 Cancel = 1;
					 //printf("Cancel:%d res:%f\n",Cancel,res);

				}



				

                            	if(Cancel!=1){
					NewOldLHoodPerPop += log(oldres);
				}else{
					NewOldLHoodPerPop = -300;
				}


				//printf("Pop:%d oldres:%e NewOldLHoodPerPop:%f\n",Pop,oldres,NewOldLHoodPerPop); //getc(stdin);
				 
       				//gsl_rng_free (r);
	} //Pop loop

	OldPost = NewOldLHoodPerPop;
	

	//////////// Adjusting for prior /////////////
		for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		  ParamChoice2 = ParamNum[PC2];
		      
		  NoGo = 0;
		
		  
		
	          if((PC2==1)&&(kHtg>1000)) NoGo = 1;
		  if((ParamChoice2==5)&&(ratioFix==1)) NoGo = 1;
		  if((ParamChoice2==6)&&(deltaFix==1)) NoGo = 1;
		  if((ParamChoice2==9)&&(mFix==1)) NoGo = 1;
			
		  


	          

		  //if((DataSet==15)&&(PC==1)) NoGo = 1;
		  if(NoGo!=1){
		 		
	               float tmp2;
		        switch(PriorDist[ParamChoice2]){
			case 1 :
			  if(FxdPars[ParamChoice2]>MaxParam[ParamChoice2]){
			    
			    Stop = 1;
			    NewPost = 0;

			  } else {
			    tmp2 =  (gsl_ran_flat_pdf(FxdPars[ParamChoice2],0,PriorVar[ParamChoice2]));
				


			  }
			  break;
			case 2 :tmp2 = gsl_ran_gaussian_pdf(log(FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // a slightly bogus kind of normal, really sort log normal

			
				break;

			case 3 :  tmp2 = (gsl_ran_lognormal_pdf(FxdPars[ParamChoice2],PriorMean[ParamChoice2],PriorVar[ParamChoice2]));

					 
					 break;
	

			 case 4 :tmp2 = gsl_ran_gaussian_pdf((FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // more like a real normal...

				
				 break;
	                 default : printf("No prior distribution specified.  Bailing.\n"); getc(stdin); exit(1); break;
		      	} 
			if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
				Cancel = 1;
				//printf("Cancel:%d tmp2:%f\n",Cancel,tmp2);
				tmp2 = 1e-300;
			} else{
				//printf("Inside switch OldPost:%f\n",OldPost);
			    	OldPost += log(tmp2);
				//printf("No cancel Paramchoice2:%d FxdPars[ParamChoice2]:%f log(tmp2):%f NewPost:%f\n",ParamChoice2,FxdPars[ParamChoice2],log(tmp2),OldPost);
			}
			if((isnan(OldPost)!=0)||(isinf(OldPost)!=0)){
				Cancel = 1;
				//printf("ParamChoice2:%d PriorMean:%f FxdPars[ParamChoice2]:%f Cancel:%d OldPost:%f\n",ParamChoice2,PriorMean[ParamChoice2],(FxdPars[ParamChoice2]),Cancel,OldPost);
				OldPost = -1000;
								
			}
			



		      } //ParamChoice2!=4, we're not doing sigma
	
		    } //PC2

		    //printf("OldPost:%f\n",OldPost); getc(stdin);
	
			

/////////////////  End of Initial Parameter Set //////////////////////


////////////  Start of Iterations ////////////////////

Itn = 1;

 Stop = 0;
 int Thin = 0;
 OldLHoodPerPop = 0.0;

 while(Itn<MaxItn){


		

		if((Verbose==1)&&(Itn>1)){exit(1);}



		/////////// Draw a complete set in PCA space (NO ACTUALLY JUST DRAW ONE)

	//printf("just about to draw a complete set...\n");
	int PCAParCnt = 0;
  	for(i=0;i<ParCnt2;i++){
			

		
		NoGo = 0;
			
		if((i==1)&&(kHtg>=1e5)){NoGo = 1;}
		if((i==5)&&(ratioFix==1)){NoGo = 1;}
	        if((i==6)&&(deltaFix==1)){NoGo=1;}
		if((i==7)&&(mFix==1)){NoGo = 1;}


		
		if(NoGo!=1){
			 

			PCA[PCAParCnt] = gsl_ran_gaussian (r2, SD[PCAParCnt]); 
			//printf("Itns PCAParCnt:%d PCA:%f Prob:%f \n",PCAParCnt,PCA[PCAParCnt],gsl_ran_gaussian_pdf(PCA[PCAParCnt],SD[PCAParCnt]));
			//getc(stdin);
			


		///////////// Convert complete set to original axes
		

		int PCAParCnt2 = 0;
 		for(PC=1;PC<=ParamCount[DataSet];PC++){

		 
			
		  NoGo = 0;
		  ParamChoice = ParamNum[PC];

		  

		  if((ParamChoice==1)&&(kHtg>=1e5)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = kHtg; }
		  if((ParamChoice==5)&&(ratioFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = ratioFixVal; }

		  if((ParamChoice==6)&&(deltaFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = deltaFixVal; }

		  if((ParamChoice==9)&&(mFix==1)) {NoGo = 1; FxdPars[ParamChoice] = CurrentPars[ParamChoice] = ParamVals.FxdPars[ParamChoice] = mFixVal; }


		  
		 
		  
		  if(NoGo!=1){

				ParamChoice = ParamNum[PC];

				for (j=0;j<ParCnt2; j++) Holder[j] = Coefficients[PCAParCnt2][j]; //Holder just holds the jth column of Coefficients ALL AT ONCE
    

				   double temp2;


				   /*
				   int ii;
				   for(ii=0;ii<ParCnt2;ii++)
					printf("ii:%d PCA:%f\n",ii,PCA[ii]);
				   getc(stdin);
				   */

				   //printf("ParamChoice:%d\n",ParamChoice);

				   double dummyDot = DotProduct(ParCnt2,Holder,PCA);

				   
				   temp = exp(dummyDot*Scale[PCAParCnt2]+Center[PCAParCnt2]);


				  //printf("ParamChoice:%d dummyDot:%f Scale:%f Center:%f temp:%e\n",ParamChoice,dummyDot,Scale[PCAParCnt2],Center[PCAParCnt2],temp);
					
				   

		 
				  
				  
				   if(temp<0) temp = 1e-10;


				    ParamVals.FxdPars[ParamChoice] = temp;
				    FxdPars[ParamChoice] = temp;

				    
				    if(FxdPars[ParamChoice]>MaxParam[ParamChoice]){
			    			
						//printf("ParamChoice:%d FxdPars[ParamChoice]:%f MaxParam[ParamChoice]:%f \n",ParamChoice,FxdPars[ParamChoice],MaxParam[ParamChoice]);  getc(stdin);
			    			Stop = 1;
			    			NewPost = 0;
			  	    } 
				    

				   // printf("Itns ParamChoice:%d temp:%f\n",ParamChoice,temp); 
				    PCAParCnt2++;
				  	  					  
				    
		  } //NoGo 
		} //PC
		//getc(stdin);



		if(Stop<1){
			CrashOut = 0;
			NewOldLHoodPerPop = 0.0; 
			LHoodPerPop = 0.0;
       		        AvgLHood = 0.0;

			
			for(Pop=1;Pop<=NumPops[DataSet];Pop++){
				

					Cancel = 0;
					ParamVals.FxdPars[0] = Pop; //the Plot Number
			
					
       					dim = GetDim(ParamVals.FxdPars[11],Pop)/ParamVals.FxdPars[8];
					
					//Monte function G()
      					//First time: do it for proposed parameter values
					
					


					gsl_monte_function G = { &LHoodRK45, dim, &ParamVals };
					gsl_monte_miser_state *s = gsl_monte_miser_alloc (dim);


         				gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &res, &err);
				      
					
					//Now set ParamVals.FxdPars back to the OLD set of parameter values...
					for(PC=1;PC<=ParamCount[DataSet];PC++){
						
						ParamChoice = ParamNum[PC];
						ParamVals.FxdPars[ParamChoice] = CurrentPars[ParamChoice]; //Go back and re-calculate the likelihood again for the current parameters
						//printf("ParamChoice:%d ParamVals.FxdPars:%f\n",ParamChoice,ParamVals.FxdPars[ParamChoice]);

					}
					//getc(stdin);

					
					
				       //and re-calculate the OLD likelihood, to avoid winner's curse
         			       gsl_monte_miser_integrate (&G, xl, xu, dim, calls, r2, s, &oldres, &err); 
				       gsl_monte_miser_free(s);
					

					
					
					//Now set ParamVals.FxdPars back to the proposed set of parameters				
				      	for(PC=1;PC<=ParamCount[DataSet];PC++){
						
						ParamChoice = ParamNum[PC];
						ParamVals.FxdPars[ParamChoice] = FxdPars[ParamChoice]; //CurrentPars holds the old parameters, so go back to using ParamVals.FxdPars to  hold the proposed parameters

					}
					
					
							

				
				
				NewOldLHoodPerPop += log(oldres);

				
				if(res<=0){
					 Cancel = 1;
					// printf("Cancel:%d res:%f\n",Cancel,res);
				}
				if((isnan(res)!=0)||(isinf(res)!=0)){
					 Cancel = 1;
					// printf("Cancel:%d res:%f\n",Cancel,res);

				}



				

                            	if(Cancel!=1){
					LHoodPerPop += log(res);
				}else{
					LHoodPerPop = -1000;
					//printf("LHood Crash...\n");
				}

				//printf("Itn:%f Pop:%d Cancel:%d res:%e LHood:%e OldPost:%f\n",Itn,Pop,Cancel,res,LHoodPerPop,OldPost);
				//getc(stdin);

				 
       				//gsl_rng_free (r);
			} //Pop loop



		   NewPost = LHoodPerPop;
	      //  printf("LHoodPerPop:%f\n",LHoodPerPop); getc(stdin);
		
		//printf("just before if(Itn>1) Itn:%f NewPost:%f OldPost:%f\n",Itn,NewPost,OldPost); getc(stdin);
		
		//if(Itn>1){ //THIS IS WRONG, BECAUSE IT ASSUMES THAT OldPost HAS NO VALUE WHEN Itn=1, BUT IN FACT IT DOES
			
				 OldPost += NewOldLHoodPerPop - OldLHoodPerPop;
				 OldLHoodPerPop = NewOldLHoodPerPop;

	
				// if(Itn>2101){("Itn:%f OldPost:%f NewOldLHoodPerPop:%f OldLHoodPerPop:%f OldPost:%f\n",Itn,NewOldLHoodPerPop,OldLHoodPerPop,OldPost); getc(stdin);}
				
			 
		//}else{ WRONG
		//	OldLHoodPerPop = NewOldLHoodPerPop; WRONG

			 
		//} WRONG

		//printf("just after if(Itn>1) OldPost:%f\n",OldPost);


		double dummyNew, dummyOld;

		// Almost sure this is obsolete 28 Feb 2015
		/*
		if(InitFlag[ParamChoice]>1){
			Inflate = InitScale;
		}else{
			Inflate = 1;
		}
		*/
		
	
		/////////////// Adjusting for proposal ////////////////////

		NewPropAdj = log(gsl_ran_gaussian_pdf(PCA[PCAParCnt],SD[PCAParCnt])); 
		if(Itn>1)  //WRONG?? SHOULD PROBABLY BE IF(ACCEPTCOUNT>0)
			VarScale = 1.0; 
		OldPropAdj = log(gsl_ran_gaussian_pdf(OldPCA[PCAParCnt],VarScale*SD[PCAParCnt]));

		//if(Itn>2101) printf("PCAParCnt:%d PCA:%f OldPCA:%f SD:%f NewPropAdj:%f OldPropAdj:%f\n",PCAParCnt,PCA[PCAParCnt],OldPCA[PCAParCnt],SD[PCAParCnt],NewPropAdj,OldPropAdj);

		/////////////// Done adjusting for Proposal 



 
		//////////// Adjusting for prior /////////////
		for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		  ParamChoice2 = ParamNum[PC2];
		      
		  NoGo = 0;
		
		  
		
	          if((PC2==1)&&(kHtg>1000)) NoGo = 1;
		  if((PC2==5)&&(ratioFix==1)) NoGo = 1;
		  if((PC2==6)&&(deltaFix==1)) NoGo = 1;
		  if((PC2==7)&&(mFix==1)) NoGo = 1;
		  

	          

		  //if((DataSet==15)&&(PC==1)) NoGo = 1;
		  if(NoGo!=1){
		 		
	               float tmp2;
		        switch(PriorDist[ParamChoice2]){
			case 1 :
			  if(FxdPars[ParamChoice2]>MaxParam[ParamChoice2]){
			    
			    Stop = 1;
			    NewPost = 0;
			    
			  } else {
			    tmp2 =  (gsl_ran_flat_pdf(FxdPars[ParamChoice2],0,PriorVar[ParamChoice2]));
				


			  }
			  break;
			case 2 :tmp2 = gsl_ran_gaussian_pdf(log(FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // a slightly bogus kind of normal, really sort log normal

			
				break;

			case 3 :  tmp2 = (gsl_ran_lognormal_pdf(FxdPars[ParamChoice2],PriorMean[ParamChoice2],PriorVar[ParamChoice2]));

					 
					 break;
	

			 case 4 :tmp2 = gsl_ran_gaussian_pdf((FxdPars[ParamChoice2])-PriorMean[ParamChoice2],PriorVar[ParamChoice2]); // more like a real normal...

				
				 break;
	                 default : printf("No prior distribution specified.  Bailing.\n"); getc(stdin); exit(1); break;
		      	} //End of switch
			if((isnan(tmp2)!=0)||(isinf(tmp2)!=0)){
				Cancel = 1;
				//printf("Cancel:%d tmp2:%f\n",Cancel,tmp2);
				tmp2 = 1e-300;
			} else{

				//if(Itn>168) printf("NewPost before adding tmp2 NewPost:%f log(tmp2):%f \n",NewPost,log(tmp2));
			    	NewPost += log(tmp2);
				//printf("log(tmp2):%f NewPost:%f\n",log(tmp2),NewPost);
			}
			if((isnan(NewPost)!=0)||(isinf(NewPost)!=0)){
				Cancel = 1;
				//printf("ParamChoice2:%d PriorMean:%f log(FxdPars[ParamChoice2]):%f Cancel:%d NewPost:%f\n",ParamChoice2,PriorMean[ParamChoice2],log(FxdPars[ParamChoice2]),Cancel,NewPost);
				NewPost = -300;
				
				
			}
			



		      } //ParamChoice2!=4, we're not doing sigma
	
		    } //PC2
	
	
	/*
	if(fmod(Itn,1.0)==0){

		printf("New Itn:%f k:%f nubar:%f mu:%f sig:%f ratio:%f delta:%f m:%f gamma:%f LH:%f OldLH:%f NewPost:%f OldPost:%f\n",Itn,FxdPars[1],FxdPars[2],FxdPars[3],FxdPars[4],FxdPars[5],FxdPars[6],FxdPars[9],FxdPars[20],LHoodPerPop,NewOldLHoodPerPop,NewPost-NewPropAdj,OldPost-OldPropAdj);
		printf("Old Itn:%f k:%f nubar:%f mu:%f sig:%f ratio:%f delta:%f m:%f gamma:%f LH:%f OldLH:%f NewPost:%f OldPost:%f\n",Itn,CurrentPars[1],CurrentPars[2],CurrentPars[3],CurrentPars[4],CurrentPars[5],CurrentPars[6],CurrentPars[9],CurrentPars[20],LHoodPerPop,NewOldLHoodPerPop,NewPost-NewPropAdj,OldPost-OldPropAdj);
		//getc(stdin);

		//if(Itn>2101) getc(stdin);

		//if(NewPost>LHoodPerPop) getc(stdin);
		//if(OldPost>NewOldLHoodPerPop) getc(stdin);


		printf("Itn:%f Cancel:%d Stop:%d ParamChoice:%d FxdPars:%f LH:%f NewPost:%e OldPost:%e\n",Itn,Cancel,Stop,ParamChoice,FxdPars[ParamChoice],LHoodPerPop,NewPost,OldPost);
	       //if(FxdPars[2]>25) getc(stdin);


	}
	*/
	
	
		if(Itn>0){ 

		  //printf("Cancel:%d\n",Cancel);
		  if(Cancel!=1){//Cancel if prior gives bad calc or lhood gives bad calc
		  

			if(OldPost>-900){
				Criterion = exp(NewPost - NewPropAdj -(OldPost-OldPropAdj));
			}else{
				Criterion = 100.0;
			}
			//printf("NewPost-NewPropAdj:%f OldPost-OldPropAdj:%f Criterion:%f\n",NewPost-NewPropAdj,OldPost-OldPropAdj,Criterion); getc(stdin);
			


			if((isnan(Criterion)!=0)||(isinf(Criterion)!=0)){
				Criterion = -1;
			}
			
			temp = gsl_rng_uniform(r2);


			
			
			/*
			if(fmod(Itn,1)==0){
				printf("%lf ",Itn);
				for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
				 ParamChoice2 = ParamNum[PC2];
				 printf("%f ",FxdPars[ParamChoice2]);
				}
				if(sigma>0){
					printf("Itn:%f LH:%f NewOldLH:%f NewPost:%e OldPost:%e Criterion:%e\n",Itn,LHoodPerPop,NewOldLHoodPerPop,NewPost,OldPost,Criterion);
				}else{
					printf("Itn:%f LH:%f OldLH:%f NewPost:%e OldPost:%e Criterion:%e\n",Itn,LHoodPerPop,OldLHoodPerPop,NewPost,OldPost,Criterion);
				}

				//getc(stdin);

			}
			*/
			
				
			
			

			if((Criterion>1)||(temp<Criterion)){

			 //printf("Accept...\n"); getc(stdin);
		       

			  //Acceptance means changing the current parameters to the proposed parameters
			  for(PC2=1;PC2<=ParamCount[DataSet];PC2++){  //Isn't this whole loop pointless?  No, this is where FxdPars gets reset to CurrentPars.  Soon, ONE element of FxdPars will change, BUT ONLY ONE 
		  		 
		                ParamChoice = ParamNum[PC2];
				CurrentPars[ParamChoice] = FxdPars[ParamChoice];
			     
			                
		  	  }



			  OldPCA[PCAParCnt] = PCA[PCAParCnt];
				
			  
			  


			 
			
			 // printf("Inside Itn:%f LH:%f OldLH:%f NewPost:%e OldPost:%e Criterion:%e\n",Itn,LHoodPerPop,OldLHoodPerPop,NewPost,OldPost,Criterion);
			  //if(Itn>2101) getc(stdin);

			  OldPost = NewPost;
			  OldLHoodPerPop = LHoodPerPop;
			  
			 // if(InitFlag[ParamChoice]==1) InitFlag[ParamChoice] = 0; //You have moved on from the initial value of this parameter - OBSOLETE as of 28 Feb 2015
			   
			  
			  //AcceptCount++;
			  //Accept[ParamChoice]++;

			   
			} else { //Criterion < 1 and temp > Criterion
				OldLHoodPerPop = NewOldLHoodPerPop;
				
			}
		 } //if Cancel

		      
		 //Reset FxdPars to the current parameters, which may or may not have been updated	         	
		  for(PC2=1;PC2<=ParamCount[DataSet];PC2++){  //Isn't this whole loop pointless?  No, this is where FxdPars gets reset to CurrentPars.  Soon, ONE element of FxdPars will change, BUT ONLY ONE 
		  		 
		                ParamChoice2 = ParamNum[PC2];
			     
			        FxdPars[ParamChoice2] = CurrentPars[ParamChoice2];
				ParamVals.FxdPars[ParamChoice2] = CurrentPars[ParamChoice2];
				
				
			         
		  }


   
		  if(Thin==ThinStop){

			int prnt = 0;
			 strcpy(testbuff,test[0]);
		  	 fp2 = fopen(testbuff,"a");
		        
		        fprintf(fp2,"%f  ",Itn);
			 if(prnt==1) printf("%f  ",Itn);

			 for(PC2=1;PC2<=ParamCount[DataSet];PC2++){
		  		 
		                 ParamChoice2 = ParamNum[PC2];
			     
			         
				 NoGo = 0;
                             if((ParamChoice2==1)&&(kHtg>1000)) NoGo = 1;
			     if((ParamChoice2==3)&&(mu<=0.0)) NoGo = 1;
		             if((ParamChoice2==4)&&(sigma==0.0)) NoGo = 1;
		             if((ParamChoice2==5)&&(ratioFix==1)) NoGo = 1;
			     if((ParamChoice2==6)&&(deltaFix==1)) NoGo = 1;
			     if((ParamChoice2==9)&&(mFix==1)) NoGo = 1;

				 

				 
				 fprintf(fp2,"%e  ",FxdPars[ParamChoice2]);
				 if(prnt==1) printf("%e  ",FxdPars[ParamChoice2]);

			 }
			 
			 fprintf(fp2,"%e %e\n",OldLHoodPerPop,OldPost);
			 if(prnt==1) printf("%e %e\n",OldLHoodPerPop,OldPost);
			

			 fclose(fp2); 
			 
			 Thin = 0;
		 }//Thin == ThinStop
			  
			
	} //else{ //Itn == 1
			
			//OldPost = NewPost;
						
       //}

		 
		 //getc(stdin);

   
		          	  

		 Itn++;
		 Thin++;
		 ParamItn[ParamChoice]++;
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
		

		 // } // PC!=4, meaning not sigma

	
	         } //Stop 

			Stop = 0;
			PCAParCnt++;
		} //PCA loop?
 	}  //Huh...




	







} //Itn










	



//gsl_rng_free (r);


 return 0;

}


int Debug = 0;
extern double PointwiseLH[30][100];
extern double PointwiseVar[30][100];
extern double PointwiseMean[30][100];


double LogBetaBnmlLhood(float Inf, float Total, float p,float gamma, float delta){

	double n, x;
        double a, b;
        float LognChoosek;
	double LHood;
	float FractInf = p;
	FILE *fp;



	/*
	if(p==1.0)
        	FractInf = 1 - (1e-10);
	
	if(p==0)
		FractInf = 1e-10;
        */


	//LognChoosek = gsl_sf_lnchoose(n,k);

	//LHood = LognChoosek + k*log(FractInf) + (n-k)*log(1-FractInf);

     			
 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"Just after variable declarations, in LogLhoodBetaBnml. Inf:%f Total:%f p:%f gamma:%f delta:%f\n",Inf,Total,p,gamma,delta);
     fclose(fp);
   }

       n = Total;
       x = Inf;
	
	
	if(p>1) p = 1;
	if(p<0) p = 0;
	if(isnan(p)==1){
		//printf("nan!\n");
 		p = 0;
	}
	if(x>n) x = n;
       //if(x>1e100) x = 1e10;
	//if(n>1e100) n = 1e10;
	if(isnan(x)==1) x = 0;
       if(isnan(n)==1) n = 1;



        a = p*exp(gamma) + delta; 
        b = (1-p)*exp(gamma) + delta;


 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"Just after a and b calc, in LogLhoodBetaBnml. a:%f b:%f \n",a,b);
     fclose(fp);
   }





	LHood = gsl_sf_lnbeta(x+a,n-x+b) - gsl_sf_lnbeta(a,b) - gsl_sf_lnbeta(x+1,n-x+1) - log(n+1);





 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"Just after x:%f n:%f Inf:%f Total:%f\n",x,n,Inf,Total);
     fclose(fp);
   }


	return LHood;

}

double LhoodFuncBetaBnml(int MaxRlzns, int MaxPlots, int MaxWk,  float ***Data, float ***Model, double *lhoodVec,float gamma, float delta){


	FILE *fp,*fp2, *fp3;
	int i, Plot, Rlzn;
	double LHood = 0;
	double dummy;

//int Plot, i;


	//printf("Inside LhoodFunc...gamma:%f\n",gamma);
	//getc(stdin);
	
/*
fp = fopen("DataOutMoreau.dat","w");
for(Plot=1;Plot<=3;Plot++)
	for(i=1;i<=10;i++)
		if(Data[1][Plot][i]>=0)
			fprintf(fp,"%d\t%f\n",i,((float) Data[1][Plot][i]/Data[2][Plot][i]));
fclose(fp);
exit(1);
*/

	if(Verbose==1){
		fp2 = fopen("PlotOSnumuFxdXC.dat","a");
		fp = fopen("LHOSnumuFxdXC.dat","a");
		//fp = fopen("LHoodJunk.dat","a");
		//fp2 = fopen("PlotModelJunk.dat","a");
	}
	/*
	if(DIC==1){
		fp3 = fopen("G8O4HRlzns.dat","a");
	}
	*/


     
	
 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"After variable declarations in LhoodFuncBetaBnml. \n");
     fclose(fp);
   }

//printf("MaxPlots:%d\n",MaxPlots); getc(stdin);
	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
		for(Plot=MaxPlots;Plot<=MaxPlots;Plot++){
			for(i=1;i<=MaxWk;i++){
			 float dummy = ((float) Data[1][Plot][i])/((float) Data[2][Plot][i]);
			/*if(Plot==8){
				printf("here i:%d Plot:%d S:%f I:%f Model:%f dataFI:%f gamma:%f delta:%f\n",i,Plot,Data[1][Plot][i],Data[2][Plot][i],Model[Rlzn][Plot][i],dummy,gamma,delta);
				getc(stdin);
			}*/

					if((Verbose==1)&&(Data[1][Plot][i]>=-1)){
					
						//fprintf(fp2,"%d\t%f\n",i,Data[2][Plot][i],Model[Rlzn][Plot][i]);
						float dummy = Data[1][Plot][i]/Data[2][Plot][i];
						if(dummy==1) dummy = -1;
						fprintf(fp2,"%d\t%f\t%f\n",i,dummy,Model[Rlzn][Plot][i]);
                                                //fprintf(fp2,"%d\t%f\t%f\n",i,Data[1][Plot][i],Model[Rlzn][Plot][i]);
        
                                                
					
						//fprintf(fp2,"%d\t%f\t%f\t%f\t%f\t%f\t%f\n",i,Data[1][Plot][i],Data[2][Plot][i],Model[Rlzn][Plot][i],dummy,gamma,delta);
					}
		
					if(DIC==1){ //What is the point of this?

						float dummy = Data[1][Plot][i]/Data[2][Plot][i];
						if(dummy==1) dummy = -1;
						
					}

						


					if(Data[1][Plot][i]>-1){
						//dummy = LogBnmlLhood(Data[1][Plot][i],Data[2][Plot][i],Model[Rlzn][Plot][i]);
	
 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"Just before LogBetaBnmlLhood call in LhoodFuncBetaBnml. \n");
     fclose(fp);
   }

						dummy = LogBetaBnmlLhood(Data[1][Plot][i],Data[2][Plot][i],Model[Rlzn][Plot][i],gamma,delta);
						if(WAIC==1){  //Apparently this is dead code...WAIC is no longer calculated this way
							 if(PointwiseLH[Plot][i]<=-1e10){	
							 	PointwiseLH[Plot][i] = dummy;									
							 } else{
								PointwiseLH[Plot][i] += dummy;
							 }
								
			
							 if(PointwiseVar[Plot][i]<=-1e10){	
							 	PointwiseVar[Plot][i] = dummy*dummy;
								PointwiseMean[Plot][i] = dummy;										
							 } else{
								PointwiseVar[Plot][i] += dummy*dummy;
								PointwiseMean[Plot][i] += dummy;
							 }

							 
							 //printf("Plot:%d i:%d dummy:%f PointwiseVar:%f PointwiseMean:%f \n",Plot,i,dummy,PointwiseVar[Plot][i],PointwiseMean[Plot][i]); 
							 //getc(stdin);
							 
						}
						
 if(Debug==1){
     fp = fopen("DebugOut7.dat","at");
     fprintf(fp,"Just after LogBetaBnmlLhood call in LhoodFuncBetaBnml. \n");
     fclose(fp);
   }

					


						LHood += dummy;


					}
			}
		}
		if(Verbose==1)
			fprintf(fp,"%f\n",lhoodVec[Rlzn]);


	}

	if(Verbose==1){
		fclose(fp);
		fclose(fp2);
	}
	
	//if(DIC==1) fclose(fp3);

	if(LHood<-700) LHood = -700;

	lhoodVec[1] = LHood; //this is actually what the calling function uses...
	//getc(stdin);

	return LHood; //this by comparison is pretty much irrelevant...

}
double LogBnmlLhood(float Inf, float Total, float p){

	unsigned int n, k;
        float LognChoosek;
	double LHood;
	float FractInf = p;
	FILE *fp;

	//fp = fopen("GMBnmlOut.dat","at");

	if(p==1.0)
        	FractInf = 1 - (1e-10);
	
	if(p==0)
		FractInf = 1e-10;

	n = round(Total);
	k = round(Inf);

	LognChoosek = gsl_sf_lnchoose(n,k);

	LHood = LognChoosek + k*log(FractInf) + (n-k)*log(1-FractInf);
	
	//fprintf(fp,"%f\t%f\t%f\t%f\n",LognChoosek,FractInf,k*log(FractInf),(n-k)*log(1-FractInf));
	//fprintf(fp,"%f\t%f\t%f\n",LognChoosek,FractInf,Inf/Total);
	//fclose(fp);
	return LHood;

}

double LhoodFuncBnml(int MaxRlzns, int MaxPlots, int MaxWk,  float ***Data, float ***Model, double *lhoodVec){


	FILE *fp,*fp2;
	int i, Plot, Rlzn;
	double LHood = 0;
	double dummy;

	//fp = fopen("PlotModelOutBnml.dat","at");
	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
		for(Plot=MaxPlots;Plot<=MaxPlots;Plot++){
			for(i=1;i<=MaxWk;i++){
					if(Data[1][Plot][i]>-1){
	//Model[Rlzn][Plot][i] = 0.5;
						dummy = LogBnmlLhood(Data[1][Plot][i],Data[2][Plot][i],Model[Rlzn][Plot][i]);
						//dummy = LogBnmlLhood(Data[1][Plot][i],Data[2][Plot][i],0.5);
						LHood += dummy;
//fprintf(fp,"%d\t%f\n",i,Model[Rlzn][Plot][i]);
//printf("%d\t%f\n",i,Model[Rlzn][Plot][i]);
//getc(stdin);
					}
			}
		}
	}
//fclose(fp);

//fp2 = fopen("LHoodOutBnml.dat","at");
//fprintf(fp2,"%f\n",LHood);
//fclose(fp2);
	 lhoodVec[1] = LHood; //this is actually what the calling function uses...

	return LHood; //this by comparison is pretty much irrelevant...

}

double SSECalc(int MaxRlzns, int Plot, int MaxWk,  float ***Data, float ***Model, float ***SSE,int WhichData){


	FILE *fp,*fp2,*fp3;
	int i, Rlzn;


/*
fp = fopen("DataOut.dat","w");
for(Plot=1;Plot<=5;Plot++)
	for(i=1;i<=10;i++)
		if(Data[WhichData][Plot][i]>=0)
			fprintf(fp,"%d\t%f\n",i,Data[WhichData][Plot][i]);
fclose(fp);
exit(1);

*/

	if(Verbose==1)
		fp = fopen("PlotModelOutGauss.dat","at");
	for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
		//for(Plot=MaxPlots;Plot<=MaxPlots;Plot++){
			for(i=1;i<=MaxWk;i++){
				SSE[Rlzn][Plot][i] = -1;
					if(Data[WhichData][Plot][i]>-1){
						SSE[Rlzn][Plot][i] = (Model[Rlzn][Plot][i]-Data[WhichData][Plot][i])*(Model[Rlzn][Plot][i]-Data[WhichData][Plot][i]); 
							if(Verbose==1)
								fprintf(fp,"%d\t%f\n",i,Model[Rlzn][Plot][i]);
//						//printf("%d\t%d\t%d\t%f\t%f\n",Rlzn,Plot,i,Model[Rlzn][Plot][i],Data[WhichData][Plot][i]);
					}
			}
		//}  //Plot
	}

	if(Verbose==1)
		fclose(fp);
	//fclose(fp2);
	//exit(1);

	return 0;

}

double LHoodFunc2(int MaxRlzns, int Plot, int MaxWk, float ***SSE, double tau, double *lhoodVec)
{

 int i;
 double pi = 3.141592654;

 double arg2 = (tau*sqrt(2*pi));
 double LHood = 1; int Num;
 double AvgLHood = 0; 
 double logLHood = 0;
 double TOOLOW = 1e-299;
 double VarLHood = 0;
 //int Plot, 
 int Rlzn;
 double arg1, arg3; 
 double junk;
 int PlotJunk;
 FILE *fp,*fp2;

//fp = fopen("avglhoodout.dat","wt");
//fclose(fp);

	if(Verbose==1)
		fp = fopen("LHoodOutGauss.dat","a");


 for(Rlzn=1;Rlzn<=MaxRlzns;Rlzn++){
	 Num = 0; 

	 LHood = 1; logLHood = 0.0;
	 //for(Plot=MaxPlots;Plot<=MaxPlots;Plot++){
		 for(i=1;i<=MaxWk;i++){
			   if(SSE[Rlzn][Plot][i]>-1){

				   //arg1 = (Model[Rlzn][Plot][i]-Data[Plot][i])*(Model[Rlzn][Plot][i]-Data[Plot][i]); 

				   arg1 = SSE[Rlzn][Plot][i];
				   //SSE += arg1;
				   arg3 = (2*tau*tau);
				   LHood *= exp(-arg1/arg3)/arg2;

				   logLHood += (-arg1/arg3) - log(arg2);
				   Num++;
			   } //Data not missing
			
		 } //i loop
		 
	 //} //Plt
	 
	 	lhoodVec[Rlzn] = logLHood;
		if(Verbose==1)
			fprintf(fp,"%f\n",lhoodVec[Rlzn]);
		AvgLHood += LHood;
		VarLHood += LHood*LHood;


 } //Rlzns
 
 if(Verbose==1)
	fclose(fp);
 //exit(1);
 


 AvgLHood /= MaxRlzns;
 VarLHood /= MaxRlzns;
 VarLHood -= AvgLHood*AvgLHood;


  if(AvgLHood<=0) AvgLHood = 1e-299;


  return -log(AvgLHood); 
 

} //LhoodFunc2




#define MIN 0.000
#define MIN2 0.000

int swplots5(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval, int *SData,int 
*FractDead)
{

		 int Extra = 0;
		 int Plot, MaxPlot = 5;
		 int Wk;

		 *SData = 0;
		 *FractDead = 1;
		 *MaxWk = 7;
		 *DataInterval = 7;
                 int DayTemp = 7-2;


		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	     //FractI[1] = 0.25; //0.0885;
	     FractI[1] = 0.0885;
	      Data[1][1][1] = 0.05;
	     Data[1][1][2] = 0.20; Data[1][1][3] = 0.50; Data[1][1][4] = 0.65;
	     Split[1] = 1; MaxDays[1] = DayTemp*7 + Extra;

		 S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	     FractI[2] = 0.0407;
	     Data[1][2][1] = MIN2;
	     Data[1][2][2] = 0.01; Data[1][2][3] = 0.04; Data[1][2][4] = 0.13; Data[1][2][5] = 0.14;
	     Split[2] = 0; MaxDays[2] = (DayTemp+1)*7 + Extra;

		 S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.0233;
	     Data[1][3][1] = MIN2;
	     Data[1][3][2] = 0.03; Data[1][3][3] = 0.06; Data[1][3][4] = 0.17;
	     Split[3] = 0; MaxDays[3] = (DayTemp)*7 + Extra;

		 S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     FractI[4] = 0.2854;
	     Data[1][4][1] = 0.05;
	     Data[1][4][2] = 0.25; Data[1][4][3] = 0.50; Data[1][4][4] = 0.55;
	     Split[4] = 1; MaxDays[4] = (DayTemp)*7 + Extra;

		 S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0943;
	     Data[1][5][1] = 0.02; Data[1][5][2] = 0.03;
	     Data[1][5][3] = 0.02; Data[1][5][4] = 0.07;Data[1][5][5] = 0.15; Data[1][5][6] = 0.17;
	     Split[5] = 0; MaxDays[5] = (DayTemp+2)*7 + Extra;

		 for(Plot=1;Plot<=MaxPlot;Plot++)
			for(Wk=1;Wk<=(*MaxWk);Wk++)
				Data[2][Plot][Wk] = -2;


		 return MaxPlot;
}


//This next version actually has 8 plots 
/* DATA FOR THE NEXT FUNCTION HAVE BEEN CORRECTED twice, 
second time is 17 January 2015 when numbers were checked with DataThief on original figure, and plot 5 1984, Plot 5 1985, and Plot 16 1985 were added back in */ 
/*But as of 9 August 2014, we are going back to EMs hatched in lab for plot 1...*/ /* As of 18 August 2014, we are going back to the values observed in field for plot 1*/
int swplots8Bnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval, int *SData,int 
*FractDead)
{


		 int Extra = 0;
		 int Plot, MaxPlot = 8;
		 int Wk;

		 *SData = 0;
		 *FractDead = 1;


		 *MaxWk = 7;
		 *DataInterval = 7;
                 int DayTemp = 7-2; //don't quite know why...

		//Plot 1 '83
		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	        //FractI[1] = 0.0885; //This is from EMs hatched in lab
		 FractI[1] = 0.169; //This is the value in the first week
		//FractI[1][2] = 0.286; //Don't need this here, it's in main()
		
		Data[1][1][1] = round(0.043*141);
	       Data[1][1][2] = round(0.202*104); 
		Data[1][1][3] = round(0.496*124); 
		Data[1][1][4] = round(0.691*92);

	       Data[2][1][1] = 141;
	       Data[2][1][2] = 104; 
		Data[2][1][3] = 124; 
		Data[2][1][4] = 92;
	       Split[1] = 1; MaxDays[1] = 5*7 + Extra; //It's 5 instead of 4 because Split[1] = 1


		//Plot 3 '83
		S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	       FractI[2] = 0.0167; //0.0407; NOTE 0.0407 MUST BE FROM EMS HATCHED IN LAB, 0.017 IS OBSERVED VALUE
	       Data[1][2][1] = 0;
	       Data[1][2][2] = round(0.0*158); 
		Data[1][2][3] = round(0.00606*145); 
		Data[1][2][4] = round(0.0333*151); 
		Data[1][2][5] = round(0.1159*122);
		Data[1][2][6] = round(0.1364*120); //Sample size is just the reported minimum

		Data[2][2][1] = 163;
	       Data[2][2][2] = 158; 
		Data[2][2][3] = 145;
		Data[2][2][4] = 151; 
		Data[2][2][5] = 122;
		Data[2][2][6] = 120;
		
	       Split[2] = 0; MaxDays[2] = 6*7 + Extra;


	     
	     //Plot 5 '83
            S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.014; //0.026; //0.0233; NOTE 0.0233 MUST BE FROM EMs HATCHING IN LAB, 0.026 IS OBSERVED 
	     Data[1][3][1] = 0;
	     Data[1][3][2] = 0;
	     Data[1][3][3] = round(0.0359*132); 
	     Data[1][3][4] = round(0.074*128);  
            Data[1][3][5] = round(0.177*111);

	     Data[2][3][1] = 150;
	     Data[2][3][2] = 120;  //just the minimum reported sample size
	     Data[2][3][3] = 132; 
	     Data[2][3][4] = 128; 
            Data[2][3][5] = 111;
	     Split[3] = 0; 
	     MaxDays[3] = 5*7 + Extra;

	     //Plot 10 '84
	     S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     //FractI[4] = 0.2854;
	     FractI[4] = 0.1149;
	    // FractI[4][2] = 0.134;
	     Data[1][4][1] = round(0.0427*466);
	     Data[1][4][2] = round(0.2741*451); 
	     Data[1][4][3] = round(0.5078*391); 
            Data[1][4][4] = round(0.5447*181);

	     Data[2][4][1] = 466;
	     Data[2][4][2] = 451; 
	     Data[2][4][3] = 391; 
	     Data[2][4][4] = 181;
	     Split[4] = 1; MaxDays[4] = 5*7 + Extra;  //it's 5*7 instead of 4*7 because Split = 1

	     //Plot 1 '85
	     S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0561; //0.0943; NOTE: 0.0943 must be the initial fraction as estimated from hatching larvae from egg masses

	     Data[1][5][1] = round(0.00564*1083); 
	     Data[1][5][2] = round(0.0182*660);
	     Data[1][5][3] = round(0.0315*607); 
	     Data[1][5][4] = round(0.0184*558);
	     Data[1][5][5] = round(0.0746*498); 
	     Data[1][5][6] = round(0.147*48);
	     Data[1][5][7] = round(0.182*48); //sample size is a guess, based on previous week
		 
            Data[2][5][1] = 1083; 
	     Data[2][5][2] = 660;
	     Data[2][5][3] = 607; 
            Data[2][5][4] = 558;
            Data[2][5][5] = 498; 
	     Data[2][5][6] = 48;
	     Data[2][5][7] = 48; //sample size is a guess

	     Split[5] = 0; MaxDays[5] = (DayTemp+2)*7 + Extra;

	     //Plot 5 '85
	     S0[6] = 2355.0 * 192.0/(1e4 * 1.4);
	     FractI[6] = 0.0681; //From EMs in lab?: 0.1034;
	     Data[1][6][1] = round(0.0482*861); 
	     Data[1][6][2] = round(0.0272*851);
	     Data[1][6][3] = round(0.029*845); 
	     Data[1][6][4] = round(0.1073*840); //Sample size is a guess, rounding down from the minimum of the other 3, 17 January 2015
		 
	     Data[2][6][1] = 861; 
	     Data[2][6][2] = 851;
	     Data[2][6][3] = 845; 
	     Data[2][6][4] = 840;
	     Split[6] = 0; //1; Why was this ever 1? 
	     MaxDays[6] = 4*7 + Extra;

	
	     //Plot 5 '84     	     	  
	     S0[7] = 427.0 * 204/(1e4 * 1.4);  //EM size here is the average of EM sizes for plot 5 '83, and plot 5 '84, because I don't know where the original got to
	     FractI[7] = 0.00934;
	     Data[1][7][1] = round(0.0154*450); 
	     Data[1][7][2] = round(0.0206*450);
	     Data[1][7][3] = round(0.0343*450); 
	     Data[1][7][4] = round(0.0214*450);
	     Data[1][7][5] = round(0.016*450); 
	     Data[1][7][6] = round(0.0586*450); 
	     Data[1][7][7] = round(0.0388*450);  
		 
	     //Sample sizes are minimum sample sizes given in the paper, because I don't know where the actual sample sizes are
	     Data[2][7][1] = 450; 
	     Data[2][7][2] = 450;
	     Data[2][7][3] = 450; 
	     Data[2][7][4] = 450;
	     Data[2][7][5] = 450;
	     Data[2][7][6] = 450;
	     Data[2][7][7] = 450;
	     Split[7] = 0; 
	     MaxDays[7] = 7*7 + Extra;

	     
	     
	   //  *FractI = 0.2869;
	   //  S0 = 827.0 * 146.0/(1e4 * 1.4);
	   //  *maxt = 49;
	   //  maxData = 6;
	   //  Data[0] = -1.0; Data[1] = -1.0; Data[2] = 0.08; Data[3] = 0.02;
	   //  Data[4] = 0.07; Data[5] = 0.15;
	   //  for(i=0;i<maxData;i++){
	   //   S[i] = 450;
	   //   I[i] = round(Data[i]*S[i]);
	   //  }
	   //  	     S[0] = 1; I[0] = 0;
	   //  S[1] = 609; I[1] = 62; S[2] = 844; I[2] = 19; S[3] = 1; I[3] = -1;
	   //  S[4] = 842; I[4] = 121;
	   //


	     
	       S0[8] = 827.0 * 146.0/(1e4 * 1.4);
   	   	FractI[8] = 0.0739;  //
	   	Data[1][8][1] = round(0.0755*500);  //value is missing above?  But not in published paper...500 seems conservative, given the other values
	   	Data[1][8][2] = round(0.0674*842);
		Data[1][8][3] = round(0.1605*450); //Sample size is the minimum, from the paper

		 
	   	Data[2][8][1] = 500; 
	   	Data[2][8][2] = 842;
		Data[2][8][3] = 450;

	   	Split[8] = 0; MaxDays[8] = 3*7 + Extra;


	     
		 return MaxPlot;
}



/* DATA FOR THE NEXT FUNCTION HAVE BEEN CORRECTED once */ /*But as of 9 August 2014, we are going back to EMs hatched in lab for plot 1...*/ /* As of 18 August 2014, we are going back to the values observed in field for plot 1*/

int swplots5Bnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval, int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 int Plot, MaxPlot = 5;
		 int Wk;

		 *SData = 0;
		 *FractDead = 1;


		 *MaxWk = 6;
		 *DataInterval = 7;
                 int DayTemp = 7-2; //don't quite know why...

		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	    //FractI[1] = 0.0885; //This is from EMs hatched in lab
		FractI[1] = 0.179;
		//FractI[1][2] = 0.286; //Don't need this here, it's in main()
		
		 Data[1][1][1] = round(0.05*141);
	     Data[1][1][2] = round(0.20*104); Data[1][1][3] = round(0.5*124); Data[1][1][4] = round(0.65*92);
	      Data[2][1][1] = 141;
	     Data[2][1][2] = 104; Data[2][1][3] = 124; Data[2][1][4] = 92;
	     Split[1] = 1; MaxDays[1] = DayTemp*7 + Extra;

		 S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	     FractI[2] = 0.017; //0.0407; NOTE 0.0407 MUST BE FROM EMS HATCHED IN LAB, 0.017 IS OBSERVED VALUE
	     Data[1][2][1] = 0;
	     Data[1][2][2] = round(0.01*158); Data[1][2][3] = round(0.04*145); Data[1][2][4] = round(0.13*151); 
		 Data[1][2][5] = round(0.14*122);
		 Data[2][2][1] = 163;
	     Data[2][2][2] = 158; Data[2][2][3] = 145; Data[2][2][4] = 151; Data[2][2][5] = 122;
		

	     Split[2] = 0; MaxDays[2] = (DayTemp+1)*7 + Extra;

		 S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.026; //0.0233; NOTE 0.0233 MUST BE FROM EMs HATCHING IN LAB, 0.026 IS OBSERVED 
		 Data[1][3][1] = 0;
	     Data[1][3][2] = round(0.03*132); Data[1][3][3] = round(0.06*128); Data[1][3][4] = round(0.17*111);
	     Data[2][3][1] = 150;
	     Data[2][3][2] = 132; Data[2][3][3] = 128; Data[2][3][4] = 111;
	     Split[3] = 0; MaxDays[3] = (DayTemp)*7 + Extra;

	
		 S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     //FractI[4] = 0.2854;
	     FractI[4] = 0.119;
	    // FractI[4][2] = 0.134;
	     Data[1][4][1] = round(0.05*466);
	     Data[1][4][2] = round(0.25*451); Data[1][4][3] = round(0.50*391); Data[1][4][4] = round(0.55*181);

		 Data[2][4][1] = 466;
	     Data[2][4][2] = 451; Data[2][4][3] = 391; Data[2][4][4] = 181;
	     Split[4] = 1; MaxDays[4] = (DayTemp)*7 + Extra;

		 S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0569; //0.0943; NOTE: 0.0943 must be the initial fraction as estimated from hatching larvae from egg masses

	     Data[1][5][1] = round(0.02*1083); Data[1][5][2] = round(0.03*660);
	     Data[1][5][3] = round(0.02*607); Data[1][5][4] = round(0.07*558);Data[1][5][5] = round(0.15*498); 
		 Data[1][5][6] = round(0.17*48);
		 Data[2][5][1] = 1083; Data[2][5][2] = 660;
	     Data[2][5][3] = 607; Data[2][5][4] = 558;Data[2][5][5] = 498; 
		 Data[2][5][6] = 48;
	     Split[5] = 0; MaxDays[5] = (DayTemp+2)*7 + Extra;

	     
		 return MaxPlot;
}


/*this next one is really swplots6Bnml, because it has 6 plots */ /* DATA ARE NOT CORRECTED */
/*
int swplots5Bnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval, int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 int Plot, MaxPlot = 6; //5;
		 int Wk;

		 *SData = 0;
		 *FractDead = 1;


		 *MaxWk = 6;
		 *DataInterval = 7;
                 int DayTemp = 7-2; //why 7-2??
		
	     //S[0] = 186; I[0] = 8; 
           // S[1] = 546; I[1] = 94; 
          //  S[2] = 1; I[2] = -1;
	  //   S[3] = 144; I[3] = 22; 
          //  S[4] = 116; I[4] = 24;  
          //  S[5] = 139; I[5] = 70;
	  //   S[6] = 110; I[6] = 57;
		

		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	     FractI[1] = 0.0885;
		 //FractI[1] = 0.25;
	     Data[1][1][1] = round(0.05*141);
	     Data[1][1][2] = round(0.20*104); 
	     Data[1][1][3] = round(0.5*124); 
	     Data[1][1][4] = round(0.65*92);
	     Data[2][1][1] = 141;
	     Data[2][1][2] = 104; 
	     Data[2][1][3] = 124; 
            Data[2][1][4] = 92;
	     Split[1] = 1; MaxDays[1] = DayTemp*7 + Extra;  //DayTemp is numbers of weeks in data (4), plus one extra week to avoid numerical errors

	     S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	     FractI[2] = 0.0407;
	     Data[1][2][1] = 0;
	     Data[1][2][2] = round(0.01*158); 
            Data[1][2][3] = round(0.04*145); 
            Data[1][2][4] = round(0.13*151); 
	     Data[1][2][5] = round(0.14*122);
	     Data[2][2][1] = 163;
	     Data[2][2][2] = 158; 
            Data[2][2][3] = 145; 
            Data[2][2][4] = 151; 
            Data[2][2][5] = 122;
		

	     Split[2] = 0; MaxDays[2] = (DayTemp+1)*7 + Extra; //5 weeks in data plus 1

	     S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.0233;
	     Data[1][3][1] = 0;
	     Data[1][3][2] = round(0.03*132); 
	     Data[1][3][3] = round(0.06*128); 
	     Data[1][3][4] = round(0.17*111);
	     Data[2][3][1] = 150;
	     Data[2][3][2] = 132; 
            Data[2][3][3] = 128; 
            Data[2][3][4] = 111;
	     Split[3] = 0; MaxDays[3] = (DayTemp)*7 + Extra;

	
	     S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     FractI[4] = 0.2854;
	     Data[1][4][1] = round(0.05*466);
	     Data[1][4][2] = round(0.25*451); 
            Data[1][4][3] = round(0.50*391); 
            Data[1][4][4] = round(0.55*181);

	     Data[2][4][1] = 466;
	     Data[2][4][2] = 451; 
            Data[2][4][3] = 391; 
            Data[2][4][4] = 181;
	     Split[4] = 1; MaxDays[4] = (DayTemp)*7 + Extra;

		 S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0943;
	     Data[1][5][1] = round(0.02*1083); 
            Data[1][5][2] = round(0.03*660);
	     Data[1][5][3] = round(0.02*607); 
            Data[1][5][4] = round(0.07*558);
            Data[1][5][5] = round(0.15*498); 
	     Data[1][5][6] = round(0.17*48);
	     Data[2][5][1] = 1083; 
            Data[2][5][2] = 660;
	     Data[2][5][3] = 607; 
            Data[2][5][4] = 558;
            Data[2][5][5] = 498; 
	     Data[2][5][6] = 48;
	     Split[5] = 0; MaxDays[5] = (DayTemp+2)*7 + Extra;

	     
	     //     case 4 : *S0 = 1720.0 * 216.0/(1e4 * 1.4);
	    // *FractI = 0.0446;
	    //   maxData = 10;
	   //  *maxt = 77;
	   //  S[0] = 1; I[0] = 0;
	   //  S[1] = 88; I[1] = 1; S[2] = 146; I[2] = 9; S[3] = 126; I[3] = 6;
	  //   S[4] = 88; I[4] = 5; S[5] = 147; I[5] = 5;
	  //   S[6] = 136; I[6] = 11; S[7] = 153; I[7] = 24; S[8] = 122; I[8] = 27;
	  //   S[9] = 113; I[9] = 28;  //I think 8 and 9 are not in the published paper, so I'm starting at 7 and going backwards 
	     

	     
	   //  S0[6] = 1720.0 * 216.0/(1e4 * 1.4);
	   //  FractI[6] = 0.0446;
	   //  Data[1][6][1] = round(0.0*147); 
	   //  Data[1][6][2] = round(0.07*136);
	   //  Data[1][6][3] = round(0.02*153); 
		 
	    // Data[2][6][1] = 147; 
	    // Data[2][6][2] = 136;
	    // Data[2][6][3] = 153; 
	    // Split[6] = 0; MaxDays[6] = (DayTemp-1)*7 + Extra;
	    

	  
	     
	  //   *S0 = 2356.0 * 192.0/(1e4 * 1.4);
	  //   *FractI = 0.1034;
	  //   *maxt = 50;
	  //   maxData = 6;
	  //   Data[0] = -1.0; Data[1] = 0.07; Data[2] = 0.05; Data[3] = 0.03;
	  //   Data[4] = 0.05; Data[5] = 0.12;
	  //   S[0] = 604; S[1] = 843; S[2] = 861; S[3] = 851;
	  //   S[4] = 845;
	  //   

	   //This next one was 7, but we're making it 6 in order to only have 6 plots
	     S0[6] = 2355.0 * 192.0/(1e4 * 1.4);
	     FractI[6] = 0.1034;
	     Data[1][6][1] = round(0.035*861); 
	     Data[1][6][2] = round(0.045*851);
	     Data[1][6][3] = round(0.11*845); 
		 
	     Data[2][6][1] = 861; 
	     Data[2][6][2] = 851;
	     Data[2][6][3] = 845; 
	     Split[6] = 1; 
	     MaxDays[6] = (DayTemp-1)*7 + Extra;

	     
	     
	   //  *FractI = 0.2869;
	   //  S0 = 827.0 * 146.0/(1e4 * 1.4);
	   //  *maxt = 49;
	   //  maxData = 6;
	   //  Data[0] = -1.0; Data[1] = -1.0; Data[2] = 0.08; Data[3] = 0.02;
	   //  Data[4] = 0.07; Data[5] = 0.15;
	   //  for(i=0;i<maxData;i++){
	   //   S[i] = 450;
	   //   I[i] = round(Data[i]*S[i]);
	   //  }
	   //  	     S[0] = 1; I[0] = 0;
	   //  S[1] = 609; I[1] = 62; S[2] = 844; I[2] = 19; S[3] = 1; I[3] = -1;
	   //  S[4] = 842; I[4] = 121;
	   //


	     
	   //  S0[8] = 827.0 * 146.0/(1e4 * 1.4);
   	   //  FractI[8] = 0.08;  //the value above looks wrong, wrong, wrong...
	   //  Data[1][8][1] = round(0.075*500);  //value is missing above?  But not in published paper...500 seems conservative, given the other values
	   //  Data[1][8][2] = round(0.14*842);

		 
	   //  Data[2][8][1] = 500; 
	   //  Data[2][8][2] = 842;

	   //  Split[8] = 0; MaxDays[8] = (DayTemp-2)*7 + Extra;
	   //	

	     
		 return MaxPlot;
}
*/


/*THIS NEXT ONE IS ACTUALLY swplots8Bnml, because it has all 8 plots */  /* DATA ARE NOT CORRECTED */
/*
int swplots5Bnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval, int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 int Plot, MaxPlot = 5; //5;
		 int Wk;

		 *SData = 0;
		 *FractDead = 1;


		 *MaxWk = 6;
		 *DataInterval = 7;
                 int DayTemp = 7-2; //why 7-2??
		
	     //S[0] = 186; I[0] = 8; 
           // S[1] = 546; I[1] = 94; 
           // S[2] = 1; I[2] = -1;
	    // S[3] = 144; I[3] = 22; 
           // S[4] = 116; I[4] = 24;  
          //  S[5] = 139; I[5] = 70;
	   //  S[6] = 110; I[6] = 57;
		

		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	     FractI[1] = 0.0885;
		 //FractI[1] = 0.25;
	     Data[1][1][1] = round(0.05*141);
	     Data[1][1][2] = round(0.20*104); 
	     Data[1][1][3] = round(0.5*124); 
	     Data[1][1][4] = round(0.65*92);
	     Data[2][1][1] = 141;
	     Data[2][1][2] = 104; 
	     Data[2][1][3] = 124; 
            Data[2][1][4] = 92;
	     Split[1] = 1; MaxDays[1] = DayTemp*7 + Extra;  //DayTemp is numbers of weeks in data (4), plus one extra week to avoid numerical errors

	     S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	     FractI[2] = 0.0407;
	     Data[1][2][1] = 0;
	     Data[1][2][2] = round(0.01*158); 
            Data[1][2][3] = round(0.04*145); 
            Data[1][2][4] = round(0.13*151); 
	     Data[1][2][5] = round(0.14*122);
	     Data[2][2][1] = 163;
	     Data[2][2][2] = 158; 
            Data[2][2][3] = 145; 
            Data[2][2][4] = 151; 
            Data[2][2][5] = 122;
		

	     Split[2] = 0; MaxDays[2] = (DayTemp+1)*7 + Extra; //5 weeks in data plus 1

	     S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.0233;
	     Data[1][3][1] = 0;
	     Data[1][3][2] = round(0.03*132); 
	     Data[1][3][3] = round(0.06*128); 
	     Data[1][3][4] = round(0.17*111);
	     Data[2][3][1] = 150;
	     Data[2][3][2] = 132; 
            Data[2][3][3] = 128; 
            Data[2][3][4] = 111;
	     Split[3] = 0; MaxDays[3] = (DayTemp)*7 + Extra;

	
	     S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     FractI[4] = 0.2854;
	     Data[1][4][1] = round(0.05*466);
	     Data[1][4][2] = round(0.25*451); 
            Data[1][4][3] = round(0.50*391); 
            Data[1][4][4] = round(0.55*181);

	     Data[2][4][1] = 466;
	     Data[2][4][2] = 451; 
            Data[2][4][3] = 391; 
            Data[2][4][4] = 181;
	     Split[4] = 1; MaxDays[4] = (DayTemp)*7 + Extra;

		 S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0943;
	     Data[1][5][1] = round(0.02*1083); 
            Data[1][5][2] = round(0.03*660);
	     Data[1][5][3] = round(0.02*607); 
            Data[1][5][4] = round(0.07*558);
            Data[1][5][5] = round(0.15*498); 
	     Data[1][5][6] = round(0.17*48);
	     Data[2][5][1] = 1083; 
            Data[2][5][2] = 660;
	     Data[2][5][3] = 607; 
            Data[2][5][4] = 558;
            Data[2][5][5] = 498; 
	     Data[2][5][6] = 48;
	     Split[5] = 0; MaxDays[5] = (DayTemp+2)*7 + Extra;

	     
	     //     case 4 : *S0 = 1720.0 * 216.0/(1e4 * 1.4);
	    // *FractI = 0.0446;
	    //   maxData = 10;
	   //  *maxt = 77;
	   //  S[0] = 1; I[0] = 0;
	   //  S[1] = 88; I[1] = 1; S[2] = 146; I[2] = 9; S[3] = 126; I[3] = 6;
	  //   S[4] = 88; I[4] = 5; S[5] = 147; I[5] = 5;
	 //    S[6] = 136; I[6] = 11; S[7] = 153; I[7] = 24; S[8] = 122; I[8] = 27;
	 //    S[9] = 113; I[9] = 28;  //I think 8 and 9 are not in the published paper, so I'm starting at 7 and going backwards 

	     

	     S0[6] = 1720.0 * 216.0/(1e4 * 1.4);
	     FractI[6] = 0.0446;
	     Data[1][6][1] = round(0.0*147); 
	     Data[1][6][2] = round(0.07*136);
	     Data[1][6][3] = round(0.02*153); 
		 
	     Data[2][6][1] = 147; 
	     Data[2][6][2] = 136;
	     Data[2][6][3] = 153; 
	     Split[6] = 0; MaxDays[6] = (DayTemp-1)*7 + Extra;

	  
	     
	   //  *S0 = 2356.0 * 192.0/(1e4 * 1.4);
	  //   *FractI = 0.1034;
	  //   *maxt = 50;
	  //   maxData = 6;
	  //   Data[0] = -1.0; Data[1] = 0.07; Data[2] = 0.05; Data[3] = 0.03;
	  //   Data[4] = 0.05; Data[5] = 0.12;
	  //   S[0] = 604; S[1] = 843; S[2] = 861; S[3] = 851;
	  //   S[4] = 845;
	     

	     S0[7] = 2356.0 * 192.0/(1e4 * 1.4);
	     FractI[7] = 0.1034;
	     Data[1][7][1] = round(0.035*861); 
	     Data[1][7][2] = round(0.045*851);
	     Data[1][7][3] = round(0.11*845); 
		 
	     Data[2][7][1] = 861; 
	     Data[2][7][2] = 851;
	     Data[2][7][3] = 845; 
	     Split[7] = 1; MaxDays[7] = (DayTemp-1)*7 + Extra;

	     
	     
     //    *FractI = 0.2869;
    //     S0 = 827.0 * 146.0/(1e4 * 1.4);
	//     *maxt = 49;
	 //    maxData = 6;
	 //    Data[0] = -1.0; Data[1] = -1.0; Data[2] = 0.08; Data[3] = 0.02;
	 //    Data[4] = 0.07; Data[5] = 0.15;
	 //    for(i=0;i<maxData;i++){
	 //     S[i] = 450;
	 //     I[i] = round(Data[i]*S[i]);
	 //    }
	//     	     S[0] = 1; I[0] = 0;
	//     S[1] = 609; I[1] = 62; S[2] = 844; I[2] = 19; S[3] = 1; I[3] = -1;
	//     S[4] = 842; I[4] = 121;
	   

	     S0[8] = 827.0 * 146.0/(1e4 * 1.4);
   	     FractI[8] = 0.08;  //the value above looks wrong, wrong, wrong...
	     Data[1][8][1] = round(0.075*500);  //value is missing above?  But not in published paper...500 seems conservative, given the other values
	     Data[1][8][2] = round(0.14*842);

		 
	     Data[2][8][1] = 500; 
	     Data[2][8][2] = 842;

	     Split[8] = 0; MaxDays[8] = (DayTemp-2)*7 + Extra;

	     
		 return MaxPlot;
}

*/


int swplots8(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 int Plot, MaxPlot = 8;
		 int Wk;
		 *SData = 0;
		 *FractDead = 1;
		 
		 *MaxWk = 6;
		 *DataInterval = 7;

		 S0[1] = 12550.0 * 185.0/(1e4 * 1.4);
	     FractI[1] = 0.0885;
	     Data[1][1][1] = 0.05;
	     Data[1][1][2] = 0.20; Data[1][1][3] = 0.50; Data[1][1][4] = 0.65;
	     Split[1] = 1; MaxDays[1] = (7-2)*7 + Extra;

		 S0[2] = 314.0 * 304.0/(1e4 * 1.4);
	     FractI[2] = 0.0407;
	     Data[1][2][1] = MIN2;
	     Data[1][2][2] = 0.01; Data[1][2][3] = 0.04; Data[1][2][4] = 0.13; Data[1][2][5] = 0.14;
	     Split[2] = 0; MaxDays[2] = (8-2)*7 + Extra;

		 S0[3] = 1525.0 * 216.0/(1e4 * 1.4);
	     FractI[3] = 0.0233;
	     Data[1][3][1] = MIN2;
	     Data[1][3][2] = 0.03; Data[1][3][3] = 0.06; Data[1][3][4] = 0.17;
	     Split[3] = 0; MaxDays[3] = (7-2)*7 + Extra;

		 S0[4] = 4100.0 * 213.0/(1e4 * 1.4);
	     FractI[4] = 0.2854;
	     Data[1][4][1] = 0.05;
	     Data[1][4][2] = 0.25; Data[1][4][3] = 0.50; Data[1][4][4] = 0.55;
	     Split[4] = 1; MaxDays[4] = (7-2)*7 + Extra;

		 S0[5] = 153.0 * 315.0/(1e4 * 1.4);
	     FractI[5] = 0.0943;
	     Data[1][5][1] = 0.02; Data[1][5][2] = 0.03;
	     Data[1][5][3] = 0.02; Data[1][5][4] = 0.07;Data[1][5][5] = 0.15; Data[1][5][6] = 0.17;
	     Split[5] = 0; MaxDays[5] = (9-2)*7 + Extra;
		 
         S0[6] = 1720.0 * 216.0/(1e4 * 1.4);
	     FractI[6] = 0.0446;
	     Data[1][6][1] = MIN2; Data[1][6][2] = 0.07; Data[1][6][3] = 0.02;
	     Split[6] = 0; MaxDays[6] = (6-2)*7 + Extra;
    
         S0[7] = 2356.0 * 192.0/(1e4 * 1.4);
	     FractI[7] = 0.1034;
	     Data[1][7][1] = 0.03;
	     Data[1][7][2] = 0.05; Data[1][7][3] = 0.12;
	     Split[7] = 1; MaxDays[7] = (6-2)*7 + Extra;

		 S0[8] = 827.0 * 146.0/(1e4 * 1.4);
	     FractI[8] = 0.2869;
	     Data[1][8][1] = 0.02;
	     Data[1][8][2] = 0.07; Data[1][8][3] = 0.15;
		 Split[8] = 0; MaxDays[8] = (6-2)*7 + Extra;

		 for(Plot=1;Plot<=MaxPlot;Plot++)
			for(Wk=1;Wk<=(*MaxWk);Wk++){
				Data[2][Plot][Wk] = -2;
			}

		 return MaxPlot;

}

int Mandarte(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData, int *FractDead)
{
		 double Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 3;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 1;
		 *DataInterval = 56;

		 S0[1] = 40;
	     FractI[1] = 0.27;
	     Data[1][1][1] = 0.40;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 100;
	     FractI[2] = 0.27;
	     Data[1][2][1] = 0.50;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 400;
	     FractI[3] = 0.50;
	     Data[1][3][1] = 0.70;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		  for(Plot=1;Plot<=MaxPlot;Plot++)
			for(Wk=1;Wk<=(*MaxWk);Wk++)
				Data[2][Plot][Wk] = -2;

		 return MaxPlot;

}

/*
int ShepherdControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 3;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

		 S0[1] = 180;  //C7, aka "Plot 4"
	     FractI[1] = 0.01;
	     Data[1][1][1] = 0.00; Data[1][1][2] = 0.0; Data[1][1][3] = 0.0; Data[1][1][4] = 0.0;
		 Data[1][1][5] = 0.0; Data[1][1][6] = 0.06; Data[1][1][7] = 0.10; Data[1][1][8] = 0.55;
	
	 Data[2][1][1] = 140; Data[2][1][2] = 135; Data[2][1][3] = 150;
		 Data[2][1][4] = 220.0; Data[2][1][5]  = 160.0; Data[2][1][6] = 150.0; Data[2][1][7] = 80;
		 Data[2][1][8] = 40;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 18;  //C8, aka "plot 16"
	     FractI[2] = 0.01;
	     Data[1][2][1] = 0.0; Data[1][2][2] = 0.00; Data[1][2][3] = 0.00; Data[1][2][4] = 0.0;
		 Data[1][2][5] = 0.024; Data[1][2][6] = 0.0; Data[1][2][7] = 0.0; Data[1][2][8] = 0.0;

		 Data[2][2][1] = 10; Data[2][2][2] = 16; Data[2][2][3] = 8;
		 Data[2][2][4] = 9; Data[2][2][5]  = 15.0; Data[2][2][6] = 18; Data[2][2][7] = 30;
		 Data[2][2][8] = 12;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 80;  //C9, aka "plot 19"
	     FractI[3] = 0.015;

	     Data[1][3][1] = 0.00; Data[1][3][2] = 0.0; Data[1][3][3] = 0.015; Data[1][3][4] = 0.005;
		 Data[1][3][5] = 0.022; Data[1][3][6] = 0.017; Data[1][3][7] = 0.138; Data[1][3][8] = 0.474;

		 Data[2][3][1] = 85; Data[2][3][2] = 45; Data[2][3][3] = 85;
		 Data[2][3][4] = 110; Data[2][3][5]  = 55; 
		 Data[2][3][6] = 80; Data[2][3][7] = 50; Data[2][3][8] = 45;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;


		 return MaxPlot;

}
*/

/*

int ShepherdControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 3;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

//	     Data[1][1][1] = 0.00; Data[1][1][2] = 0.0; Data[1][1][3] = 0.0; Data[1][1][4] = 0.0;
//		 Data[1][1][5] = 0.0; Data[1][1][6] = 0.06; Data[1][1][7] = 0.10; Data[1][1][8] = 0.55;

		 S0[1] = 180;  //C7
	     FractI[1] = 0.01;
	     Data[1][1][1] = 0.00; Data[1][1][2] = 0.0; Data[1][1][3] = 0.0; Data[1][1][4] = 0.0;
		 Data[1][1][5] = 0.0; Data[1][1][6] = 0.06; Data[1][1][7] = 0.10; Data[1][1][8] = 0.55;
		 Data[2][1][1] = 140; Data[2][1][2] = 135; Data[2][1][3] = 150;
		 Data[2][1][4] = 220.0; Data[2][1][5]  = 160.0; Data[2][1][6] = 150.0; Data[2][1][7] = 80;
		 Data[2][1][8] = 40;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;



		 S0[2] = 18;  //C8
	     FractI[2] = 0.01;
	     Data[1][2][1] = 0.0; Data[1][2][2] = 0.00; Data[1][2][3] = 0.02; Data[1][2][4] = 0.0;
		 Data[1][2][5] = 0.024; Data[1][2][6] = 0.0; Data[1][2][7] = 0.0; Data[1][2][8] = 0.0;
		 Data[2][2][1] = 10; Data[2][2][2] = 16; Data[2][2][3] = 8;
		 Data[2][2][4] = 9; Data[2][2][5]  = 15.0; Data[2][2][6] = 18; Data[2][2][7] = 30;
		 Data[2][2][8] = 12;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

//	     Data[1][3][1] = 0.00; Data[1][3][2] = 0.0; Data[1][3][3] = 0.015; Data[1][3][4] = 0.005;
//		 Data[1][3][5] = 0.022; Data[1][3][6] = 0.017; Data[1][3][7] = 0.138; Data[1][3][8] = 0.474;

		 S0[3] = 80;  //C9
	     FractI[3] = 0.01;
	     Data[1][3][1] = 0.00; Data[1][3][2] = 0.0; Data[1][3][3] = 0.015; Data[1][3][4] = 0.005;
		 Data[1][3][5] = 0.022; Data[1][3][6] = 0.017; Data[1][3][7] = 0.138; Data[1][3][8] = 0.474;
		 Data[2][3][1] = 85; Data[2][3][2] = 45; Data[2][3][3] = 85;
		 Data[2][3][4] = 110; Data[2][3][5]  = 55; 
		 Data[2][3][6] = 80; Data[2][3][7] = 50; Data[2][3][8] = 45;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;


		 return MaxPlot;

}

*/



int ShepherdAll(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 9;
		 int Plot, MaxPlot = 9;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 9;
		 *DataInterval = 7;
		 float InitFICon = -1;



		 S0[1] = 172.5;  //C7
	     FractI[1] = InitFICon;
Data[1][1][1] = 0.0; Data[1][1][2] = 0.0; Data[1][1][3] = 0.0;
		 Data[1][1][4] = 0.0; Data[1][1][5] = 0.00; Data[1][1][6] = 0.006; Data[1][1][7] = 0.102; Data[1][1][8] = 0.536;
		 
Data[2][1][1] = 150; Data[2][1][2] = 140;
		 Data[2][1][3] = 220.0; Data[2][1][4]  = 140.0; Data[2][1][5] = 100.0; Data[2][1][6] = 120;
		 Data[2][1][7] = 90; Data[2][1][8] = 60;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;



		 S0[2] = 13;  //C8
	     FractI[2] = InitFICon;
Data[1][2][1] = 0.00; Data[1][2][2] = 0.0; Data[1][2][3] = 0.0;
		 Data[1][2][4] = 0.0; Data[1][2][5] = 0.024; Data[1][2][6] = 0.0; Data[1][2][7] = 0.0; Data[1][2][8] = 0.0;

Data[2][2][1] = 10; 
Data[2][2][2] = 16;
		 Data[2][2][3] = 9; Data[2][2][4]  = 9; Data[2][2][5] = 15; Data[2][2][6] = 18;
		 Data[2][2][7] = 30; Data[2][2][8] = 17;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 70;  //C9
	     FractI[3] = InitFICon;
Data[1][3][1] = 0.0; Data[1][3][2] = 0.00; Data[1][3][3] = 0.015;
		 Data[1][3][4] = 0.005; Data[1][3][5] = 0.022; Data[1][3][6] = 0.017; Data[1][3][7] = 0.138; Data[1][3][8] = 0.474;
		 
Data[2][3][1] = 80; Data[2][3][2] = 40;
		 Data[2][3][3] = 85; Data[2][3][4]  = 110; 
		 Data[2][3][5] = 65; Data[2][3][6] = 75; Data[2][3][7] = 60; Data[2][3][8] = 55;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 240;  //G5
	     FractI[4] = 0.149;
Data[1][4][1] = 0.233;
		 Data[1][4][2] = 0.291; Data[1][4][3] = 0.625; Data[1][4][4] = 0.524; 
Data[1][4][5] = 0.585;
Data[1][4][6] = 0.849;

Data[2][4][1] = 420;
		 Data[2][4][2] = 360; Data[2][4][3]  = 300; 
		 Data[2][4][4] = 100; Data[2][4][5] = 50; Data[2][4][6] = 20;
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

		 S0[5] = 251.25;  //G6
	     FractI[5] = 0.176;
Data[1][5][1] = 0.205;
		 Data[1][5][2] = 0.183; Data[1][5][3] = 0.244; Data[1][5][4] = 0.548; 
Data[1][5][5] = 0.835;
Data[1][5][6] = 0.831;

Data[2][5][1] = 430;
		 Data[2][5][2] = 320; Data[2][5][3]  = 210; 
		 Data[2][5][4] = 170; Data[2][5][5] = 100; Data[2][5][6] = 100;
	     Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

		 S0[6] = 29.5;  //A1
	     FractI[6] = 0.0112;
Data[1][6][1] = 0.034; Data[1][6][2] = 0.309;
		 Data[1][6][3] = 0.083; Data[1][6][4] = 0.556; Data[1][6][5] = 0.327; 
Data[1][6][6] = 0.63;
Data[1][6][7] = 0.972;

Data[2][6][1] = 26; Data[2][6][2] = 25;
		 Data[2][6][3] = 24; Data[2][6][4]  = 18; 
		 Data[2][6][5] = 12; Data[2][6][6] = 8; Data[2][6][7] = 10;
	     Split[6] = 0; MaxDays[6] = Weeks*7 + Extra;

		 S0[7] = 221.25;  //A2
	     FractI[7] = 0.0131;
Data[1][7][1] = 0.084; Data[1][7][2] = 0.275;
		 Data[1][7][3] = 0.176; Data[1][7][4] = 0.987; Data[1][7][5] = 0.942; 
Data[1][7][6] = 0.969;
Data[1][7][7] = -2;

Data[2][7][1] = 190; Data[2][7][2] = 280;
		 Data[2][7][3] = 260; Data[2][7][4]  = 150; 
		 Data[2][7][5] = 140; Data[2][7][6] = 30; Data[2][7][7] = 10;
	     Split[7] = 0; MaxDays[7] = Weeks*7 + Extra;

		 S0[8] = 88.75;  //A3
	     FractI[8] = 0.04;
Data[1][8][1] = 0.0; Data[1][8][2] = 0.006;
		 Data[1][8][3] = 0.162; Data[1][8][4] = 0.72; Data[1][8][5] = 0.932; 
Data[1][8][6] = 0.952;
Data[1][8][7] = 1.0;

Data[2][8][1] = 55; Data[2][8][2] = 140;
		 Data[2][8][3] = 160; Data[2][8][4]  = 80; 
		 Data[2][8][5] = 70; Data[2][8][6] = 50; Data[2][8][7] = 20;
	     Split[8] = 0; MaxDays[8] = Weeks*7 + Extra;

		 S0[9] = 77.5;  //A4
	     FractI[9] = 0.032;
Data[1][9][1] = 0.025; Data[1][9][2] = 0.041;
		 Data[1][9][3] = 0.066; Data[1][9][4] = 0.234; Data[1][9][5] = 0.25; 
Data[1][9][6] = 0.556;
Data[1][9][7] = 0.772;

Data[2][9][1] = 40; Data[2][9][2] = 100;
		 Data[2][9][3] = 40; Data[2][9][4]  = 30; 
		 Data[2][9][5] = 110; Data[2][9][6] = 95; Data[2][9][7] = 80;
	     Split[9] = 0; MaxDays[9] = Weeks*7 + Extra;


		 return MaxPlot;

}

int ShepherdSpray(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
                 int i;
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 6;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;
		 float InitFICon = -1;

		 S0[1] = 420; //240;  //G5
	     FractI[1] = 0.149;
Data[1][1][1] = 0.233;
		 Data[1][1][2] = 0.291; Data[1][1][3] = 0.625; Data[1][1][4] = 0.524; 
Data[1][1][5] = 0.585;
Data[1][1][6] = 0.849;

Data[2][1][1] = 420;
		 Data[2][1][2] = 360; Data[2][1][3]  = 300; 
		 Data[2][1][4] = 100; Data[2][1][5] = 50; Data[2][1][6] = 20;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 430; //251.25;  //G6
	     FractI[2] = 0.176;
Data[1][2][1] = 0.205;
		 Data[1][2][2] = 0.183; Data[1][2][3] = 0.244; Data[1][2][4] = 0.548; 
Data[1][2][5] = 0.835;
Data[1][2][6] = 0.831;

Data[2][2][1] = 430;
		 Data[2][2][2] = 320; Data[2][2][3]  = 210; 
		 Data[2][2][4] = 170; Data[2][2][5] = 100; Data[2][2][6] = 100;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 25; //29.5;  //A1
	     FractI[3] = 0.0112;
Data[1][3][1] = 0.034; Data[1][3][2] = 0.309;
		 Data[1][3][3] = 0.083; Data[1][3][4] = 0.556; Data[1][3][5] = 0.327; 
Data[1][3][6] = 0.63;
Data[1][3][7] = 0.972;

Data[2][3][1] = 26; Data[2][3][2] = 25;
		 Data[2][3][3] = 24; Data[2][3][4]  = 18; 
		 Data[2][3][5] = 12; Data[2][3][6] = 8; Data[2][3][7] = 10;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 280; //221.25;  //A2
	     FractI[4] = 0.0131;
Data[1][4][1] = 0.084; Data[1][4][2] = 0.275;
		 Data[1][4][3] = 0.176; Data[1][4][4] = 0.987; Data[1][4][5] = 0.942; 
Data[1][4][6] = 0.969;
Data[1][4][7] = -2;

Data[2][4][1] = 190; Data[2][4][2] = 280;
		 Data[2][4][3] = 260; Data[2][4][4]  = 150; 
		 Data[2][4][5] = 140; Data[2][4][6] = 30; Data[2][4][7] = 10;
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

		 S0[5] = 140; //88.75;  //A3
	     FractI[5] = 0.04;
Data[1][5][1] = 0.0; Data[1][5][2] = 0.006;
		 Data[1][5][3] = 0.162; Data[1][5][4] = 0.72; Data[1][5][5] = 0.932; 
Data[1][5][6] = 0.952;
Data[1][5][7] = 1.0;

Data[2][5][1] = 55; Data[2][5][2] = 140;
		 Data[2][5][3] = 160; Data[2][5][4]  = 80; 
		 Data[2][5][5] = 70; Data[2][5][6] = 50; Data[2][5][7] = 20;
	     Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

/*
		 S0[6] = 100; //77.5;  //A4
	     FractI[6] = 0.032;
Data[1][6][1] = 0.025; Data[1][6][2] = 0.041;
		 Data[1][6][3] = 0.066; Data[1][6][4] = 0.234; Data[1][6][5] = 0.25; 
Data[1][6][6] = 0.556;
Data[1][6][7] = 0.772;

Data[2][6][1] = 40; Data[2][6][2] = 100;
		 Data[2][6][3] = 40; Data[2][6][4]  = 30; 
		 Data[2][6][5] = 110; Data[2][6][6] = 95; Data[2][6][7] = 80;
	     Split[6] = 0; MaxDays[6] = Weeks*7 + Extra;
*/


/*
		for(i=1;i<=6;i++)
		    S0[i] = Data[2][i][2];
*/

                 //S0[3] = Data[2][3][1];


		 return MaxPlot;

}
int ShepherdAerial(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
                 int i;
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 4;
		 int Wk;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;
		 float InitFICon = -1; //1e-6;

		 S0[1] = 29.5;  //A1
	     FractI[1] = 0.0112;
Data[1][1][1] = 0.034; Data[1][1][2] = 0.309;
		 Data[1][1][3] = 0.083; Data[1][1][4] = 0.556; Data[1][1][5] = 0.327; 
Data[1][1][6] = 0.63;
Data[1][1][7] = 0.972;

Data[2][1][1] = 26; Data[2][1][2] = 25;
		 Data[2][1][3] = 24; Data[2][1][4]  = 18; 
		 Data[2][1][5] = 12; Data[2][1][6] = 8; Data[2][1][7] = 10;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 221.25;  //A2
	     FractI[2] = 0.0131;
Data[1][2][1] = 0.084; Data[1][2][2] = 0.275;
		 Data[1][2][3] = 0.176; Data[1][2][4] = 0.987; Data[1][2][5] = 0.942; 
Data[1][2][6] = 0.969;
Data[1][2][7] = -2;

Data[2][2][1] = 190; Data[2][2][2] = 280;
		 Data[2][2][3] = 260; Data[2][2][4]  = 150; 
		 Data[2][2][5] = 140; Data[2][2][6] = 30; Data[2][2][7] = 10;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 88.75;  //A3
	     FractI[3] = 0.04;
Data[1][3][1] = 0.0; Data[1][3][2] = 0.006;
		 Data[1][3][3] = 0.162; Data[1][3][4] = 0.72; Data[1][3][5] = 0.932; 
Data[1][3][6] = 0.952;
Data[1][3][7] = 1.0;

Data[2][3][1] = 55; Data[2][3][2] = 140;
		 Data[2][3][3] = 160; Data[2][3][4]  = 80; 
		 Data[2][3][5] = 70; Data[2][3][6] = 50; Data[2][3][7] = 20;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 40.0; //77.5;  //A4
	     FractI[4] = 0.032;
Data[1][4][1] = 0.025; Data[1][4][2] = 0.041;
		 Data[1][4][3] = 0.066; Data[1][4][4] = 0.234; Data[1][4][5] = 0.25; 
Data[1][4][6] = 0.556;
Data[1][4][7] = 0.772;

Data[2][4][1] = 40; Data[2][4][2] = 100;
		 Data[2][4][3] = 40; Data[2][4][4]  = 30; 
		 Data[2][4][5] = 110; Data[2][4][6] = 95; Data[2][4][7] = 80;
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;


/*
		for(i=1;i<=6;i++)
		    S0[i] = Data[2][i][2];
*/

                 //S0[3] = Data[2][3][1];


		 return MaxPlot;

}


int ShepherdControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 3;
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

             		 S0[1] = 172.5;  //C7
	     FractI[1] = InitFICon;
Data[1][1][1] = 0.0; Data[1][1][2] = 0.0; Data[1][1][3] = 0.0;
		 Data[1][1][4] = 0.0; Data[1][1][5] = 0.00; Data[1][1][6] = 0.006; Data[1][1][7] = 0.102; Data[1][1][8] = 0.536;
		 
Data[2][1][1] = 150; Data[2][1][2] = 140;
		 Data[2][1][3] = 220.0; Data[2][1][4]  = 140.0; Data[2][1][5] = 100.0; Data[2][1][6] = 120;
		 Data[2][1][7] = 90; Data[2][1][8] = 60;
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;



		 S0[2] = 13;  //C8
	     FractI[2] = InitFICon;
Data[1][2][1] = 0.00; Data[1][2][2] = 0.0; Data[1][2][3] = 0.0;
		 Data[1][2][4] = 0.0; Data[1][2][5] = 0.024; Data[1][2][6] = 0.0; Data[1][2][7] = 0.0; Data[1][2][8] = 0.0;

Data[2][2][1] = 10; 
Data[2][2][2] = 16;
		 Data[2][2][3] = 9; Data[2][2][4]  = 9; Data[2][2][5] = 15; Data[2][2][6] = 18;
		 Data[2][2][7] = 30; Data[2][2][8] = 17;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 70;  //C9
	     FractI[3] = InitFICon;
Data[1][3][1] = 0.0; Data[1][3][2] = 0.00; Data[1][3][3] = 0.015;
		 Data[1][3][4] = 0.005; Data[1][3][5] = 0.022; Data[1][3][6] = 0.017; Data[1][3][7] = 0.138; Data[1][3][8] = 0.474;
		 
Data[2][3][1] = 80; Data[2][3][2] = 40;
		 Data[2][3][3] = 85; Data[2][3][4]  = 110; 
		 Data[2][3][5] = 65; Data[2][3][6] = 75; Data[2][3][7] = 60; Data[2][3][8] = 55;
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 return MaxPlot;

}
/*
int MoreauControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 4;
		 int Plot, MaxPlot = 3;
		 int Wk;
                 double InitFICon = 1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

             		 S0[1] = 172.5;  //01_C1, avg of first 2 wks
	     FractI[1] = InitFICon;
Data[1][1][1] = 0.0; Data[1][1][2] = 0.0; Data[1][1][3] = 0.2;
		 Data[1][1][4] = 0.419; 
		 
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 226.23;  //02_C1, avg of first 3 wks
	     FractI[2] = InitFICon;
Data[1][2][1] = 0.00; Data[1][2][2] = 0.0; Data[1][2][3] = 0.0;
		 Data[1][2][4] = 0.14; 

	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 150.0;  //02_C2
	     FractI[3] = InitFICon;
Data[1][3][1] = 0.0; Data[1][3][2] = 0.00; Data[1][3][3] = 0.25;
		 Data[1][3][4] = 0.393939; 
		 
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 return MaxPlot;

}
*/
int MoreauControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 5;
		 double Weeks = 5;
		 int Plot, MaxPlot = 3;
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

             		 S0[1] = 172.5;  //01_C1, avg of first 2 wks
	     FractI[1] = InitFICon;
Data[1][1][1] = 0.0; Data[1][1][2] = 0.2; Data[1][1][3] = 0.419;
		 
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 226.23;  //02_C1, avg of first 3 wks
	     FractI[2] = InitFICon;
Data[1][2][1] = 0.00; Data[1][2][2] = 0.0; Data[1][2][3] = 0.14;

	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 150.0;  //02_C2
	     FractI[3] = InitFICon;
Data[1][3][1] = 0.0; Data[1][3][2] = 0.25; Data[1][3][3] = 0.393939;
		 
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 return MaxPlot;

}

int OtvosControl(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 4;
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

             		 S0[1] = 197.5;  //C1
	     FractI[1] = InitFICon;
Data[1][1][1] = -1; Data[1][1][2] = -1; Data[1][1][3] = 0.0;
		 Data[1][1][4] = -1; Data[1][1][5] = 0.043; 
Data[1][1][6] = 0.007; Data[1][1][7] = 0.091; 
		 
Data[2][1][1] = 211.25; Data[2][1][2] = 239.66;
		 Data[2][1][3] = 209.57; Data[2][1][4]  = 121.22; 
Data[2][1][5] = 25.07; Data[2][1][6] = 19.07;
		 Data[2][1][7] = -1; 
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;



		 S0[2] = 136.88;  //C2
	     FractI[2] = InitFICon;
Data[1][2][1] = -1; Data[1][2][2] = -1; Data[1][2][3] = 0.0;
		 Data[1][2][4] = 0.007; Data[1][2][5] = 0.013; 
Data[1][2][6] = 0.014; 

Data[2][2][1] = 168.23; 
Data[2][2][2] = 137.54;
		 Data[2][2][3] = 130.36; Data[2][2][4]  = 98.29; 
Data[2][2][5] = 53.74; Data[2][2][6] = 28.66;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 360.59;  //C3
	     FractI[3] = InitFICon;
Data[1][3][1] = -1; Data[1][3][2] = -1; Data[1][3][3] = 0.0;
		 Data[1][3][4] = -1; Data[1][3][5] = 0.103; 
Data[1][3][6] = 0.11; Data[1][3][7] = 0.434; 
		 
Data[2][3][1] = 272.67; Data[2][3][2] = 214.20;
		 Data[2][3][3] = 120.97; Data[2][3][4]  = 103.21; 
		 Data[2][3][5] = 26.77; Data[2][3][6] = 16.57; 
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 81.66;  //C4
	     FractI[4] = 100.0*InitFICon;
Data[1][4][1] = -1; Data[1][4][2] = -1; Data[1][4][3] = -1;
		 Data[1][4][4] = -1; Data[1][4][5] = 0.013; 
Data[1][4][6] = 0.049; Data[1][4][7] = 0.234; 
		 
Data[2][4][1] = 55.86; Data[2][4][2] = 81.66;
		 Data[2][4][3] = 51.20; Data[2][4][4]  = 37.51; 
		 Data[2][4][5] = 27.53; Data[2][4][6] = 29.69; 
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;



		 return MaxPlot;

}

int OtvosConBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 3; //used to be 4, then I nuked plot C2
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
S0[1] = 216.14;  //C1 - avg first 3 weeks
FractI[1] = InitFICon;
Data[1][1][1] = -1; 
Data[1][1][2] = -1; 
Data[1][1][3] = 0;  //out of 148
Data[1][1][4] = -1; 
Data[1][1][5] = 7;  //out of 161 
Data[1][1][6] = 1; //out of 148 
Data[1][1][7] = 16; //out of 176
		 
Data[2][1][1] = 1; 
Data[2][1][2] = 1;
Data[2][1][3] = 148; 
Data[2][1][4]  = 1; 
Data[2][1][5] = 161; 
Data[2][1][6] = 148;
Data[2][1][7] = 176; 

Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;


//S0[2] = 360.59;  //C3 - 1st wk
S0[2] = 282.49;  //C3 - avg of 1st 3 wks
FractI[2] = InitFICon;
Data[1][2][1] = -1; 
Data[1][2][2] = -1; 
Data[1][2][3] = 0; //out of 150
Data[1][2][4] = -1; 
Data[1][2][5] = 16; //out of 155
Data[1][2][6] = 17; //out of 155
Data[1][2][7] = 56; //out of 129 
		 
Data[2][2][1] = 1; 
Data[2][2][2] = 1;
Data[2][2][3] = 150; 
Data[2][2][4]  = 1; 
Data[2][2][5] = 155; 
Data[2][2][6] = 155; 
Data[2][2][7] = 129; 
Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

//S0[4] = 81.66;  //C4
S0[3] = 62.91;  //C4
FractI[3] = InitFICon;
Data[1][3][1] = -1; 
Data[1][3][2] = -1; 
Data[1][3][3] = -1;
Data[1][3][4] = -1; 
Data[1][3][5] = 2; //out of 150 
Data[1][3][6] = 7; //out of 144 
Data[1][3][7] = 33; //out of 141 
		 
Data[2][3][1] = 1; 
Data[2][3][2] = 1;
Data[2][3][3] = 1; 
Data[2][3][4]  = 1; 
Data[2][3][5] = 150; 
Data[2][3][6] = 144; 
Data[2][3][7] = 141; 
Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;



		 return MaxPlot;

}


int OtvosShepConBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 5; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
S0[1] = 216.14;  //C1
FractI[1] = 7.0/161.0;
Data[1][1][1] = 1; //out of 148 
Data[1][1][2] = 16; //out of 176
		 
Data[2][1][1] = 148;
Data[2][1][2] = 176; 

Split[1] = 0; MaxDays[1] = 2*7 + Extra; //3 is 3 weeks

//S0[2] = 360.59;  //C3 - 1st wk
S0[2] = 282.49;  //C3 - avg of 1st 3 wks
FractI[2] = 16/155;
Data[1][2][1] = 17; //out of 155
Data[1][2][2] = 56; //out of 129 
		 
Data[2][2][1] = 155; 
Data[2][2][2] = 129; 
Split[2] = 0; MaxDays[2] = 2*7 + Extra;

S0[3] = 68.76;  //C4
FractI[3] = 2/150;
Data[1][3][1] = 7; //out of 144 
Data[1][3][2] = 33; //out of 141 

Data[2][3][1] = 144; 
Data[2][3][2] = 141; 
Split[3] = 0; MaxDays[3] = 2*7 + Extra;

S0[4] = 172.5;  //C7 - Shepherd
//S0[4] = 220.0;  //C7 - Shepherd - this is the max density for whole season, which happened in week 4 - FOR S2!

FractI[4] = 1/174; 
Data[1][4][1] = 21; //out of 206 
Data[1][4][2] = 118; //out of 220 
		 
Data[2][4][1] = 206; 
Data[2][4][2] = 220; 
Split[4] = 0; MaxDays[4] = 2*7 + Extra;

S0[5] = 70;  //C9 - Shepherd
//S0[5] = 85;  //C9 - Shepherd - this is just for week 3 - APPARENTLY UNUSED?
//S0[5] = 110.0;  //C9 - Shepherd - this is the max density for whole season, which happened in week 4 - FOR S2!
FractI[5] = 3/199;
Data[1][5][1] = 1; //out of 193
Data[1][5][2] = 4; //out of 178 
Data[1][5][3] = 3; //out of 179
Data[1][5][4] = 26; //out of 188 
Data[1][5][5] = 45; //out of 95 
		 
Data[2][5][1] = 193; 
Data[2][5][2]  = 178; 
Data[2][5][3] = 179; 
Data[2][5][4] = 188; 
Data[2][5][5] = 95; 
Split[5] = 0; MaxDays[5] = 5*7 + Extra;

		 return MaxPlot;

}


//oN AUGUST 11 2015, THE OTVOS DATA WERE REVISED SO THAT THE INITIAL HOST DENSITY IS THE DENSITY THAT WAS RECORDED ON THE SAME DAY AS THE FRACTION INFECTED, AS OPPOSED TO THE INITIAL HOST DENSITY BEING THE VALUE PRE-SPRAY
//FOR SOME OF THE CONTROL PLOTS, THE FRACTION INFECTED ONLY GOES ABOVE ZERO VERY LATE IN THE SEASON, BY WHICH TIME THE HOST DENSITY HAS FALLEN STRONGLY
int OtvosAllBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 7; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 6;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
//S0[1] = 216.14;  //C1 - first 3 weeks
//S0[1] = 204.37; //SOME KIND OF AVERAGE?
S0[1] = 25.07; //FIRST WEEK WHEN FRACTION INFECTED WAS ABOVE 0
FractI[1] = 7.0/161.0;
Data[1][1][1] = 1; //out of 148 
Data[1][1][2] = 16; //out of 176
		 
Data[2][1][1] = 148;
Data[2][1][2] = 176; 

Split[1] = 0; MaxDays[1] = 3*7 + Extra; //3 is 3 weeks

//S0[2] = 360.59;  //C3 - 1st wk
//S0[2] = 282.49;  //C3 - avg of 1st 3 wks
//S0[2] = 316.63; //SOME KIND OF AVERAGE? 
S0[2] = 26.77; //FIRST WEEK WHEN FI > 0
FractI[2] = 16.0/155.0;
Data[1][2][1] = 17; //out of 155
Data[1][2][2] = 56; //out of 129 
		 
Data[2][2][1] = 155; 
Data[2][2][2] = 129; 
Split[2] = 0; MaxDays[2] = 3*7 + Extra;

//S0[3] = 68.76;  //C4
//S0[3] = 55.86;
S0[3] = 27.53; //FIRST WEEK WHEN FI > 0
 

FractI[3] = 2.0/150.0;
Data[1][3][1] = 7; //out of 144 
Data[1][3][2] = 33; //out of 141 

Data[2][3][1] = 144; 
Data[2][3][2] = 141; 
Split[3] = 0; MaxDays[3] = 3*7 + Extra;

S0[4] = 172.65;  //T1 - Otvos - avg 1st 2 weeks
S0[4] = 185.73; //First week when FI > 0
FractI[4] = 18.0/194.0; 
Data[1][4][1] = 3; //out of 197 
Data[1][4][2] = 80; //out of 123 
Data[1][4][3] = 61; //out of 89 
Data[1][4][4] = 125; //out of 142 
Data[1][4][5] = 75; //out of 103 
/*
194	18
197	3
123	80
89	61
142	125
103	75
*/
		 
Data[2][4][1] = 197; 
Data[2][4][2] = 123;
Data[2][4][3] = 89; 
Data[2][4][4] = 142;
Data[2][4][5] = 103; 
Split[4] = 0; MaxDays[4] = 6*7 + Extra;

//S0[5] = 124.08;  //T2 - Otvos - avg 1st 2 weeks
S0[5] = 133.61;  //FIRST WEEK WHEN FI > 0

FractI[5] = 16.0/178.0;
Data[1][5][1] = 10; //out of 169
Data[1][5][2] = 66; //out of 102 
Data[1][5][3] = 87; //out of 103
Data[1][5][4] = 80; //out of 104 
Data[1][5][5] = 121; //out of 133 

/*
178	16
169	10
102	66
103	87
104	80
133	121
*/
		 
Data[2][5][1] = 169; 
Data[2][5][2]  = 102; 
Data[2][5][3] = 103; 
Data[2][5][4] = 104; 
Data[2][5][5] = 133; 
Split[5] = 0; MaxDays[5] = 6*7 + Extra;

//S0[6] = 232.03;  //T3 - Otvos - avg 1st 2 weeks
S0[6] = 128.96; //FIRST WEEK WHEN FI > 0
FractI[6] = 28.0/191.0;
Data[1][6][1] = 53; //out of 172
Data[1][6][2] = 56; //out of 64 
Data[1][6][3] = 130; //out of 130
Data[1][6][4] = 58; //out of 60 
Data[1][6][5] = 58; //out of 110 

/*
191	28
172	53
64	56
130	130
60	58
110	58
*/
		 
Data[2][6][1] = 172; 
Data[2][6][2]  = 64; 
Data[2][6][3] = 130; 
Data[2][6][4] = 60; 
Data[2][6][5] = 58; 
Split[6] = 0; MaxDays[6] = 6*7 + Extra;

//S0[7] = 34.27;  //T4 - Otvos - avg 1st 2 weeks
S0[7] = 26.76;  //FIRST WEEK WHEN FI > 0 
FractI[7] = 31.0/111.0;
Data[1][7][1] = 23; //out of 103
Data[1][7][2] = 20; //out of 150 
Data[1][7][3] = 57; //out of 112
Data[1][7][4] = 128; //out of 131 
Data[1][7][5] = 78; //out of 88 
Data[1][7][6] = 12; //out of 20 

/*
111	31
103	23
150	20
112	57
131	128
88	78
20	12
*/
		 
Data[2][7][1] = 103; 
Data[2][7][2]  = 150; 
Data[2][7][3] = 112; 
Data[2][7][4] = 131; 
Data[2][7][5] = 88; 
Data[2][7][6] = 20; 
Split[7] = 0; MaxDays[7] = 7*7 + Extra;


return MaxPlot;

}


 
int OldOtvosAllBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)   //OBSOLETE AS OF 11 AUGUST 2015
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 7; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 6;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
//S0[1] = 216.14;  //C1 - first 3 weeks
S0[1] = 204.37;
FractI[1] = 7.0/161.0;
Data[1][1][1] = 1; //out of 148 
Data[1][1][2] = 16; //out of 176
		 
Data[2][1][1] = 148;
Data[2][1][2] = 176; 

Split[1] = 0; MaxDays[1] = 3*7 + Extra; //3 is 3 weeks

//S0[2] = 360.59;  //C3 - 1st wk
//S0[2] = 282.49;  //C3 - avg of 1st 3 wks
S0[2] = 316.63;
FractI[2] = 16.0/155.0;
Data[1][2][1] = 17; //out of 155
Data[1][2][2] = 56; //out of 129 
		 
Data[2][2][1] = 155; 
Data[2][2][2] = 129; 
Split[2] = 0; MaxDays[2] = 3*7 + Extra;

S0[3] = 68.76;  //C4
S0[3] = 55.86;

FractI[3] = 2.0/150.0;
Data[1][3][1] = 7; //out of 144 
Data[1][3][2] = 33; //out of 141 

Data[2][3][1] = 144; 
Data[2][3][2] = 141; 
Split[3] = 0; MaxDays[3] = 3*7 + Extra;

S0[4] = 172.65;  //T1 - Otvos - avg 1st 2 weeks
FractI[4] = 18.0/194.0; 
Data[1][4][1] = 3; //out of 197 
Data[1][4][2] = 80; //out of 123 
Data[1][4][3] = 61; //out of 89 
Data[1][4][4] = 125; //out of 142 
Data[1][4][5] = 75; //out of 103 
/*
194	18
197	3
123	80
89	61
142	125
103	75
*/
		 
Data[2][4][1] = 197; 
Data[2][4][2] = 123;
Data[2][4][3] = 89; 
Data[2][4][4] = 142;
Data[2][4][5] = 103; 
Split[4] = 0; MaxDays[4] = 6*7 + Extra;

S0[5] = 124.08;  //T2 - Otvos - avg 1st 2 weeks
FractI[5] = 16.0/178.0;
Data[1][5][1] = 10; //out of 169
Data[1][5][2] = 66; //out of 102 
Data[1][5][3] = 87; //out of 103
Data[1][5][4] = 80; //out of 104 
Data[1][5][5] = 121; //out of 133 

/*
178	16
169	10
102	66
103	87
104	80
133	121
*/
		 
Data[2][5][1] = 169; 
Data[2][5][2]  = 102; 
Data[2][5][3] = 103; 
Data[2][5][4] = 104; 
Data[2][5][5] = 133; 
Split[5] = 0; MaxDays[5] = 6*7 + Extra;

S0[6] = 232.03;  //T3 - Otvos - avg 1st 2 weeks
FractI[6] = 28.0/191.0;
Data[1][6][1] = 53; //out of 172
Data[1][6][2] = 56; //out of 64 
Data[1][6][3] = 130; //out of 130
Data[1][6][4] = 58; //out of 60 
Data[1][6][5] = 58; //out of 110 

/*
191	28
172	53
64	56
130	130
60	58
110	58
*/
		 
Data[2][6][1] = 172; 
Data[2][6][2]  = 64; 
Data[2][6][3] = 130; 
Data[2][6][4] = 60; 
Data[2][6][5] = 58; 
Split[6] = 0; MaxDays[6] = 6*7 + Extra;

S0[7] = 34.27;  //T4 - Otvos - avg 1st 2 weeks
FractI[7] = 31.0/111.0;
Data[1][7][1] = 23; //out of 103
Data[1][7][2] = 20; //out of 150 
Data[1][7][3] = 57; //out of 112
Data[1][7][4] = 128; //out of 131 
Data[1][7][5] = 78; //out of 88 
Data[1][7][6] = 12; //out of 20 

/*
111	31
103	23
150	20
112	57
131	128
88	78
20	12
*/
		 
Data[2][7][1] = 103; 
Data[2][7][2]  = 150; 
Data[2][7][3] = 112; 
Data[2][7][4] = 131; 
Data[2][7][5] = 88; 
Data[2][7][6] = 20; 
Split[7] = 0; MaxDays[7] = 7*7 + Extra;


return MaxPlot;

}


int OtvosShepAllBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 7; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
//S0[1] = 216.14;  //C1 - first 3 weeks
S0[1] = 204.37;
FractI[1] = 7.0/161.0;
Data[1][1][1] = 1; //out of 148 
Data[1][1][2] = 16; //out of 176
		 
Data[2][1][1] = 148;
Data[2][1][2] = 176; 

Split[1] = 0; MaxDays[1] = 3*7 + Extra; //3 is 3 weeks

//S0[2] = 360.59;  //C3 - 1st wk
//S0[2] = 282.49;  //C3 - avg of 1st 3 wks
S0[2] = 316.63;
FractI[2] = 16.0/155.0;
Data[1][2][1] = 17; //out of 155
Data[1][2][2] = 56; //out of 129 
		 
Data[2][2][1] = 155; 
Data[2][2][2] = 129; 
Split[2] = 0; MaxDays[2] = 3*7 + Extra;

S0[3] = 68.76;  //C4
S0[3] = 55.86;

FractI[3] = 2.0/150.0;
Data[1][3][1] = 7; //out of 144 
Data[1][3][2] = 33; //out of 141 

Data[2][3][1] = 144; 
Data[2][3][2] = 141; 
Split[3] = 0; MaxDays[3] = 3*7 + Extra;

S0[4] = 172.65;  //T1 - Otvos - avg 1st 2 weeks
FractI[4] = 18.0/194.0; 
Data[1][4][1] = 3; //out of 197 
Data[1][4][2] = 80; //out of 123 
Data[1][4][3] = 61; //out of 89 
Data[1][4][4] = 125; //out of 142 
Data[1][4][5] = 75; //out of 103 
/*
194	18
197	3
123	80
89	61
142	125
103	75
*/
		 
Data[2][4][1] = 197; 
Data[2][4][2] = 123;
Data[2][4][3] = 89; 
Data[2][4][4] = 142;
Data[2][4][5] = 103; 
Split[4] = 0; MaxDays[4] = 6*7 + Extra;

S0[5] = 124.08;  //T2 - Otvos - avg 1st 2 weeks
FractI[5] = 16.0/178.0;
Data[1][5][1] = 10; //out of 169
Data[1][5][2] = 66; //out of 102 
Data[1][5][3] = 87; //out of 103
Data[1][5][4] = 80; //out of 104 
Data[1][5][5] = 121; //out of 133 

/*
178	16
169	10
102	66
103	87
104	80
133	121
*/
		 
Data[2][5][1] = 169; 
Data[2][5][2]  = 102; 
Data[2][5][3] = 103; 
Data[2][5][4] = 104; 
Data[2][5][5] = 133; 
Split[5] = 0; MaxDays[5] = 6*7 + Extra;

S0[6] = 232.03;  //T3 - Otvos - avg 1st 2 weeks
FractI[6] = 28.0/191.0;
Data[1][6][1] = 53; //out of 172
Data[1][6][2] = 56; //out of 64 
Data[1][6][3] = 130; //out of 130
Data[1][6][4] = 58; //out of 60 
Data[1][6][5] = 58; //out of 110 

/*
191	28
172	53
64	56
130	130
60	58
110	58
*/
		 
Data[2][6][1] = 172; 
Data[2][6][2]  = 64; 
Data[2][6][3] = 130; 
Data[2][6][4] = 60; 
Data[2][6][5] = 58; 
Split[6] = 0; MaxDays[6] = 6*7 + Extra;

S0[7] = 34.27;  //T4 - Otvos - avg 1st 2 weeks
FractI[7] = 31.0/111.0;
Data[1][7][1] = 23; //out of 103
Data[1][7][2] = 20; //out of 150 
Data[1][7][3] = 57; //out of 112
Data[1][7][4] = 128; //out of 131 
Data[1][7][5] = 78; //out of 88 
Data[1][7][6] = 12; //out of 20 

/*
111	31
103	23
150	20
112	57
131	128
88	78
20	12
*/
		 
Data[2][7][1] = 103; 
Data[2][7][2]  = 150; 
Data[2][7][3] = 112; 
Data[2][7][4] = 131; 
Data[2][7][5] = 88; 
Data[2][7][6] = 20; 
Split[7] = 0; MaxDays[7] = 7*7 + Extra;

S0[8] = 172.5;  //C7 - Shepherd
FractI[8] = 1.0/174.0; 
Data[1][8][1] = 21; //out of 206 
Data[1][8][2] = 118; //out of 220 
		 
Data[2][8][1] = 206; 
Data[2][8][2] = 220; 
Split[8] = 0; MaxDays[8] = 2*7 + Extra;

//S0[5] = 70;  //C9 - Shepherd
S0[9] = 85;  //C9 - Shepherd - this is just for week 4
FractI[9] = 3.0/199.0;
Data[1][9][1] = 1; //out of 193
Data[1][9][2] = 4; //out of 178 
Data[1][9][3] = 3; //out of 179
Data[1][9][4] = 26; //out of 188 
Data[1][9][5] = 45; //out of 95 
		 
Data[2][9][1] = 193; 
Data[2][9][2]  = 178; 
Data[2][9][3] = 179; 
Data[2][9][4] = 188; 
Data[2][9][5] = 95; 
Split[9] = 0; MaxDays[9] = 5*7 + Extra;


S0[10] = 420; //240;  //G5
FractI[10] = 11.0/74.0;
/*
No. Insects	%NPV	No NPV
65	1.5	1
74	14.9	11
176	23.3	41
175	29.1	51
136	62.5	85
170	52.4	89
182	58.5	106
146	84.9	124
*/

Data[1][10][1] = 41; //out of 176
Data[1][10][2] = 51; //out of 175 
Data[1][10][3] = 85; //out of 136
Data[1][10][4] = 89; //out of 170 
Data[1][10][5] = 106; //out of 182
Data[1][10][6] = 124; //out of 146

Data[2][10][1] = 176; 
Data[2][10][2]  = 175; 
Data[2][10][3] = 136; 
Data[2][10][4] = 170; 
Data[2][10][5] = 182; 
Data[2][10][6] = 146; 
Split[10] = 0; MaxDays[10] = 7*7 + Extra;

S0[11] = 430; //251.25;  //G6 Shepherd
FractI[11] = 16.0/91.0;
/*
70	1.4	1
91	17.6	16
205	20.5	42
202	18.3	37
127	24.4	31
146	54.8	80
164	83.5	137
183	83.1	152
*/

//Data[1][11][1] = 16; //out of 91 //Change in plans...
Data[1][11][1] = 42; //out of 205 
Data[1][11][2] = 37; //out of 202
Data[1][11][3] = 31; //out of 127 
Data[1][11][4] = 80; //out of 146
Data[1][11][5] = 137; //out of 164
Data[1][11][6] = 152; //out of 183

//Data[2][11][1] = 91; 
Data[2][11][1]  = 205; 
Data[2][11][2] = 202; 
Data[2][11][3] = 127; 
Data[2][11][4] = 146; 
Data[2][11][5] = 164;
Data[2][11][6] = 183; 
Split[11] = 0; MaxDays[11] = 7*7 + Extra;

S0[12] = 25; //29.5;  //A1 Shepherd
FractI[12] = 11.0/98.0;
/*
98	11.2	11
58	3.4	2
136	30.9	42
180	8.3	15
142	55.6	79
113	32.7	37
108	63.0	68
71	97.2	69
*/

Data[1][12][1] = 2; //out of 58
Data[1][12][2] = 42; //out of 136 
Data[1][12][3] = 15; //out of 180
Data[1][12][4] = 79; //out of 142 
Data[1][12][5] = 37; //out of 113
Data[1][12][6] = 68; //out of 108
Data[1][12][7] = 69; //out of 71

Data[2][12][1] = 58; 
Data[2][12][2]  = 136; 
Data[2][12][3] = 180; 
Data[2][12][4] = 142; 
Data[2][12][5] = 113; 
Data[2][12][6] = 108;
Data[2][12][7] = 71; 
Split[12] = 0; MaxDays[12] = 8*7 + Extra;

S0[13] = 280; //221.25;  //A2 - Shepherd
FractI[13] = 13.0/99.0;
/*
99	13.1	13
83	8.4	7
138	27.5	38
148	17.6	26
76	98.7	75
173	94.2	163
64	96.9	62
*/

Data[1][13][1] = 7; //out of 83
Data[1][13][2] = 38; //out of 138 
Data[1][13][3] = 26; //out of 148
Data[1][13][4] = 75; //out of 76 
Data[1][13][5] = 163; //out of 173
Data[1][13][6] = 62; //out of 64

Data[2][13][1] = 83; 
Data[2][13][2]  = 138; 
Data[2][13][3] = 148; 
Data[2][13][4] = 76; 
Data[2][13][5] = 173; 
Data[2][13][6] = 64;
Split[13] = 0; MaxDays[13] = 7*7 + Extra;

S0[14] = 140; //88.75;  //A3
FractI[14] = 4.0/100.0;
/*
100	4.0	4
91	0.0	0
180	0.6	1
142	16.2	23
143	72.0	103
110	93.2	103
124	95.2	118
55	100.0	55
*/
Data[1][14][1] = 0; //out of 91
Data[1][14][2] = 1; //out of 180 
Data[1][14][3] = 23; //out of 142
Data[1][14][4] = 103; //out of 143 
Data[1][14][5] = 103; //out of 110
Data[1][14][6] = 118; //out of 124
Data[1][14][7] = 55; //out of 55

Data[2][14][1] = 91; 
Data[2][14][2]  = 180; 
Data[2][14][3] = 142; 
Data[2][14][4] = 143; 
Data[2][14][5] = 110; 
Data[2][14][6] = 124;
Data[2][14][7] = 55;
Split[14] = 0; MaxDays[14] = 8*7 + Extra;

S0[15] = 77.5;  //A4 - Shepherd
FractI[15] = 3.0/94.0;
/*
94	3.2	3
120	2.5	3
197	4.1	8
182	6.6	12
171	23.4	40
168	25.0	42
162	55.6	90
114	77.2	88
*/

Data[1][15][1] = 3; //out of 120
Data[1][15][2] = 8; //out of 197 
Data[1][15][3] = 12; //out of 182
Data[1][15][4] = 40; //out of 171 
Data[1][15][5] = 42; //out of 168
Data[1][15][6] = 90; //out of 162
Data[1][15][7] = 88; //out of 114

Data[2][15][1] = 120; 
Data[2][15][2]  = 197; 
Data[2][15][3] = 182; 
Data[2][15][4] = 171; 
Data[2][15][5] = 168; 
Data[2][15][6] = 162;
Data[2][15][7] = 114;
Split[15] = 0; MaxDays[15] = 8*7 + Extra;
return MaxPlot;

}


int OtvosShepTmtBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 7; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 8;
		 *DataInterval = 7;

S0[1] = 172.65;  //T1 - Otvos - avg 1st 2 weeks
FractI[1] = 18.0/194.0; 
Data[1][1][1] = 3; //out of 197 
Data[1][1][2] = 80; //out of 123 
Data[1][1][3] = 61; //out of 89 
Data[1][1][4] = 125; //out of 142 
Data[1][1][5] = 75; //out of 103 
/*
194	18
197	3
123	80
89	61
142	125
103	75
*/
		 
Data[2][1][1] = 197; 
Data[2][1][2] = 123;
Data[2][1][3] = 89; 
Data[2][1][4] = 142;
Data[2][1][5] = 103; 
Split[1] = 0; MaxDays[1] = 6*7 + Extra;

S0[2] = 124.08;  //T2 - Otvos - avg 1st 2 weeks
FractI[2] = 16.0/178.0;
Data[1][2][1] = 10; //out of 169
Data[1][2][2] = 66; //out of 102 
Data[1][2][3] = 87; //out of 103
Data[1][2][4] = 80; //out of 104 
Data[1][2][5] = 121; //out of 133 

/*
178	16
169	10
102	66
103	87
104	80
133	121
*/
		 
Data[2][2][1] = 169; 
Data[2][2][2]  = 102; 
Data[2][2][3] = 103; 
Data[2][2][4] = 104; 
Data[2][2][5] = 133; 
Split[2] = 0; MaxDays[2] = 6*7 + Extra;

S0[3] = 232.03;  //T3 - Otvos - avg 1st 2 weeks
FractI[3] = 28.0/191.0;
Data[1][3][1] = 53; //out of 172
Data[1][3][2] = 56; //out of 64 
Data[1][3][3] = 130; //out of 130
Data[1][3][4] = 58; //out of 60 
Data[1][3][5] = 58; //out of 110 

/*
191	28
172	53
64	56
130	130
60	58
110	58
*/
		 
Data[2][3][1] = 172; 
Data[2][3][2]  = 64; 
Data[2][3][3] = 130; 
Data[2][3][4] = 60; 
Data[2][3][5] = 58; 
Split[3] = 0; MaxDays[3] = 6*7 + Extra;

S0[4] = 34.27;  //T4 - Otvos - avg 1st 2 weeks
FractI[4] = 31.0/111.0;
Data[1][4][1] = 23; //out of 103
Data[1][4][2] = 20; //out of 150 
Data[1][4][3] = 57; //out of 112
Data[1][4][4] = 128; //out of 131 
Data[1][4][5] = 78; //out of 88 
Data[1][4][6] = 12; //out of 20 

/*
111	31
103	23
150	20
112	57
131	128
88	78
20	12
*/
		 
Data[2][4][1] = 103; 
Data[2][4][2]  = 150; 
Data[2][4][3] = 112; 
Data[2][4][4] = 131; 
Data[2][4][5] = 88; 
Data[2][4][6] = 20; 
Split[4] = 0; MaxDays[4] = 7*7 + Extra;


S0[5] = 420; //240;  //G5 Shepherd
FractI[5] = 11.0/74.0;
/*
No. Insects	%NPV	No NPV
65	1.5	1
74	14.9	11
176	23.3	41
175	29.1	51
136	62.5	85
170	52.4	89
182	58.5	106
146	84.9	124
*/

Data[1][5][1] = 41; //out of 176
Data[1][5][2] = 51; //out of 175 
Data[1][5][3] = 85; //out of 136
Data[1][5][4] = 89; //out of 170 
Data[1][5][5] = 106; //out of 182
Data[1][5][6] = 124; //out of 146

Data[2][5][1] = 176; 
Data[2][5][2]  = 175; 
Data[2][5][3] = 136; 
Data[2][5][4] = 170; 
Data[2][5][5] = 182; 
Data[2][5][6] = 146; 
Split[5] = 0; MaxDays[5] = 7*7 + Extra;

S0[6] = 430; //251.25;  //G6 Shepherd
FractI[6] = 16.0/91.0;
/*
70	1.4	1
91	17.6	16
205	20.5	42
202	18.3	37
127	24.4	31
146	54.8	80
164	83.5	137
183	83.1	152
*/

//Data[1][6][1] = 16; //out of 91 //Change in plans...
Data[1][6][1] = 42; //out of 205 
Data[1][6][2] = 37; //out of 202
Data[1][6][3] = 31; //out of 127 
Data[1][6][4] = 80; //out of 146
Data[1][6][5] = 137; //out of 164
Data[1][6][6] = 152; //out of 183

//Data[2][11][1] = 91; 
Data[2][6][1]  = 205; 
Data[2][6][2] = 202; 
Data[2][6][3] = 127; 
Data[2][6][4] = 146; 
Data[2][6][5] = 164;
Data[2][6][6] = 183; 
Split[6] = 0; MaxDays[6] = 7*7 + Extra;

S0[7] = 25; //29.5;  //A1 Shepherd
FractI[7] = 11.0/98.0;
/*
98	11.2	11
58	3.4	2
136	30.9	42
180	8.3	15
142	55.6	79
113	32.7	37
108	63.0	68
71	97.2	69
*/

Data[1][7][1] = 2; //out of 58
Data[1][7][2] = 42; //out of 136 
Data[1][7][3] = 15; //out of 180
Data[1][7][4] = 79; //out of 142 
Data[1][7][5] = 37; //out of 113
Data[1][7][6] = 68; //out of 108
Data[1][7][7] = 69; //out of 71

Data[2][7][1] = 58; 
Data[2][7][2]  = 136; 
Data[2][7][3] = 180; 
Data[2][7][4] = 142; 
Data[2][7][5] = 113; 
Data[2][7][6] = 108;
Data[2][7][7] = 71; 
Split[7] = 0; MaxDays[7] = 8*7 + Extra;

S0[8] = 280; //221.25;  //A2 - Shepherd
FractI[8] = 13.0/99.0;
/*
99	13.1	13
83	8.4	7
138	27.5	38
148	17.6	26
76	98.7	75
173	94.2	163
64	96.9	62
*/

Data[1][8][1] = 7; //out of 83
Data[1][8][2] = 38; //out of 138 
Data[1][8][3] = 26; //out of 148
Data[1][8][4] = 75; //out of 76 
Data[1][8][5] = 163; //out of 173
Data[1][8][6] = 62; //out of 64

Data[2][8][1] = 83; 
Data[2][8][2]  = 138; 
Data[2][8][3] = 148; 
Data[2][8][4] = 76; 
Data[2][8][5] = 173; 
Data[2][8][6] = 64;
Split[8] = 0; MaxDays[8] = 7*7 + Extra;

S0[9] = 140; //88.75;  //A3
FractI[9] = 4.0/100.0;
/*
100	4.0	4
91	0.0	0
180	0.6	1
142	16.2	23
143	72.0	103
110	93.2	103
124	95.2	118
55	100.0	55
*/
Data[1][9][1] = 0; //out of 91
Data[1][9][2] = 1; //out of 180 
Data[1][9][3] = 23; //out of 142
Data[1][9][4] = 103; //out of 143 
Data[1][9][5] = 103; //out of 110
Data[1][9][6] = 118; //out of 124
Data[1][9][7] = 55; //out of 55

Data[2][9][1] = 91; 
Data[2][9][2]  = 180; 
Data[2][9][3] = 142; 
Data[2][9][4] = 143; 
Data[2][9][5] = 110; 
Data[2][9][6] = 124;
Data[2][9][7] = 55;
Split[9] = 0; MaxDays[9] = 8*7 + Extra;

S0[10] = 77.5;  //A4 - Shepherd
FractI[10] = 3.0/94.0;
/*
94	3.2	3
120	2.5	3
197	4.1	8
182	6.6	12
171	23.4	40
168	25.0	42
162	55.6	90
114	77.2	88
*/

Data[1][10][1] = 3; //out of 120
Data[1][10][2] = 8; //out of 197 
Data[1][10][3] = 12; //out of 182
Data[1][10][4] = 40; //out of 171 
Data[1][10][5] = 42; //out of 168
Data[1][10][6] = 90; //out of 162
Data[1][10][7] = 88; //out of 114

Data[2][10][1] = 120; 
Data[2][10][2]  = 197; 
Data[2][10][3] = 182; 
Data[2][10][4] = 171; 
Data[2][10][5] = 168; 
Data[2][10][6] = 162;
Data[2][10][7] = 114;
Split[10] = 0; MaxDays[10] = 8*7 + Extra;
return MaxPlot;

}




int KPBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 7; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 6;
		 *DataInterval = 7;

/* Cub Creek
Day	Virus	N
0	1	32
1	2	32
8	4	52
22	3	25
29	2	32
36	11	48
44	40	49
Covars[1][1] = 150*0.062*0.5; 
*/


S0[1] = 150.0*0.062*0.5;
FractI[1] = 1.0/32.0;

Data[1][1][1] = 2; //out of 32 
Data[1][1][2] = 4; //out of 52
Data[1][1][3] = 3; //out of 25;
Data[1][1][4] = 2; //out of 25;
Data[1][1][5] = 11; //out of 48;
Data[1][1][6] = 40; //out of 49;
		 
Data[2][1][1] = 32;
Data[2][1][2] = 52; 
Data[2][1][3] = 25;
Data[2][1][4] = 25;
Data[2][1][5] = 48;
Data[2][1][6] = 49;

Split[1] = 0; MaxDays[1] = 7*7 + Extra; //3 is 3 weeks


/* Klipchuk
Day	Virus	N
20	8	102
27	2	45
28	3	25
34	0	24
35	1	14
43	33	46

Covars[1][1] = 100*2.79*0.5; 

*/

S0[2] = 100.0*2.79*0.5;
FractI[2] = 8.0/102.0;
Data[1][2][1] = 2; //out of 45
Data[1][2][2] = 3; //out of 25
Data[1][2][3] = 0; //out of 24
Data[1][2][4] = 1; //out of 14;
Data[1][2][5] = 33; //out of 46 
		 
Data[2][2][1] = 45; 
Data[2][2][2] = 25; 
Data[2][2][3] = 24;
Data[2][2][4] = 14; 
Data[2][2][5] = 46; 

Split[2] = 0; MaxDays[2] = 6*7 + Extra;

/* 8 Mile 
Day	Virus	N
0	44	60
6	57	72
14	18	42
21	15	20
28	34	58

Covars[1][1] = 150*0.062*0.5;


*/

S0[3] = 150.0*0.062*0.5; 

FractI[3] = 44.0/60.0;

Data[1][3][1] = 57; //out of 72 
Data[1][3][2] = 18; //out of 42
Data[1][3][3] = 15; //out of 20;
Data[1][3][4] = 34; //out of 58; 

Data[2][3][1] = 72; 
Data[2][3][2] = 42; 
Data[2][3][3] = 20;
Data[2][3][4] = 58;
Split[3] = 0; MaxDays[3] = 5*7 + Extra;

return MaxPlot;

}





int oldOtvosShepConBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 5; //used to be 4, then I nuked plot C2, and added 2 Shepherd plots
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
S0[1] = 216.14;  //C1
FractI[1] = InitFICon;
Data[1][1][1] = -1; 
Data[1][1][2] = -1; 
Data[1][1][3] = 0;  //out of 148
Data[1][1][4] = -1; 
Data[1][1][5] = 7;  //out of 161 
Data[1][1][6] = 1; //out of 148 
Data[1][1][7] = 16; //out of 176
		 
Data[2][1][1] = 1; 
Data[2][1][2] = 1;
Data[2][1][3] = 148; 
Data[2][1][4]  = 1; 
Data[2][1][5] = 161; 
Data[2][1][6] = 148;
Data[2][1][7] = 176; 

Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;


//S0[2] = 360.59;  //C3 - 1st wk
S0[2] = 282.49;  //C3 - avg of 1st 3 wks
FractI[2] = InitFICon;
Data[1][2][1] = -1; 
Data[1][2][2] = -1; 
Data[1][2][3] = 0; //out of 150
Data[1][2][4] = -1; 
Data[1][2][5] = 16; //out of 155
Data[1][2][6] = 17; //out of 155
Data[1][2][7] = 56; //out of 129 
		 
Data[2][2][1] = 1; 
Data[2][2][2] = 1;
Data[2][2][3] = 150; 
Data[2][2][4]  = 1; 
Data[2][2][5] = 155; 
Data[2][2][6] = 155; 
Data[2][2][7] = 129; 
Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

S0[3] = 68.76;  //C4
FractI[3] = InitFICon;
Data[1][3][1] = -1; 
Data[1][3][2] = -1; 
Data[1][3][3] = -1;
Data[1][3][4] = -1; 
Data[1][3][5] = 2; //out of 150 
Data[1][3][6] = 7; //out of 144 
Data[1][3][7] = 33; //out of 141 
		 
Data[2][3][1] = 1; 
Data[2][3][2] = 1;
Data[2][3][3] = 1; 
Data[2][3][4]  = 1; 
Data[2][3][5] = 150; 
Data[2][3][6] = 144; 
Data[2][3][7] = 141; 
Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

S0[4] = 172.5;  //C7 - Shepherd
FractI[4] = InitFICon;
Data[1][4][1] = -1; 
Data[1][4][2] = 0; 
Data[1][4][3] = 0;
Data[1][4][4] = 0; 
Data[1][4][5] = 1; //out of 174 
Data[1][4][6] = 21; //out of 206 
Data[1][4][7] = 118; //out of 220 
		 
Data[2][4][1] = 1; 
Data[2][4][2] = 212;
Data[2][4][3] = 194; 
Data[2][4][4]  = 174; 
Data[2][4][5] = 174; 
Data[2][4][6] = 206; 
Data[2][4][7] = 220; 
Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

//S0[5] = 70;  //C9 - Shepherd
//S0[5] = 85;  //C9 - Shepherd - this is just for week 4
FractI[5] = InitFICon;
Data[1][5][1] = -1; 
Data[1][5][2] = 3; //out of 199
Data[1][5][3] = 1; //out of 193
Data[1][5][4] = 4; //out of 178 
Data[1][5][5] = 3; //out of 179
Data[1][5][6] = 26; //out of 188 
Data[1][5][7] = 45; //out of 95 
		 
Data[2][5][1] = 1; 
Data[2][5][2] = 199;
Data[2][5][3] = 193; 
Data[2][5][4]  = 178; 
Data[2][5][5] = 179; 
Data[2][5][6] = 188; 
Data[2][5][7] = 95; 
Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

		 return MaxPlot;

}


/* No idea what this next one is, apparently it's obsolet? */
int OtvosShepConBnml2(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 3;
		 int Plot, MaxPlot = 6; //used to be 5, but C2 has been added back in then I nuked plot C2
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
S0[1] = 216.14;  //C1
FractI[1] = 7.0/161.0;
Data[1][1][1] = 1; //out of 148 
Data[1][1][2] = 16; //out of 176
		 
Data[2][1][1] = 148;
Data[2][1][2] = 176; 

Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

//S0[2] = 360.59;  //C3 - 1st wk
S0[2] = 282.49;  //C3 - avg of 1st 3 wks
FractI[2] = 16.0/155.0;
Data[1][2][1] = 17; //out of 155
Data[1][2][2] = 56; //out of 129 
		 
Data[2][2][1] = 155; 
Data[2][2][2] = 129; 
Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

S0[3] = 68.76;  //C4
FractI[3] = 2.0/150.0;
Data[1][3][1] = 7; //out of 144 
Data[1][3][2] = 33; //out of 141 
		 
Data[2][3][1] = 144; 
Data[2][3][2] = 141; 
Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

S0[4] = 145.16;  //C2 - Otvos
FractI[4] = 1.0/143.0;
Data[1][4][1] = 2; //out of 150 
Data[1][4][2] = 2; //out of 143 
		 
Data[2][4][1] = 150; 
Data[2][4][2] = 143; 
Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

S0[5] = 172.5;  //C7 - Shepherd
FractI[5] = 1.0/174.0;
Data[1][5][1] = 21; //out of 206 
Data[1][5][2] = 118; //out of 220 

Data[2][5][1] = 206; 
Data[2][5][2] = 220; 
Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

//S0[6] = 70.0;  //C9 - Shepherd
S0[6] = 85.0;  //C9 - Shepherd
FractI[6] = 3.0/199.0;
Data[1][6][1] = 1; //out of 193
Data[1][6][2] = 4; //out of 178 
Data[1][6][3] = 3; //out of 179
Data[1][6][4] = 26; //out of 188 
Data[1][6][5] = 45; //out of 95 
		 
Data[2][6][1] = 193; 
Data[2][6][2]  = 178; 
Data[2][6][3] = 179; 
Data[2][6][4] = 188; 
Data[2][6][5] = 95; 
Split[6] = 0; MaxDays[6] = 6*7 + Extra;

		 return MaxPlot;

}

int oldOtvosShepConBnml2(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int 
*FractDead)
{
		 int Extra = 0;
		 double Weeks = 8;
		 int Plot, MaxPlot = 6; //used to be 5, but C2 has been added back in then I nuked plot C2
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

//S0[1] = 197.5;  //C1 - first week
S0[1] = 216.14;  //C1
FractI[1] = InitFICon;
Data[1][1][1] = -1; 
Data[1][1][2] = -1; 
Data[1][1][3] = 0;  //out of 148
Data[1][1][4] = -1; 
Data[1][1][5] = 7;  //out of 161 
Data[1][1][6] = 1; //out of 148 
Data[1][1][7] = 16; //out of 176
		 
Data[2][1][1] = 1; 
Data[2][1][2] = 1;
Data[2][1][3] = 148; 
Data[2][1][4]  = 1; 
Data[2][1][5] = 161; 
Data[2][1][6] = 148;
Data[2][1][7] = 176; 

Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

//S0[2] = 360.59;  //C3 - 1st wk
S0[2] = 282.49;  //C3 - avg of 1st 3 wks
FractI[2] = InitFICon;
Data[1][2][1] = -1; 
Data[1][2][2] = -1; 
Data[1][2][3] = 0; //out of 150
Data[1][2][4] = -1; 
Data[1][2][5] = 16; //out of 155
Data[1][2][6] = 17; //out of 155
Data[1][2][7] = 56; //out of 129 
		 
Data[2][2][1] = 1; 
Data[2][2][2] = 1;
Data[2][2][3] = 150; 
Data[2][2][4]  = 1; 
Data[2][2][5] = 155; 
Data[2][2][6] = 155; 
Data[2][2][7] = 129; 
Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

S0[3] = 68.76;  //C4
FractI[3] = InitFICon;
Data[1][3][1] = -1; 
Data[1][3][2] = -1; 
Data[1][3][3] = -1;
Data[1][3][4] = -1; 
Data[1][3][5] = 2; //out of 150 
Data[1][3][6] = 7; //out of 144 
Data[1][3][7] = 33; //out of 141 
		 
Data[2][3][1] = 1; 
Data[2][3][2] = 1;
Data[2][3][3] = 1; 
Data[2][3][4]  = 1; 
Data[2][3][5] = 150; 
Data[2][3][6] = 144; 
Data[2][3][7] = 141; 
Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

S0[4] = 145.16;  //C2 - Otvos
FractI[4] = InitFICon;
Data[1][4][1] = -1; 
Data[1][4][2] = -1; //
Data[1][4][3] = 0; //out of 200
Data[1][4][4] = -1; // 
Data[1][4][5] = 1; //out of 143
Data[1][4][6] = 2; //out of 150 
Data[1][4][7] = 2; //out of 143 
		 
Data[2][4][1] = 1; 
Data[2][4][2] = 1;
Data[2][4][3] = 200; 
Data[2][4][4]  = 1; 
Data[2][4][5] = 143; 
Data[2][4][6] = 150; 
Data[2][4][7] = 143; 
Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

S0[5] = 172.5;  //C7 - Shepherd
FractI[5] = InitFICon;
Data[1][5][1] = -1; 
Data[1][5][2] = 0; 
Data[1][5][3] = 0;
Data[1][5][4] = 0; 
Data[1][5][5] = 1; //out of 174 
Data[1][5][6] = 21; //out of 206 
Data[1][5][7] = 118; //out of 220 

Data[2][5][1] = 1; 
Data[2][5][2] = 212;
Data[2][5][3] = 194; 
Data[2][5][4]  = 174; 
Data[2][5][5] = 174; 
Data[2][5][6] = 206; 
Data[2][5][7] = 220; 
Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

//S0[6] = 70;  //C9 - Shepherd
S0[6] = 85;  //C9 - Shepherd
FractI[6] = InitFICon;
Data[1][6][1] = -1; 
Data[1][6][2] = 3; //out of 199
Data[1][6][3] = 1; //out of 193
Data[1][6][4] = 4; //out of 178 
Data[1][6][5] = 3; //out of 179
Data[1][6][6] = 26; //out of 188 
Data[1][6][7] = 45; //out of 95 
		 
Data[2][6][1] = 1; 
Data[2][6][2] = 199;
Data[2][6][3] = 193; 
Data[2][6][4]  = 178; 
Data[2][6][5] = 179; 
Data[2][6][6] = 188; 
Data[2][6][7] = 95; 
Split[6] = 0; MaxDays[6] = Weeks*7 + Extra;

		 return MaxPlot;

}


int MoreauConBnml(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 0;
		 double Weeks = 5;
		 int Plot, MaxPlot = 3; 
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 4;
		 *DataInterval = 7;

S0[1] = 558.329;  //C1 - avg of 1st 2 wks
//S0[1] = 910.917;  //C1
FractI[1] = InitFICon;
Data[1][1][1] = 0; 
Data[1][1][2] = 0; 
Data[1][1][3] = 10;
Data[1][1][4] = 21; 
		 
Data[2][1][1] = 50; 
Data[2][1][2] = 50;
Data[2][1][3] = 50; 
Data[2][1][4]  = 50; 

Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;


S0[2] = 282.49;  //avg of 1st 2 wks
//S0[2] = 167.757; //1st wk only
FractI[2] = InitFICon;
Data[1][2][1] = 0; 
Data[1][2][2] = 0; 
Data[1][2][3] = 0; //out of 150
Data[1][2][4] = 7; 
		 
Data[2][2][1] = 50; 
Data[2][2][2] = 50;
Data[2][2][3] = 50; 
Data[2][2][4]  = 50; 

Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

S0[3] = 429.362; //second week only...have to check with Chris Lucarotti about this...  
//S0[3] = 149.988;  
FractI[3] = InitFICon;
Data[1][3][1] =  0; 
Data[1][3][2] =  0; 
Data[1][3][3] =  13;
Data[1][3][4] =  20; 
		 
Data[2][3][1] = 50; 
Data[2][3][2] = 50;
Data[2][3][3] = 50; 
Data[2][3][4]  = 50; 

Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 return MaxPlot;

}


int OtvosAll(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
int i;
		 int Plot, MaxPlot = 8;
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;
		 FILE *fp;

             		 S0[1] = 197.5;  //C1
	     FractI[1] = InitFICon;
Data[1][1][1] = -1; Data[1][1][2] = -1; Data[1][1][3] = 0.0;
		 Data[1][1][4] = -1; Data[1][1][5] = 0.043; 
Data[1][1][6] = 0.007; Data[1][1][7] = 0.091; 
		 
Data[2][1][1] = 211.25; Data[2][1][2] = 239.66;
		 Data[2][1][3] = 209.57; Data[2][1][4]  = 121.22; 
Data[2][1][5] = 25.07; Data[2][1][6] = 19.07;
		 Data[2][1][7] = -1; 
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;



		 S0[2] = 136.88;  //C2
	     FractI[2] = InitFICon;
Data[1][2][1] = -1; Data[1][2][2] = -1; Data[1][2][3] = 0.0;
		 Data[1][2][4] = 0.007; Data[1][2][5] = 0.013; 
Data[1][2][6] = 0.014; Data[1][2][7] = -1;

Data[2][2][1] = 168.23; 
Data[2][2][2] = 137.54;
		 Data[2][2][3] = 130.36; Data[2][2][4]  = 98.29; 
Data[2][2][5] = 53.74; Data[2][2][6] = 28.66;
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 360.59;  //C3
	     FractI[3] = InitFICon;
Data[1][3][1] = -1; Data[1][3][2] = -1; Data[1][3][3] = 0.0;
		 Data[1][3][4] = -1; Data[1][3][5] = 0.103; 
Data[1][3][6] = 0.11; Data[1][3][7] = 0.434; 
		 
Data[2][3][1] = 272.67; Data[2][3][2] = 214.20;
		 Data[2][3][3] = 120.97; Data[2][3][4]  = 103.21; 
		 Data[2][3][5] = 26.77; Data[2][3][6] = 16.57; 
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 81.66;  //C4
	     FractI[4] = 100.0*InitFICon;
Data[1][4][1] = -1; Data[1][4][2] = -1; Data[1][4][3] = -1;
		 Data[1][4][4] = -1; Data[1][4][5] = 0.013; 
Data[1][4][6] = 0.049; Data[1][4][7] = 0.234; 
		 
Data[2][4][1] = 55.86; Data[2][4][2] = 81.66;
		 Data[2][4][3] = 51.20; Data[2][4][4]  = 37.51; 
		 Data[2][4][5] = 27.53; Data[2][4][6] = 29.69; 
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;

		 S0[5] = 194.34;  //T1 ("Plot 2")
	     FractI[5] = 0.093;
Data[1][5][1] = -1; Data[1][5][2] = -1; Data[1][5][3] = 0.015;
		 Data[1][5][4] = 0.65; Data[1][5][5] = 0.685; 
Data[1][5][6] = 0.88; Data[1][5][7] = 0.728; 

Data[2][5][1] = 150.96; Data[2][5][2] = 185.73;
		 Data[2][5][3] = 103.39; Data[2][5][4]  = 82.31; 
		 Data[2][5][5] = 41.31; Data[2][5][6] = 6.72; 
	     Split[5] = 0; MaxDays[5] = Weeks*7 + Extra;

		 S0[6] = 145.76;  //T2 ("Plot 3")
	     FractI[6] = 0.09;
Data[1][6][1] = -1; Data[1][6][2] = -1; Data[1][6][3] = 0.059;
		 Data[1][6][4] = 0.647; Data[1][6][5] = 0.845; 
Data[1][6][6] = 0.769; Data[1][6][7] = 0.910; 
		 
Data[2][6][1] = 102.39; Data[2][6][2] = 133.61;
		 Data[2][6][3] = 99.68; Data[2][6][4]  = 70.93; 
		 Data[2][6][5] = 13.05; Data[2][6][6] = 2.77; 
	     Split[6] = 0; MaxDays[6] = Weeks*7 + Extra;

		 S0[7] = 302.02;  //T3 ("Plot 4")
	     FractI[7] = 0.147;
Data[1][7][1] = -1; Data[1][7][2] = -1; Data[1][7][3] = 0.308;
		 Data[1][7][4] = 0.875; Data[1][7][5] = 1.0; 
Data[1][7][6] = 0.967; Data[1][7][7] = 0.527; 
		 
Data[2][7][1] = 162.04; Data[2][7][2] = 128.96;
		 Data[2][7][3] = 130.90; Data[2][7][4]  = 96.86; 
		 Data[2][7][5] = 5.84; Data[2][7][6] = 0.96; 
	     Split[7] = 0; MaxDays[7] = Weeks*7 + Extra;

		 S0[8] = 41.77;  //T4 ("Plot 5")
	     FractI[8] = (0.279+0.223)/2.0;
Data[1][8][1] = -1; Data[1][8][2] = -1; Data[1][8][3] = 0.133;
		 Data[1][8][4] = 0.509; Data[1][8][5] = 0.977; 
Data[1][8][6] = 0.886; Data[1][8][7] = 0.60; 
		 
Data[2][8][1] = 26.76; Data[2][8][2] = 20.70;
		 Data[2][8][3] = 18.21; Data[2][8][4]  = 15.51; 
		 Data[2][8][5] = 6.71; Data[2][8][6] = 2.00; 
	     Split[8] = 0; MaxDays[8] = Weeks*7 + Extra;

fp = fopen("FIOut.dat","wt");
for(i=1;i<=8;i++)
	fprintf(fp,"%e\n",FractI[i]);
fclose(fp);


		 return MaxPlot;

}

int OtvosSpray(float *S0, float *FractI,float ***Data,int *Split, int *MaxDays, int *MaxWk, int *DataInterval,int *SData,int *FractDead)
{
		 int Extra = 1;
		 double Weeks = 8;
		 int Plot, MaxPlot = 4;
		 int Wk;
                 double InitFICon = -1; //1e-6;
		 *SData = 0;
		 *FractDead = 0;
		 
		 *MaxWk = 7;
		 *DataInterval = 7;

		 S0[1] = 194.34;  //T1 ("Plot 2")
	     FractI[1] = 0.093;
Data[1][1][1] = -1; Data[1][1][2] = -1; Data[1][1][3] = 0.015;
		 Data[1][1][4] = 0.65; Data[1][1][5] = 0.685; 

Data[1][1][6] = 0.88; Data[1][1][7] = 0.728; 
		 
Data[2][1][1] = 150.96; Data[2][1][2] = 185.73;
		 Data[2][1][3] = 103.39; Data[2][1][4]  = 82.31; 
		 Data[2][1][5] = 41.31; Data[2][1][6] = 6.72; 
	     Split[1] = 0; MaxDays[1] = Weeks*7 + Extra;

		 S0[2] = 145.76;  //T2 ("Plot 3")
	     FractI[2] = 0.09;
Data[1][2][1] = -1; Data[1][2][2] = -1; Data[1][2][3] = 0.059;
		 Data[1][2][4] = 0.647; Data[1][2][5] = 0.845; 
Data[1][2][6] = 0.769; Data[1][2][7] = 0.910; 
		 
Data[2][2][1] = 102.39; Data[2][2][2] = 133.61;
		 Data[2][2][3] = 99.68; Data[2][2][4]  = 70.93; 
		 Data[2][2][5] = 13.05; Data[2][2][6] = 2.77; 
	     Split[2] = 0; MaxDays[2] = Weeks*7 + Extra;

		 S0[3] = 302.02;  //T3 ("Plot 4")
	     FractI[3] = 0.147;
Data[1][3][1] = -1; Data[1][3][2] = -1; Data[1][3][3] = 0.308;
		 Data[1][3][4] = 0.875; Data[1][3][5] = 1.0; 
Data[1][3][6] = 0.967; Data[1][3][7] = 0.527; 
		 
Data[2][3][1] = 162.04; Data[2][3][2] = 128.96;
		 Data[2][3][3] = 130.90; Data[2][3][4]  = 96.86; 
		 Data[2][3][5] = 5.84; Data[2][3][6] = 0.96; 
	     Split[3] = 0; MaxDays[3] = Weeks*7 + Extra;

		 S0[4] = 41.77;  //T4 ("Plot 5")
	     FractI[4] = (0.279+0.223)/2.0;
Data[1][4][1] = -1; Data[1][4][2] = -1; Data[1][4][3] = 0.133;
		 Data[1][4][4] = 0.509; Data[1][4][5] = 0.977; 
Data[1][4][6] = 0.886; Data[1][4][7] = 0.60; 
		 
Data[2][4][1] = 26.76; Data[2][4][2] = 20.70;
		 Data[2][4][3] = 18.21; Data[2][4][4]  = 15.51; 
		 Data[2][4][5] = 6.71; Data[2][4][6] = 2.00; 
	     Split[4] = 0; MaxDays[4] = Weeks*7 + Extra;



		 return MaxPlot;

}


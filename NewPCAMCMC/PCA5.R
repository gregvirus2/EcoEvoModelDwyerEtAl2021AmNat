

#Data = read.table("AllGLSA8O41YA.dat");

#Data = read.table("AllGLSA8O1YC.dat");
#Data = read.table("AllGLSA8O42YC.dat");
#Data = read.table("AllGLSA8O0YC.dat");
#Data = read.table("AllGLSA8P0YC.dat");
#Data = read.table("AllGLSA8O5YC.dat");
#Data = read.table("AllGLSA8O42ZC.dat");
#Data = read.table("AllGLSA8O15ZC.dat");

#Data = read.table("AllGLSA8O15ZD.dat");

KeepFract = 0.5; #0.25 for GLSA8O1YC, 0.1 for 8O42YC, 0.65 for 8P0YC, 
#0.5 for 8O42ZC, 8O15ZC #0.5 for 8O15ZD

Index = which(Data[,9]!=0) #Huh??

Data = Data[Index,1:10]

ColNum = ncol(Data); #likelihoods/posteriors are in last TWO columns

yy = order(Data[,ColNum],decreasing=TRUE); #ordering the likelihoods

zz = cbind(Data)[yy,] #finishing off the likelihoods

mFixed = length(unique(Data[,7]));
deltaFixed = length(unique(Data[,6]));

ratioFixed = length(unique(Data[,5])); 

if(Data[1,1]>1e5){ #no k parameter
Data1 = (log(zz[,2])); #already on a log scale
Data2 = log(zz[,3]); #already on a log scale
Data3 = log(zz[,4]);
Data4 = log(zz[,5]);
Data5 = log(zz[,6]);
Data6 = log(zz[,7]);
Data7 = log(zz[,8]);


	 

	if(mFixed==1){ 

		if(deltaFixed==1){ #maybe wrong?  Maybe no Data1?
			if(ratioFixed==1){
				Data8 = cbind(Data1,Data2,Data3,Data7);
			}else{
				Data8 = cbind(Data1,Data2,Data3,Data4,Data7);
			}
		}else{

		     if(ratioFixed==1){
				Data8 = cbind(Data1,Data2,Data3,Data5,Data7);
			}else{				
				Data8 = cbind(Data1,Data2,Data3,Data4,Data5,Data7);
			}
		}
		

	} else{
		
		if(ratioFixed==1){
			Data8 = cbind(Data1,Data2,Data3,Data5,Data6,Data7);
		}else{
			Data8 = cbind(Data1,Data2,Data3,Data4,Data5,Data6,Data7);
		}
		
	}


}else{


Data1 = (log(zz[,1])); #already on a log scale
Data2 = log(zz[,2]); #already on a log scale
Data3 = log(zz[,3]);
Data4 = log(zz[,4]);
Data5 = log(zz[,5]);
Data6 = log(zz[,6]);
Data7 = log(zz[,7]);
Data8 = log(zz[,8]);


	if(mFixed==1){ 
	
		if(deltaFixed==1){
			cat("m fixed and delta fixed\n");
			Data8 = cbind(Data1,Data2,Data3,Data4,Data5,Data8);
		}else{
			Data8 = cbind(Data1,Data2,Data3,Data4,Data5,Data6,Data8);
		}

	} else{
		
		Data8 = cbind(Data1,Data2,Data3,Data4,Data5,Data6,Data7,Data8);

	}


}

#Data8 = Data8[1:nrow(zz),];

if(0){
Best = zz[1,10];
Add = -20;
BestAdd = Best + Add;

Index1 = which(zz[,10]<BestAdd);
Last = max(zz[Index1,10])
Index2 = which(zz[,10]==Last);
KeepFract = Index2/nrow(zz);
}




End = KeepFract*nrow(zz);  #CHANGE THIS LINE TO CHANGE HOW MANY PARAMETER SETS YOU KEEP
Data8b = Data8[1:End,];

par(mfrow=c(1,1))
par(ask="TRUE")
for(i in 1:ncol(Data8b))
	hist((Data8b[,i]),breaks=20);

x = prcomp(Data8b,scale=TRUE); #PCA starts here
x$rotation  
x$scale
x$center


#Now a piece of Dave's code, showing the correlations

Inverse = solve(x$rotation);
PC_Matrix = matrix(,nrow(Data8b),ncol(Data8b));
Data8B = as.matrix(Data8b);  #THIS LINE LOOKS POINTLESS
for(i in 1:nrow(Data8b)){

	PC_Matrix[i,] = Inverse %*% ((Data8B[i,]-(x$center))/x$scale);
	


}

par(ask=TRUE); 

plot(data.frame(PC_Matrix),cex=0.01)






#A is all, B is 0.5, C is 0.25, D is 0.1, E is 0.05, F is 0.4
#GLS8O17XB is KeepFract = 0.5
#GLS8O19ZA is KeepFract = 1.0
#GLS8O25XD has KeepFract = 0.75, because otherwise you get huge nubars that stop the code
#GLS8O27ZA has KeepFract = 1.0
#GLS8O29QC has KeepFract = 1.0
#GLS8O17QC has KeepFract = 0.5
#GLS8P23ZB has KeepFract = 1.0, also GLS8O23ZC and GLS8O29ZC
#OLO11F has KeepFract = 1.0
#GzLS8O29QC has KeepFract = 1.0
#OzLSO16B has KeepFract = 0.25 (probably messed up)
#OzLSO31F has KeepFract = 1.0 (TOTALLY messed up - delta and m SHOULD have had zero variances)
#OzLSO16ZB has KeepFract = 1.0
#OzLSP11B has KeepFract = 1.0
#OzLSP11H has KeepFract = 0.25
#OLSP23ZA has KeepFract = 1.0, and it's from OzLSP23ZE
#OLSP23Z2E has KeepFract = 0.05, also from OzLSP23ZE
#ORLSP11F has KeepFract = 1.0
#ORLSP33D has KeepFract = 1.0
#ORLSP23ZC has KeepFract = 1.0
#ORLSO23ZC "         "      "
#ORLSP11ZD has KeepFract = 1.0

#GLSA8O23ZA has KeepFract = 0.5;
#GLSA8P11ZA has KeepFract = 1.0;
#GLSA8P11ZB has KeepFract = 0.5;
#GLSA8P34ZA has KeepFract = 0.5;
#GLSA8P1ZA has KeepFract = 0.01;
#GLSA8P35XA has KeepFract = 0.25;
#GLSA8P1RA has KeepFract = 1.0;
write(x$rotation,file="GLSA8O15ZDRotations.txt",ncolumns=length(Data8B[1,]));
write(x$scale,file="GLSA8O15ZDScale.txt",ncolumns=length(Data8B[1,]));
write(x$center,file="GLSA8O15ZDCenter.txt",ncolumns=length(Data8B[1,]));
write(x$sd,file="GLSA8O15ZDsd.txt",ncolumns=length(Data8B[1,]));


pD = -2*(mean(zz[,9]) - max(zz[,9]))
Dbar = -2*(mean(zz[,9]));
AltpD = 2.0*var(zz[,9]);
DIC = -2*max(zz[,9]) + 2*pD
AltDIC = Dbar + AltpD
cat("max:",max(zz[,9]),"mean:",mean(zz[,9]),"pD:",pD,"DIC:",DIC,"AltpD:",AltpD,"AltDIC:",AltDIC,"\n");

par(ask=FALSE);

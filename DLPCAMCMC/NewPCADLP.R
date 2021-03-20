#Data = read.table("AllLSEG3B28O28YA.dat");

#Data = read.table("AllLSNEG48O28YA.dat");

#Data = read.table("AllLSE8O42YC.dat");
#Data = read.table("AllLSNE8O42YC.dat");

#Data = read.table("AllDLSE8O42ZC.dat");
#Data = read.table("AllD4DLSE8O42ZC.dat");
#Data = read.table("AllD7CLSE8O42ZC.dat");
#Data = read.table("AllD2CLSNE8O42ZC.dat");
#Data = read.table("AllD2CLSNE8P0ZC.dat");
#Data = read.table("AllD6CLSNE8P0ZC.dat");
#Data = read.table("AllD7CLSNE8O15ZC.dat");
#Data = read.table("AllE8CLSE8O42ZC.dat");
#Data = read.table("AllE7CLSE8O42ZC.dat");
#Data = read.table("AllG8CLSE8O42ZC.dat");
#Data = read.table("AllH8CLSE8O15ZC.dat");
#Data = read.table("AllH9CLSE8O15ZC.dat");
#Data = read.table("AllI7CLSNE8O42ZC.dat");
#Data = read.table("AllJ9CLSE8O42ZC.dat");
#Data = read.table("AllI9CLSE8O42ZC.dat");
#Data = read.table("AllJ7CLSNE8O42ZC.dat");

#Data = read.table("AllL8CLSE8O42ZC.dat");
#Data = read.table("AllL7CNE8O42ZC.dat");

#Data = read.table("AllL6CLSE8O42ZCFirst.dat");
#Data = read.table("AllL6CLSE8O42ZCThird.dat");
#Data = read.table("AllL6CLSE8O42ZCFourth.dat");
#Data = read.table("AllL6CLSE8O42ZC.dat");

#Data = read.table("AllM0CLSE8O42ZC.dat");
#Data = read.table("AllN0CLSE8O42ZC.dat");

#Data = read.table("AllN1CLSNE8O42ZC.dat");
#Data = read.table("AllN3CLSE8O42ZC.dat");
#Data = read.table("AllN5CLSE8O42ZC.dat"); #evoln with stoch, but stoch doesn't help 
#Data = read.table("AllN6CLSE8O15ZC.dat");
#Data = read.table("AllN7CLSE8O15ZC.dat");
#Data = read.table("AllN8CLSNE8O15ZC.dat");
#Data = read.table("AllN2e4CLSE8O42ZC.dat");

#Data = read.table("AllN2e4CLSE8O42ZC.dat");
#Data = read.table("AllN2e6CLSE8O42ZC.dat");
#Data = read.table("AllN2e55CLSE8O42ZC.dat");
#Data = read.table("AllN2e56CLSE8O42ZC.dat");
#Data = read.table("AllN2e7CLSE8O42ZC.dat");
#Data = read.table("AllN8e4CLSE8O42ZC.dat");
#Data = read.table("AllN1e5CLSNE8O42ZC.dat");
#Data = read.table("AllN1e6CLSNE8O42ZC.dat");
#Data = read.table("AllN2e5CLSE8O42ZC.dat");
#Data = read.table("AllN4e5CLSE8O42ZC.dat");
#Data = read.table("AllN4e6CLSE8O42ZC.dat");
#Data = read.table("AllN4e56CLSE8O42ZC.dat");
#Data = read.table("AllN1e3CLSNE8O42ZC.dat");
#Data = read.table("AllN4CLSE8O42ZC.dat");
#Data = read.table("AllN4e5CLSE8O42ZC.dat");
#Data = read.table("AllN2BCLSE8O42ZC.dat");

#Data = read.table("AllN7CLSE8O42ZC.dat");
Data = read.table("AllN6BCLSE8O15ZC.dat");




Evoln = 1;
if(Evoln==1){
	 if(Data[1,5]!=0){
	    Num = 4;
	 }else{ 
	    Num = 3;
	 }
}else{
	Num = 4;
}

KeepFract = 1.0; #8O42YC uses 0.25, #NE8O42YC uses 1.0, E8O42ZC uses 0.5(3), NE8O42ZC uses 1.0
#D3CLSNE8O42ZC could use 0.75, same for D4CLSE8O42ZC
#D6CLSNE8O42ZC uses 0.25; D2CLSNE8O42ZC uses 0.5
#D6CLSNE8P0ZC uses 0.15 D7CLNE8O15ZC uses 0.5; D8CE8O42ZC uses 0.1
#G7CLSE8O42ZC uses 0.75; G8CLSE8O42ZC uses 0.25
#H8CLSE8O15ZC uses 0.5; H9CLSE8O15ZC uses 1.0
#I9CLSE8O42ZC uses 0.5; J9 and J7 use 0.5
#J8CLSE8O42ZC uses 0.5; 
#L8CLSE8O42ZC uses 0.5; L7CNE8O42ZC uses 0.5 (note notation error!)
#L8CLSE8O42ZC2 uses 1.0, because 0.5 was a disaster
#L9CLSE8O42ZC uses 0.75, L6CLSE8O42ZC uses 0.75
#L0CLSE8042ZC is L6 but uses 1.0; L1CLSE8042ZC uses 0.25
#N0CLSE8O42ZC uses 0.5; N1CLSNE8O42ZC uses 1.0;
#N2CLSE8O42ZC is N0, but uses 1.0
#N4, N6, N7 use 1.0


ColNum = ncol(Data); #likelihoods/posteriors are in last TWO columns
yy = order(Data[,ColNum],decreasing=TRUE); #ordering the posteriors
zz = cbind(Data)[yy,] #finishing off the ordering

if(Evoln==1){
  Data1 = log(zz[,1]); 
  Data3 = log(zz[,3]);
  Data5 = log(zz[,5])
  Data6 = log(zz[,6]);
}else{
  Data1 = log(zz[,1]); 
  Data2 = log(zz[,2]);
  Data3 = log(zz[,3]);
  Data4 = log(zz[,4]);
}

#Data3 = log(zz[,3]);
#Data4 = log(zz[,4]);
#else{
 # Data3 = log(zz[,3])
#}
if(Evoln<1){
	#TotData = cbind(Data1,Data3,Data4,Data5,Data6);
	#TotData = cbind(Data1,Data2,Data3,Data4);
	#TotData = cbind(Data1,Data3,Data6);
	TotData = cbind(Data1,Data2,Data3,Data4);
}else{

  if(Data[1,5]!=0){
    TotData = cbind(Data1,Data3,Data5,Data6);
  }else{
    TotData = cbind(Data1,Data3,Data6);
  }
  

}

End = KeepFract*nrow(zz);  #CHANGE THIS LINE TO CHANGE HOW MANY PARAMETER SETS YOU KEEP
TotDatab = TotData[1:End,];

par(ask=FALSE);

#par(mai=c(0.5,0.5,0.5,0.5));

#par(mfrow=c(Num,1)); 
par(mfrow=c(1,1));
par(ask="TRUE");
for(i in 1:Num){
  
  #hist(log(Data[,i]),breaks=50);
  hist((TotDatab[,i]),breaks=100);
  
}

x = prcomp(TotDatab,scale=TRUE); #PCA starts here
x$rotation  
x$scale
x$center


#Now a piece of Dave's code, showing the correlations

Inverse = solve(x$rotation);
PC_Matrix = matrix(,nrow(TotDatab),ncol(TotDatab));
Data8B = as.matrix(TotDatab);  #THIS LINE LOOKS POINTLESS, but it is used below
for(i in 1:nrow(TotDatab)){

	PC_Matrix[i,] = Inverse %*% ((TotDatab[i,]-(x$center))/x$scale);

}

par(ask=TRUE); 

plot(data.frame(PC_Matrix),cex=0.01)


write(x$rotation,file="N6BCLSE8O15ZCRotations.txt",ncolumns=length(Data8B[1,]));
write(x$scale,file="N6BCLSE8O15ZCScale.txt",ncolumns=length(Data8B[1,]));
write(x$center,file="N6BCLSE8O15ZCCenter.txt",ncolumns=length(Data8B[1,]));
write(x$sd,file="N6BCLSE8O15ZCsd.txt",ncolumns=length(Data8B[1,]))




#DataTemp = read.table("AllOR4P11ZD.dat");


#DataTemp = read.table("AllMEG3B28O28YA.dat");
#DataTemp = read.table("AllMEG3B28O28YA.dat");
#DataTemp = read.table("AllMNEG3A8O28YA.dat");

#DataTemp = read.table("AlltestD7CE8O42ZC.dat");
#DataTemp = read.table("AllD6CNE8P0ZC.dat");
#DataTemp = read.table("AlltestD7CE8O42ZC.dat");
#DataTemp = read.table("AllD6CNE8O42ZC.dat");
#DataTemp = read.table("AllD7CNE8O15ZC.dat");
#DataTemp = read.table("AllD7CE8O15ZC.dat");
#DataTemp = read.table("AllD8CE8O42ZC.dat");
#DataTemp = read.table("AllH8CE8O42ZC.dat");
#DataTemp = read.table("AllH9CE8O15ZC.dat");
#DataTemp = read.table("AllH7CE8O42ZC.dat");
#DataTemp = read.table("AllI9BCE8O42ZC.dat");
#DataTemp = read.table("AllI7CNE8O42ZC.dat");
#DataTemp = read.table("AllJ9CE8O42ZC.dat");
#DataTemp = read.table("AllJ7CNE8O42ZC.dat");

#DataTemp = read.table("AllL8CE8O42ZC.dat");
#DataTemp = read.table("AllL7CNE8O42ZC.dat");
#DataTemp = read.table("AllL9CE8O42ZC.dat");
#DataTemp = read.table("AllL6CE8O42ZC.dat");
#DataTemp = read.table("AllL0CE8O42ZC.dat");
#DataTemp = read.table("AllL1CE8O42ZC.dat");
#DataTemp = read.table("AllM0CE8O42ZC.dat");

#DataTemp = read.table("AllN1CNE8O42ZC.dat");
#DataTemp = read.table("AllN2CE8O42ZC.dat");

#DataTemp = read.table("AllN0CE8O42ZC.dat");
#DataTemp = read.table("AllN3CE8O42ZC.dat");
#DataTemp = read.table("AllN4CE8O42ZC.dat");
#DataTemp = read.table("AllN6CE8O15ZC.dat");
#DataTemp = read.table("AllN7CE8O15ZC.dat");
#DataTemp = read.table("AllN8CNE8O15ZC.dat");
#DataTemp = read.table("AllN2e4CE8O42ZC.dat");
#DataTemp = read.table("AllN2e6CE8O42ZC.dat");
#DataTemp = read.table("AllN2e56CE8O42ZC.dat");
#DataTemp = read.table("AllN2e7CE8O42ZC.dat");
#DataTemp = read.table("AllN2e55CE8O42ZC.dat");
#DataTemp = read.table("AllN2e4CE8O42ZC.dat");
#DataTemp = read.table("AllN2e56CE8O42ZC.dat");
#DataTemp = read.table("AllN1CNE8O42ZC.dat");
#DataTemp = read.table("AllN8e4CNE8O42ZC.dat");
#DataTemp = read.table("AllN1e6CNE8O42ZC.dat");
#DataTemp = read.table("AllN1e5CNE8O42ZC.dat");
#DataTemp = read.table("AllN2e5CE8O42ZC.dat");
#DataTemp = read.table("AllN2e7CE8O42ZC.dat");
#DataTemp = read.table("AllN1e6CNE8O42ZC.dat");
#DataTemp = read.table("AllN8e4CNE8O42ZC.dat");
#DataTemp = read.table("AllN4e5CE8O42ZC.dat");
#DataTemp = read.table("AllN1e3CNE8O42ZC.dat");
#DataTemp = read.table("AllN4e5CE8O42ZC.dat");
#DataTemp = read.table("AllN1e3CNE8O42ZC.dat");
#DataTemp = read.table("AllN4e5CE8O42ZC.dat");
#DataTemp = read.table("AllN2BCE8O42ZC.dat");
#DataTemp = read.table("AllN8CE8O15ZC.dat");
DataTemp = read.table("AllN6BCE8O15ZC.dat");



Evoln = 1;

if(Evoln==1){
  #DataAll = cbind(log(DataTemp[,1]),log(DataTemp[,2]),log(DataTemp[,3]),log(DataTemp[,4]),log(DataTemp[,6]));

  #DataAll = cbind(log(DataTemp[,2]),log(DataTemp[,3]),log(DataTemp[,4]),log(DataTemp[,5]),log(DataTemp[,7]));
  if(DataTemp[1,6]!=0){
    
    DataAll = cbind(log(DataTemp[,2]),log(DataTemp[,4]),log(DataTemp[,6]),log(DataTemp[,7]));
  }else{
    DataAll = cbind(log(DataTemp[,2]),log(DataTemp[,4]),log(DataTemp[,7]));
  }
  
  
  #DataAll = cbind(log(DataTemp[,2]),log(DataTemp[,4]),log(DataTemp[,7]));
  
  
  }else{
  DataAll = cbind(log(DataTemp[,2]),log(DataTemp[,3]),log(DataTemp[,4]),log(DataTemp[,5]));
}

require(coda);

NumFiles = 1;
xLast = numeric();
Diff = numeric();
String = 1;
ColLength = ncol(DataTemp);
for(i in 1:(nrow(DataTemp)-1)){
#for(i in 1:1e2){
	if(DataTemp[i,1] > DataTemp[i+1,1]){
		cat("String:",String," i:",i," i+1:",i+1,"\n");
		xLast[String] = i;
		String = String + 1;
		NumFiles = NumFiles + 1;
		cat("i:",i,"NumFiles:",NumFiles,"xLast:",xLast[String-1],"\n");
	}
	if(0){
	for(j in 1:(ColLength-3)){
		if(DataTemp[i,j]<=1e-5){
			DataTemp[i,j] = 1e-5;
		}
		if(is.na(DataTemp[i,j])) cat("found an na...\n");
		if(is.infinite(DataTemp[i,j])) cat("found an inf...\n");
		if(is.nan(DataTemp[i,j])) cat("found an inf...\n");


	}
	}
	
}

xLast[NumFiles] = nrow(DataTemp);
for(i in 1:(NumFiles-1)){
	Diff[i] = xLast[i+1]-xLast[i];
	cat("i:",i," xLast:",xLast[i]," xLast:",xLast[i+1]," Diff:",Diff,"\n");
}

Last = min(Diff);
cat("NumFiles:",NumFiles,"\n");

#xLast = 1e6*c(1:10);


thinVal = 1;



x1 = mcmc(DataAll[1:Last,],thin=thinVal);
Index = which(x1==-Inf);
x1[Index] = -13;
#x1b = DataAll[seq(1,Last,thinVal),];


x1b = DataAll[seq(1,Last,thinVal),];
temp1 = seq(1,Last,thinVal);
cat("temp1:",length(temp1),"\n");
dummy = xLast[1]+1;
x2 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
Index = which(x2==-Inf);
x2[Index] = -13;
temp2 = seq(dummy,dummy+Last-1,thinVal);
cat("temp2:",length(temp2),"\n");

x2b = DataAll[seq(dummy,(dummy+Last-1),thinVal),];
dummy = xLast[2]+1;

if(NumFiles>2){
x3 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
Index = which(x3==-Inf);
x3[Index] = -13;
temp3 = seq(dummy,dummy+Last-1,thinVal);
cat("temp3:",length(temp3),"\n");
x3b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
dummy = xLast[3]+1;
}
if(NumFiles>3){
	cat("just before mcmc...\n");
	x4 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x4==-Inf);
	x4[Index] = -13;
	x4b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[4]+1;
}
if(NumFiles>4){
	x5 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x5==-Inf);
	x5[Index] = -13;
	x5b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[5]+1;
}
if(NumFiles>5){
	x6 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x6==-Inf);
	x6[Index] = -13;
	x6b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[6]+1;
}
if(NumFiles>6){
	x7 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x7==-Inf);
	x7[Index] = -13;
	x7b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[7]+1;
}
if(NumFiles>7){
	x8 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x8==-Inf);
	x8[Index] = -13;
	x8b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[8]+1;
}
if(NumFiles>8){
	x9 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x9==-Inf);
	x9[Index] = -13;
	x9b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[9]+1;
}
if(NumFiles>9){
	x10 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	Index = which(x10==-Inf);
	x10[Index] = -13;
	x10b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[10]+1;
}

if(NumFiles>10){
	x11 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x11b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[11]+1;
}

if(NumFiles>11){
	x12 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x12b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[12]+1;
}

if(NumFiles>12){
	x13 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x13b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[13]+1;
}

if(NumFiles>13){
	x14 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x14b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[14]+1;
}

if(NumFiles>14){
	x15 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x15b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[15]+1;
}

if(NumFiles>15){
	x16 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x16b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[16]+1;
}

if(NumFiles>16){
	x17 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x17b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[17]+1;
}

if(NumFiles>17){
	x18 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x18b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[18]+1;
}

if(NumFiles>18){
	#cat("19 Files...\n");
	x19 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x19b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[19]+1;
}

if(NumFiles>19){
	x20 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x20b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[20]+1;
}

if(NumFiles>20){
	x21 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x21b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[21]+1;
}

if(NumFiles>21){
	x22 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x22b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[22]+1;
}


if(NumFiles>22){
	x23 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x23b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[23]+1;
}

if(NumFiles>23){
	x24 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x24b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[24]+1;
}

if(NumFiles>24){
	x25 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	x25b = DataAll[seq(dummy,dummy+Last-1,thinVal),];
	dummy = xLast[25]+1;
}



if(NumFiles==2) xAll = mcmc.list(x1,x2);
if(NumFiles==3) xAll = mcmc.list(x1,x2,x3);
if(NumFiles==4) xAll = mcmc.list(x1,x2,x3,x4);
if(NumFiles==5) xAll = mcmc.list(x1,x2,x3,x4,x5);
if(NumFiles==6) xAll = mcmc.list(x1,x2,x3,x4,x5,x6);
if(NumFiles==7) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7);
if(NumFiles==8) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8);
if(NumFiles==9) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9);
if(NumFiles==10) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);
if(NumFiles==11) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11);
if(NumFiles==12) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12);
if(NumFiles==13) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13);
if(NumFiles==16) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16);

if(NumFiles==19) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19);


if(NumFiles==20){
	 xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20);
	 #xAll = mcmc.list(x3,x4,x5);
}


if(NumFiles==21) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21);



if(NumFiles==25) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25);

#xAll = mcmc.list(x1,x2,x4,x5,x6,x7,x8,x9,x10); #x3 is messed up for 5O1C, x10 for 5P2A, x10 for 5O18XA)

#xAll = mcmc.list(x2,x3,x4,x5,x6,x8,x9,x10); #x6 and x7 messed up in OAO1C

#xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x9,x10);

#gelman.diag(xAll);

if(0){ #why is this next bit here?
par(ask="TRUE");
x1c = mcmc(x1b);
x2c = mcmc(x2b);
if(NumFiles>2) x3c = mcmc(x3b);
if(NumFiles>3) x4c = mcmc(x4b);
if(NumFiles>4) x5c = mcmc(x5b);
if(NumFiles>5) x6c = mcmc(x6b);
if(NumFiles>6) x7c = mcmc(x7b);
if(NumFiles==3) xAllc = mcmc.list(x1c,x2c,x3c);
if(NumFiles==4) xAllc = mcmc.list(x1c,x2c,x3c,x4c);
if(NumFiles==5) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c);
if(NumFiles==6) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c);
if(NumFiles==7) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c);
}
print(gelman.diag(xAll));
diagOut = gelman.diag(xAll);
par(mai=c(0.25,0.25,0.25,0.25));
par(ask="TRUE");
plot(xAll);
summary(xAll)

if(0){
par(mai=c(0.5,0.5,0.5,0.5));
par(mfrow=c(1,1))
hist(DataTemp[,8])
}



#test = H9CE8O15ZC$statistics[,1:2]
#N8e5CNE8O15ZC = summary(xAll);
#test = N4e5CE8O42ZC$statistics[,1:2]
N4e5CE8O42ZC = summary(xAll);




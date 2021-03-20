
#Last - last possible data point over all chains, based on shortest chain
StartFract = 0.5; #What fraction to throw out, so 0.1 means throw out first 10%
thinStart = round(StartFract*Last);


Target = 10; #target value for how LITTTLE to thin by, sort of. Last/Target is how much to jump by (small Target, big jump)
thinVal = round(StartFract*Last/Target); #how much to jump by when thinning


x1 = mcmc(DataAll[1:Last,],thin=thinVal);
#x1b = DataAll[seq(1,Last,thinVal),];

cat("x1 thinStart:",thinStart,"Last:",Last,"Diff:",Last-thinStart,"\n");

######### x1b etc have the thinned versions
temp1 = seq(thinStart,Last,thinVal);
temp1 = temp1[1:(length(temp1))]
#x1b = DataAll[seq(from=thinStart,to=Last,by=thinVal),];
x1b = DataAll[temp1,];

cat("temp1:",length(temp1),"\n");

dummy = xLast[1]+1;
x2 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
temp2 = seq(dummy+thinStart,dummy+Last-1,thinVal);
cat("temp2:",length(temp2),"\n");
cat("x2 Start:",dummy+thinStart,"dummy+Last-1:",dummy+Last-1,"Diff:",Last-1-thinStart,"\n");
x2b = DataAll[seq(dummy+thinStart,(dummy+Last-1),thinVal),];
dummy = xLast[2]+1;
if(NumFiles>2){
x3 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
temp3 = seq(dummy+thinStart,dummy+Last-1,thinVal);
cat("temp3:",length(temp2),"\n");
x3b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
dummy = xLast[3]+1;
}
if(NumFiles>3){
	x4 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp4 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp4:",length(temp4),"\n");

	x4b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[4]+1;
}
if(NumFiles>4){
	x5 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp5 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp5:",length(temp5),"\n");
	x5b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[5]+1;
}
if(NumFiles>5){
	x6 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp6 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp6:",length(temp6),"\n");
	x6b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[6]+1;
}
if(NumFiles>6){
	x7 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp7 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp7:",length(temp7),"\n");
	x7b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[7]+1;
}
if(NumFiles>7){
	x8 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp8 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp8:",length(temp8),"\n");
	x8b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[8]+1;
}
if(NumFiles>8){
	x9 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp9 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp9:",length(temp9),"\n");
	x9b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[9]+1;
}
if(NumFiles>9){
	x10 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp10 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp10:",length(temp10),"\n");
	x10b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[10]+1;
}

if(NumFiles>10){
	x11 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp11 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp11:",length(temp11),"\n");
	x11b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[11]+1;
}

if(NumFiles>11){
	x12 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp12 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp12:",length(temp12),"\n");
	x12b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[12]+1;
}

if(NumFiles>12){
	x13 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp13 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp13:",length(temp13),"\n");
	x13b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[13]+1;
}


if(NumFiles>13){
	x14 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp14 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp14:",length(temp14),"\n");
	x14b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[14]+1;
}

if(NumFiles>14){
	x15 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp15 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp15:",length(temp15),"\n");
	x15b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[15]+1;
}


if(NumFiles>15){
	x16 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp16 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp16:",length(temp16),"\n");
	x16b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[16]+1;
}

if(NumFiles>16){
	x17 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp17 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp17:",length(temp16),"\n");
	x17b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[17]+1;
}

if(NumFiles>17){
	x18 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp18 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp15:",length(temp18),"\n");
	x18b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[18]+1;
}

if(NumFiles>18){
	x19 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp19 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp19:",length(temp19),"\n");
	x19b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[19]+1;
}


if(NumFiles>19){
	x20 = mcmc(DataAll[(dummy):(dummy+Last-1),],thin=thinVal);
	temp20 = seq(dummy+thinStart,dummy+Last-1,thinVal);
	cat("temp20:",length(temp20),"\n");
	x20b = DataAll[seq(dummy+thinStart,dummy+Last-1,thinVal),];
	dummy = xLast[20]+1;
}

cat("Length x1:",length(x1b),"\n");

if(NumFiles==2) xAll = mcmc.list(x1,x2);
if(NumFiles==3) xAll = mcmc.list(x1,x2,x3);
if(NumFiles==4) xAll = mcmc.list(x1,x2,x3,x4);
if(NumFiles==5) xAll = mcmc.list(x1,x2,x3,x4,x5);
if(NumFiles==6) xAll = mcmc.list(x1,x2,x3,x4,x5,x6);
if(NumFiles==7) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7);
if(NumFiles==8) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8);
if(NumFiles==9) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9);
if(NumFiles==10) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);
if(NumFiles==10) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10);

xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19);

if(NumFiles==19) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19);
if(NumFiles==20) xAll = mcmc.list(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20);


if(0){ #this chunk of code seems to have some problem...
par(ask="TRUE");
x1c = mcmc(x1b);
x2c = mcmc(x2b);
if(NumFiles>2) x3c = mcmc(x3b);
if(NumFiles>3) x4c = mcmc(x4b);
if(NumFiles>4) x5c = mcmc(x5b);
if(NumFiles>5) x6c = mcmc(x6b);
if(NumFiles>6) x7c = mcmc(x7b);
if(NumFiles>7) x8c = mcmc(x8b);
if(NumFiles>8) x9c = mcmc(x9b);
if(NumFiles>9) x10c = mcmc(x10b);
if(NumFiles>10) x11c = mcmc(x11b);
if(NumFiles>11) x12c = mcmc(x12b);
if(NumFiles>12) x13c = mcmc(x13b);
if(NumFiles>13) x14c = mcmc(x14b);
if(NumFiles>14) x15c = mcmc(x15b);
if(NumFiles>15) x16c = mcmc(x16b);
if(NumFiles>16) x17c = mcmc(x17b);
if(NumFiles>17) x18c = mcmc(x18b);
if(NumFiles>18) x19c = mcmc(x19b);
if(NumFiles>19) x20c = mcmc(x20b);

if(NumFiles==2) xAllc = mcmc.list(x1c,x2c);
if(NumFiles==3) xAllc = mcmc.list(x1c,x2c,x3c);
if(NumFiles==4) xAllc = mcmc.list(x1c,x2c,x3c,x4c);
if(NumFiles==5) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c);
if(NumFiles==6) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c);
if(NumFiles==7) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c);
if(NumFiles==8) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c,x8c);
if(NumFiles==9) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c,x8c,x9c);
if(NumFiles==10) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c,x8c,x9c,x10c);

if(NumFiles==19) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c,x8c,x9c,x11c,x12,x13c,x14c,x15c,x16c,x17c,x18c,x19c);
if(NumFiles==20) xAllc = mcmc.list(x1c,x2c,x3c,x4c,x5c,x6c,x7c,x8c,x9c,x11c,x12c,x13c,x14c,x15c,x16c,x17c,x18c,x19c,x20c);


par(mai=c(0.25,0.25,0.25,0.25));
plot(xAllc);
print(gelman.diag(xAllc));

}

if(0){  #ACFs are screwing up for ME8O42YC, dunno why


if(Evoln==1){
	Params = 5;
}else{	
	Params = 4;
}
par(mfrow=c(1,Params));
par(ask="TRUE");
for(i in 1:Params){
	acf(x1c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}

for(i in 1:Params){
	acf(x2c[,i]);
  #cat("i=",i,"\n");
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x3c[,i]);
	abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x4c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x5c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}

for(i in 1:Params){
	acf(x6c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x7c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x8c[,i]);
	abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x9c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x10c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}

for(i in 1:Params){
	acf(x11c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x12c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x13c[,i]);
	abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x14c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x15c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}

for(i in 1:Params){
	acf(x16c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x17c[,i]);
      abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x18c[,i]);
	abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x19c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}
for(i in 1:Params){
	acf(x20c[,i]);
       abline(h=0.1);
	abline(h=0.2,col="RED");
	abline(h=-0.1);
	abline(h=-0.2,col="RED");
}

}


Index = which(DataAll[,1]<=0);
print(length(Index)/length(DataAll[,1]));


#need to allow for variable number of files
#rm(temp);

if(NumFiles==10)
	temp = rbind(x1b,x2b,x3b,x4b,x5b,x6b,x7b,x8b,x9b,x10b);

if(NumFiles==14)
	temp = rbind(x1b,x2b,x3b,x4b,x5b,x6b,x7b,x8b,x9b,x10b,x11b,x12b,x13b,x14b);

if(NumFiles==19)
  temp = rbind(x1b,x2b,x3b,x4b,x5b,x6b,x7b,x8b,x9b,x10b,x11b,x12b,x13b,x14b,x15b,x16b,x17b,x18b,x19b);

if(NumFiles==20)
	temp = rbind(x1b,x2b,x3b,x4b,x5b,x6b,x7b,x8b,x9b,x10b,x11b,x12b,x13b,x14b,x15b,x16b,x17b,x18b,x19b,x20b);

if(Evoln==1){
	NumberColumns = 3;
	#xOut = cbind(exp(temp[,1]),exp(temp[,2]),exp(temp[,3]),exp(temp[,4]),exp(temp[,5])); 
	
	#xOut = cbind(exp(temp[,1]),exp(temp[,2]),exp(temp[,3]),exp(temp[,4]),exp(temp[,5])); 
	xOut = cbind(exp(temp[,1]),exp(temp[,2]),exp(temp[,3])); 
	
	#write(x=t(xOut),file="N6BCE8O15ZCtmpC.in",sep="\t",ncolumns=NumberColumns);#You MUST transpose xOut! 
	write(x=t(xOut),file="N6BCE8O15ZC.in",sep="\t",ncolumns=NumberColumns);#You MUST transpose xOut! 
	
}else{ #No evoln
  NumberColumns = 4;
  
	xOut = cbind(exp(temp[,1]),exp(temp[,2]),exp(temp[,3]),exp(temp[,4])); 
	 
		#write(x=t(xOut),file="MNEG48O28YA.in",sep="\t",ncolumns=NumberColumns);#You MUST transpose xOut! 
	write(x=t(xOut),file="N1e3CNE8O42ZC.in",sep="\t",ncolumns=NumberColumns);#You MUST transpose xOut! 
}

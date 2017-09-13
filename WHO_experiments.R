
install.packages("energy")
install.packages("acepack")
install.packages("rPython")
install.packages("hexbin")
install.packages("ggplot2") 

library("hexbin")
library("ggplot2")
library("energy")
library("acepack")
library("rPython") 

python.load("HC_estimator.py")
python.exec("import numpy as np")

##############################  Hypercontractivity ############################## 

myHC <- function(x,y,bw){
  python.assign('x',x)
  python.assign('y',y)
  python.assign('bandwidth',bw)
  python.exec("res = HC(x,y,bandwidth)")
  val <- python.get("res")
  return(val)
}


##############################  MIC   ############################## 
# download mine.jar from http://www.exploredata.net/Downloads/MINE-Application

mymine <- function(x,y){
  xx=cbind(x,y)
  write("x,y",file="mytest.csv")
  write(t(xx),sep=",",file="mytest.csv",ncol=2,append=T)
  command <- 'java -jar MINE.jar "mytest.csv" -allPairs'
  system(command)
  res=scan("mytest.csv,allpairs,cv=0.0,B=n^0.6,Results.csv",what="",sep=",")
  val=as.numeric(res[11])
  return(val)
}


##############################  Scatterplots   ############################## 

myWHO_country_scatterplots <- function(var1,var2,bw){
	
	val.cor=val.dcor=val.mine=val.maxCor=val.HC=val.HCR=val.var1=val.var2=rep(NA,1) 
	set.seed(1)
	
	mydata = read.csv("WHO.csv")
	
	vx <- mydata[,var1]
	vy <- mydata[,var2]
	
	x <- vx[!is.na(vx+vy)]
	y <- vy[!is.na(vx+vy)]
	
	# Compute correlation estimates
	ind<- 1
	val.cor[ind]=(cor(x,y))^2            # Calculate the correlation
	val.dcor[ind]=dcor(x,y)              # Calculate dcor
	val.mine[ind]=mymine(x,y)            # Calculate mic
	a=ace(x,y)
	val.maxCor[ind]=cor(a$tx,a$ty)^2    # Calculate maxCorr
	val.HC[ind]=myHC(x,y,bw)
	val.HCR[ind]=myHC(y,x,bw)
	
	# Plot
	plot.new()
	title<-paste('WHO_figures/bw',bw,'_var1_',toString(var1),'_var2_',toString(var2),'.pdf',sep="")
	pdf(title)
	
	figtitle<-paste('HC:',round(val.HC, digits = 2), ' Cor:', round(val.cor, digits = 2),' dCor:',round(val.dcor, digits = 2),' mCor:',round(val.maxCor, digits = 2),' MIC:',round(val.mine, digits = 2), ' HCR:',round(val.HCR, digits = 2))
	plot(x,y,type="p", main=figtitle, xlab = colnames(mydata)[var1], ylab = colnames(mydata)[var2],cex.lab=1.5, cex.axis=1.50, cex.main=1.3, cex.sub=1.50)     
	dev.off()
	
	country <- mydata[,1]
	country <- country[!is.na(vx+vy)]
	
	## Remove 1 country
	lx <- length(x)
	x1 = x[c(1:(which.max(x)-1),(which.max(x)+1):lx)]
	y1 = y[c(1:(which.max(x)-1),(which.max(x)+1):lx)]
	
	ct1 <- country[which.max(x)]
	
	country = country[c(1:(which.max(x)-1),(which.max(x)+1):lx)]
	x = x1
	y = y1
	
	# Compute correlation estimates
	ind<- 1
	val.cor[ind]=(cor(x,y))^2            # Calculate the correlation
	val.dcor[ind]=dcor(x,y)              # Calculate dcor
	val.mine[ind]=mymine(x,y)            # Calculate mic
	a=ace(x,y)
	val.maxCor[ind]=cor(a$tx,a$ty)^2    # Calculate maxCorr
	val.HC[ind]=myHC(x,y,bw)
	val.HCR[ind]=myHC(y,x,bw)
	
	# Plot
	title<-paste('WHO_figures/bw',bw,'_var1_',toString(var1),'_var2_',toString(var2),'_without1country.pdf',sep="")
	plot.new()
	pdf(title)
	
	figtitle<-paste('HC:',round(val.HC, digits = 2), ' Cor:', round(val.cor, digits = 2),' dCor:',round(val.dcor, digits = 2),' mCor:',round(val.maxCor, digits = 2),' MIC:',round(val.mine, digits = 2), ' HCR:',round(val.HCR, digits = 2))
	plot(x,y,type="p", main=figtitle, xlab = colnames(mydata)[var1], ylab = colnames(mydata)[var2],cex.lab=1.5, cex.axis=1.50, cex.main=1.3, cex.sub=1.50)     
	dev.off()
	
	
	## Remove one more country
	lx1 <- length(x1)
	x2 = x1[c(1:(which.max(x1)-1),(which.max(x1)+1):lx1)]
	y2 = y1[c(1:(which.max(x1)-1),(which.max(x1)+1):lx1)]
	
	x = x2
	y = y2
	
	ct2 <- country[which.max(x1)]
	
	# Compute correlation estimates
	ind<- 1
	val.cor[ind]=(cor(x,y))^2            # Calculate the correlation
	val.dcor[ind]=dcor(x,y)              # Calculate dcor
	val.mine[ind]=mymine(x,y)            # Calculate mic
	a=ace(x,y)
	val.maxCor[ind]=cor(a$tx,a$ty)^2    # Calculate maxCorr
	val.HC[ind]=myHC(x,y,bw)
	val.HCR[ind]=myHC(y,x,bw)
	
	
	# Plot
	title<-paste('WHO_figures/bw',bw,'_var1_',toString(var1),'_var2_',toString(var2),'_without2countries.pdf',sep="")
	plot.new()
	pdf(title)
	
	figtitle<-paste('HC:',round(val.HC, digits = 2),' Cor:', round(val.cor, digits = 2),' dCor:',round(val.dcor, digits = 2),' mCor:',round(val.maxCor, digits = 2),' MIC:',round(val.mine, digits = 2), ' HCR:',round(val.HCR, digits = 2))
	plot(x,y,type="p", main=figtitle, xlab = colnames(mydata)[var1], ylab = colnames(mydata)[var2],cex.lab=1.5, cex.axis=1.50, cex.main=1.3, cex.sub=1.50)     
	dev.off()
	
	return(paste(ct1,ct2))
}

###################################### MAIN 
# Compute Hypercontractivity and other correlation measures for 1600 pairs  
# To run estimators, set Run_Estimator = True (takes ~ 40 min to run)
# Otherwise, saved estimates WHO_cor_1600.RData will be loaded.

Run_Estimator = FALSE

if(Run_Estimator){
	run_num = 1
	set.seed(1)
	mydata = read.csv("WHO.csv")
	ind=0
	
	num1 = 40
	num2 = 40
	
	total=num1*num2
	val.cor=val.dcor=val.mine=val.maxCor=val.HC=val.HCR=val.var1=val.var2=rep(NA,total) 
	
	for(var1 in 170+2*c(1:num1)){      #(172:2:250, 173:2:251)
		for(var2 in 171+2*c(1:num2)){
			vx <- mydata[,var1]
			vy <- mydata[,var2]
			x <- vx[!is.na(vx+vy)]
			y <- vy[!is.na(vx+vy)]
			ind <- ind+1
	   		val.cor[ind]=(cor(x,y))^2            # Calculate the correlation
	      	val.dcor[ind]=dcor(x,y)              # Calculate dcor
	      	val.mine[ind]=mymine(x,y)            # Calculate mic
	     	 a=ace(x,y)
	      	 val.maxCor[ind]=cor(a$tx,a$ty)^2    # Calculate maxCorr
	     	 val.HC[ind]=myHC(x,y,1.06)
	      	 val.HCR[ind]=myHC(y,x,1.06)
			 val.var1[ind]=var1
			 val.var2[ind]=var2		    	  
	  	    
		}
	}
	
	save(val.cor, val.dcor,val.mine, val.maxCor,val.HC,val.HCR,val.var1,val.var2, file=paste("WHO_vals_run",run_num,".RData",sep="")) 								
	load(paste("WHO_vals_run",run_num,".RData",sep="")) 
}else{
	mydata = read.csv("WHO.csv")
	load("WHO_1600.RData") # saved HC and correlation coefficients for 1600 pairs	
}


#########  Scatter plots 
# Fig 9-A (Cor vs HC)

plot.new()
par(mfrow = c(1,1), cex = 0.45)
df = data.frame(val.HC,val.cor,z=c(0.05,0.95,0.8,0.61),w=c(0.01,0.83,0.13,0.1),v=c(0.13,0.1,0.1,0.1))
df2 = data.frame(val.HC,val.cor,z=c(0.91,0.99,0.99,0.94),w=c(0.42,0.95,0.05,0.05),v=c(0.13,-0.1,0.1,0.1))

ggplot(df,aes(x=val.HC,y=val.cor))   + geom_point(color="blue") +xlab("Hypercontractivity")+ylab("Pearson") + ylim(0, 1) + xlim(0, 1) + coord_fixed()+ theme(axis.text=element_text(size=20),text = element_text(size=20))  + annotate("text", x=c(0.05,0.95,0.8,0.61), y=c(0.01,0.83,0.13,0.1)+c(0.18,0.15,0.15,0.15), label= c("B","C","E","F"), size=13)  + geom_segment(aes(x=z, y=w+v, xend=z, yend=w), data=df, arrow = arrow(length = unit(0.03, "npc")),size=1) + annotate("text", x=c(0.91,0.99,0.99,0.94), y=c(0.42,0.95,0.05,0.05)+c(0.18,-0.15,0.15,0.15), label= c("G","H","  I","J"), size=13)  + geom_segment(aes(x=z, y=w+v, xend=z, yend=w), data=df2, arrow = arrow(length = unit(0.03, "npc")),size=1)

ggsave("WHO_figures/WHO_Cor_vs_HC.pdf")

dev.off()



# Fig 9-D (MIC vs HC)
plot.new()
par(mfrow = c(1,1), cex = 0.45)
df = data.frame(val.HC,val.mine,z=c(0.05,0.95,0.8,0.61),w=c(0.18,0.98,0.19,0.28),v=c(0.1,0.1,0.1,0.1))
df2 = data.frame(val.HC,val.mine,z=c(0.91,0.99,0.99,0.94),w=c(0.29,0.21,0.26,0.21),v=c(0.1,0.1,-0.1,0.1))
ggplot(df,aes(x=val.HC,y=val.mine))   + geom_point(color="blue") +xlab("Hypercontractivity")+ylab("MIC") + ylim(0, 1) + xlim(0, 1) + coord_fixed()+ theme(axis.text=element_text(size=20), text = element_text(size=20)) + geom_segment(aes(x=z, y=w-v, xend=z, yend=w), data=df, arrow = arrow(length = unit(0.03, "npc")),size=1) + annotate("text", x=c(0.05,0.95,0.8,0.61), y=c(0.18,0.98,0.19,0.28)-c(0.15,0.15,0.15,0.15), label= c("B","C","E","F"), size=13) + geom_segment(aes(x=z, y=w-v, xend=z, yend=w), data=df2, arrow = arrow(length = unit(0.03, "npc")),size=1) + annotate("text", x=c(0.91,0.99,0.99,0.94), y=c(0.29,0.21,0.26,0.21)-c(0.15,0.15,-0.15,0.15), label= c("G "," H"," I","J"), size=13)  

ggsave("WHO_figures/WHO_MIC_vs_HC.pdf")

dev.off()


###### Find indices of samples for which correlation is small but Hypercontractivity is large
for(ind0 in c(1:1600)){
	if (val.cor[ind0]<0.5 && val.HC[ind0]>0.6){
		print(c(val.var1[ind0],val.var2[ind0]))
		print(c(colnames(mydata)[val.var1[ind0]],colnames(mydata)[val.var2[ind0]]))
		#myWHO_country_scatterplots(val.var1[ind],val.var2[ind],1.06)		
	}
}

# Figure 9- B,C,E,F,G,H,I,J   (Scatter plots of countries) and Figure 10 - 15
myWHO_country_scatterplots(178,217,1.06) # Bad teeth per child vs. Democracy score
myWHO_country_scatterplots(186,221,1.06) # CO2 vs. energy use
myWHO_country_scatterplots(172,249,1.06) # aid vs. income growth
myWHO_country_scatterplots(176,241,1.06) # arms exports vs. health expenditure total
myWHO_country_scatterplots(244,221,1.06) # hydroE vs. Electricity consumption
myWHO_country_scatterplots(210,241,1.06) # colon male deaths vs. health expenditure
myWHO_country_scatterplots(184,237,1.06) # broadband subscribers vs. health expenditure private
myWHO_country_scatterplots(232,231,1.06) # investment outflow vs. inflow






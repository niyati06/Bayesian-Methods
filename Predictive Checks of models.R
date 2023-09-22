df=read.csv("paleo_dat.csv")
summary(df)
df_quant=df[!is.na(df$Temperature.C),]
df_quant=df_quant %>% group_by(Paleocoordinate.Age) %>% dplyr::mutate(ID = cur_group_id())
df_quant=df_quant[order(df_quant$ID),]
df_quant

#Creating Y and X
Y_quant=cbind(df_quant$ID,df_quant$Temperature.C)
colnames(Y_quant)=c("ID","Y")
X_lat=scale(df_quant$Paleo.Lat,TRUE,TRUE)
X2_lat=scale((df_quant$Paleo.Lat)^2,TRUE,TRUE)
X_quant_scale=cbind(df_quant$ID, X_lat,X2_lat)
X_quant=cbind(df_quant$ID,df_quant$Paleo.Lat,(df_quant$Paleo.Lat)^2)
colnames(X_quant)=c("ID","X","X^2")


#Second Model has lower DIC and WAIC so that is the better model for the fit 
D=samples2[[1]][,1:3]


#Diagnostics Checks and Bayesian P value
bayesian_p_val=rep(0,3)
D0   <- c(   mean(Y_quant[,2]),   max(Y_quant[,2])-min(Y_quant[,2]), 
             sd(Y_quant[,2]))
Dnames <- c("Mean Y", "Range Y","SD Y")

par(mfrow=c(3,2))
for(j in 1:3){
  plot(density(D[,j]),xlim=range(c(D0[j],D[,j])),
       xlab="D",ylab="Posterior probability",
       main=Dnames[j])
  abline(v=D0[j],col=2,lwd=2)
  legend("topright",c("Model 2","True Data"),lty=1,col=1:2,bty="n")
  
  bayesian_p_val[j] <- mean(D[,j]>D0[j]) 
}

names(bayesian_p_val)=Dnames
bayesian_p_val


#Plotting curves with uncertainty
sum=summary(samples2)
samp=samples2[[1]]
beta0=samp[,4:12]
beta1=samp[,13:21]
beta2=samp[,22:30]



X_quant=data.frame(X_quant)
Y_quant=data.frame(Y_quant)
pred=seq(-70,70,length=1133)
pred_scale=scale(pred,TRUE,TRUE)

par(mar=c(2,2,3,4),mfrow=c(3,3))
for (i in 1:9){
  X=X_quant[X_quant$ID==i,]
  Y=Y_quant[Y_quant$ID==i,]
  
  # Plot the posterior of the mean alpha1+age[j]*alpha2
  
  fit <- NULL
  for(j in 1:length(pred)){
    fit <- cbind(fit,beta0[,i]+pred_scale[j]*beta1[,i]+(pred_scale[j]^2)*beta2[,i])
  }
  q  <- apply(fit,2,quantile,c(0.025,0.5,0.975))
  plot(X[,2],Y[,2],xlab="PaleoLatitude",ylab="Temperature",
       cex.lab=1.5,cex.axis=1.5,xlim=c(-70,70),ylim=c(0,50),main=paste("g(X)",i))
  lines(pred,q[1,],lty=2,col='red')
  lines(pred,q[2,],lty=1,col='green')
  lines(pred,q[3,],lty=2,col='red')
}

#Plotting MAT
par(mfrow=c(1,1))
plot(NA,NA,xlab='Paleo Latitude',ylab='MAT',xlim=c(-70,70),ylim=c(0,50),main='Comparative MAT across time slice')

for (i in 1:9){
  X=X_quant[X_quant$ID==i,]
  Y=Y_quant[Y_quant$ID==i,]
  
  # Plot the posterior of the mean alpha1+age[j]*alpha2
  
  fit <- NULL
  for(j in 1:length(pred)){
    fit <- cbind(fit,mean(beta0[,i])+pred_scale[j]*mean(beta1[,i])+(pred_scale[j]^2)*mean(beta2[,i]))
  }
  q=apply(fit,2,quantile,c(.025,0.5,0.75))
  lines(pred,q[1,],lty=2,col=i)
  lines(pred,q[2,],lty=1,col=i)
  lines(pred,q[3,],lty=2,col=i)
}
legend("topleft",c('g1(x)','g2(x)','g3(x)','g4(x)','g5(x)','g6(x)','g7(x)','g8(x)','g9(x)'),col=1:9,lty=1,bty='n')

 
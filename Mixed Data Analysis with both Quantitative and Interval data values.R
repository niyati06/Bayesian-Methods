#Part2

#Use model 1 from before 
#Create data
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

df_mixed=df[!is.na(df$Min.Temp),]
df_mixed=df_mixed %>% group_by(Paleocoordinate.Age) %>% dplyr::mutate(ID = cur_group_id())
df_mixed=df_mixed[order(df_mixed$ID),]
Y_mixed=cbind(df_mixed$ID, df_mixed$Min.Temp,df_mixed$Max.Temp)
X_lat=scale(df_mixed$Paleo.Lat)
X2_lat=scale((df_mixed$Paleo.Lat)^2)
X_mixed_scale=cbind(df_mixed$ID,X_lat,X2_lat)
X_mixed=cbind(df_mixed$ID,df_mixed$Paleo.Lat,(df_mixed$Paleo.Lat)^2)
Y_mixed=data.frame(Y_mixed)
X_mixed=data.frame(X_mixed)
one=rep(1,2120)

model_string=textConnection("model{
                            
                            #likelihood for quantitative observations
                            for (i in 1:n1){
                                Y1[i,2]~dnorm(mn[i],taue[i])
                                taue[i]=1/sig2[i]
                                log(sig2[i])=alpha0[gk[i]] + alpha1[gk[i]]*X[i,2]+alpha2[gk[i]]*X[i,3]
                                mn[i]=beta0[gk[i]] + beta1[gk[i]]*X[i,2] + beta2[gk[i]]*X[i,3]
                            }
                              
                            #likelihood for interval observations
                            for (j in 1:n2){
                                one[j]~dbern(p[j])
                                logit(p[j])=pnorm(Y2[j,3],m[j],taue[j])-pnorm(Y2[j,2],m[j],tauf[j])
                                tauf[j]=1/sig2f[j]
                                log(sig2f[j])=alpha0[gk[j]] + alpha1[gk[j]]*X2[j,2]+alpha2[gk[j]]*X2[j,3]
                                m[j]=beta0[gk2[j]] + beta1[gk2[j]]*X2[j,2] + beta2[gk2[j]]*X2[j,3]
                              }
                            
                           #Priors
                           for (j in 1:9){
                            beta0[j]~dnorm(0,0.001)
                           beta1[j]~dnorm(0,taub)
                           beta2[j]~dnorm(0,taub)
                           alpha0[j]~dnorm(0,0.001)
                           alpha1[j]~dnorm(0,taua)
                           alpha2[j]~dnorm(0,taua)
                            }
                           taua~dgamma(0.1,0.1)
                           taub~dgamma(0.1,0.1)
                           
                           #Diagnostic Checks
                           for (i in 1:n1){
                           Ypred[i]~dnorm(mn[i],taue[i])
                           }
                           
                           D[1]=mean(Ypred[])
                           D[2]=max(Ypred[])-min(Ypred[])
                           D[3]=sd(Ypred[])
                           
                      
                      }")

#compiling model 
model=jags.model(model_string,data=list(one=one,n1=nrow(Y_quant),n2=nrow(Y_mixed),Y1=Y_quant,Y2=Y_mixed,X=X_quant_scale,X2=X_mixed_scale,gk=df_quant$ID,gk2=df_mixed$ID),n.chain=2,quiet=TRUE)
#burn in for 5000 samples
update(model, 5000, progress.bar="none")
params  <- c("beta0","beta1","beta2","D")

#Generate post burn out samples
samples1 <- coda.samples(model, 
                         variable.names=params, 
                         n.iter=5000, progress.bar="none")

#summarise ouputs
summary(samples1)
D=samples1[[1]][,1:3]

#Bayesian P Values
bayesian_p_val=rep(0,3)
D0   <- c(   mean(Y_quant[,2]),   max(Y_quant[,2])-min(Y_quant[,2]), 
             sd(Y_quant[,2]))
Dnames <- c("MeanY", "Range Y","SD Y")

par(mfrow=c(3,2))
for(j in 1:3){
  plot(density(D[,j]),xlim=range(c(D0[j],D[,j])),
       xlab="D",ylab="Posterior probability",
       main=Dnames[j])
  abline(v=D0[j],col=2,lwd=2)
  legend("topright",c("Model Output","True Data"),lty=1,col=1:2,bty="n")
  
  bayesian_p_val[j] <- mean(D[,j]>D0[j]) 
}

names(bayesian_p_val)=Dnames
bayesian_p_val

#plots
sum=summary(samples1)
samp=rbind(samples1[[1]],samples1[[2]])
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
  lines(pred,q[1,],lty=2,col=i,lwd=2)
  lines(pred,q[2,],lty=1,col=i,lwd=2)
  lines(pred,q[3,],lty=2,col=i,lwd=2)
}
legend("topleft",c('g1(x)','g2(x)','g3(x)','g4(x)','g5(x)','g6(x)','g7(x)','g8(x)','g9(x)'),col=1:9,lty=1,bty='n')


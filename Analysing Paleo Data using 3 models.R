#Exam2
df=read.csv("paleo_dat.csv")
summary(df)
df_quant=df[!is.na(df$Temperature.C),]
df_quant=df_quant %>% group_by(Paleocoordinate.Age) %>% dplyr::mutate(ID = cur_group_id())
df_quant=df_quant[order(df_quant$ID),]
df_quant

#Analysing Data
plot(subset(df_quant,ID=='2')$Paleo.Lat,subset(df_quant,ID=='2')$Temperature.C)
plot(subset(df_quant,ID=='3')$Paleo.Lat,subset(df_quant,ID=='3')$Temperature.C)
plot(subset(df_quant,ID=='1')$Paleo.Lat,subset(df_quant,ID=='1')$Temperature.C)

 
#Creating Y and X
Y_quant=cbind(df_quant$ID,df_quant$Temperature.C)
colnames(Y_quant)=c("ID","Y")
X_lat=scale(df_quant$Paleo.Lat,TRUE,TRUE)
X2_lat=scale((df_quant$Paleo.Lat)^2,TRUE,TRUE)
X_quant_scale=cbind(df_quant$ID, X_lat,X2_lat)
X_quant=cbind(df_quant$ID,df_quant$Paleo.Lat,(df_quant$Paleo.Lat)^2)
colnames(X_quant)=c("ID","X","X^2")



library(rjags)
#Varying slopes with normal uninformed priors
#creating model specification
model_string1=textConnection("model{
                            
                           #Likelihood
                            for (i in 1:n){
                            Y[i,2]~dnorm(mn[i],taue[gk[i]])
                            mn[i]=beta0[gk[i]] + beta1[gk[i]]*X[i,2] + beta2[gk[i]]*X[i,3]
                            }
                              
                           
                           #Priors
                           for (j in 1:9){
                           beta0[j]~dnorm(0,0.001)
                           beta1[j]~dnorm(0,0.001)
                           beta2[j]~dnorm(0,0.001)
                           taue[j]~dgamma(0.1,0.1)
                           sigmae[j]=1/sqrt(taue[j])
                           }
                           
                           
                           # WAIC calculations
                          for(i in 1:n){
                              like[i]    <- dnorm(Y[i,2],mu[i],taue[gk[i]])
                              mu[i] <- beta0[gk[i]] + beta1[gk[i]]*X[i,2] + beta2[gk[i]]*X[i,3]
                          }
                          
                          #Diagnostic Checks
                          for (i in 1:n){
                          Y2[i]~dnorm(mn[i],taue[gk[i]])
                          }
                          
                            D[1]<- mean(Y2[])
                            D[2]<- max(Y2[])-min(Y2[])
                            D[3]<-sd(Y2[])
                          
                          
                          
}")


#compiling model 
model=jags.model(model_string1,data=list(n=nrow(X_quant),Y=Y_quant,X=X_quant_scale,gk=df_quant$ID),n.chain=2,quiet=TRUE)
#burn in for 5000 samples
update(model, 5000, progress.bar="none")
params  <- c("beta0","beta1","beta2","sigmae",'D')

#Generate post burn out samples
samples1 <- coda.samples(model, 
                        variable.names=params, 
                        n.iter=10000, progress.bar="none")

 #summarise ouputs
summary(samples1)

#Calculate ESS
effectiveSize(samples1)
#Calculate Gelman Diagnostic
gelman.diag(samples1)

# Compute DIC
dic1   <- dic.samples(model,n.iter=10000,progress.bar="none")

# Compute WAIC
waic1   <- coda.samples(model, 
                        variable.names=c("like"), 
                        n.iter=10000, progress.bar="none")
like1   <- waic1[[1]]
fbar1   <- colMeans(like1)
P1      <- sum(apply(log(like1),2,var))
WAIC1   <- -2*sum(log(fbar1))+2*P1



#Varying slope with normal informed priors and heteroscadastic variance
#creating model specification
model_string2=textConnection("model{
                            
                           #Likelihood
                            for (i in 1:n){
                            Y[i,2]~dnorm(mn[i],taue[i])
                            taue[i]=1/sig2[i]
                            log(sig2[i])=alpha0[gk[i]] + alpha1[gk[i]]*X[i,2]+alpha2[gk[i]]*X[i,3]
                            mn[i]=beta0[gk[i]] + beta1[gk[i]]*X[i,2] + beta2[gk[i]]*X[i,3]
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
                           
                           taub~dgamma(0.1,0.1)
                           taua~dgamma(0.1,0.1)
                           
                           
                           # WAIC calculations
                          for(i in 1:n){
                              like[i]    <- dnorm(Y[i,2],mn[i],taue[i])
                          }
                               
                            #Diagnostic Checks
                          for (i in 1:n){
                          Y2[i]~dnorm(mn[i],taue[i])
                          }
                          
                            D[1]<- mean(Y2[])
                            D[2]<- max(Y2[])-min(Y2[])
                            D[3]<-sd(Y2[])
                          
}")

#compiling model 
model2=jags.model(model_string2,data=list(n=nrow(X_quant_scale),Y=Y_quant,X=X_quant_scale,gk=df_quant$ID),n.chain=2,quiet=TRUE)

#burn in for 5000 samples
update(model2, 10000, progress.bar="none")
params  <- c("beta0","beta1","beta2","D")

#Generate post burn out samples
samples2 <- coda.samples(model2, 
                         variable.names=params, 
                         n.iter=25000, progress.bar="none")

#summarise ouputs
summary(samples2)


#Calculate ESS
effectiveSize(samples2)
#Calculate Gelman Diagnostic
gelman.diag(samples2)

# Compute DIC
dic2  <- dic.samples(model2,n.iter=10000,progress.bar="none")

# Compute WAIC
waic2   <- coda.samples(model2, 
                        variable.names=c("like"), 
                        n.iter=10000, progress.bar="none")
like2   <- waic2[[1]]
fbar2   <- colMeans(like2)
P2      <- sum(apply(log(like2),2,var))
WAIC2   <- -2*sum(log(fbar2))+2*P2

#varying regression coeffcient with LASSO priors and constant variance
model_string3=textConnection("model{
                            
                           #Likelihood
                            for (i in 1:n){
                            Y[i,2]~dnorm(mn[i],taue)
                            mn[i]=beta0 + beta1[gk[i]]*X[i,2] + beta2[gk[i]]*X[i,3]
                            }
                              
                           
                           #Priors
                           beta0~dnorm(0,0.001)
                           for (j in 1:9){
                           beta1[j]~ddexp(0,taue*taub)
                           beta2[j]~ddexp(0,taue*taub)
                           }
                           taue~dgamma(0.1,0.1)
                           taub~dgamma(0.1,0.1)
                           
                           # WAIC calculations
                          for(i in 1:n){
                              like[i]    <- dnorm(Y[i,2],mn[i],taue)
                          }
                          
                            
                           
}")

#compiling model 
model3=jags.model(model_string3,data=list(n=nrow(X_quant_scale),Y=Y_quant,X=X_quant_scale,gk=df_quant$ID),n.chain=2,quiet=TRUE)
#burn in for 5000 samples
update(model3, 5000, progress.bar="none")
params  <- c("beta0","beta1","beta2")

#Generate post burn out samples
samples3 <- coda.samples(model3, 
                         variable.names=params, 
                         n.iter=10000, progress.bar="none")

#summarise ouputs
summary(samples3)


#Calculate ESS
effectiveSize(samples3)
#Calculate Gelman Diagnostic
gelman.diag(samples3)

# Compute DIC
dic3    <- dic.samples(model3,n.iter=10000,progress.bar="none")

# Compute WAIC
waic3   <- coda.samples(model3, 
                        variable.names=c("like"), 
                        n.iter=10000, progress.bar="none")
like3   <- waic3[[1]]
fbar3   <- colMeans(like3)
P3      <- sum(apply(log(like3),2,var))
WAIC3   <- -2*sum(log(fbar3))+2*P3


#Comparing DIC and WAIC for best model
dic1
dic2
dic3

WAIC1
WAIC2
WAIC3




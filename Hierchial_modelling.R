df=read.csv("train_80.csv")
unique(df$lat)
unique(df$lon)
df$latlon<-paste0(as.character(df$lat),", ", as.character(df$lon))
df1=df %>% group_by(latlon) %>% dplyr::mutate(ID = cur_group_id())
df1 <- df1[, -which(names(df1) == "latlon")]
df1=df1[,-1]
df1 <- select(df1, ID, everything())


# Load the ggplot2 library
library(ggplot2)

# Set the y-column to plot against
y_col <- "contest.tmp2m.14d__tmp2m"

# Create an empty list to store the plots
plots <- list()

# Loop through each column (except the x-column) and create a plot
for (col in names(df1)[-27]) {
  # Create a plot and add it to the list
  plots[[col]] <- ggplot(df, aes(x = !!sym(col), y = !!sym(y_col))) +
    geom_point() +
    labs(title = paste(col, "vs.", y_col))
}

# Arrange and print all the plots
library(gridExtra)
grid.arrange(grobs = plots)


#Defining Y and X
Y=as.data.frame(cbind(df1$ID,df1$contest.tmp2m.14d__tmp2m))
X=as.data.frame(df1[, -which(names(df1) == "contest.tmp2m.14d__tmp2m")])
colnames(Y)=c("ID","Temperature")
n=nrow(X)
p=ncol(X)


#SPline Variables: daynum, nmme.tmp2m.56w__gfdlflorb, 
#contest.rhum.sig995.14d__rhum,nmme0.prate.34w__ccsm30,contest.wind.vwnd.250.14d__wind.vwnd.250,contest.wind.h850.14d__wind.hgt.850
#wind.vwnd.250.2010.3,wind.vwnd.250.2010.13,wind.uwnd.250.2010.18,sst.2010.2,wind.hgt.500.2010.2,wind.uwnd.925.2010.8,wind.uwnd.925.2010.12
#wind.hgt.10.2010.7,wind.vwnd.925.2010.3,

#Linear Variables: contest.slp.14d__slp,contest.wind.vwnd.925.14d__wind.vwnd.925,contest.pres.sfc.gauss.14d__pres
#nmme.tmp2m.34w__gfdlflorb,nmme.tmp2m.34w__nasa,contest.prwtr.eatm.14d__prwtr,contest.wind.h500.14d__wind.hgt.500
#Elevation

#Setting up spline function for Spline variables

X_spline = df1[, names(df1) %in% c("daynum", "nmme.tmp2m.56w__gfdlflorb", "contest.rhum.sig995.14d__rhum", "nmme0.prate.34w__ccsm30", "contest.wind.vwnd.250.14d__wind.vwnd.250", "contest.wind.h850.14d__wind.hgt.850", "wind.vwnd.250.2010.3", "wind.vwnd.250.2010.13", "wind.uwnd.250.2010.18", "sst.2010.2", "wind.hgt.500.2010.2", "wind.uwnd.925.2010.8", "wind.uwnd.925.2010.12", "wind.hgt.10.2010.7", "wind.vwnd.925.2010.3")]
X_lin = df1[,names(df1) %in% c("lat","lon","contest.slp.14d__slp","contest.wind.vwnd.925.14d__wind.vwnd.925","contest.pres.sfc.gauss.14d__pres",
                               "nmme.tmp2m.34w__gfdlflorb","nmme.tmp2m.34w__nasa","contest.prwtr.eatm.14d__prwtr","contest.wind.h500.14d__wind.hgt.500",
                               "elevation__elevation")]

X_spline=cbind(df1$ID,X_spline)
colnames(X_spline)[1]='ID'

X_lin=cbind(df1$ID,X_lin)
colnames(X_lin)[1]='ID'

dim(X_spline)
dim(X_lin)

library(splines)

#Linear Hierchial 
set.seed(1)

X_scale=as.data.frame(cbind(X$ID,scale(X[,2:26])))
colnames(Y)=c("ID","Temperature")
library(rjags)


modelstring_heir1=textConnection("model{
                                 
                                 #Data Layer
                                 for (i in 1:n){
                                 Y[i,2]~dnorm(mn[i],taue[i])
                                 taue[i]=1/sig2[i]
                                 log(sig2[i])=alpha0[gk[i]]+inprod(gamma[gk[i],],X[i,2:26])
                                 mn[i]= beta0[gk[i]]+ inprod(beta[gk[i], ],X[i,2:26])
                                 }
                                 
                                 #Process layer
                                 for (j in 1:26){
                                 beta0[j]~dnorm(0.1,0.1)
                                  alpha0[j]~dnorm(0.1,0.1)
                                 }
                                 
                                 for (j in 1:26){
                                  for (k in 1:25){
                                    beta[j,k]~dnorm(0,taub)
                                    gamma[j,k]~dnorm(0,taug)
                                  }}
                                 
                                 
                                 #Prior Layer
                                 taub~dgamma(0.1,0.1)
                                 taug~dgamma(0.1,0.1)
                                
}")

#compiling model 
model1=jags.model(modelstring_heir1,data=list(n=nrow(X_scale),Y=Y,X=X_scale,gk=Y$ID),n.chain=2,quiet=TRUE)
#burn in for 5000 samples
update(model1, 5000, progress.bar="none")
params  <- c("beta0","beta")

#Generate post burn out samples
samples1 <- coda.samples(model1, 
                         variable.names=params, 
                         n.iter=10000, progress.bar="none")

#summarise ouputs
summary(samples1)

#Calculate ESS
effectiveSize(samples1)
#Calculate Gelman Diagnostic
gelman.diag(samples1)

# Compute DIC
dic2  <- dic.samples(model1,n.iter=10000,progress.bar="none")
dic2



samp=rbind(samples1[[1]],samples1[[2]])
beta1=samp[,1:26]
beta2=samp[,27:52]
beta3=samp[,53:78]
beta4=samp[,79:104]
beta5=samp[,105:130]
beta6=samp[,131:156]
beta7=samp[,157:182]
beta8=samp[,183:208]
beta9=samp[,209:234]
beta10=samp[,235:260]
beta11=samp[,261:286]
beta12=samp[,287:312]
beta13=samp[,313:338]
beta14=samp[,339:364]
beta15=samp[,365:390]
beta16=samp[,391:416]
beta17=samp[,417:442]
beta18=samp[,443:468]
beta19=samp[,469:494]
beta20=samp[,495:520]
beta21=samp[,521:546]
beta22=samp[,547:572]
beta23=samp[,573:598]
beta24=samp[,599:624]
beta25=samp[,625:650]
beta0=samp[,651:676]

write.csv(beta1, "beta1.csv")
write.csv(beta2, "beta2.csv")
write.csv(beta3, "beta3.csv")
write.csv(beta4, "beta4.csv")
write.csv(beta5, "beta5.csv")
write.csv(beta6, "beta6.csv")
write.csv(beta7, "beta7.csv")
write.csv(beta8, "beta8.csv")
write.csv(beta9, "beta9.csv")
write.csv(beta10, "beta10.csv")
write.csv(beta11, "beta11.csv")
write.csv(beta12, "beta12.csv")
write.csv(beta13, "beta13.csv")
write.csv(beta14, "beta14.csv")
write.csv(beta15, "beta15.csv")
write.csv(beta16, "beta16.csv")
write.csv(beta17, "beta17.csv")
write.csv(beta18, "beta18.csv")
write.csv(beta19, "beta19.csv")
write.csv(beta20, "beta20.csv")
write.csv(beta21, "beta21.csv")
write.csv(beta22, "beta22.csv")
write.csv(beta23, "beta23.csv")
write.csv(beta23, "beta24.csv")
write.csv(beta24, "beta25.csv")
write.csv(beta0, "beta0.csv")



#Predictions
testdf=read.csv("test.csv")
testdf$latlon<-paste0(as.character(testdf$lat),", ", as.character(testdf$lon))
test1=testdf %>% group_by(latlon) %>% dplyr::mutate(ID = cur_group_id())
test1 <- test1[, -which(names(test1) == "latlon")]
test1=test1[,-1]
test1 <- select(test1, ID, everything())

Y_test=as.data.frame(cbind(test1$ID,test1$contest.tmp2m.14d__tmp2m))
X_test=as.data.frame(test1[, -which(names(test1) == "contest.tmp2m.14d__tmp2m")])
colnames(Y_test)=c("ID","Temperature")


X_scale_test=as.data.frame(cbind(X_test$ID,scale(X_test[,2:26])))
n_test=nrow(X_test)

#Fitting
prediction=matrix(0,n_test,3)
for (i in 1:n_test){
  g=X_scale_test$V1[i]
  fit=beta0[,g]+beta1[,g]*X_scale_test[i,2]+beta2[,g]*X_scale_test[i,3]+beta3[,g]*X_scale_test[i,4]+beta4[,g]*X_scale_test[i,5]+beta5[,g]*X_scale_test[i,6]+
    beta6[,g]*X_scale_test[i,7]+beta7[,g]*X_scale_test[i,8]+beta8[,g]*X_scale_test[i,9]+beta9[,g]*X_scale_test[i,10]+
    beta10[,g]*X_scale_test[i,11]+beta11[,g]*X_scale_test[i,12]+beta12[,g]*X_scale_test[i,13]+beta13[,g]*X_scale_test[i,14]+
    beta14[,g]*X_scale_test[i,15]+beta15[,g]*X_scale_test[i,16]+beta16[,g]*X_scale_test[i,17]+beta17[,g]*X_scale_test[i,18]+
    beta18[,g]*X_scale_test[i,19]+beta19[,g]*X_scale_test[i,20]+beta20[,g]*X_scale_test[i,21]+beta21[,g]*X_scale_test[i,22]+beta22[,g]*X_scale_test[i,23]+
    beta23[,g]*X_scale_test[i,24]+beta24[,g]*X_scale_test[i,25]+beta25[,g]*X_scale_test[i,26]
  q=c(quantile(fit,0.025),quantile(fit,0.5),quantile(fit,0.975))
  prediction[i,]=q
}

plot(X_scale_test$daynum,prediction[,2])

#RMSE
MSE=sum((prediction[,2]-Y_test[,2])^2)/n_test
MSE

RMSE=sqrt(MSE)
RMSE

MAE=mae(prediction[,2],Y_test[,2])
MAE

#R^2
RSS=sum((Y_test[,2]-prediction[,2])^2)
TSS=sum((Y_test[,2]-mean(Y_test[,2]))^2)
R2=1-(RSS/TSS)
R2

plot(X_test[X_test$ID==1,]$daynum,Y_test[Y_test$ID==1,][,2],main="Temperature v/s Days at Location 1",xlab="Number of Days",ylab="Temperature")
lines(X_test[X_test$ID==1,]$daynum,prediction[1:134,1],col=2,lwd=2)
lines(X_test[X_test$ID==1,]$daynum,prediction[1:134,2],col=3,lwd=2)
lines(X_test[X_test$ID==1,]$daynum,prediction[1:134,3],col=4,lwd=2)
legend("bottomright",c('True Data','2.5% CI','50% CI','97.5% CI'),col=1:4,lty=1,bty='n')


plot(X_test$lat,Y_test[,2],main="Temperature v/s Latitude",xlab="Latitude",ylab="Temperature")
lines(X_test$lat,prediction[,1],col=2,lwd=2)
lines(X_test$lat,prediction[,2],col=3,lwd=2)
lines(X_test$lat,prediction[,3],col=4,lwd=2)
legend("bottomleft",c('True Data','2.5% CI','50% CI','97.5% CI'),col=1:4,lty=1,bty='n')



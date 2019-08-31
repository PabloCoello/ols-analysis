#Analisis econometrico previo de los datos
if(!require(R6)){install.packages("R6")} 
if(!require(faraway)){install.packages("faraway")} 
if(!require(lmtest)){install.packages("lmtest")} 
if(!require(readxl)){install.packages("readxl")} 
if(!require(leaps)){install.packages("leaps")} 
if(!require(mctest)){install.packages("mctest")} 
if(!require(EnvStats)){install.packages("EnvStats")} 
if(!require(MASS)){install.packages("MASS")} 
library(R6)
library(faraway)
library(lmtest)
library(readxl)
library(leaps)
library(mctest)
library(EnvStats)
library(MASS)

ols_analysis <- R6Class("ols_analysis",
                       public=list(
                         initialize=function(){},
                         get_model = function(fitdata, regresando, fit, coef){
                           fitdata=fitdata[,c(regresando,"Cranes","Tugs","meters","storages","wages")]
                           regresion= function(formula, data){
                             lm(formula,
                                data,
                                na.action = na.exclude)
                           }
                           data.leaps=fitdata[,c(regresando,"Cranes","Tugs","meters","storages","wages")]
                           test=mctest(data.leaps[,-1],
                                       data.frame(data.leaps[,1]))
                           if (sum(test$odiags[,2])>3){
                             return("MULTICOLINEALIDAD")
                           } else{
                           fit.leaps<- regsubsets(data.leaps[,1]~
                                                    Cranes+
                                                    Tugs+
                                                    meters+
                                                    storages+
                                                    wages,
                                                  data=data.leaps,
                                                  method = "exhaustive")
                           sum=summary(fit.leaps)
                           regresores=c("intercept","Cranes","Tugs","meters","storages","wages")
                           for (i in 2:length(data.leaps)-1){
                             n=which(sum$which[i,])[-1]
                             assign(paste("mod",i,sep=""),paste(regresando,"~",paste(regresores[n],sep="",collapse = "+"),sep=""))
                             
                           }
                           modelos=c("mod1","mod2","mod3","mod4","mod5","fit")
                           
                           #TEST BOXCOX
                           for (i in 1:length(modelos)){
                             assign(paste(modelos[i],"box",sep=""),boxcox(object =as.formula(get(modelos[i])),
                                                                   data = fitdata
                             ))
                             assign(paste(modelos[i],"lambda",sep=""),get(paste(modelos[i],"box",sep=""))$x[which.max(get(paste(modelos[i],"box",sep=""))$y)])
                           
                           }
                           
                           
                           for (i in 1:length(modelos)){
                             
                             if (get(paste(modelos[i],"lambda",sep="")) < (-0.75)){
                               assign(paste(modelos[i]),regresion(formula=get(modelos[i]),data=1/fitdata))
                             } else if (get(paste(modelos[i],"lambda",sep="")) < (-0.25) & 
                                      (-0.75) < get(paste(modelos[i],"lambda",sep=""))){
                               assign(paste(modelos[i]),regresion(formula=get(modelos[i]),data=1/sqrt(fitdata)))
                             } else if (get(paste(modelos[i],"lambda",sep="")) < (0.25) &
                                      (-0.25) < get(paste(modelos[i],"lambda",sep=""))){
                               assign(paste(modelos[i]),regresion(formula=get(modelos[i]),data=log(fitdata)))
                             } else if (get(paste(modelos[i],"lambda",sep="")) < (0.75) &
                                      (0.25) < get(paste(modelos[i],"lambda",sep=""))){
                               assign(paste(modelos[i]),regresion(formula=get(modelos[i]),data=sqrt(fitdata)))
                             } else {
                               assign(paste(modelos[i]),regresion(formula=get(modelos[i]),data=fitdata))
                             }
                             
                             assign(paste("ak",modelos[i],sep="_"),round(AIC(get(modelos[i])),3)) # Criterio Akaike
                             assign(paste("bic",modelos[i],sep="_"),round(BIC(get(modelos[i])),3)) # Criterio de información bayesiano
                             assign(paste("r2ajustado",modelos[i],sep="_"),round(summary(get(modelos[i]))$adj.r.squared,3)) # R cuadrado ajustado
                             assign(paste("bp",modelos[i],sep="_"),round(bptest(get(modelos[i]))$p.value[1],5)) # heterocedasticidad: Breusch-Pagan Test
                             assign(paste("shap",modelos[i],sep="_"),round(shapiro.test(rstandard(get(modelos[i])))$p.value,5)) # normalidad de los residuos: test shapiro
                             assign(paste("ks",modelos[i],sep="_"),round(ks.test(rstandard(get(modelos[i])),pnorm)$p.value,3)) # normalidad de los residuos: test kolgomorov
                           }
                           akaike=cbind(mget(paste("ak",modelos,sep="_")))
                           infbayes=cbind(mget(paste("bic",modelos,sep="_")))
                           r2ajustado=cbind(mget(paste("r2ajustado",modelos,sep="_")))
                           breuschpagan=cbind(mget(paste("bp",modelos,sep="_")))
                           shapiro=cbind(mget(paste("shap",modelos,sep="_")))
                           kolgomorov=cbind(mget(paste("ks",modelos,sep="_")))
                           lambda=cbind(mget(paste(modelos,"lambda",sep="")))
                           
                           (result=cbind.data.frame(modelos,akaike,infbayes,r2ajustado,breuschpagan,shapiro, kolgomorov,lambda))
                           if (any(result[,5]>coef & result[,6]>coef & result[,7]>coef)){
                             res=result[which(result[,5]>coef & result[,6]>coef & result[,7]>coef),]
                           } else if (any(result[,5]>coef & result[,7]>coef)){
                             res=result[which(result[,5]>coef & result[,7]>coef),]
                           } else if (any(result[,5]>coef)){
                             res=result[which(result[,5]>coef),]
                           } else if ((any(result[,5]>coef | result[,6]>coef | result[,7]>coef))){
                             res=result[which(result[,5]>coef | result[,6]>coef | result[,7]>coef),]
                           } else{
                             res=result
                           }
                           res=as.matrix(res[(which.max(as.numeric(res[,4]))),])
                           especification=get(res[[1]])
                           return(list(especification, res, result, mod1,mod2,mod3,mod4,mod5,fit ))
                           }
                         }
                         ),
                         private = list(
                         )
                       )
                         


set.seed(12345)
library(MASS)

#get distr of NB parameter distribution in tmp2Rna

n <- 500 #sample size
parEstim <- lapply(sample(x=1:dim(tmp2Rna)[1],size=n,replace=F), function(x){
  fitdistr(x=round(as.numeric(tmp2Rna[x,])),densfun='negative binomial')
})

sizeEst <- sapply(1:n, function(x){
  parEstim[[x]]$estimate[1]
})
summary(sizeEst)

muEst <- sapply(1:n, function(x){
  parEstim[[x]]$estimate[2]
})
summary(muEst)

##simulation framework####

#what is the interesting effect size
#the estimator would relate to a 1% change of mutation load
# so eg. exp(0.01)=1.01005, so an increase by 50% would be 1.01005^50=1.648708
# i.e. comparing an indivdual with 0% AKT1mut and 100 reads for gene X to one with 50% AKTmut, the latter
# has 164 reads of gene X
#power as fx of effect size#####
for(i in c(0.005,0.01,0.015,0.02,0.025,0.03)){
  power60 <- lapply(1:100, function(x){
    nMut <- 20
    nNotMut <- nMut*2 #not mutated samples have a 2:1 ratio to mutated ones
    n <- nMut+nNotMut
    #decide on distribution of Mutations, here we'll use uniform
    AKT <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    PIK3CA <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    TP53 <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    GeneX <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    age <- round(rnorm(n=n, mean=68.8, sd=9.4))
    df <- cbind(AKT,PIK3CA,TP53,GeneX,age)
    beta <- c(i,rnorm(1,sd=0.01),rnorm(1,sd=0.01),rnorm(1,sd=0.001),rnorm(1,sd=0.0001))
    y <- round(exp(beta%*%t(df))) # at this point we just round y to get integers, this also
    #adds random noise
    out <- tryCatch({summary(glm.nb(formula=y[1,]~AKT+PIK3CA+TP53+GeneX+age))},
                    error=function(e){e})
    return(out)
  })
  pval <- sapply(1:length(power60), function(x){
    if(class(power60[[x]])[1] == "summary.negbin"){
      return(power60[[x]]$coefficients['AKT','Pr(>|z|)'])
    }else{
      return(NA)
    }
  })
  print(i)
  print(exp(i)^50)
  print(sum(pval<=0.05,na.rm=T))
  print(sum(pval*22000<=0.05,na.rm=T)/length(pval))
  cat('##############################\n')
}

#sample size, we expect ~31% samples to be mutatted for AKT1 


#power as fx of sample size only, effect size = 2####
for(i in c(10,20,40,80,160)){
  power60 <- lapply(1:22000, function(x){
    nMut <- i
    nNotMut <- nMut*2 #not mutated samples have a 2:1 ratio to mutated ones
    n <- nMut+nNotMut
    #decide on distribution of Mutations, here we'll use uniform
    AKT <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    PIK3CA <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    TP53 <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    GeneX <- c(rep(0,nNotMut),runif(n=nMut,min=0,max=100))
    age <- round(rnorm(n=n, mean=68.8, sd=9.4))
    df <- cbind(AKT,PIK3CA,TP53,GeneX,age)
    beta <- c(0.015,rnorm(1,sd=0.01),rnorm(1,sd=0.01),rnorm(1,sd=0.001),rnorm(1,sd=0.0001))
    y <- round(exp(beta%*%t(df))) # at this point we just round y to get integers, this also
    #adds random noise
    out <- tryCatch({summary(glm.nb(formula=y[1,]~AKT+PIK3CA+TP53+GeneX+age))},
                    error=function(e){e})
    return(out)
  })
  pval <- sapply(1:length(power60), function(x){
    if(class(power60[[x]])[1] == "summary.negbin"){
      return(power60[[x]]$coefficients['AKT','Pr(>|z|)'])
    }else{
      return(NA)
    }
  })
  print(i*3)
#   print(exp(i)^50)
  print(sum(pval<=0.05,na.rm=T))
  print(sum(p.adjust(pval,method='BH') <= 0.05,na.rm=T))
  print(sum(p.adjust(pval,method='bonf') <= 0.05, na.rm=T))
  print(sum(pval*22000 <= 0.05,na.rm=T)/length(pval))
  cat('##############################\n')
}




pval <- sapply(1:length(power63), function(x){
  if(class(power63[[x]])[1] == "summary.negbin"){
    return(power63[[x]]$coefficients['AKT','Pr(>|z|)'])
  }else{
    return(NA)
  }
})
hist(pval)
hist(pval*22000)
sum(pval*22000<=0.05,na.rm=T)/length(pval)

#beta 0.01 -> 0.002


# n <- 63
# #decide on distribution of Mutations, here we'll use uniform
# AKT <- runif(n=n,min=0,max=100)
# PIK3CA <- runif(n=n,min=0,max=100)
# TP53 <- runif(n=n,min=0,max=100)
# age <- round(rnorm(n=n, mean=68.8, sd=9.4))
# df <- cbind(AKT,PIK3CA,TP53,age)
# beta <- c(0.05,0.000000001,-0.000000001,0)
# y <- exp(beta%*%t(df))
# summary(glm.nb(formula=round(y[1,])~AKT+PIK3CA+TP53+age))$coef['AKT','Pr(>|z|)']
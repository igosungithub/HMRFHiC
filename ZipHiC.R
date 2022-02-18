
library(parallel)
#2D Neighbours for the sufficient statistic in the potts model
Neighbours<-function(potts_data,N){
  
  mydata1=matrix(0,(N+1),N)
  mydata1[1:N,1:N]=potts_data
  neighbour_down=matrix(0,(N+1),N)
  neighbour_down[2:(N+1),1:N]=mydata1[1:N,1:N]
  
  b1=which(mydata1[1:(N+1),]==neighbour_down[1:(N+1),], arr.ind = T)
  b11=matrix(0,N,N)
  b11[b1]<-1
  
  
  neighbour_up=b11[-1,]
  a1=rep(0,N)
  a11=rbind(neighbour_up,a1)
  
  mydata2=matrix(0,N,(N+1))
  mydata2[1:N,2:(N+1)]=potts_data
  
  neighbour_left=matrix(0,N,(N+1))
  neighbour_left[1:N,1:N]=mydata2[1:N,2:(N+1)]
  c1=which(mydata2[,1:(N+1)]==neighbour_left[,1:(N+1)], arr.ind = T)
  c11=matrix(0,N,N)
  c11[c1]<-1
  
  neighbour_right=c11[,-1]
  d1=rep(0,N)
  d11=cbind(neighbour_right,d1)
  d11[d1]<-1
  
  
  
  Neighbours_total=a11+b11+c11+d11
  return(Neighbours_total)
}


#2D Neighbours for the proposed new component in the sufficient statistic of the potts model.
Neighbours_proposed<-function(potts_data, N, proposed_value){
  mydata1=matrix(0,(N+1),N)
  mydata1[2:(N+1),1:N]=potts_data
  neighbour_down=matrix(0,(N+1),N)
  neighbour_down[1:N,1:N]=proposed_value[1:N,1:N]
  
  a1=which(mydata1[1:(N+1),]==neighbour_down[1:(N+1),], arr.ind = T)
  a11=matrix(0,N,N)
  a11[a1]<-1
  
  
  mydata2=matrix(0,(N+1),N)
  mydata2[1:N,1:N]=potts_data
  neighbour_up=matrix(0,(N+1),N)
  neighbour_up[2:(N+1),1:N]=proposed_value[1:N,1:N]
  b1=which(mydata2[1:(N+1),]==neighbour_up[1:(N+1),], arr.ind = T)
  
  b11=matrix(0,N,N)
  b11[b1]<-1
  neighbour_up=b11[-1,]
  n1=rep(0,N)
  b111=rbind(neighbour_up,n1)
  
  mydata3=matrix(0,N,(N+1))
  mydata3[1:N,2:(N+1)]=potts_data
  
  neighbour_left=matrix(0,N,(N+1))
  neighbour_left[1:N,1:N]=proposed_value[1:N,1:N]
  c1=which(mydata3[,1:(N+1)]==neighbour_left[,1:(N+1)], arr.ind = T)
  c11=matrix(0,N,N)
  c11[c1]<-1
  
  mydata4=matrix(0,N,(N+1))
  mydata4[1:N,1:N]=potts_data
  neighbour_right=matrix(0,N,(N+1))
  neighbour_right[1:N,2:(N+1)]=proposed_value[1:N,1:N]
  d1=which(mydata4[,1:(N+1)]==neighbour_right[,1:(N+1)], arr.ind = T)
  
  d11=matrix(0,N,N)
  d11[d1]<-1
  neighbour_right=d11[,-1]
  p1=rep(0,N)
  d111=cbind(neighbour_right,p1)
  
  
  
  Neighbours_total=a11+b111+c11+d111
  return(Neighbours_total)
}



gamma_prior=rbeta(1,10,5)



expand<-function(ss,N){
  psp=list()
  df=ss
  #df<-do.call(cbind, ss)
  diff1<-length(df[[1]])-length(df[[length(df)]])
  add<-rep(0,diff1)
  bind<-c(df[[length(df)]],add)
  df[[length(df)]]<-bind
  
  for (i in 1:length(df)) {
    psp[[i]]=matrix(df[[i]],N,N)
    
  }
  return(psp)
}



run_metropolis_MCMC_betas <- function(N,startvalue_betas1,startvalue_betas2,startvalue_betas3, gamma_prior,iterations,x1,x2,x3,x4,y,theta){
  
  ####for the lambda of the first component#########
  pred1 <- function(param1,z,x1,x2,x3,x4){
    a1 = param1[1]
    b1 = param1[2]
    c1 = param1[3]
    d1 = param1[4]
    e1 = param1[5]
    
    x11 = x1[z==1]
    x12 = x2[z==1]
    x13 = x3[z==1]
    x14 = x4[z==1]
    
    pred1 = a1 + b1*x11 + c1*x12 + d1*x13 + e1*x14
    return(pred1)
  }
  
  ####for the lambda of the second component#########
  pred2 <- function(param2,z,x1,x2,x3,x4){
    a2 = param2[1]
    b2 = param2[2]
    c2 = param2[3]
    d2 = param2[4]
    e2 = param2[5]
    
    x21 = x1[z==2]
    x22 = x2[z==2]
    x23 = x3[z==2]
    x24 = x4[z==2]
    
    pred2 = a2 + b2*x21 + c2*x22 + d2*x23 + e2*x24
    return(pred2)  
  }
  
  ####for the lambda of the third component#########
  pred3 <- function(param3,z,x1,x2,x3,x4){
    a3 = param3[1]
    b3 = param3[2]
    c3 = param3[3]
    d3 = param3[4]
    e3 = param3[5]
    
    x31 = x1[z==3]
    x32 = x2[z==3]
    x33 = x3[z==3]
    x34 = x4[z==3]
    
    pred3 = a3 + b3*x31 + c3*x32 + d3*x33 + e3*x34
    return(pred3)  
  }
  
  
  ####define the model and likelihood functions  for component 1########
  likelihood1 <- function(pred1,param1,z,y,x1,x2,x3,x4,theta){
    y1=y[z==1]
    singlelikelihoods11 =  ifelse(y1==0,((theta+(1-theta))*(exp(-pred1(param1,z,x1,x2,x3,x4)))),(dpois(y1, lambda = exp(pred1(param1,z,x1,x2,x3,x4)), log = F)*(1-theta)))
    s1=singlelikelihoods11+1
    sumll1 = sum(log(s1))
    #cat("s1=",s1,"\n")
    return(sumll1)  
  }
  
  
  #####define the model and likelihood functions for component 2#######
  likelihood2 <- function(pred2,param2,z,y,x1,x2,x3,x4){
    y2 = y[z==2]
    singlelikelihoods2 = dpois(y2, lambda = exp(pred2(param2,z,x1,x2,x3,x4)), log = T)
    sumll2 = sum(singlelikelihoods2)
    return(sumll2)  
  }
  
  ######define the model and likelihood functions for component 3########
  likelihood3 <- function(pred3,param3,z,y,x1,x2,x3,x4){
    y3 = y[z==3]
    singlelikelihoods3 = dpois(y3, lambda = exp(pred3(param3,z,x1,x2,x3,x4)), log = T)
    sumll3 = sum(singlelikelihoods3)
    return(sumll3)  
  }
  
  
  prior1 <- function(param1,y,x1,x2,x3,x4,z){
    a1 = param1[1]
    b1 = param1[2]
    c1 = param1[3]
    d1 = param1[4]
    e1 = param1[5]
    
    x1=x1[z==1]
    x2=x2[z==1]
    x3=x3[z==1]
    x4=x4[z==1]
    y1=y[z==1]
    
    
    inversesdy=rgamma(1,(10-1)/2,10)
    sdy=1/inversesdy
    meany1=rnorm(1,2,sd(y1)/2000)
    
    inversesd1=rgamma(1,(length(x1)-1)/2,(sum((x1-mean(x1))^2)/2))
    sd1=1/inversesd1
    mean1=rnorm(1,mean(x1),sd(x1)/length(x1))
    
    inversesd12=rgamma(1,(length(x2)-1)/2,(sum((x2-mean(x2))^2)/2))
    sd12=1/inversesd12
    mean12=rnorm(1,mean(x2),sd(x2)/length(x2))
    
    inversesd13=rgamma(1,(length(x3)-1)/2,(sum((x3-mean(x3))^2)/2))
    sd13=1/inversesd13
    mean13=rnorm(1,mean(x3),sd(x3)/length(x3))
    
    inversesd14=rgamma(1,(length(x4)-1)/2,(sum((x4-mean(x4))^2)/2))
    sd14=1/inversesd14
    mean14=rnorm(1,mean(x4),sd(x4)/length(x4))
    
    a1prior = dnorm(a1,meany1,sdy, log=T)
    b1prior = dnorm(b1,mean1, sd1, log = T) 
    c1prior = dnorm(c1,mean12, sd12, log = T) 
    d1prior = dnorm(d1,mean13, sd13, log = T)
    e1prior = dnorm(e1,mean14, sd14, log = T)
    
    #cat("mean13=",mean13,",","sd13=",sd13,"\n")
    return(a1prior+b1prior+c1prior+d1prior+e1prior)
  }
  
  ######Prior distribution for component 2######
  prior2 <- function(param2,y,x1,x2,x3,x4,z){
    a2 = param2[1]
    b2 = param2[2]
    c2 = param2[3]
    d2 = param2[4]
    e2 = param2[5]
    
    x1=x1[z==2]
    x2=x2[z==2]
    x3=x3[z==2]
    x4=x4[z==2]
    y2=y[z==2]
    
    #inversesdy2=rgamma(1,(length(y2)-1)/2,(sum((y2-mean(y2))^2)/2))
    inversesdy2=rgamma(1,(200-1)/2,2000)
    sdy2=1/inversesdy2
    #meany2=rnorm(1,mean(y2),sdy2/length(y2))
    meany2=rnorm(1,100,sd(y2)/5000)
    
    inversesdx1=rgamma(1,(length(x1)-1)/2,(sum((x1-mean(x1))^2)/2))
    sdx1=1/inversesdx1
    meanx1=rnorm(1,mean(x1),sd(x1)/length(x1))
    
    inversesdx2=rgamma(1,(length(x2)-1)/2,(sum((x2-mean(x2))^2)/2))
    sdx2=1/inversesdx2
    meanx2=rnorm(1,mean(x2),sdx2/length(x2))
    
    inversesdx3=rgamma(1,(length(x3)-1)/2,(sum((x3-mean(x3))^2)/2))
    sdx3=1/inversesdx3
    meanx3=rnorm(1,mean(x3),sd(x3)/length(x3))
    
    inversesdx4=rgamma(1,(length(x4)-1)/2,(sum((x4-mean(x4))^2)/2))
    sdx4=1/inversesdx4
    meanx4=rnorm(1,mean(x4),sd(x4)/length(x4))
    
    a2prior = dnorm(a2,meany2,sdy2, log=T)
    b2prior = dnorm(b2, mean = meanx1, sd = sdx1, log = T)
    c2prior = dnorm(c2, mean = meanx2, sd = sdx2, log = T)
    d2prior = dnorm(d2, mean = meanx3, sd = sdx3, log = T)
    e2prior = dnorm(e2, mean = meanx4, sd = sdx4, log = T)
    
    #cat("meany2=",meany2,",","sdy2=",sdy2,"\n")
    
    return(a2prior+b2prior+c2prior+d2prior+e2prior)
  }
  
  
  
  ########Prior distribution for component 3########
  prior3 <- function(param3,y,x1,x2,x3,x4,z){
    a3 = param3[1]
    b3 = param3[2]
    c3 = param3[3]
    d3 = param3[4]
    e3 = param3[5]
    
    x1=x1[z==3]
    x2=x2[z==3]
    x3=x3[z==3]
    x4=x4[z==3]
    y3=y[z==3]
    
    #inversesdy2=rgamma(1,(length(y2)-1)/2,(sum((y2-mean(y2))^2)/2))
    inversesdy3=rgamma(1,(50000-1)/2,5000)
    sdy3=1/inversesdy3
    #meany2=rnorm(1,mean(y2),sdy2/length(y2))
    meany3=rnorm(1,100,sd(y3)/7000)
    
    inversesdx1=rgamma(1,(length(x1)-1)/2,(sum((x1-mean(x1))^2)/2))
    sdx1=1/inversesdx1
    meanx1=rnorm(1,mean(x1),sd(x1)/length(x1))
    
    inversesdx2=rgamma(1,(length(x2)-1)/2,(sum((x2-mean(x2))^2)/2))
    sdx2=1/inversesdx2
    meanx2=rnorm(1,mean(x2),sdx2/length(x2))
    
    inversesdx3=rgamma(1,(length(x3)-1)/2,(sum((x3-mean(x3))^2)/2))
    sdx3=1/inversesdx3
    meanx3=rnorm(1,mean(x3),sd(x3)/length(x3))
    
    inversesdx4=rgamma(1,(length(x4)-1)/2,(sum((x4-mean(x4))^2)/2))
    sdx4=1/inversesdx4
    meanx4=rnorm(1,mean(x4),sd(x4)/length(x4))
    
    
    
    a3prior = dnorm(a3,meany3,sdy3, log=T)
    b3prior = dnorm(b3, mean = meanx1, sd= sdx1, log = T)
    c3prior = dnorm(c3, meanx2, sdx2, log = T)
    d3prior = dnorm(d3, mean = meanx3, sd = sdx3, log = T)
    e3prior = dnorm(e3, mean = meanx4, sd = sdx4, log = T)
    
    sump3=a3prior+b3prior+ c3prior + d3prior + e3prior
    
    #cat("sump3=",sump3,"\n")
    
    return(sump3)
  }
  
  
  #####Posterior distribution for component 1######
  posterior1 <- function(pred1,param1,z,y,x1,x2,x3,x4,theta){
    return (likelihood1(pred1,param1,z,y,x1,x2,x3,x4,theta) + prior1(param1,y,x1,x2,x3,x4,z))
  }
  
  #####Posterior distribution for component 2######
  posterior2 <- function(pred2,param2,z,y,x1,x2,x3,x4){
    return (likelihood2(pred2,param2,z,y,x1,x2,x3,x4) + prior2(param2,y,x1,x2,x3,x4,z))
  }
  
  #####Posterior distribution for component 3#####
  posterior3 <- function(pred3,param3,z,y,x1,x2,x3,x4){
    l3=likelihood3(pred3,param3,z,y,x1,x2,x3,x4)
    p3=prior3(param3,y,x1,x2,x3,x4,z)
    sumposterior3= l3 + p3
    #cat("l3=",l3,",","p3=",p3,"\n")
    return (sumposterior3)
    
  }
  
  
  ########Proposal function for the 3 components################
  
  proposalfunction1 <- function(param1){
    return(rnorm(5,mean = param1, sd= c(0.07,0.005,0.005,0.005,0.01)))
  }
  
  proposalfunction2 <- function(param2){
    return(rnorm(5,mean = param2, sd= c(0.07,0.005,0.005,0.05,0.01)))
  }
  
  proposalfunction3 <- function(param3){
    return(rnorm(5,mean = param3, sd= c(0.7,0.2,0.07,0.07,0.07)))
  }
  
  
  
  #############To calculate the Mixture model for the z update###########
  #############To calculate the Mixture model for the z update###########
  pz_123=function(z,sum_neighbours,y,pred1,pred2,pred3,chain,chain1,chain2,chain3,x1,x2,x3,x4,theta, iterations){
    for (iter in 1:iterations) {
      y11=y[z==1]
      y22=y[z==2]
      y33=y[z==3]
      
      a=exp(chain[iter]*sum_neighbours[z==1])*ifelse(y11==0,(theta+(1-theta))*(exp(-pred1(chain1[iter,],z,x1,x2,x3,x4))),(dpois(y11,exp(pred1(chain1[iter,],z,x1,x2,x3,x4)),log=F)*(1-theta)))
      b=exp(chain[iter]*sum_neighbours[z==2])*dpois(y22,exp(pred2(chain2[iter,],z,x1,x2,x3,x4)),log=F)
      c=exp(chain[iter]*sum_neighbours[z==3])*dpois(y33,exp(pred3(chain3[iter,],z,x1,x2,x3,x4)),log=F)
      
      psum1=t(sapply(1:NCOL(y11),function(i) a))
      psum2=t(sapply(1:NCOL(y22),function(i) b))
      psum3=t(sapply(1:NCOL(y33),function(i) c))
      prob_all=matrix(c(psum1,psum2,psum3),nrow=N,ncol=N)
      prob_sum=t(sapply(1:NROW(z), function(i) prob_all[i,][order(z[i,])]))
      return(log(prob_sum+1))
    }
  }
  
  ##########proposal density component 1############
  proposaldensity1 <- function(param1){
    a1 = param1[1]
    b1 = param1[2]
    c1 = param1[3]
    d1 = param1[4]
    e1 = param1[5]
    
    
    aprodensity1 = dgamma(a1,1, 0.1, log=T)
    bprodensity1 = dnorm(b1, mean = 0.1, sd = 0.01, log=T)
    cprodensity1 = dnorm(c1, 0.2, 0.01,log = T)
    dprodensity1 = dnorm(d1, mean=0, sd=0.01,log = T)
    eprodensity1 = dnorm(e1, mean = 0.1, sd = 0.01, log = T)
    
    return(aprodensity1+bprodensity1 +cprodensity1+dprodensity1+eprodensity1)
  }
  
  ########proposal density component 2###############
  proposaldensity2 <- function(param2){
    a2 = param2[1]
    b2 = param2[2]
    c2 = param2[3]
    d2 = param2[4]
    e2 = param2[5]
    
    
    aprodensity2 = dnorm(a2, 3,0.2, log=T)
    bprodensity2 = dnorm(b2, mean = 2, sd = 0.2, log=T)
    cprodensity2 = dnorm(c2, 0.4, 0.1, log = T)
    dprodensity2 = dnorm(d2, mean=1, sd=1, log = T)
    eprodensity2 = dnorm(e2, mean = 1, sd = 0.1, log = T)
    
    
    return(aprodensity2+bprodensity2+cprodensity2+dprodensity2+eprodensity2)
  }
  
  
  ##proposal density component 3#########
  proposaldensity3 <- function(param3){
    a3 = param3[1]
    b3 = param3[2]
    c3 = param3[3]
    d3 = param3[4]
    e3 = param3[5]
    
    
    aprodensity3 = dnorm(a3, mean = 70, sd = 0.2, log=T)
    bprodensity3 = dnorm(b3, mean = 2, sd = 0.5, log=T)
    cprodensity3 = dnorm(c3, mean = 0.8, sd = 0.5, log=T)
    dprodensity3 = dnorm(d3, mean = 1, sd = 0.2, log=T)
    eprodensity3 = dnorm(e3, mean = 2, sd = 0.10, log = T)
    
    return(aprodensity3+bprodensity3+ cprodensity3 + dprodensity3 + eprodensity3)
  }
  
  #########likelihood to generate the simulated data for the ABC step.###########
  likelihood_gamma = function(x,pair_neighbours_DA_x1) {
    a=exp(x*(pair_neighbours_DA_x1))/sum(exp(x*(pair_neighbours_DA_x1)))
    potts_DA<-sapply(x,function(i) a)
    
    return(matrix(potts_DA,N,N))
  }
  
  #########The initial interaction parameter is simulated from the Uniform distribution.############
  proposalfunction = function(){
    rbeta(1,10,5)
  }
  
  
  chain1 = array(dim = c(iterations+1,5))
  chain2 = array(dim = c(iterations+1,5))
  chain3 = array(dim = c(iterations+1,5))
  
  chain1[1,] = startvalue_betas1
  chain2[1,] = startvalue_betas2
  chain3[1,] = startvalue_betas3
  chain = rep(0,iterations)
  chain[1] = gamma_prior
  
  theta=rep(0,iterations)
  theta[1]<-thetap              #extra zero probability parameter in ZIP
  
  n1=rep(0,iterations)
  n2=rep(0,iterations)
  n0=rep(0,iterations)
  
  
  z=list()
  z[[1]]<-matrix(0,N,N)
  
  theta=rep(0,iterations)
  theta[1]<-thetap              #extra zero probability parameter in ZIP
  n1=list()
  n0=list()
  #cat("y=",y , "\n")
  for(i in 1:N){
    for(j in 1:N){
      if(y[i,j]==0){
        z[[1]][i,j]<-1
      }else{
        z[[1]][i,j]<-sample(1:3,1,replace=T)
      }
    }
  }
  #z[[1]]=matrix(sample(1:3,N*N,replace = T),N,N)
  
  for (iter in 1:iterations){
    #cat("iter=", iter, "\n")
    z[[iter+1]]=matrix(0,N,N)
    Pros1<-matrix(0,N,N)
    Pros2<-matrix(0,N,N)
    
    for(i in 1:N){
      for(j in 1:N){
        if(y[i,j]==0){
          Pros1[i,j]<-1
        }else{
          Pros1[i,j]<-sample(1:3,1,replace=T)
        }
      }
    }
    
    
    for(i in 1:N){
      for(j in 1:N){
        if(y[i,j]==0){
          Pros2[i,j]<-1
        }else{
          Pros2[i,j]<-sample(1:3,1,replace=T)
        }
      }
    }
    
    Pros1=matrix(sample(1:3,N*N,replace = T),N,N)
    Pros2=matrix(sample(1:3,N*N,replace = T),N,N)
    
    for (i in 1:N) {
      for (j in 1:N) {
        if(y[i,j]>0){
          while (z[[iter]][i,j]==Pros1[i,j]) {
            Pros1[i,j]=sample(1:3,1)
          }
          while (Pros2[i,j]==z[[iter]][i,j] | Pros2[i,j]==Pros1[i,j]) {
            Pros2[i,j]=sample(1:3,1)
          }
        }
      }
    }
    
    sum1=Neighbours(z[[iter]],N) ######sum of the initial z
    sum2=Neighbours_proposed(z[[iter]],N,Pros1) ######sum of the first proposed z
    sum3=Neighbours_proposed(z[[iter]],N,Pros2) ######sum of the second proposed z
    
    P_initials=pz_123(z[[iter]],sum1,y,pred1,pred2,pred3,chain,chain1,chain2,chain3,x1,x2,x3,x4,theta,iterations) ##### full Conditional using Initial z
    P_proposed1=pz_123(Pros1,sum2,y,pred1,pred2,pred3,chain,chain1,chain2,chain3,x1,x2,x3,x4,theta,iterations) ##### full Conditional using first proposed z
    P_proposed2=pz_123(Pros2,sum3,y,pred1,pred2,pred3,chain,chain1,chain2,chain3,x1,x2,x3,x4,theta,iterations) ##### full Conditional using second proposed z
    
    #cat("P_initials=",P_initials,",","P_proposed1",P_proposed1,",","P_proposed2",P_proposed2,"\n")
    log_P_initials = log(P_initials+1)
    log_P_proposed1 = log(P_proposed1+1)
    log_P_proposed2 = log(P_proposed2+1)
    
    psum=log_P_initials + log_P_proposed1 + log_P_proposed2 ######sum of the full conditional
    
    
    probab1=exp(log_P_initials)/exp(psum) #####Conditional probability using the initial z
    probab2=exp(log_P_proposed1)/exp(psum) #####Conditional probability using the first proposed z
    probab3=exp(log_P_proposed2)/exp(psum) #####Conditional probability using the second proposed z
    
    U=matrix(replicate(N,(runif(N))),N,N)
    
    for (i in 1:N) {
      for (j in 1:N) {
        #cat("probab1[i,j]=",probab1[i,j],",","probab2[i,j]",probab2[i,j],",","probab3[i,j]",probab3[i,j],"\n")
        
        
        if(y[i,j]==0){
          z[[iter+1]][i,j]<-1
        }else{
          z[[iter+1]][i,j]<-sample(1:3,1,prob = c(probab1[i,j],probab2[i,j],probab3[i,j]))
        }
        
      }
    }
    
    #    cat("z[[1]]=",z[[1]],"\n")
    
    ##########update of the Betas###########
    proposal1 = proposalfunction1(chain1[iter,])
    proposal2 = proposalfunction2(chain2[iter,])
    proposal3 = proposalfunction3(chain3[iter,])
    a1=posterior1(pred1,proposal1,z[[iter]],y,x1,x2,x3,x4,theta)
    b1=posterior1(pred1,chain1[iter,],z[[iter]],y,x1,x2,x3,x4,theta)
    a2=posterior2(pred2,proposal2,z[[iter]],y,x1,x2,x3,x4)
    b2=posterior2(pred2,chain2[iter,],z[[iter]],y,x1,x2,x3,x4)
    a3=posterior3(pred3,proposal3,z[[iter]],y,x1,x2,x3,x4)
    b3=posterior3(pred3,chain3[iter,],z[[iter]],y,x1,x2,x3,x4)
    probab11 = (a1 - b1-proposaldensity1(chain1[iter,])+proposaldensity1(proposal1))
    #cat("a1=",a1,",","b1=",b1,"\n")
    #cat("a1=",a1,",","b1=",b1,",","probab1=",probab1,"\n")
    probab22 = (a2 -b2 -proposaldensity2(chain2[iter,])+proposaldensity2(proposal2))
    #cat("a2=",a2,",","b2=",b2,",","probab2=",probab2,"\n")
    probab33 = (a3 - b3-proposaldensity3(chain3[iter,])+proposaldensity3(proposal3))
    #cat("a3=",a3,",","b3=",b3,",","probab3=",probab3,"\n")
    
    
    if (log(runif(1)) < probab11){
      chain1[iter+1,] = proposal1
    }else{
      chain1[iter+1,] = chain1[iter,]
    }
    if (log(runif(1)) < probab22){
      chain2[iter+1,] = proposal2
    }else{
      chain2[iter+1,] = chain2[iter,]
    }
    if (log(runif(1)) < probab33){
      chain3[iter+1,] = proposal3
    }else{
      chain3[iter+1,] = chain3[iter,]
    }
    
    
    #######update the theta in ZIP######
    n1[[iter]]     <- sum(z[[iter]]==1)     
    n0[[iter]]     <- length(y)-n1[[iter]]
    theta[iter+1]  <- rbeta(1,n1[[iter]],n0[[iter]])
    
    
    ###########update of gamma##########
    #simulated data using the initial interaction parameter
    x_data =likelihood_gamma(gamma_prior,Neighbours(z[[iter]],N))
    dist_data1= abs(y - x_data)  #difference between the original data and the simulated data.
    dist_data=as.vector(dist_data1)
    epsilon = 2
    #epsilon = round(quantile(dist_data,0.01))  
    proposal = proposalfunction()
    for (n in 1:(N*N)) {
      if(dist_data[n]<as.numeric(epsilon)){
        chain[iter+1] = proposal
      }else{
        chain[iter+1] = chain[iter]
      }
    }
  }
  return(list(chain1,chain2,chain3,chain,z,theta))
  
}


		     
N=100
sim_data=sample_data
y_dat=sim_data$interaction
x111=log(abs(sim_data$end.j.-sim_data$start.i.)+1)
x222=log(sim_data$GC+1)
x333=log(sim_data$Tes+1)
x444=log(sim_data$Acc+1)


#####convert into a matrix		     
]as=split(x111, ceiling(seq_along(x111)/(N*N)))
ab=split(x222, ceiling(seq_along(x222)/(N*N)))
ac=split(x333, ceiling(seq_along(x333)/(N*N)))
ad=split(x444, ceiling(seq_along(x444)/(N*N)))
yy=split(y_dat, ceiling(seq_along(y_dat)/(N*N)))

ass=expand(as,N)
abb=expand(ab,N)
acc=expand(ac,N)
add=expand(ad,N)
y1=expand(yy,N)

a=ass[c(1)]
b=abb[c(1)]
c=acc[c(1)]
d=add[c(1)]
y_sim1=y1[c(1)]

startvalue_betas1 = c(0,0,0,0,0)
startvalue_betas2 = c(12,-1,0,0,0)
startvalue_betas3 = c(1,1,0,1,1)
startvalue_betas1=as.data.frame(startvalue_betas1)
startvalue_betas2=as.data.frame(startvalue_betas2)
startvalue_betas3=as.data.frame(startvalue_betas3)

		     #####run the chain###########
chain_betas = mcmapply(run_metropolis_MCMC_betas,N,startvalue_betas1,startvalue_betas2,startvalue_betas3, gamma_prior, 20000, a, b,c,d,y_sim1, thetap=0.5, SIMPLIFY = F, mc.cores=1)




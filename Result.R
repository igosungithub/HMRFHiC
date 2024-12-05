sample_data=as.data.frame(sample_data)

burnin=(iterations/2)

intercept_1=mean(chain_betas1[[1]][["chains"]][[1]][burnin:iterations,1])
distance_1=mean(chain_betas1[[1]][["chains"]][[1]][burnin:iterations,2])
intercept_2=mean(chain_betas1[[1]][["chains"]][[2]][burnin:iterations,1])
distance_2=mean(chain_betas1[[1]][["chains"]][[2]][burnin:iterations,2])
intercept_3=mean(chain_betas1[[1]][["chains"]][[3]][burnin:iterations,1])
distance_3=mean(chain_betas1[[1]][["chains"]][[3]][burnin:iterations,2])
GC_1=mean(chain_betas1[[1]][["chains"]][[1]][burnin:iterations,3])
GC_2=mean(chain_betas1[[1]][["chains"]][[2]][burnin:iterations,3])
GC_3=mean(chain_betas1[[1]][["chains"]][[3]][burnin:iterations,3])
TES_1=mean(chain_betas1[[1]][["chains"]][[1]][burnin:iterations,4])
TES_2=mean(chain_betas1[[1]][["chains"]][[2]][burnin:iterations,4])
TES_3=mean(chain_betas1[[1]][["chains"]][[3]][burnin:iterations,4])
ACC_1=mean(chain_betas1[[1]][["chains"]][[1]][burnin:iterations,5])
ACC_2=mean(chain_betas1[[1]][["chains"]][[2]][burnin:iterations,5])
ACC_3=mean(chain_betas1[[1]][["chains"]][[3]][burnin:iterations,5])
theta=mean(chain_betas1[[1]][["theta"]][burnin:iterations])


lambda_1<-((intercept_1) +(distance_1)*log(x111+1) + (GC_1)*log(x222+1) + (TES_1)*log(x333+1) + (ACC_1)*log(x444+1))
lambda_2<-((intercept_2) +(distance_2)*log(x111+1) + (GC_2)*log(x222+1) + (TES_2)*log(x333+1) + (ACC_2)*log(x444+1))
lambda_3<-((intercept_3) +(distance_3)*log(x111+1) + (GC_3)*log(x222+1) + (TES_3)*log(x333+1) + (ACC_3)*log(x444+1))

zip_prob <- function(y, lambda, theta) {
  ifelse(
    y == 0,
    log(theta + (1 - theta) * exp(-lambda)),  # Zero inflation for y = 0
    log(1 - theta) + dpois(y, lambda = lambda, log = TRUE)  # Poisson for y > 0
  )
}


prob1=zip_prob(y_dat,exp(lambda_1),theta)
prob2=dpois(y_dat,exp(lambda_2),log=T)
prob3=dpois(y_dat,exp(lambda_3),log=T)

prob1=exp(prob1)
prob2=exp(prob2)
prob3=exp(prob3)

prob11=prob1/(prob1+prob2+prob3)
prob22=prob2/(prob1+prob2+prob3)
prob33=prob3/(prob1+prob2+prob3)


sample_data$prob1<-as.vector(prob11)
sample_data$prob2<-as.vector(prob22)
sample_data$prob3<-as.vector(prob33)


#####the probability can be changed either to a lower value or higher value to determine the certainty of the interacting pair 
sample_data$Z1<-rep(0,nrow(sample_data))
sample_data$Z1[sample_data$prob1>=0.5]<-1
sample_data$Z2<-rep(0,nrow(sample_data))
sample_data$Z2[sample_data$prob2>=0.5]<-1
sample_data$Z3<-rep(0,nrow(sample_data))
sample_data$Z3[sample_data$prob3>=0.5]<-1

mysample=cbind(sample_data$start.i., sample_data$end.j., sample_data$distance, sample_data$interaction)
mysample=as.data.frame(mysample)

burnin=10000
iterations=20000

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


lambda_1<-((intercept_1) +(distance_1)*(x111) + (GC_1)*(x222) + (TES_1)*(x333) + (ACC_1)*(x444))
lambda_2<-((intercept_2) +(distance_2)*(x111) + (GC_2)*(x222) + (TES_2)*(x333) + (ACC_2)*(x444))
lambda_3<-((intercept_3) +(distance_3)*(x111) + (GC_3)*(x222) + (TES_3)*(x333) + (ACC_3)*(x444))

prob1=dpois(y_dat,exp(lambda_1),log=T)
prob2=dpois(y_dat,exp(lambda_2),log=T)
prob3=dpois(y_dat,exp(lambda_3),log=T)


prob11=prob1/(prob1+prob2+prob3)
prob22=prob2/(prob1+prob2+prob3)
prob33=prob3/(prob1+prob2+prob3)


mysample$prob1<-as.vector(prob11)
mysample$prob2<-as.vector(prob22)
mysample$prob3<-as.vector(prob33)


#####the probability can be changed either to a lower value or higher value to determine the certainty of the interacting pair 
mysample$Z1<-rep(0,nrow(mysample))
mysample$Z1[mysample$prob1>=0.5]<-1
mysample$Z2<-rep(0,nrow(mysample))
mysample$Z2[mysample$prob2>=0.5]<-1
mysample$Z3<-rep(0,nrow(mysample))
mysample$Z3[mysample$prob3>=0.5]<-1

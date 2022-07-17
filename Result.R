
myresult=cbind(sim_data$start.i., sim_data$end.j., sim_data$interaction)
myresult=as.data.frame(myresult)
names(myresult)=c("start.i.","end.j.","interaction")

intercept_1=mean(chain_betas[[1]][[1]][burnin:iterations,1][])
distance_1=mean(chain_betas[[1]][[1]][burnin:iterations,2])
intercept_2=mean(chain_betas[[1]][[2]][burnin:iterations,1])
distance_2=mean(chain_betas[[1]][[2]][burnin:iterations,2])
intercept_3=mean(chain_betas[[1]][[3]][burnin:iterations,1])
distance_3=mean(chain_betas[[1]][[3]][burnin:iterations,2])
GC_1=mean(chain_betas[[1]][[1]][burnin:iterations,3])
GC_2=mean(chain_betas[[1]][[2]][burnin:iterations,3])
GC_3=mean(chain_betas[[1]][[3]][burnin:iterations,3])
TES_1=mean(chain_betas[[1]][[1]][burnin:iterations,4])
TES_2=mean(chain_betas[[1]][[2]][burnin:iterations,4])
TES_3=mean(chain_betas[[1]][[3]][burnin:iterations,4])
ACC_1=mean(chain_betas[[1]][[1]][burnin:iterations,5])
ACC_2=mean(chain_betas[[1]][[2]][burnin:iterations,5])
ACC_3=mean(chain_betas[[1]][[3]][burnin:iterations,5])


y_datt=sim_data$interaction

x1111=log(abs(sim_data$end.j.-sim_data$start.i.)+1)

x2222=log(sim_data$GC+1)

x3333=log(sim_data$Tes+1)

x4444=log(sim_data$Acc+1)


lambda_1<-exp(intercept_1 +(distance_1*x111) + (GC_1*x222) + (TES_1*x333) + (ACC_1*x444))
lambda_2<-exp(intercept_2 +(distance_2*x111) + (GC_2*x222) + (TES_2*x333) + (ACC_2*x444))
lambda_3<-exp(intercept_3 +(distance_3*x111) + (GC_3*x222) + (TES_3*x333) + (ACC_3*x444))


prob1=dpois(y_dat,lambda_1)
prob2=dpois(y_dat,lambda_2)
prob3=dpois(y_dat,lambda_3)


prob11=prob1/(prob1+prob2+prob3)
prob22=prob2/(prob1+prob2+prob3)
prob33=prob3/(prob1+prob2+prob3)

myresult$prob2<-as.vector(prob22)
myresult$prob2[is.na(myresult$prob2)]<-0
myresult$prob3<-as.vector(prob33)
myresult$prob3[is.na(myresult$prob3)]<-0

myresult$Z2<-rep(0,nrow(myresult))
myresult$Z2[myresult$prob2>=0.9]<-1  #####the probability can be changed either to a lower value or higher value to determine the certainty of the interacting pair 

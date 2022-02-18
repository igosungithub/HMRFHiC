#########The burnin and iteration is user-specified###########
burnin<-10000
iterations<-20000


########noise component############
par(mfrow=c(1,5))

        plot(chain_betas1[[1]][[1]][,1], type="l", ylab = "Estimated values", xlab = "iterations", main = "Intercept (noise)", las=1)
        plot(chain_betas1[[1]][[1]][,2], type="l", ylab = "Estimated values", xlab = "iterations", main = "Distance (noise)", las=1)
        plot(chain_betas1[[1]][[1]][,3], type="l", ylab = "Estimated values", xlab = "iterations", main = "GC (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[1]][,4], type="l", ylab = "Estimated values", xlab = "iterations", main = "Tes (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[1]][,5], type="l", ylab = "Estimated values", xlab = "iterations", main = "Acc (noise)",col="red",, las=1)


########signal component############
par(mfrow=c(1,5))

        plot(chain_betas1[[1]][[2]][,1], type="l", ylab = "Estimated values", xlab = "iterations", main = "Intercept (noise)", las=1)
        plot(chain_betas1[[1]][[2]][,2], type="l", ylab = "Estimated values", xlab = "iterations", main = "Distance (noise)", las=1)
        plot(chain_betas1[[1]][[2]][,3], type="l", ylab = "Estimated values", xlab = "iterations", main = "GC (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[2]][,4], type="l", ylab = "Estimated values", xlab = "iterations", main = "Tes (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[2]][,5], type="l", ylab = "Estimated values", xlab = "iterations", main = "Acc (noise)",col="red",, las=1)
       

########false signal component############
par(mfrow=c(1,5))

        plot(chain_betas1[[1]][[3]][,1], type="l", ylab = "Estimated values", xlab = "iterations", main = "Intercept (noise)", las=1)
        plot(chain_betas1[[1]][[3]][,2], type="l", ylab = "Estimated values", xlab = "iterations", main = "Distance (noise)", las=1)
        plot(chain_betas1[[1]][[3]][,3], type="l", ylab = "Estimated values", xlab = "iterations", main = "GC (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[3]][,4], type="l", ylab = "Estimated values", xlab = "iterations", main = "Tes (noise)",col="red", las=1)
        plot(chain_betas1[[1]][[3]][,5], type="l", ylab = "Estimated values", xlab = "iterations", main = "Acc (noise)",col="red",, las=1)
        

######### Note that present version of the code may require users input, e.g certain parts of the priors, the proposals,etc., but new updated version in the next couple of weeks will require very minimal input from users########### 
library(parallel)
#2D Neighbours for the sufficient statistic in the potts model
# Combined function for 2D Neighbours in the Potts model
Neighbours_combined <- function(potts_data, N, proposed_value = NULL){
  is_proposed <- !is.null(proposed_value)
  
  # Function to create shifted matrices
  shift_matrix <- function(data, shift_direction, N) {
    # Initialize the shifted matrix with the same dimension as data
    shifted = matrix(0, nrow = N, ncol = N)
    
    if (shift_direction == "down") {
      shifted[2:N, ] = data[1:(N-1), ]
    } else if (shift_direction == "up") {
      shifted[1:(N-1), ] = data[2:N, ]
    } else if (shift_direction == "right") {
      shifted[, 2:N] = data[, 1:(N-1)]
    } else if (shift_direction == "left") {
      shifted[, 1:(N-1)] = data[, 2:N]
    }
    
    return(shifted)
  }
  
  
  # Initialize matrices
  mydata1 = shift_matrix(potts_data, if (is_proposed) "up" else "down", N)
  compare_matrix1 = if (is_proposed) proposed_value else mydata1
  
  # Compute neighbour relationships
  compute_neighbours <- function(data1, compare_matrix, N) {
    comparison_result = (data1 == compare_matrix)
    result_indices = which(comparison_result, arr.ind = TRUE)
    result_matrix = matrix(0, N, N)
    result_matrix[result_indices] <- 1
    return(result_matrix)
  }
  
  neighbour1 = compute_neighbours(mydata1, compare_matrix1, N)
  neighbour2 = compute_neighbours(shift_matrix(potts_data, if (is_proposed) "down" else "up", N), if (is_proposed) proposed_value else neighbour1, N)
  neighbour3 = compute_neighbours(shift_matrix(potts_data, if (is_proposed) "left" else "right", N), if (is_proposed) proposed_value else neighbour1, N)
  neighbour4 = compute_neighbours(shift_matrix(potts_data, if (is_proposed) "right" else "left", N), if (is_proposed) proposed_value else neighbour1, N)
  
  # Calculating the total neighbours
  Neighbours_total = neighbour1 + neighbour2 + neighbour3 + neighbour4
  return(Neighbours_total)
}

#lambda equation
pred_combined <- function(params, z, x_vars, component,N){
  # Extract parameters
  a = params[1]
  b = params[2]
  c = params[3]
  d = params[4]
  e = params[5]
  
  # Create a logical matrix for subsetting
  logical_mask = z == component
  
  # Subset each matrix in x_vars based on the component
  x1_sub1 = x_vars[[1]][[1]]
  x2_sub1 = x_vars[[2]][[1]]
  x3_sub1 = x_vars[[3]][[1]]
  x4_sub1 = x_vars[[4]][[1]]
  
  x1_sub = x1_sub1[logical_mask]
  x2_sub = x2_sub1[logical_mask]
  x3_sub = x3_sub1[logical_mask]
  x4_sub = x4_sub1[logical_mask]
  
  pred=NULL
  # Compute the pred based on the component
  if (component == 1) {
    pred = a + b * x1_sub + c * x2_sub + d * x3_sub + e * x4_sub
  } else if (component == 2) {
    pred = a + b * x1_sub + c * x2_sub + d * x3_sub + e * x4_sub
  } else {
    pred = a + b * x1_sub + c * x2_sub + d * x3_sub + e * x4_sub
  }
  return(pred)
}


# Combined likelihood function for different components
likelihood_combined <- function(pred_combined, params, z, y, x_vars, component, theta,N){
  # Subset the data based on the component
  yc = y[z == component]
  # Calculate the likelihood based on the component
  if (component == 1) {
    singlelikelihoods = ifelse(yc == 0, 
                               ((theta + (1 - theta)) * (exp(-pred_combined(params, z, x_vars, component, N)))), 
                               ((1 - theta) * dpois(yc, lambda=exp(pred_combined(params, z, x_vars, component, N)), log=F)))
    s = singlelikelihoods+1
    sumll = sum(log(s))
    
  } else {
    singlelikelihoods = dpois(yc, lambda = exp(pred_combined(params, z, x_vars,component,N)), log = TRUE)
    #print(singlelikelihoods)
    sumll = sum(singlelikelihoods)
  }
  return(sumll)
}

# Combined prior function for different components
prior_combined <- function(params, component, y, x_vars, z) {
    # Extract parameters
    a = params[1]
    b = params[2]
    c = params[3]
    d = params[4]
    e = params[5]
    
    # Define mean and standard deviation values for each component
    if (component == 1) {
        meany = 3
        meanx1 = 8
        meanx2 = 0.3
        meanx3 = 1
        meanx4 = 4
        sdy = 0.1
        sdx1 = 0.1
        sdx2 = 0.1
        sdx3 = 0.2
        sdx4 = 0.5
    } else if (component == 3) {
        meany = 200
        meanx1 = 10
        meanx2 = 0.3
        meanx3 = 1.6
        meanx4 = 4.5
        sdy = 100
        sdx1 = 0.6
        sdx2 = 0.1
        sdx3 = 0.2
        sdx4 = 0.2
        
    } else if (component == 2) {  # Corrected this line
        meany = 800
        meanx1 = 8
        meanx2 = 0.3
        meanx3 = 1
        meanx4 = 4
        sdy = 0.1
        sdx1 = 0.1
        sdx2 = 0.1
        sdx3 = 0.2
        sdx4 = 0.2
    }
    
    # Calculate priors for each parameter
    a_prior = dnorm(a, meany, sdy, log = TRUE)
    b_prior = dnorm(b, meanx1, sdx1, log = TRUE)
    c_prior = dnorm(c, meanx2, sdx2, log = TRUE)
    d_prior = dnorm(d, meanx3, sdx3, log = TRUE)
    e_prior = dnorm(e, meanx4, sdx4, log = TRUE)
    
    # Return the sum of log priors
    return(a_prior + b_prior + c_prior + d_prior + e_prior)
}


  

# Combined posterior function for different components
posterior_combined <- function(pred_combined, params, z, y, x_vars, component, theta,N){
  likelihood = likelihood_combined(pred_combined, params, z, y, x_vars, component, theta,N)
  prior = prior_combined(params, component, y, x_vars, z)
  #print(likelihood)
  return(likelihood + prior)
}

# Combined proposal function for different components
proposal_function_combined <- function(params, component){
    # Define standard deviations for each component
    sd_values = list(
        component1 = c(0.007, 0.05, 0.005, 0.05, 0.01),
        component2 = c(7, 0.5, 0.5, 0.5, 0.1),
        component3 = c(0.7, 0.02, 0.07, 0.07, 0.07)
    )
    
    # Select the appropriate standard deviations
    selected_sd = sd_values[[paste0("component", component)]]
    
    # Generate and return the new proposal
    new_proposal = rnorm(5, mean = params, sd = selected_sd)
    return(new_proposal)
}


# Combined proposal density function for different components
proposaldensity_combined <- function(params, component){
  # Define mean and standard deviation values for each component
  densities = list(
    component1 = list(means = c(5, 0.15, 0.2, 0, 0.1), sds = c(0.1, 0.01, 0.01, 0.01, 0.1)),
    component2 = list(means = c(300, 2, 0.4, 0.5, 1), sds = c(0.5, 0.7, 0.1, 0.9, 0.1)),
    component3 = list(means = c(700, 2, 0.8, 1, 2), sds = c(0.2, 0.5, 0.5, 0.2, 0.1))
  )
  
  # Select the appropriate mean and standard deviation values
  selected_densities = densities[[paste0("component", component)]]
  
  # Calculate the proposal density for each parameter
  proposaldensity = mapply(dnorm, params, selected_densities$means, selected_densities$sds, MoreArgs = list(log = TRUE))
  
  # Sum and return the log proposal densities
  sum_proposaldensity = sum(proposaldensity)
  return(sum_proposaldensity)
}



likelihood_gamma = function(x,pair_neighbours_DA_x1) {
  a=exp(x*(pair_neighbours_DA_x1))/sum(exp(x*(pair_neighbours_DA_x1)))
  potts_DA<-sapply(x,function(i) a)
  
  return(matrix(potts_DA,N,N))
}


# Proposal function for interaction parameter
proposalfunction <- function() {
  rbeta(1, 10, 5)  # Assuming component 1's settings are relevant here
}

# Prior value for the interaction parameter
gamma_prior_value <- function() {
  rbeta(1, 10, 5)
}


expand <- function(ss, N) {
  lapply(ss, function(x) {
    length_diff = max(N^2 - length(x), 0)
    matrix(c(x, rep(0, length_diff)), N, N)
  })
}





pz_123 <- function(z, sum_neighbours, y, pred_combined,chains, chain_gamma, x_vars, theta, N, iter) {
  # Initialize a matrix to store probabilities for all z values
  prob_sum <- matrix(0, N, N)
  
  # Loop over components
  for (comp in 1:3) {
    # Extract y values for the current component based on z
    y_comp <- y[z == comp]
    
    # Extract the parameters for the current iteration and component
    params_comp_current_iter <- chain_gamma[iter]
    
    # Here, ensure sum_neighbours is correctly used for each component
    # Assuming sum_neighbours is a matrix where each entry corresponds to the sum of neighbours for each position
    sum_neighbours_comp <- sum_neighbours[z == comp]
    chain1=chains[[comp]][iter,]
    
    # Calculate component-specific probabilities
    if (comp == 1) {
      # For component 1, calculate probabilities with special handling for y_comp == 0
      pp = ifelse(y_comp == 0, 
                  (theta + (1 - theta))* (exp(-pred_combined(chain1, z, x_vars, comp, N))), 
                  (1 - theta) * dpois(y_comp, lambda = exp(pred_combined(chain1, z, x_vars, comp, N)), log = F))
      
      prob_comp <- exp(params_comp_current_iter * sum_neighbours_comp) * pp
    } else {
      #print(chain1)
      # For components 2 and 3, directly calculate log-probabilities
      ps = dpois(y_comp, lambda = exp(pred_combined(chain1, z, x_vars, comp, N)), log = F)
      prob_comp <- exp(params_comp_current_iter * sum_neighbours_comp) * ps
    }
    # Assign calculated log-probabilities back to the appropriate locations in prob_sum
    # This line ensures that each cell in prob_sum corresponding to 'comp' gets updated
    prob_sum[z == comp] <-(prob_comp)  # Convert back from log-probabilities to probabilities for aggregation
  }
  # After loop completion, return the log of the sum of prob_sum to avoid log(0)
  return(log(prob_sum+1))
}

		     

run_metropolis_MCMC_betas <- function(N, gamma_prior, iterations, x_vars, y, theta_start) {
  # Check and convert y and x_vars to lists if they are single matrices
  #y<-y_sim1[[1]]
  #print(dim(y))
  chains = list(matrix(NA, nrow = (iterations + 1), ncol = 5), 
                matrix(NA, nrow = (iterations + 1), ncol = 5), 
                matrix(NA, nrow = (iterations + 1), ncol = 5))
  
  chains[[1]][1, ] = runif(5)
  chains[[2]][1, ] = runif(5)
  chains[[3]][1, ] = runif(5)
  
  chain_gamma = rep(NA, iterations)
  chain_gamma[1] = gamma_prior
  
  theta = rep(NA, iterations)
  theta[1] = theta_start
  z = list(matrix(0, N, N))
  # Initialize z based on y
  for (i in 1:N) {
    for (j in 1:N) {
      z[[1]][i, j] = if (y[i, j] == 0) 1 else sample(1:3, 1, replace = TRUE)
    }
  }
  
  
  for (iter in 1:iterations) {
    print(paste("Iteration:", iter))
    # Update z
    
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
    
    Pros1 = matrix(sample(1:3, N*N, replace = TRUE), N, N)
    Pros2 = matrix(sample(1:3, N*N, replace = TRUE), N, N)
    
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
    
    sum1 = Neighbours_combined(z[[iter]], N)
    sum2 = Neighbours_combined(z[[iter]], N, Pros1)
    sum3 = Neighbours_combined(z[[iter]], N, Pros2)
    
    
    P_initials = pz_123(z[[iter]], sum1, y, pred_combined,chains, chain_gamma, x_vars, theta[iter],N,iter)
    P_proposed1 = pz_123(Pros1, sum2, y, pred_combined, chains, chain_gamma, x_vars, theta[iter],N, iter)
    P_proposed2 = pz_123(Pros2, sum3, y, pred_combined, chains, chain_gamma, x_vars, theta[iter],N, iter)
    
    
    log_P_initials =  (P_initials )
    log_P_proposed1 = (P_proposed1)
    log_P_proposed2 = (P_proposed2)
    
    psum = (log_P_initials + log_P_proposed1 + log_P_proposed2)
    
    
    probab1 = exp(log_P_initials) /  exp(psum)
    probab2 = exp(log_P_proposed1) / exp(psum)
    probab3 = exp(log_P_proposed2) / exp(psum)
    
    
    for (i in 1:N) {
      for (j in 1:N) {
        z[[iter + 1]][i, j] <- ifelse(y[i, j] == 0, 1, sample(1:3, 1, prob = c(probab1[i,j],probab2[i,j],probab3[i,j])))
      }
    }
    
    
    # Update chains
    for (comp in 1:3) {  
      proposal = proposal_function_combined(chains[[comp]][iter, ], comp)
      
      posterior_current = posterior_combined(pred_combined, chains[[comp]][iter, ], z[[iter + 1]], y, x_vars, comp, theta[iter],N)
      posterior_proposal = posterior_combined(pred_combined, proposal, z[[iter + 1]], y, x_vars, comp, theta[iter],N)
      
      proposaldensity_current = proposaldensity_combined(chains[[comp]][iter, ], comp)
      proposaldensity_proposal = proposaldensity_combined(proposal, comp)
      
      #print(posterior_proposal)
      #print(posterior_current)
      #print(proposaldensity_current)
      #print(proposaldensity_proposal)
      
      probab = (posterior_proposal - posterior_current - proposaldensity_current + proposaldensity_proposal)
      
      
      #print(probab)
      
      # Accept or reject the proposal based on its specific probab value
      if (runif(1) < probab) {
        chains[[comp]][iter + 1, ] = proposal
      } else {
        chains[[comp]][iter + 1, ] = chains[[comp]][iter, ]
      }
    }
    
    
    
    # Update theta and chain_gamma
    update_theta <- function(z_current, y) {
      n1 <- sum(z_current == 1)
      n0 <- length(y) - n1
      return(rbeta(1, n1+1, n0 + 1))
    }
    
    update_gamma <- function(gamma_current, y, z_current, N, proposalfunction) {
      # Generate synthetic data using Poisson likelihood with current gamma parameter
      x_data <- rpois(length(y), lambda = likelihood_gamma(gamma_current, Neighbours_combined(z_current, N))) # Simulate data from Poisson distribution
      
      # Calculate mean for both observed and synthetic data
      mean_y = mean(y)
      mean_x_data = mean(x_data)
      
      # Use mean difference as distance metric for Poisson data
      dist_data <- abs(mean_y - mean_x_data)
      
      proposal <- proposalfunction() # Generate a proposal
      epsilon=quantile(c(mean_y,mean_x_data),0.01)
      
      if (epsilon < dist_data) {
        return(proposal)
      } else {
        return(gamma_current)
      }
    }
    
    
    theta[iter + 1] = update_theta(z[[iter]], y)
    chain_gamma[iter + 1] = update_gamma(chain_gamma[iter], y, z[[iter]], N, gamma_prior)
  }
  
  return(list(
    chains=chains,
    gamma=chain_gamma,
    theta=theta
  ))
}

		     

#########load the data######
#sample_data <- read.csv("~/Download/sample_data.csv")

sim_data=sample_data

y_dat=sim_data$interaction

x111=log(abs(sim_data$end.j.-sim_data$start.i.)+1)

x222=log(sim_data$GC+1)

x333=log(sim_data$Tes+1)

x444=log(sim_data$Acc+1)


N=floor(sqrt(nrow(sim_data)))      ##### that is the square-root


######covert the data to a symmetric matrix##########
as=split(x111, ceiling(seq_along(x111)/(N*N)))
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

thetap=0.6

gamma_prior=rbeta(1,10,5)
iterations=20000

 #####run the chain###########
chain_betas1 = mcmapply(
  FUN = run_metropolis_MCMC_betas,
  N = rep(N, 1),  # Assuming 'times' is the number of times you want to run the function
  gamma_prior = rep(gamma_prior, 1),
  iterations = rep(iterations, 1),
  MoreArgs = list(
    x_vars = list(x1 = a, x2 = b, x3 = c, x4 = d),
    y = y_sim1[[1]],
    theta_start = thetap
  ),
  SIMPLIFY = FALSE,
  mc.cores = 1  # Or set to the desired number of cores
)
		    

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
pred_combined <- function(params, z, x_vars, component, N) {
  # Extract parameters
  a = params[1]
  b = params[2]
  c = params[3]
  d = params[4]
  e = params[5]
  
  # Create logical mask and subset matrices
  logical_mask = z == component
  x1_sub = x_vars[[1]][[1]][logical_mask]
  x2_sub = x_vars[[2]][[1]][logical_mask]
  x3_sub = x_vars[[3]][[1]][logical_mask]
  x4_sub = x_vars[[4]][[1]][logical_mask]
  
  # Compute pred for any component
  pred = a + b * x1_sub + c * x2_sub + d * x3_sub + e * x4_sub
  return(pred)
}


# Combined likelihood function for different components
likelihood_combined <- function(pred_combined, params, z, y, x_vars, component, theta,N){
  # Subset the data based on the component
  yc = y[z == component]
  # Calculate the likelihood based on the component
  if (component == 1) {
    singlelikelihoods = ifelse(yc == 0, 
                               log(theta + (1 - theta) * exp(-pred_combined(params, z, x_vars, component, N))), 
                               log(1 - theta) + dpois(yc, lambda=exp(pred_combined(params, z, x_vars, component, N)), log=T))
    #s = singlelikelihoods
    sumll = sum(singlelikelihoods)
    
  } else {
    singlelikelihoods = dpois(yc, lambda = exp(pred_combined(params, z, x_vars,component,N)), log = TRUE)
    #print(singlelikelihoods)
    sumll = sum(singlelikelihoods)
  }
  return(sumll)
}
					       
# Combined prior function for different components

prior_combined <- function(params, component, y, x_vars, z, use_data_priors = TRUE, user_fixed_priors) {
  # Extract parameters from the 'params' vector
  a = params[1]
  b = params[2]
  c = params[3]
  d = params[4]
  e = params[5]
  
  # Subset data based on component if using data-driven priors
  if (use_data_priors) {
    logical_mask = z == component
    x1 = x_vars[[1]][[1]][logical_mask]
    x2 = x_vars[[2]][[1]][logical_mask]
    x3 = x_vars[[3]][[1]][logical_mask]
    x4 = x_vars[[4]][[1]][logical_mask]
    yc = y[logical_mask]
    
    # Data-driven priors: Calculate means and standard deviations from the data
    # Different settings for meany based on the component
    if (component == 1) {
      inversesdy = rgamma(1,10,100)
      sdy = 1 / inversesdy
      meany = rnorm(1, 5, sdy/200)  # Use mean of 2 for component 1
     
    } else if (component == 2) {
      inversesdy = rgamma(1,500,2000000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 700, sdy/10)  # Use mean of 1 for other components
      
    } else{
      inversesdy = rgamma(1,500,100)
      sdy = 1 / inversesdy
      meany = rnorm(1, 500, sdy/100)
     
    }
    
    inversesdx1 = rgamma(1, (length(x1)-1)/2, sum((x1 - mean(x1))^2)/2)
    sdx1 = 1 / inversesdx1
    
    meanx1 = rnorm(1, mean(x1), sdx1)
    
    inversesdx2 = rgamma(1, (length(x2)-1)/2, sum((x2 - mean(x2))^2)/2)
    sdx2 = 1 / inversesdx2
    
    meanx2 = rnorm(1, mean(x2), sdx2)
    
    inversesdx3 = rgamma(1, (length(x3)-1)/2, sum((x3 - mean(x3))^2)/2)
    sdx3 = 1 / inversesdx3
    
    meanx3 = rnorm(1, mean(x3), sdx3)
    
    inversesdx4 = rgamma(1, (length(x4)-1)/2, sum((x4 - mean(x4))^2)/2)
    sdx4 = 1 / inversesdx4
    
    meanx4 = rnorm(1, mean(x4), sdx4)
  } else {
    # Fixed priors: Use user-provided fixed priors for the component
    if (component == 1) {
      priors <- user_fixed_priors$component1
    } else if (component == 2) {
      priors <- user_fixed_priors$component2
    } else if (component == 3) {
      priors <- user_fixed_priors$component3
    } else {
      stop("Invalid component specified!")
    }
    
    # Assign user-provided priors
    meany <- priors$meany
    meanx1 <- priors$meanx1
    meanx2 <- priors$meanx2
    meanx3 <- priors$meanx3
    meanx4 <- priors$meanx4
    sdy <- priors$sdy
    sdx1 <- priors$sdx1
    sdx2 <- priors$sdx2
    sdx3 <- priors$sdx3
    sdx4 <- priors$sdx4
  }
  
  # Calculate the log priors for each parameter based on the means and standard deviations
  a_prior = dnorm(a, meany, sdy, log = TRUE)
  b_prior = dnorm(b, meanx1, sdx1, log = TRUE)
  c_prior = dnorm(c, meanx2, sdx2, log = TRUE)
  d_prior = dnorm(d, meanx3, sdx3, log = TRUE)
  e_prior = dnorm(e, meanx4, sdx4, log = TRUE)
  
  # Return the sum of log priors for all parameters
  return(a_prior + b_prior + c_prior + d_prior + e_prior)
}


  
# Combined posterior function for different components
posterior_combined <- function(pred_combined, params, z, y, x_vars, component, theta,N,use_data_priors, user_fixed_priors){
  likelihood = likelihood_combined(pred_combined, params, z, y, x_vars, component, theta,N)
  prior = prior_combined(params, component, y, x_vars, z,use_data_priors, user_fixed_priors)
  #print(likelihood)
  return(likelihood + prior)
}




# Combined proposal function for different components
proposal_function_combined <- function(params, component) {
  # Define standard deviations for each component
  sd_values = list(
    component1 = c(0.007, 0.05, 0.005, 0.05, 0.01),
    component2 = c(7, 0.5, 0.5, 0.5, 0.1),
    component3 = c(0.7, 0.02, 0.07, 0.07, 0.07)
  )
  
  # Select the appropriate standard deviations
  selected_sd = sd_values[[paste0("component", component)]]
  
  # Ensure the length of params matches the length of selected_sd
  if (length(params) != length(selected_sd)) {
    stop("The length of params does not match the expected length for component ", component, ".")
  }
  
  # Ensure standard deviations are valid (not negative or zero)
  epsilon = 1e-6  # Small positive value to ensure positivity
  selected_sd = pmax(selected_sd, epsilon)
  
  # Generate and return the new proposal
  new_proposal = rnorm(length(params), mean = params, sd = selected_sd)
  return(new_proposal)
}



# Combined proposal density function for different components
proposaldensity_combined <- function(params, component) {
  # Define mean and standard deviation values for each component
  densities = list(
    component1 = list(means = c(5, 1, 2, 0, 1), sds = c(1000, 1000, 1000, 1000, 1000)),
    component2 = list(means = c(300, 2, 4, 5, 1), sds = c(5000, 7000, 1000, 9000, 1000)),
    component3 = list(means = c(700, 2, 8, 1, 2), sds = c(2000, 5000, 5000, 2000, 1000))
  )
  
  # Select the appropriate mean and standard deviation values
  selected_densities = densities[[paste0("component", component)]]
  means = selected_densities$means
  sds = selected_densities$sds
  
  # Ensure the length of params matches the length of means and sds
  if (length(params) != length(means) || length(params) != length(sds)) {
    stop("The length of params does not match the expected length for component ", component, ".")
  }
  
  # Ensure standard deviations are valid (not negative or zero)
  epsilon = 1e-6  # Small positive value to ensure positivity
  sds = pmax(sds, epsilon)
  
  # Calculate the proposal density for each parameter
  proposaldensity = mapply(dnorm, params, means, sds, MoreArgs = list(log = TRUE))
  
  # Sum and return the log proposal densities
  sum_proposaldensity = sum(proposaldensity)
  return(sum_proposaldensity)
}

likelihood_gamma <- function(x, pair_neighbours_DA_x1, N) {
  # Step 1: Calculate max value for numerical stability
  max_val <- max(x * pair_neighbours_DA_x1)
  
  # Step 2: Calculate adjusted exponents
  exponent_diff <- x * pair_neighbours_DA_x1 - max_val
  exp_values <- exp(exponent_diff)
  
  # Step 3: Sum the adjusted exponentiated values
  sum_exp_values <- sum(exp_values)
  
  # Step 4: Calculate normalized probabilities
  a <- exp_values / sum_exp_values
  
  # Step 5: Create Potts DA matrix using the normalized probabilities
  potts_DA <- matrix(a, nrow = N, ncol = N)
  
  return(potts_DA)
}



# Proposal function for interaction parameter
proposalfunction <- function() {
  rbeta(1, 10, 5)  # Assuming component 1's settings are relevant here
}

# Prior value for the interaction parameter
gamma_prior_value <- function() {
  rbeta(1, 10, 5)
}


pz_123 <- function(z, sum_neighbours, y, pred_combined, chains, chain_gamma, x_vars, theta, N, iter) {
  # Initialize a matrix to store probabilities for all z values
  prob_sum <- matrix(0, nrow = N, ncol = N)
  
  # Define a small value for clipping to avoid log(0)
  epsilon <- 1e-6
  
  # Loop over components
  for (comp in 1:3) {
    # Extract y values for the current component based on z
    y_comp <- y[z == comp]
    
    # Extract the parameters for the current iteration and component
    params_comp_current_iter <- chain_gamma[iter]
    
    # Ensure sum_neighbours is correctly used for each component
    sum_neighbours_comp <- sum_neighbours[z == comp]
    chain1 <- chains[[comp]][iter, ]
    
    # Calculate component-specific probabilities
    if (comp == 1) {
      # For component 1, calculate probabilities with special handling for y_comp == 0
      pp <- ifelse(y_comp == 0, 
                   (theta + (1 - theta) * exp(-pred_combined(chain1, z, x_vars, comp, N))), 
                   (1 - theta) * dpois(y_comp, lambda = exp(pred_combined(chain1, z, x_vars, comp, N)), log = FALSE))
      
      prob_comp <- exp(params_comp_current_iter * sum_neighbours_comp) * pp
    } else {
      # For components 2 and 3, directly calculate log-probabilities
      ps <- dpois(y_comp, lambda = exp(pred_combined(chain1, z, x_vars, comp, N)), log = FALSE)
      prob_comp <- exp(params_comp_current_iter * sum_neighbours_comp) * ps
    }
    
    # Ensure that dimensions are consistent during assignment
    if (length(prob_comp) == length(prob_sum[z == comp])) {
      prob_sum[z == comp] <- prob_comp
    } else {
      stop("Dimension mismatch: prob_comp and prob_sum[z == comp] have different lengths.")
    }
  }
  
  # After loop completion, clip small values in prob_sum to avoid log(0)
  prob_sum <- pmax(prob_sum, epsilon)
  
  # Return the log of the probabilities
  return(log(prob_sum))
}



run_metropolis_MCMC_betas <- function(N, gamma_prior, iterations, x_vars, y, theta_start,use_data_priors, user_fixed_priors,epsilon=NULL, epsilon_quantile = 0.01, distance_metric = "manhattan", adaptive_epsilon = TRUE) {
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
  
  # Initialize epsilon history
  epsilon_history <- numeric()
  
  
  theta = rep(NA, iterations)
  theta[1] = theta_start
  z = list(matrix(0, N, N))
  # Initialize z based on y
  for (i in 1:N) {
    for (j in 1:N) {
      z[[1]][i, j] = if (y[i, j] == 0) 1 else sample(1:3, 1, replace = TRUE)
    }
  }
  
  # Initialize acceptance counters and sd_values
  acceptance_counts = list(0,0,0)
  sd_values = list(
    component1 = c(0.7, 0.5, 0.5, 0.5, 0.1),
    component2 = c(0.7, 0.5, 0.5, 0.5, 0.1),
    component3 = c(0.7, 0.2, 0.7, 0.7, 0.7)
  )
  
  target_acceptance_rate = 0.3
  adaptation_interval = 50
  adaptation_scaling = 1.1
  
  
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
      # Current proposal function with adaptive sd_values
      proposal = rnorm(5, mean = chains[[comp]][iter, ], sd = sd_values[[paste0("component", comp)]])
      
      posterior_current = posterior_combined(pred_combined, chains[[comp]][iter, ], z[[iter + 1]], y, x_vars, comp, theta[iter],N,use_data_priors, user_fixed_priors)
      posterior_proposal = posterior_combined(pred_combined, proposal, z[[iter + 1]], y, x_vars, comp, theta[iter],N,use_data_priors, user_fixed_priors)
      
      proposaldensity_current = proposaldensity_combined(chains[[comp]][iter, ], comp)
      proposaldensity_proposal = proposaldensity_combined(proposal, comp)

      
      probab = posterior_proposal - posterior_current + proposaldensity_current - proposaldensity_proposal
      
      
      # Accept or reject the proposal based on its specific probab value
      if (log(runif(1)) < probab) {
        chains[[comp]][iter + 1, ] = proposal
        acceptance_counts[[comp]] = acceptance_counts[[comp]] + 1
      } else {
        chains[[comp]][iter + 1, ] = chains[[comp]][iter, ]
      }
    }
    
    
    # Adaptive tuning of sd_values
    if (iter %% adaptation_interval == 0 && iter <= (iterations / 2)) { # Only adapt during burn-in phase
      for (comp in 1:3) {
        acceptance_rate = acceptance_counts[[comp]] / adaptation_interval
        
        # Adjust sd_values based on acceptance rate
        if (acceptance_rate < target_acceptance_rate) {
          # Increase sd to encourage larger proposals
          sd_values[[paste0("component", comp)]] = sd_values[[paste0("component", comp)]] / adaptation_scaling
        } else if (acceptance_rate > target_acceptance_rate) {
          # Decrease sd to encourage more acceptance
          sd_values[[paste0("component", comp)]] = sd_values[[paste0("component", comp)]] * adaptation_scaling
        }
        
        # Reset acceptance count for the next adaptation interval
        acceptance_counts[[comp]] = 0
      }
    }
    
    
    # Update theta and chain_gamma
    update_theta <- function(z_current, y) {
      n1 <- sum(z_current == 1)
      n0 <- length(y) - n1
      return(rbeta(1, n1+1, n0 + 1))
    }
    
  
    update_gamma <- function(pred_combined, gamma_current, y, z_current, N, priorfunction, proposalfunction, x_vars, params, component, epsilon = NULL, epsilon_quantile = 0.01, distance_metric = "manhattan", adaptive_epsilon = TRUE, epsilon_history = NULL, iter) {
      
      # Step 1: Simulate a new value of gamma from the proposal distribution
      gamma_proposed <- proposalfunction()  # Generate a proposal for gamma
      
      # Step 2: Compute the Poisson rate (lambda) using the pred_combined function
      # Ensure that pred_combined returns a value with the correct length based on the component
      pred_values <- pred_combined(params, z_current, x_vars, component, N)
      
      # If the component corresponds to a subset of z_current, ensure pred_values matches the subset size
      if (length(pred_values) != length(z_current[z_current == component])) {
        stop("Length of pred_values does not match the subset size of z_current for the given component.")
      }
      
      # Create a full prediction matrix for z_current based on the component
      pred_matrix <- matrix(0, nrow = N, ncol = N)
      pred_matrix[z_current == component] <- pred_values
      
      # Step 3: Introduce neighborhood dependencies using the Potts model
      neighbours <- Neighbours_combined(z_current, N)  # Compute neighbours for the Potts model
      potts_prob <- exp(gamma_proposed * neighbours)  # Neighborhood influence term
      
      # Modify lambda to introduce neighbor dependency
      lambda_modified <- exp(pred_matrix) * potts_prob  # Adjusted lambda
      
      # Step 4: Simulate synthetic data from the Poisson distribution
      synthetic_data <- rpois(length(y), lambda = lambda_modified)  # Simulate data with modified lambda
      
      # Step 5: Compute the distance metric between observed and synthetic data
      if (distance_metric == "euclidean") {
        distance <- sqrt(sum((y - synthetic_data)^2))  # Euclidean distance
      } else if (distance_metric == "manhattan") {
        distance <- sum(abs(y - synthetic_data))  # Manhattan distance
      } else {
        stop("Unsupported distance metric. Use 'euclidean' or 'manhattan'.")
      }
      
      # Step 6: Set the tolerance (epsilon) for ABC
      if (adaptive_epsilon) {
        # Initialize epsilon_history if it's NULL
        if (is.null(epsilon_history)) {
          epsilon_history <- numeric()  # Use numeric() instead of list() to maintain atomic structure
        }
        
        # Ensure epsilon_history is growing correctly by adding the new distance
        epsilon_history[iter] <- distance  # Assign the current distance to the current iteration index
        
        # Update epsilon history
        epsilon_history <- c(epsilon_history, distance)  # Concatenate new distance to the epsilon history
        
        # Calculate epsilon based on quantile from the history
        epsilon <- quantile(epsilon_history, epsilon_quantile)
      } else {
        # Set epsilon based on the current distance if no history
        if (is.null(epsilon)) {
          epsilon <- mean(y) * 0.1
        }
      }
      # Step 7: Accept or reject the proposed gamma based on the distance
      if (distance < epsilon) {
        return(list(gamma = gamma_proposed, epsilon_history = epsilon_history))  # Accept the proposed gamma
      } else {
        return(list(gamma = gamma_proposed, epsilon_history = epsilon_history))  # Reject and keep the current gamma
      }
    }
    
    # Update gamma using ABC
    result_gamma<-update_gamma(pred_combined,chain_gamma[iter], y, z[[iter]], N, gamma_prior_value,proposalfunction, x_vars, chains[[comp]][iter, ],comp, epsilon=NULL, epsilon_quantile = 0.01, distance_metric = "manhattan", adaptive_epsilon = TRUE, epsilon_history = NULL)
    
    theta[iter + 1] = update_theta(z[[iter]], y)
    chain_gamma[iter + 1] = result_gamma$gamma
    epsilon_history <- result_gamma$epsilon_history
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


					       

expand <- function(ss, N) {
  lapply(ss, function(x) {
    length_diff = max(N^2 - length(x), 0)
    matrix(c(x, rep(0, length_diff)), N, N)
  })
}


					       
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
#####run the chain###########
chain_betas1 = mcmapply(
  FUN = run_metropolis_MCMC_betas,
  N = rep(N, 1),  # Assuming 'times' is the number of times you want to run the function
  gamma_prior = rep(gamma_prior, 1),
  iterations = rep(iterations, 1),
  MoreArgs = list(
    x_vars = list(x1 = a, x2 = b, x3 = c, x4 = d),
    y = y_sim1[[1]],
    theta_start = thetap,
    use_data_priors=TRUE, 
    user_fixed_priors=0,
    epsilon = NULL,
    epsilon_quantile = 0.01, 
    distance_metric = "manhattan", 
    adaptive_epsilon = TRUE
  ),
  SIMPLIFY = FALSE,
  mc.cores = 1  # Or set to the desired number of cores
)



#########Example of how to specify the prior in the user_fixed_prior option############
user_fixed_priors <- list(
  component1 = list(
    meany = 3, meanx1 = 8, meanx2 = 0.3, meanx3 = 1, meanx4 = 4,
    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 2, sdx4 = 5
  ),
  component2 = list(
    meany = 800, meanx1 = 8, meanx2 = 0.3, meanx3 = 1, meanx4 = 4,
    sdy = 1, sdx1 = 1, sdx2 = 1, sdx3 = 2, sdx4 = 2
  ),
  component3 = list(
    meany = 200, meanx1 = 10, meanx2 = 0.3, meanx3 = 1.6, meanx4 = 4.5,
    sdy = 100, sdx1 = 6, sdx2 = 1, sdx3 = 2, sdx4 = 2
  )
)

		    

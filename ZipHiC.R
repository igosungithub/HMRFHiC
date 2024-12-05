######### Note that present version of the code may require users input, e.g certain parts of the priors, the proposals,etc., but new updated version in the next couple of weeks will require very minimal input from users########### 
library(RcppArmadillo)
library(parallel)
library(Rcpp)

#lambda equation (log(lambda)= B0 + B1*log(distance) +  B2*log(GC) + B3*log(TES)+ B4*log(ACC) as in equation 2 in the paper.
pred_combined <- function(params, z, x_vars, component, N) {
  # Extract parameters
  a <- params[1]
  b <- params[2]
  c <- params[3]
  d <- params[4]
  e <- params[5]
  
  # Create logical mask and subset matrices
  logical_mask <- z == component
  x1_sub <- x_vars[[1]][[1]][logical_mask]
  x2_sub <- x_vars[[2]][[1]][logical_mask]
  x3_sub <- x_vars[[3]][[1]][logical_mask]
  x4_sub <- x_vars[[4]][[1]][logical_mask]
  
  # Compute pred for any component
  pred <- a + b * log(x1_sub+1) + c * log(x2_sub+1) + d * log(x3_sub+1) + e * log(x4_sub+1)
  return(pred)
}


# Likelihood function for the ZIP distribution
likelihood_combined <- function(pred_combined, params, z, y, x_vars, component, theta, size, N) {
  # Subset the data based on the component
  yc = y[z == component]
  lambda=pred_combined(params, z, x_vars, component, N)
  # Calculate the likelihood based on the specified distribution
  if (component == 1) {
      singlelikelihoods = ifelse(yc == 0,
                                 log(theta + (1 - theta) * exp(-lambda)),
                                 log(1 - theta) + dpois(yc, lambda = exp(lambda), log = TRUE))
  } else if(component==2 || component==3) {
      singlelikelihoods = dpois(yc[component], lambda = exp(lambda)[component], log = TRUE)
     } 
  sumll = sum(singlelikelihoods)
  return(sumll)
}

					       
# Combined prior function for the sources of biases. wDefault is to use prior from the data, while users can set there priors (example given below).
prior_combined <- function(params, component, y, x_vars, z, use_data_priors = TRUE, user_fixed_priors) {
  # Extract parameters from the 'params' vector
  a = params[1]
  b = params[2]
  c = params[3]
  d = params[4]
  e = params[5]
  
  # Set a small positive value to avoid using zero or negative standard deviations
  epsilon = 1e-6
  
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
      inversesdy = rgamma(1,10,1000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 5, sdy/10)  # Use mean of 2 for component 1
      
    } else if (component == 2) {
      inversesdy = rgamma(1,500,2000000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 700, sdy/10)  # Use mean of 1 for other components
      
    } else{
      inversesdy = rgamma(1,10,10000)
      sdy = 1 / inversesdy
      meany = rnorm(1, 10, sdy/10)
      
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

  
# Combined posterior function for the ZIP distribution.
posterior_combined <- function(pred_combined, params, z, y, x_vars, component, theta, N,
                               use_data_priors, user_fixed_priors) {
  
  # Compute the likelihood
  likelihood = likelihood_combined(pred_combined, params, z, y, x_vars, component, theta, N)
                                   
  
  # Compute the prior for the parameters
  prior = prior_combined(params, component, y, x_vars, z, use_data_priors, user_fixed_priors)

    # Return only likelihood + prior when dist is not "NB" or "ZINB"
    return(likelihood + prior)
}

# Combined proposal density function for different components
proposaldensity_combined <- function(params, component) {
  # Define mean and standard deviation values for each component
  densities = list(
    component1 = list(means = c(5, 1, 2, 0, 1), sds = c(1000, 1000, 1000, 1000, 1000)),
    component2 = list(means = c(300, 2, 4, 5, 1), sds = c(5000, 7000, 1000, 9000, 1000)),
    component3 = list(means = c(700, 2, 8, 1, 2), sds = c(2000, 5000, 5000, 2000, 1000))
  )
  
  # Check if the component is valid
  component_key <- paste0("component", component)
  if (!component_key %in% names(densities)) {
    stop("Invalid component specified: ", component)
  }
  
  # Select the appropriate mean and standard deviation values
  selected_densities = densities[[component_key]]
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

#Likelihood for the Potts model (used in updating the interaction parameter (gamma) within the Potts model).
likelihood_gamma <- function(x, pair_neighbours_DA_x1, N) {
  # Validate input dimensions
  if (length(x) != length(pair_neighbours_DA_x1)) {
    stop("x and pair_neighbours_DA_x1 must have the same length.")
  }
  
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

####c++ codes####
######2D Neighbours for the sufficient statistic in the Potts model
# Combined function for 2D Neighbours in the Potts model

cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Helper function to create shifted matrices
NumericMatrix shift_matrix(const NumericMatrix& data, std::string shift_direction, int N) {
  NumericMatrix shifted(N, N);
  
  if (shift_direction == "down") {
    for (int i = 1; i < N; i++) {
      for (int j = 0; j < N; j++) {
        shifted(i, j) = data(i - 1, j);
      }
    }
  } else if (shift_direction == "up") {
    for (int i = 0; i < N - 1; i++) {
      for (int j = 0; j < N; j++) {
        shifted(i, j) = data(i + 1, j);
      }
    }
  } else if (shift_direction == "right") {
    for (int i = 0; i < N; i++) {
      for (int j = 1; j < N; j++) {
        shifted(i, j) = data(i, j - 1);
      }
    }
  } else if (shift_direction == "left") {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N - 1; j++) {
        shifted(i, j) = data(i, j + 1);
      }
    }
  }
  
  return shifted;
}

// Helper function to compute neighbours
NumericMatrix compute_neighbours(const NumericMatrix& data1, const NumericMatrix& compare_matrix, int N) {
  NumericMatrix result_matrix(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (data1(i, j) == compare_matrix(i, j)) {
        result_matrix(i, j) = 1;
      }
    }
  }
  return result_matrix;
}

// [[Rcpp::export]]
NumericMatrix Neighbours_combined(NumericMatrix potts_data, int N, Nullable<NumericMatrix> proposed_value = R_NilValue) {
  bool is_proposed = proposed_value.isNotNull();
  NumericMatrix compare_matrix1(N, N);
  
  // Determine the shifted matrices and comparison matrix
  NumericMatrix mydata1 = shift_matrix(potts_data, is_proposed ? "up" : "down", N);
  if (is_proposed) {
    compare_matrix1 = as<NumericMatrix>(proposed_value);
  } else {
    compare_matrix1 = mydata1;
  }
  
  // Compute neighbour relationships
  NumericMatrix neighbour1 = compute_neighbours(mydata1, compare_matrix1, N);
  NumericMatrix neighbour2 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "down" : "up", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);
  NumericMatrix neighbour3 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "left" : "right", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);
  NumericMatrix neighbour4 = compute_neighbours(shift_matrix(potts_data, is_proposed ? "right" : "left", N), is_proposed ? as<NumericMatrix>(proposed_value) : neighbour1, N);
  
  // Calculating the total neighbours
  NumericMatrix Neighbours_total(N, N);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      Neighbours_total(i, j) = neighbour1(i, j) + neighbour2(i, j) + neighbour3(i, j) + neighbour4(i, j);
    }
  }
  
  return Neighbours_total;
}
')


#######calculates probabilities for a probabilistic model, combines spatial interaction effects (neighbor dependencies) with the specified statistical distribution (ZIP).
cppFunction('

#include <Rcpp.h>
using namespace Rcpp;

// Helper function to get subset indices where z == comp
std::vector<std::pair<int, int>> get_indices(const NumericMatrix &z, int comp) {
  std::vector<std::pair<int, int>> indices;
  for (int i = 0; i < z.nrow(); i++) {
    for (int j = 0; j < z.ncol(); j++) {
      if (z(i, j) == comp) {
        indices.push_back(std::make_pair(i, j));
      }
    }
  }
  return indices;
}

// pz_123 function in C++
// [[Rcpp::export]]
NumericMatrix pz_123(NumericMatrix z, NumericMatrix sum_neighbours, NumericMatrix y, 
                     Function pred_combined, List chains, NumericVector gamma_chain, 
                     List x_vars, double theta, int N, int iter) {
  NumericMatrix prob_sum(N, N); // Matrix to store probabilities
  double epsilon = 1e-6;

  // Loop over components
  for (int comp = 1; comp <= 3; comp++) {
    // Get subset indices where z == comp
    std::vector<std::pair<int, int>> indices = get_indices(z, comp);
    
    // Extract y_comp and sum_neighbours_comp based on indices
    NumericVector y_comp(indices.size()), sum_neighbours_comp(indices.size());
    for (size_t k = 0; k < indices.size(); k++) {
      int i = indices[k].first;
      int j = indices[k].second;
      y_comp[k] = y(i, j);
      sum_neighbours_comp[k] = sum_neighbours(i, j);
    }
    
    // Get parameters for current component and iteration
    double params_comp_current_iter = gamma_chain[iter];
    NumericVector chain1 = chains[comp - 1];
    NumericVector pred_values = as<NumericVector>(pred_combined(chain1, z, x_vars, comp, N));
    NumericVector lambda = exp(pred_values); // Ensure lambda is a vector

    // Compute probabilities based on component
    NumericVector prob_comp(indices.size());
    if (comp == 1) {
      NumericVector pp(indices.size());
      for (size_t k = 0; k < indices.size(); k++) {
        pp[k] = (y_comp[k] == 0) 
          ? theta + (1 - theta) * exp(-pred_values[k]) 
          : (1 - theta) * R::dpois(y_comp[k], lambda[k], false);
      }
      prob_comp = exp(params_comp_current_iter * sum_neighbours_comp) * pp;
    } else {
      NumericVector ps(indices.size());
      for (size_t k = 0; k < indices.size(); k++) {
        ps[k] = R::dpois(y_comp[k], lambda[k], false);
      }
      prob_comp = exp(params_comp_current_iter * sum_neighbours_comp) * ps;
    }
    
    // Assign prob_comp to prob_sum at the correct indices
    if (prob_comp.size() == indices.size()) {
      for (size_t k = 0; k < indices.size(); k++) {
        int i = indices[k].first;
        int j = indices[k].second;
        prob_sum(i, j) = prob_comp[k];
      }
    } else {
      stop("Dimension mismatch: prob_comp and prob_sum[z == comp] have different lengths.");
    }
  }
  
  // Clip values in prob_sum to avoid log(0)
  for (int i = 0; i < prob_sum.nrow(); i++) {
    for (int j = 0; j < prob_sum.ncol(); j++) {
      prob_sum(i, j) = std::max(prob_sum(i, j), epsilon);
    }
  }
  
  // Take element-wise log of prob_sum and return it as a matrix
  for (int i = 0; i < prob_sum.nrow(); i++) {
    for (int j = 0; j < prob_sum.ncol(); j++) {
      prob_sum(i, j) = std::log(prob_sum(i, j));
    }
  }
  
  return prob_sum; // Return log of probabilities
}


')


########Function to run the MCMC chain##########
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;

// Helper function to update theta based on z_current and y
double update_theta(const NumericMatrix &z_current, const NumericMatrix &y) {
  int n1 = 0, n0 = 0;
  for (int i = 0; i < z_current.nrow(); i++) {
    for (int j = 0; j < z_current.ncol(); j++) {
      if (z_current(i, j) == 1) {
        n1++;
      } else {
        n0++;
      }
    }
  }
  return R::rbeta(n1 + 1, n0 + 1);
}

// [[Rcpp::export]]
List run_metropolis_MCMC_betas(int N, double gamma_prior, int iterations, 
                               List x_vars, NumericMatrix y, bool use_data_priors,
                               Nullable<List> user_fixed_priors = R_NilValue, 
                               double epsilon = 0.01, std::string distance_metric = "manhattan", Nullable<double> theta_start = R_NilValue) {
  
   // Map x1, x2, x3, and x4 to distance, GC, TES, and ACC
  List x1 = x_vars["distance"];
  List x2 = x_vars["GC"];
  List x3 = x_vars["TES"];
  List x4 = x_vars["ACC"];

  // Initialize chains and other variables
  List chains = List::create(NumericMatrix(iterations + 1, 5), 
                             NumericMatrix(iterations + 1, 5), 
                             NumericMatrix(iterations + 1, 5));
  // Initialize separate gamma chains
  NumericVector gamma_chain(iterations + 1, NA_REAL);
  gamma_chain[0] = gamma_prior; // Initialize with gamma_prior

  NumericVector theta(iterations + 1, NA_REAL);
  theta[0] = as<double>(theta_start);
  
  // Initialize acceptance counts and standard deviations for components 1, 2, and 3
  std::vector<int> acceptance_counts = {0, 0, 0};
  List sd_values = List::create(NumericVector::create(70, 50, 50, 50, 10), 
                                NumericVector::create(70, 50, 50, 50, 10), 
                                NumericVector::create(70, 20, 70, 70, 70));

  // Initialize z with zeros and setup function calls to R
  List z;
  z.push_back(NumericMatrix(N, N));
  Function Neighbours_combined("Neighbours_combined");
  Function pred_combined("pred_combined");
  Function pz_123("pz_123");
  Function gamma_prior_value("gamma_prior_value");
  Function posterior_combined("posterior_combined");
  Function proposaldensity_combined("proposaldensity_combined");
  Function likelihood_gamma("likelihood_gamma");
  
  // Initialize first z values based on y
  NumericMatrix z_current = as<NumericMatrix>(z[0]);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (y(i, j) == 0) {
        z_current(i, j) = 1;
      } else {
        z_current(i, j) = floor(R::runif(1, 4));
      }
    }
  }
  z[0] = z_current;

  for (int iter = 0; iter < iterations; iter++) {
    Rcout << "Iteration: " << iter + 1 << std::endl;
    
    // Initialize Pros matrices and fill them based on y and z_current
    NumericMatrix Pros1(N, N), Pros2(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        if (y(i, j) == 0) {
          Pros1(i, j) = 1;
          Pros2(i, j) = 1;
        } else {
          do {
            Pros1(i, j) = floor(R::runif(1, 4));
          } while (Pros1(i, j) == z_current(i, j));
          
          do {
            Pros2(i, j) = floor(R::runif(1, 4));
          } while (Pros2(i, j) == z_current(i, j) || Pros2(i, j) == Pros1(i, j));
        }
      }
    }
    
    // Call Neighbours_combined and pz_123 to compute probabilities
    NumericMatrix sum1 = Neighbours_combined(z_current, N);
    NumericMatrix sum2 = Neighbours_combined(z_current, N, Pros1);
    NumericMatrix sum3 = Neighbours_combined(z_current, N, Pros2);
    
    NumericMatrix P_initials = pz_123(z_current, sum1, y, pred_combined, chains, gamma_chain, x_vars, theta[iter], N, iter);
    NumericMatrix P_proposed1 = pz_123(Pros1, sum2, y, pred_combined, chains, gamma_chain, x_vars, theta[iter], N, iter);
    NumericMatrix P_proposed2 = pz_123(Pros2, sum3, y, pred_combined, chains, gamma_chain, x_vars, theta[iter], N, iter);
    
    // Compute probabilities for sampling next z
    NumericMatrix z_next(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double psum = exp(P_initials(i, j)) + exp(P_proposed1(i, j)) + exp(P_proposed2(i, j));
        if (y(i, j) == 0) {
          z_next(i, j) = 1;
        } else {
          double prob1 = exp(P_initials(i, j)) / psum;
          double prob2 = exp(P_proposed1(i, j)) / psum;
          double random_val = R::runif(0, 1);
          z_next(i, j) = (random_val < prob1) ? 1 : ((random_val < prob1 + prob2) ? 2 : 3);
        }
      }
    }
    z.push_back(z_next);
    
    // MCMC proposal
    for (int comp = 1; comp <= 3; comp++) {  // Adjusted to comp = 1, 2, 3
      NumericMatrix chain_matrix = chains[comp - 1];  // Adjust to match comp indexing (1,2,3 -> 0,1,2 in chains)
      NumericVector proposal(5);

      for (int j = 0; j < 5; j++) {
        proposal[j] = R::rnorm(chain_matrix(iter, j), 1.0);
      }
      
      double posterior_current = as<double>(posterior_combined(pred_combined, chain_matrix(iter, _), z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors));
      double posterior_proposal = as<double>(posterior_combined(pred_combined, proposal, z_next, y, x_vars, comp, theta[iter], N, use_data_priors, user_fixed_priors));
      
      double probab = posterior_proposal - posterior_current;
      if (log(R::runif(0, 1)) < probab) {
        for (int j = 0; j < 5; j++) {
          chain_matrix(iter + 1, j) = proposal[j];
        }
      } else {
        for (int j = 0; j < 5; j++) {
          chain_matrix(iter + 1, j) = chain_matrix(iter, j);
        }
      }
    }

    // Update theta and gamma using defined functions
    theta[iter + 1] = update_theta(z_next, y);
    
    // --- ABC for updating gamma ---
    // Step 1: Propose a new gamma
    double gamma_current = gamma_chain[iter];
    double gamma_proposed = as<double>(gamma_prior_value());

    NumericMatrix gamma_matrix(N, N); // Initialize an empty N x N matrix
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        gamma_matrix(i, j) = gamma_current; // Fill all elements with gamma_current
      }
    }

   // Step 2: Compute lambda using likelihood_gamma
    NumericMatrix neighbours = Neighbours_combined(z_current, N);
    NumericMatrix lambda_matrix = likelihood_gamma(gamma_matrix, neighbours, N);

    // Step 3: Simulate synthetic data
    NumericMatrix synthetic_data(N, N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        synthetic_data(i, j) = R::rpois(lambda_matrix(i, j));
      }
    }

    // Step 4: Compute distance metric (absolute difference of means)
      double mean_y = mean(y);
      double mean_y_data = mean(synthetic_data);
      //double distance = std::abs(mean_y - mean_y_data);


    // Step 5: Compute distance between observed and synthetic data
    double distance = 0.0;
    if (distance_metric == "manhattan") {
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          distance += std::abs(mean_y - mean_y_data);
        }
      }
    } else if (distance_metric == "euclidean") {
      for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
          distance += std::pow(abs(mean_y - mean_y_data), 2);
        }
      }
      distance = std::sqrt(distance);
    } else {
      stop("Unsupported distance metric. Use manhattan or euclidean.");
    }

    // Step 7: Accept or reject gamma based on distance
    
    if (distance < epsilon) {
      gamma_chain[iter + 1] = gamma_proposed;
    } else {
      gamma_chain[iter + 1] = gamma_current;
    }
  }
  
  return List::create(Named("chains") = chains, Named("gamma") = gamma_chain, Named("theta") = theta);
}

')


		     

#########load the data######
#sample_data <- read.csv("~/Download/sample_data.csv")

sim_data=sample_data

y_dat=sim_data$interaction

x111=abs(sim_data$end.j.-sim_data$start.i.)

x222=sim_data$GC

x333=sim_data$Tes

x444=sim_data$Acc


N=floor(sqrt(nrow(sim_data)))    


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
chain_betas1 = mcmapply(
  FUN = run_metropolis_MCMC_betas,
  N = rep(N, 1),  # 'times' is the number of times you want to run the function
  gamma_prior = rep(gamma_prior, 1),
  iterations = rep(iterations, 1),
  MoreArgs = list(
    x_vars = list(distance = a, GC = b, TES = c, ACC = d),
    y = scaled_data[[1]],
    use_data_priors=TRUE, 
    user_fixed_priors=NULL,
    epsilon = 0.01,
    distance_metric = "manhattan",
    theta_start = thetap
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

		    

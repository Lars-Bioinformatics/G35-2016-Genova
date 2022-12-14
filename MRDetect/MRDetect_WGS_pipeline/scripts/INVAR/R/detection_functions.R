#===================================================================================
#
# FILE: detection_functions.R
# USAGE: Contains functions for statistical detection of ctDNA using INVAR
#
# DESCRIPTION: Functions that are used in GLRT.R code for assigning a likelihood ratio (Methods)
#             Functions were written by E Fisher.
#
#===================================================================================


## Calculates the score statistic. M, R, AF, e are vector of length n (the number of sites).
# M = mut_sum, R = DP (reads), AF = tumour_AF, e = background_AF
calculate_score <- function(M, R, AF, e)
{
  g = AF*(1-e) + (1-AF)*e
  t = g - e  # Just a mid step to avoid computing this difference twice
  U = sum(M*t/e - (R-M)*t/(1-e))^2

  I = sum(R*t^2/e + R*(t*e)^2/(1-e))

  return(U/I)
}

## Calculates the score statistic. M, R, AF, e are vector of length n (the number of sites).
## M, R, AF, e, RL, are expected to be "flatened" vectors, each of length sum(R_i).
## The i*j element of the vector should be the jth read of the ith locus.
calculate_score_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1)
{
  g = AF*(1-e) + (1-AF)*e

  U = sum((RL_PROB_1/RL_PROB_0)*M*(g/e) + (R-M)*(1-g)/(1-e)) - sum(R)
  I = (((1-g)*RL_PROB_1 - (1-e)*RL_PROB_0) / (1-e)*RL_PROB_0)^2

  return(U^2/I)
}

## Estimate p using the dervied EM algorithm.
# M = mut_sum, R = DP, AF = tumour_AF, e = background_AF
estimate_p_EM <- function(M, R, AF, e, initial_p = 0.01, iterations)
{
  g = AF*(1-e) + (1-AF)*e
  p <- initial_p
  for(i in 1:iterations)
  {
    ## Expectation step
    Z_0 <- (1-g)*p/((1-g)*p + (1-e)*(1-p))
    Z_1 <- g*p/(g*p + e*(1-p))

    ## Maximization step
    p <- sum(M*Z_1 + (R-M)*Z_0) / sum(R)
  }
  return(p)
}

estimate_p_EM_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, iterations, initial_p = 0.01)
{
  g = AF*(1-e) + (1-AF)*e
  RL_prob_0 <- RL_PROB_0 # The probability that read of length RL is from nomral tissue. For all reads.
  RL_prob_1 <- RL_PROB_1 # The probability that read of length RL is from tumour tissue. For all reads.

  p <- initial_p

  for(i in 1:iterations)
  {
    ## Expectation step
    int_step_norm <- (1-g)*RL_prob_1*p
    Z_0 <- int_step_norm / (int_step_norm + (1-e)*RL_prob_0*(1-p))

    int_step_mut <- g*RL_prob_1*p
    Z_1 <- int_step_mut / (int_step_mut + e*RL_prob_0*(1-p))

    ## Maximization step
    p <- sum(M*Z_1 + (R-M)*Z_0) / sum(R)
  }

  return(p)
}

## Calculate the generalized likelihood ratio statistic for a sample.
calc_likelihood_ratio <- function(M, R, AF, e, iterations = 200)
{
  null_likelihood <- calc_log_likelihood(M, R, AF, e, p=0)
  p_mle <- estimate_p_EM(M, R, AF, e, iterations = iterations)
  alternative_likelihood <- calc_log_likelihood(M, R, AF, e, p=p_mle)

  return(list(LR.no_size = alternative_likelihood - null_likelihood, p_mle = p_mle, null_likelihood = null_likelihood, alternative_likelihood = alternative_likelihood))
}



calc_likelihood_ratio_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p = 0.01, iterations = 200)
{

  grid_search <- function(p_grid)
  {
    likelihood_grid <- sapply(p_grid, calc_log_likelihood_with_RL, M = M,
        R = R,
        AF = AF,
        e = e,
        RL = RL,
        RL_PROB_0 = RL_PROB_0,
        RL_PROB_1 = RL_PROB_1)

      p_mle <- p_grid[likelihood_grid == max(likelihood_grid)][1]

      return(p_mle)
  }

  null_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p=0)
  p_mle <<- estimate_p_EM_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, iterations, initial_p)
  alternative_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p=p_mle)

  if(is.na(p_mle)){
    print("debugging NA")
    print(paste("M,R,AF,e,RL,RL_PROB_0,RL_PROB_1,iterations,initial_p", M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, iterations, initial_p))
  }

  print(paste("first p_mle", p_mle))
  ## If we get a negative value, it means we did not converge to the MLE
  ## In that case, do a grid search, on the area left.
  if(alternative_likelihood < null_likelihood)
  { print("alternative likelihood < null_likelihood !!")
    print("MLE didn't converge; gonna do a grid search")

    # repeat estimate of p_mle with a lower initial p
    p_mle <<- estimate_p_EM_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, iterations, initial_p = 1e-5)
    print(paste("2nd p_mle, with lower initial p", p_mle))

    p_grid <<- seq(0, p_mle, length.out = 1000)
    p_mle <<- grid_search(p_grid)

    print(paste("p_mle from grid_search = ", p_mle))
    print(p_mle %in% p_grid)
    mle_index <- which(p_grid == p_mle)
    print(paste("mle_index = ", mle_index))

    lower_finer_p_index <- max(1, mle_index - 1)
    higher_finer_p_index <- min(length(p_grid), mle_index + 1)

    print(paste("p_grid" , head(p_grid, n = 50), tail(p_grid, n = 50)))
    print(paste("low = ", p_grid[lower_finer_p_index], "index = ", lower_finer_p_index))
    print(paste("high =", p_grid[higher_finer_p_index], "index = ", higher_finer_p_index))

    finer_grid <- seq(p_grid[lower_finer_p_index], p_grid[higher_finer_p_index], length.out = 1000)
    p_mle <- grid_search(finer_grid)
    alternative_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p=p_mle)
  }

  return(list(LR = alternative_likelihood - null_likelihood, p_mle = p_mle, null_likelihood = null_likelihood, alternative_likelihood = alternative_likelihood))
}

## Calculate log likelihood of a sample given p.
calc_log_likelihood <- function(M, R, AF, e, p)
{
  q = AF*(1-e)*p + (1-AF)*e*p + e*(1-p)
  l = sum(lchoose(R, M) + M*log(q) + (R-M)*log(1-q))

  return(l)
}

calc_log_likelihood_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1 ,p)
{
# This code will break if the tumour allele fraction is >1
  g = AF*(1-e) + (1-AF)*e
  L_0 <- (1-e)*RL_PROB_0*(1-p) + (1-g)*RL_PROB_1*p
  L_1 <- e*RL_PROB_0*(1-p) + (1-g)*RL_PROB_1*p
  l = sum(M*log(L_1) + (R-M)*log(L_0))

  return(l)
}

## generates a list with the probability to get  a fragment length L, given fragment length counts from a tissue
estimate_real_length_probability <- function(fragment_length, counts, bw_adjust = 0.03, min_length, max_length, error_tolerence = 10^-10)
{
  calc_probability <- function(frag_length)
  {
    weights <- counts / sum(counts)
    #print(paste("weights:", weights))
    den <- density(fragment_length, weights = weights, adjust = bw_adjust, from=min_length-0.5, to=max_length+0.5)
    den_function <- approxfun(den)

    result <- integrate(den_function, frag_length-0.5, frag_length+0.5, abs.tol = error_tolerence)

    result$value
  }

  lengths <- seq(min_length, max_length)
  probs <- sapply(lengths, calc_probability)

  return(data.frame(fragment_length = lengths, probability = probs))
}


calculate_likelihood_ratio_for_sample <- function(data, size_characterisation, min_length, max_length, use_size = T, smooth = 0.03, size_data.path.prefix, final_prefix, only_weigh_mutants = F)
{
  size.combined <- size_characterisation

  ## read length probabilites for mutant reads
  print("generating read length probabilities")
  fragment_length <- size.combined$size[size.combined$mut == 'TRUE']
  counts <- size.combined$count[size.combined$mut == 'TRUE']

  print("estimating READ length probabilities")
  print(paste("smooth = ", smooth))
  print(paste("min length = ", min_length))
  print(paste("max length = ", max_length))
  smooth <- as.numeric(smooth)

  print("estimating read length probability")
  probs_mut <- estimate_real_length_probability(fragment_length, counts, min_length = min_length, max_length = max_length, bw_adjust = smooth)

  ## read length probabilties of normal reads
  fragment_length <- size.combined$size[size.combined$mut == 'FALSE']
  counts <- size.combined$count[size.combined$mut == 'FALSE']

  print("default smooth on wild-type data set at 0.03")
  probs_normal<- estimate_real_length_probability(fragment_length, counts, min_length = min_length, max_length = max_length, bw_adjust = smooth)

  ## plot a bar plot of the different probabilties of read lengths
  print("plot barplots")
  probs_mut <- data.frame(mut = T, fragment_length = probs_mut$fragment_length, probability = probs_mut$probability)
  probs_normal <- data.frame(mut = F, fragment_length = probs_normal$fragment_length, probability = probs_normal$probability)

  #print(paste0("output_R/", final_prefix, ".smooth_", smooth ,".", size_data.path.prefix, ".size_probabilities.plot.rds"))
  #saveRDS(list(probs_mut = probs_mut, probs_normal = probs_normal), file = paste0("output_R/", final_prefix, ".smooth_", smooth ,".", size_data.path.prefix, ".size_probabilities.plot.rds"))

  M <- data$mutant
  R <- rep(1, length(M))
  AF <- data$tumour_AF
  e <- data$background_AF
  RL <- data$size

  RL_indeces <- RL - probs_mut$fragment_length[1] + 1  ## The rows in the probability table that we need

  print(table(RL_indeces))
  RL_PROB_0 <- probs_normal$probability[RL_indeces]
  RL_PROB_1 <- probs_mut$probability[RL_indeces]

  print(paste('length of RL_PROB', length(RL_PROB_0)))
  print(paste('length of M', length(M)))

  if (only_weigh_mutants == TRUE){
      print("only weighting mutant reads.")
      RL_PROB_0[M == 0] <- 0.1
      RL_PROB_1[M == 0] <- 0.1
  }

  if(use_size == TRUE){
    print("returning likelihood ratio with RL")
    output <- calc_likelihood_ratio_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p = 0.01, iterations = 200)
    return(output)
  } else{
    print("returning likelihood ratio without RL")
    output <- calc_likelihood_ratio(M, R, AF, e, iterations = 200)
    return(output)
  }
}

#### Functions for downsampling ---------
generate_size_probabilities <- function(size.combined, min_length = 60, max_length = 300, smooth = 0.25){
  ## read length probabilites for mutant reads
  print("generating read length probabilities")
  fragment_length <- size.combined$size[size.combined$mut == 'TRUE']
  counts <- size.combined$count[size.combined$mut == 'TRUE']

  print("estimating READ length probabilities")
  print(paste("smooth = ", smooth))
  print(paste("min length = ", min_length))
  print(paste("max length = ", max_length))
  smooth <- as.numeric(smooth)

  probs_mut <- estimate_real_length_probability(fragment_length, counts, min_length = min_length, max_length = max_length, bw_adjust = smooth)

  ## read length probabilties of normal reads
  fragment_length <- size.combined$size[size.combined$mut == 'FALSE']
  counts <- size.combined$count[size.combined$mut == 'FALSE']

  print("default smooth on wild-type data set at 0.03")
  probs_normal<- estimate_real_length_probability(fragment_length, counts, min_length = min_length, max_length = max_length, bw_adjust = smooth)

  ## plot a bar plot of the different probabilties of read lengths
  print("plot barplots")
  probs_mut <- data.frame(mut = T, fragment_length = probs_mut$fragment_length, probability = probs_mut$probability)
  probs_normal <- data.frame(mut = F, fragment_length = probs_normal$fragment_length, probability = probs_normal$probability)

  return(list(probs_mut, probs_normal))
}

calculate_likelihood_ratio_for_sample.fast <- function(data, use_size = T,
  final_prefix, only_weigh_mutants = T, probs_mut, probs_normal)
{
  ## N.B. fast version!

  M <- data$mutant
  R <- rep(1, length(M))
  AF <- data$tumour_AF
  e <- data$background_AF
  RL <- data$size

  RL_indeces <- RL - probs_mut$fragment_length[1] + 1  ## The rows in the probability table that we need

  print(table(RL_indeces))
  RL_PROB_0 <- probs_normal$probability[RL_indeces]
  RL_PROB_1 <- probs_mut$probability[RL_indeces]

  print(paste('length of RL_PROB', length(RL_PROB_0)))
  print(paste('length of M', length(M)))

  if (only_weigh_mutants == TRUE){
      print("only weighting mutant reads.")
      RL_PROB_0[M == 0] <- 0.1
      RL_PROB_1[M == 0] <- 0.1
  }

  if(use_size == TRUE & sum(M) > 0){
    print("returning likelihood ratio with RL")
    output <- calc_likelihood_ratio_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p = 0.01, iterations = 200)
    return(output)
  } else if (use_size == TRUE & sum(M) == 0){
    print("no mutant reads, no need to run GLRT")
    output <- list(0, 0, 0, 0)
  } else{
    print("returning likelihood ratio without RL")
    output <- calc_likelihood_ratio(M, R, AF, e, iterations = 200)
    return(output)
  }
}

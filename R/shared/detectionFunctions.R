#===================================================================================
#
# FILE: detection_functions.R
# USAGE: Contains functions for statistical detection of ctDNA using INVAR
#
# DESCRIPTION: Functions that are used in GLRT.R code for assigning a likelihood ratio (Methods)
#             Functions were written by E Fisher.
#
#===================================================================================

suppressPackageStartupMessages(require(assertthat))
suppressPackageStartupMessages(require(stringr))

## Calculates the score statistic. M, R, AF, e are vector of length n (the number of sites).
# M = mut_sum, R = DP (reads), AF = tumour_AF, e = background_AF
calculate_score <- function(M, R, AF, e)
{
    g = AF * (1 - e) + (1 - AF) * e
    t = g - e  # Just a mid step to avoid computing this difference twice
    U = sum(M * t / e - (R - M) * t / (1 - e)) ^ 2

    I = sum(R * t ^ 2 / e + R * (t * e) ^ 2 / (1 - e))

    return(U / I)
}

## Calculates the score statistic. M, R, AF, e are vector of length n (the number of sites).
## M, R, AF, e, RL, are expected to be "flatened" vectors, each of length sum(R_i).
## The i*j element of the vector should be the jth read of the ith locus.
calculate_score_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1)
{
    g = AF * (1 - e) + (1 - AF) * e

    U = sum((RL_PROB_1 / RL_PROB_0) * M * (g / e) + (R - M) * (1 - g) /
                (1 - e)) - sum(R)
    I = (((1 - g) * RL_PROB_1 - (1 - e) * RL_PROB_0) / (1 - e) * RL_PROB_0) ^
        2

    return(U ^ 2 / I)
}

## Estimate p using the derived EM algorithm.
# M = MUTANT, R = DP, AF = TUMOUR_AF, e = BACKGROUND_AF
estimate_p_EM <- function(M, R, AF, e, initial_p = 0.01, iterations = 200)
{
    g = AF * (1 - e) + (1 - AF) * e
    p <- initial_p
    for (i in 1:iterations)
    {
        ## Expectation step
        Z_0 <- (1 - g) * p / ((1 - g) * p + (1 - e) * (1 - p))
        Z_1 <- g * p / (g * p + e * (1 - p))

        ## Maximization step
        p <- sum(M * Z_1 + (R - M) * Z_0) / sum(R)
    }
    return(p)
}

estimate_p_EM_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p = 0.01, iterations = 200)
{
    g = AF * (1 - e) + (1 - AF) * e
    RL_prob_0 <- RL_PROB_0 # The probability that read of length RL is from nomral tissue. For all reads.
    RL_prob_1 <- RL_PROB_1 # The probability that read of length RL is from tumour tissue. For all reads.

    p <- initial_p

    for (i in 1:iterations)
    {
        ## Expectation step
        int_step_norm <- (1 - g) * RL_prob_1 * p
        Z_0 <-
            int_step_norm / (int_step_norm + (1 - e) * RL_prob_0 * (1 - p))

        int_step_mut <- g * RL_prob_1 * p
        Z_1 <- int_step_mut / (int_step_mut + e * RL_prob_0 * (1 - p))

        ## Maximization step
        p <- sum(M * Z_1 + (R - M) * Z_0) / sum(R)
    }

    return(p)
}

## Calculate the generalized likelihood ratio statistic for a sample.
calc_likelihood_ratio <- function(M, R, AF, e, iterations = 200)
{
    null_likelihood <- calc_log_likelihood(M, R, AF, e, p = 0)
    p_mle <- estimate_p_EM(M, R, AF, e, iterations = iterations)
    alternative_likelihood <-
        calc_log_likelihood(M, R, AF, e, p = p_mle)

    return(
        list(
            LR.no_size = alternative_likelihood - null_likelihood,
            p_mle = p_mle,
            null_likelihood = null_likelihood,
            alternative_likelihood = alternative_likelihood
        )
    )
}



calc_likelihood_ratio_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p = 0.01, iterations = 200)
{
    grid_search <- function(p_grid)
    {
        likelihood_grid <-
            sapply(p_grid, calc_log_likelihood_with_RL, M = M, R = R, AF = AF, e = e,
                   RL = RL, RL_PROB_0 = RL_PROB_0, RL_PROB_1 = RL_PROB_1)

        p_grid[likelihood_grid == max(likelihood_grid)][1]
    }

    null_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p = 0)
    p_mle <<- estimate_p_EM_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, initial_p, iterations)
    alternative_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p = p_mle)

    ## If we get a negative value, it means we did not converge to the MLE
    ## In that case, do a grid search, on the area left.
    if (alternative_likelihood < null_likelihood)
    {
        # repeat estimate of p_mle with a lower initial p
        p_mle <<- estimate_p_EM_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, 1e-5, iterations)

        p_grid <<- seq(0, p_mle, length.out = 1000)
        p_mle <<- grid_search(p_grid)

        mle_index <- which(p_grid == p_mle)

        lower_finer_p_index <- max(1, mle_index - 1)
        higher_finer_p_index <- min(length(p_grid), mle_index + 1)

        finer_grid <- seq(p_grid[lower_finer_p_index], p_grid[higher_finer_p_index], length.out = 1000)
        p_mle <- grid_search(finer_grid)
        alternative_likelihood <- calc_log_likelihood_with_RL(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1, p = p_mle)
    }

    list(
        LR = alternative_likelihood - null_likelihood,
        p_mle = p_mle,
        null_likelihood = null_likelihood,
        alternative_likelihood = alternative_likelihood
    )
}

## Calculate log likelihood of a sample given p.
calc_log_likelihood <- function(M, R, AF, e, p)
{
    q = AF * (1 - e) * p + (1 - AF) * e * p + e * (1 - p)
    sum(lchoose(R, M) + M * log(q) + (R - M) * log(1 - q))/length(R)
}

calc_log_likelihood_with_RL <- function(M, R, AF, e, RL, RL_PROB_0, RL_PROB_1 , p)
{
    # This code will break if the tumour allele fraction is >1
    g = AF * (1 - e) + (1 - AF) * e
    L_0 <- (1 - e) * RL_PROB_0 * (1 - p) + (1 - g) * RL_PROB_1 * p
    L_1 <- e * RL_PROB_0 * (1 - p) + (g) * RL_PROB_1 * p #changed g from (1-g) 16/03/22 EDitter
    sum(M * log(L_1) + (R - M) * log(L_0))/length(R)
}

## generates a list with the probability to get a fragment length L, given fragment length counts from a tissue
estimate_real_length_probability <- function(fragment_length, counts, bw_adjust = 0.03,
                                             min_length, max_length, error_tolerence = 1e-10)
{
    calc_probability <- function(frag_length)
    {
        probability <- 0

        # Need to guard against too few points.
        if (length(counts) > 1)
        {
            assert_that(sum(counts) > 0, msg = str_c("Sum of ", length(counts), " counts is zero. Cannot create weights."))

            # Weights take into account the TOTAL number of reads of all sizes
            weights <- counts / sum(counts)

            # KDE of fragment length function
            # but KDE only includes reads from min to max fragment length defined
            den <- density(fragment_length, weights = weights, adjust = bw_adjust,from = min_length - 0.5, to = max_length + 0.5)
            den_function <- approxfun(den)

            result <- integrate(den_function, frag_length - 0.5, frag_length + 0.5, abs.tol = error_tolerence)

            probability <- result$value
        }
        probability
    }

    lengths <- seq(min_length, max_length)
    probs <- sapply(lengths, calc_probability)

    data.frame(fragment_length = lengths, probability = probs)
}

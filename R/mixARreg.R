## Do not edit this file manually.
## It has been automatically generated from mixAR.org.

mixARreg <- function(x, y, mixARmodel, tol = 1e-6, niter = 200){
    if(length(x) != NROW(y)){ # 2020-01-17 @Davide: NROW(y) to allow for numeric y.
        stop("x and y must contain the same number of observations")
    } ## Need a first iteration to initialise all variables
    n <- length(x)
    reg_model <- lm(data.frame(x, y))
    reg_parameters <- beta_old <- reg_model$coefficients
    reg_residuals <- reg_model$residuals
    MAR_mod <- MAR_old <- fit_mixAR(reg_residuals, mixARmodel, fix = "shift")$model
    pk <- MAR_mod@order
    g <- length(pk)
    p <- max(pk)
    ufit <- numeric(n - p)
    for(k in 1:g){
        for(t in 1:(n - p)){
            ufit[t] <- ufit[t] + MAR_mod@prob[k] *
                (MAR_mod@arcoef@a[[k]] %*% reg_residuals[(t + p - 1):(t + p - pk[k])])
        }
    }
    xstar <- x ;
    xstar[(p+1):n] <- x[(p+1):n] - ufit
    diff <- 1
    i <- 1
    ## @Davide: not urgent but beware - 
    ##    (b) 'diff' is the sum of quantities, which are potentially on different
    ##         scales, whi may cause trouble.
    ## 24/04/2020 changed this, now stops if all single differences are smaller than tol

    while(diff > tol){  ## From iteration 2, basically repeats the code above
        i <- i+1
        reg_model <- lm(data.frame(xstar,y))
        reg_parameters <- reg_model$coefficients
        reg_residuals <- reg_model$residuals
        MAR_mod <- fit_mixAR(reg_residuals, mixARmodel, fix = "shift")$model
        ufit <- numeric(n - p)
        for(k in 1:g){
            for(t in 1:(n - p)){
                ufit[t] <- ufit[t] + MAR_mod@prob[k] *
                    (MAR_mod@arcoef@a[[k]] %*% reg_residuals[(t + p -1):(t + p - pk[k])])
            }
        }

       # diff <- sum(abs(reg_parameters / beta_old - 1)) +
    #                sum(abs(MAR_mod@prob / MAR_old@prob - 1)) +
    #                sum(abs(MAR_mod@scale / MAR_old@scale - 1)) +
    #                sum(abs(unlist(MAR_mod@arcoef@a) / unlist(MAR_old@arcoef@a) -1))
        diff <- c(abs(reg_parameters - beta_old),
            abs(MAR_mod@prob - MAR_old@prob),
            abs(MAR_mod@scale - MAR_old@scale),
            abs(unlist(MAR_mod@arcoef@a) - unlist(MAR_old@arcoef@a)))
        diff <- ifelse(all(diff < tol), 0, 1)
        beta_old <- reg_parameters
        MAR_old <- MAR_mod
        xstar <- x 
        xstar[(p+1):n] <- x[(p+1):n] - ufit
        if(i >= niter) 
            break
    }

    list(reg = reg_model, mixARmodel = MAR_mod, niter = i, convergence = diff)
}

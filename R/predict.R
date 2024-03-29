## Do not edit this file manually.
## It has been automatically generated from predict.org.

## extend_index <- function(ind, g){          # ind is a matrix with each multi-index in a column
##     if(length(g) == 1)   # g is either a vector of positive integers or a single one
##         g <- 1:g
##     do.call("cbind", lapply(g, function(x) rbind(ind, x)))
## }

predict_coef <- function(model, maxh){
    p <- max(model@order)
    g <- .nmix(model)
    
    next_tuple <- function(tuple){
        h <- length(tuple)
        if(h == 0)
            return(1L)
        else if(h < maxh)
            return(c(1, tuple))
        else{ # h == maxh here
            if(all(tuple == g))
                return(NA)  # no more tuples

            if(tuple[1] == g){       # drop leading g's
                wrkind <- which(tuple != g)[1]  # guaranteed at least one element
                tuple <- tuple[-(1:(wrkind-1))]
            }

            tuple[1] <- tuple[1] + 1
            return(tuple)
        }
    }

    gtoh <- g^((maxh-1):0)
    tuple2index <- function(tuple) sum( (tuple-1) * tail(gtoh, h) ) + 1

                                                  # ste tryabva library(sfsmisc) ako izpolzvam
    tuple_i2t <- integer(maxh)                   #                      sfsmisc::digitsBase()
    index2tuple <- function(ind){  # ind must be  >= 1 and < g^h
        for (i in 1:maxh) {
            tuple_i2t[i] <- ind %% g
            ind <- ind %/% g
            if(ind == 0)
                break
        }                       # add one to get values in 1,...,g, not 0,...,g-1
        rev(tuple_i2t[1:i] + 1) # guarrantees a non-zero value for the  most significant digit
    }


                           # b <- lapply(1:maxh, function(x) array(dim = c(p,rep(g,x))) )
                           # s <- lapply(1:maxh, function(x) array(dim = c(rep(g,x),x)) )
    b <- lapply(1:maxh, function(x) matrix(NA, nrow = p, ncol = g^x))
    s <- lapply(1:maxh, function(x) matrix(NA, nrow = x, ncol = g^x)) # dim's increase
                                                                      # with h for sigma
    probs <- lapply(1:maxh, function(x) rep(NA_real_, g^x))
    sStable <- lapply(1:maxh, function(x) rep(NA_real_, g^x))

    mar <- model@arcoef[] # take as matrix

    sigma <- model@scale
    prob <- model@prob

    bmat <- matrix(NA, nrow = maxh, ncol = p)
    bmat <- rbind(bmat, diag(rep(1,p)))

    baserow <- nrow(bmat) - p + 1  # 'baserow - h' gives the index for h in bmat or smat
    brows <- baserow + 1:p  # 'brows - h' gives the indices of the rows for matrix B or S

    smat <- matrix(NA, nrow = maxh, ncol = maxh)
    smat <- rbind(smat, matrix(0, nrow = p, ncol = maxh))

    probvec <- rep(1, maxh)

    curtuple <- 1
    while(length(curtuple) != 1  || !is.na(curtuple)){
        h <- length(curtuple)

        tuind <- tuple2index(curtuple)

        ## print(c(tuind, curtuple))
        ## # print(curtuple)

        bmat[baserow-h, ] <-  mar[curtuple[1],] %*% bmat[brows - h, ]
        b[[h]][ , tuind] <- bmat[baserow-h, ]

        swrk <- sigma[curtuple[1]]
        if(h > 1)
            swrk <- c( swrk, mar[curtuple[1],] %*% smat[brows - h, maxh - (h-2):0] )

        smat[baserow-h, maxh - (h-1):0] <- swrk
        s[[h]][ , tuind] <- swrk
        smat[baserow-h, -(maxh - (h-1):0)] <- 0 # the zeroes should be put when smat is
                                                # created; the current NA's are for testing.

        sStable[[h]][tuind] <- sqrt(sum(swrk^2))

        probvec[h] <- prob[ curtuple[1] ]
        if(h > 1)
            probvec[h] <- probvec[h] * probvec[h-1]

        probs[[h]][tuind] <- probvec[h]

        curtuple <- next_tuple(curtuple)
    }

    list(arcoefs = b, sigmas = s, probs = probs, sStable = sStable)
}

setGeneric("multiStep_dist",
           function(model, maxh, N, xcond, ...) standardGeneric("multiStep_dist"))

setMethod("multiStep_dist",
          c(model = "MixARGaussian", maxh = "numeric", N = "missing", xcond = "missing"),
          function(model, maxh){
              co <- predict_coef(model, maxh)

              f <- function(h, xcond, what){
                  arh <- t(co$arcoefs[[h]])
                  if(ncol(arh) != length(xcond)){
                      excess <- length(xcond) - ncol(arh)
                      if(excess > 0){
                          ## TODO: argument to suppress the message?
                          message(paste0("using the last ",  length(xcond),
                                         " values in 'xcond'"))
                          xcond <- xcond[-(1:excess)]
                      } else
                          stop("length(xcond) must be >= maximal AR order")
                  }

                  hmo <- new("MixARGaussian",
                             prob   = co$probs[[h]],
                             scale  = co$sStable[[h]],
                             shift  = as.numeric( arh %*% xcond ),
                             arcoef = arh
                             )
                  switch(what,
                         cdf = mix_cdf(hmo, xcond = xcond),
                         pdf = mix_pdf(hmo, xcond = xcond),
                         # default, mainly for testing
                         hmo
                         )
              }
              f
          })

setMethod("multiStep_dist",
          c(model = "MixARGaussian", maxh = "numeric", N = "missing", xcond = "ANY"),
                            # TODO: more computations can be done here when xcond is available
          function(model, maxh, ...){
              xcond <- xcond
              co <- predict_coef(model, maxh)
	      
              f <- function(h, what){
                  arh <- t(co$arcoefs[[h]])
                  if(ncol(arh) != length(xcond)){
                      excess <- length(xcond) - ncol(arh)
                      if(excess > 0){
                          ## TODO: argument to suppress the message?
                          message(paste0("using the last ",  length(xcond),
                                         " values in 'xcond'"))
                          xcond <- xcond[-(1:excess)]
                      } else
                          stop("length(xcond) must be >= maximal AR order")
                  }

                  hmo <- new("MixARGaussian",
                             prob   = co$probs[[h]],
                             scale  = co$sStable[[h]],
                             shift  = as.numeric(arh %*% xcond ),
                             arcoef =arh)
                  switch(what,
                         cdf = mix_cdf(hmo, xcond = xcond),
                         pdf = mix_pdf(hmo, xcond = xcond),
                         # default, mainly for testing
                         hmo
                         )
              }
              f
          })

                    # todo: drop unnecessary values from xcond (if there are more than needed?
setMethod("multiStep_dist",
          c(model = "MixAR", maxh = "numeric", N = "numeric", xcond = "numeric"),
          function(model, maxh, N, xcond){
              sampled <- replicate(N, mixAR_sim(model, n = maxh, init = xcond, nskip = 0))
              sampled <- apply(sampled, 1, sort)
              f <- function(h, what, ...){
                  if(is.character(what))
                      switch(what,
                             cdf = ecdf(sampled[,h]),
                             pdf = {
                                 d <- density(sampled[,h], ...)
                                 approxfun(x = d$x, y = d$y, yleft = 0, yright = 0)
                             },
                             location = mean(sampled[,h ]),
                             variance = var(sampled[,h ]),
                             sd       = sd(sampled[,h ]),
                             skewness = skewness(sampled[,h ]),
                             kurtosis = kurtosis(sampled[,h ]),
                             # default
                             if(exists(what, mode = "function"))
                                 (match.fun(what))(sampled[,h ])
                             else
                                 sampled
                             )
                  else if(is.function(what))
                      what(sampled[ , h])
                  else
                      stop("Argument 'what' is of incorrect type.")
              }
              f
          })

#
#
# TODO: unfinished!
#
# setMethod("mix_se",
#           signature(model="MixAR"),
#           function(model, orig, N, ...){
#               pmax <- max(model@order)
#
#               sampled <- replicate(N, mixAR_sim(model, n=length(orig), init=init, nskip=0))
#               sampled <- apply(sampled, 1, sort)
#           })

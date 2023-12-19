
 ## Define model to simulate data from (functions from mixAR package)

prob <- c(0.5,0.5)
sigma <- c(1, 2)
arco <- list(-0.5,1)


model1 <- new("MixARGaussian", prob=prob, scale=sigma, arcoef=arco)

y1 <- mixAR_sim(model1, 600, rep(0, max(model1@order)), nskip = 100)  # data

 ## Run the main function for 20000(10000 burnin) iterations for convergence of phi_k1, k=1,...,g
#help(bayes_mixAR)
xx <- bayes_mixAR(y1, model1, fix_shift=FALSE, tau=c(.15, .25), nsim=200, burnin=100)

 ## Note, more than one trial may be necessary to tune "tau"

 ##check labels switch
#help(label_switch)
xx_adj <- label_switch(xx$scale, m = 10)

 ## From xx_adj$true_perm, no sign of labels switch, so proceed

 ## Update "model" with highest density parameters

prob     <- c(density(xx$mix_weights[ , 1])$x[which.max(density(xx$mix_weights[,1])$y)],
              density(xx$mix_weights[ , 2])$x[which.max(density(xx$mix_weights[,2])$y)])

sigma    <- c(density(xx$scale[ , 1])$x[which.max(density(xx$scale[,1])$y)],
              density(xx$scale[ , 2])$x[which.max(density(xx$scale[,2])$y)])

arco       <- list(density(xx$ARcoeff$`Component_1`)$x[which.max(density(xx$ARcoeff$`Component_1`)$y)],
                density(xx$ARcoeff$`Component_2`)$x[which.max(density(xx$ARcoeff$`Component_2`)$y)])

shift    <- c(density(xx$shift[ , 1])$x[which.max(density(xx$shift[,1])$y)],
                  density(xx$shift[ , 2])$x[which.max(density(xx$shift[,2])$y)])

model1 <- new("MixARGaussian", prob=prob, scale=sigma, arcoef=arco, shift=shift)

 ## Now choose the best order pk
#help(Choose_pk)
pk_star <- Choose_pk(y1, model1, tau=c(.15,.25), pmax=5, method="NULL", nsim=100)


 ## Finally, calculate marginal loglik

#help(marg_loglik)
m_log <- marg_loglik(y1, model1, tau=c(.15, .25), nsim=100, pk_star$x[1, 3])

 ## in m_log, nsim=1000 means 1000 iterations at each step
 ## Total number of iterations is 1000*g (for AR coeff) + 3000 for remaining parameters
 ## In this case 5000

 ## This is only an explanatory example, sample size should be
 ## larger for good/reliable results.

#########################################

## Example for seasonal data

## How to build seasonal model

probS <- c(0.5, 0.5)
sigmaS <- c(5, 1)
ar1 <- list(c(0.5, -0.5), c(1.1, 0, -0.5))
ar12 <- list(0, c(-0.3, 0.1))
s <- 12

rag <- new("raggedCoefS", a=ar1, as=ar12, s=s)

modelS <- new("MixARGaussian", prob=probS, scale=sigmaS, arcoef=rag)

## Equivalently

model <- new("MixARGaussian", prob=probS, scale=sigmaS,
             arcoef= list(a=ar1, as=ar12, s=s))

yS <- mixAR_sim(model, n=500, init=rep(0,24))

fit_mixAR(yS, modelS)
fit_mixAR(yS, modelS, fix="shift")


######## fit_mixARreg

## Simulate covariates

xReg  <- data.frame(rnorm(500,7,1), rt(500,3),rnorm(500,3,2))
xReg2 <- seq(-10,10,length.out=500)
## Build mixAR part

probReg <- c(0.7, 0.3)
sigmaReg <- c(1, 5)
arReg <- list(c(-0.5, 0.5), 1.1)

modelReg <- new("MixARGaussian", prob=probReg, scale=sigmaReg, arcoef=arReg)

##Simulate from mixAR part
uReg <- mixAR_sim(modelReg, 500, c(0,0))

## Model yReg is:
## y = 10 + x1 + 3* x2 + 2 * x3 + e
## Model yReg2 is:
## y = 10 + x
## uReg is mixAR, same for both models

yReg <- xReg[,1] + 3 * xReg[,2] + 2 * xReg[,3] + uReg
yReg2 <- 10 + xReg2 + uReg
## Fit model

## Using MixARGaussian
fit_mixARreg(yReg, xReg, modelReg)

## Using EMinit - same output

EMinit <- list(prob=probReg, scale=sigmaReg, arcoef=arReg)
fit_mixARreg(yReg2, xReg2, EMinit = EMinit)


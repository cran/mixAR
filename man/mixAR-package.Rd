\name{mixAR-package}
\alias{mixAR-package}
%\alias{mixAR}
\docType{package}

\concept{mixture autoregressive model}
\concept{MAR model}
\concept{TMAR model}
\concept{mixture autoregression}

\concept{non-gaussian time series}
\concept{asymmetry in time series}
\concept{multi-modality in time series}

% attaching the package is not necessary for this file but it ensures
% that Sexpr's using functions from the package will work also in other
% Rd files, without them needing to attach the package. This is
% undocumented 'feature' and it doesn't help if a page is rendered on
% its own. It is best to use 'mixAR::' prefix but for generated examples
% is not very suitable.
\title{
  \Sexpr[stage=build,results=hide]{require("mixAR")}
  \Sexpr[stage=build,results=hide]{pd <- packageDescription("mixAR")}
  \Sexpr[stage=build,results=hide]{lb <- library(help="mixAR", character.only=TRUE)}
  \Sexpr[stage=build,results=hide]{lbinfo <- if(!is.null(lb$info[[2]])) lb$info[[2]] else "No_matches"}
  \Sexpr[stage=build,results=hide]{anypat <- function(string, pat, ...){ any(sapply(pat, function(x) grepl(x, string, ...))) }}
  \Sexpr[stage=build,results=hide]{
    lbsubset <- function(lbinfo, pat){
      if(length(lbinfo) == 0) return("No entries")
      ind.start <- which(sapply(lbinfo, function(x) anypat(x, pat)))
      if(length(ind.start) == 0) return("No entries")
      ind.ws <- which(grepl("^[[:space:]]+", lbinfo))
      res <- ind.start
      for(ind in ind.start){
	while(any((ind+1) - ind.ws == 0)){
	  ind <- ind + 1
	  res <- c(res, ind)
	}
      }
      lbinfo[ sort(res) ]
    }
  }
  \Sexpr[stage=build,results=hide]{lbpaste <- function(pat) paste("\\\\preformatted{",
    {wrk <- lbsubset(lbinfo,pat);
     wrk <- gsub("(^[^\\\\]?\%)|(([^\\\\])\%)", paste0("\\\\3", "\\\\\\\\", "\%"), wrk);
    paste(if(length(wrk) > 0) wrk else "No entries", collapse="\\n")}, "}", sep="")}
  \Sexpr[stage=build]{pd$Title}
}
\description{
  \Sexpr[stage=build]{pd$Description}
}
\details{
  % \tabular{ll}{
  %   Package:  \tab \Sexpr[stage=build]{pd$Package}\cr
  %   Type:     \tab \Sexpr[stage=build]{pd$Type}\cr
  %   Version:  \tab \Sexpr[stage=build]{pd$Version} \cr
  %   Date:     \tab \Sexpr[stage=build]{pd$Date}\cr
  %   License:  \tab \Sexpr[stage=build]{pd$License}\cr
  %   LazyLoad: \tab \Sexpr[stage=build]{pd$LazyLoad}\cr
  %   Built:    \tab \Sexpr[stage=build]{pd$Built}\cr
  % }
  % 
  % Index:
  %
  % \Sexpr[stage=build,results=rd]{paste("\\\\preformatted{", paste(if(!is.null(lb$info[[2]])) lb$info[[2]] else "No entries",collapse="\n"), "}", sep="")}


  %\strong{Index:}

%%%%   Below is an attempt to group related topics in clusters. This is for
%%%%   development purposes. It is easy to edit this, including adding and
%%%%   removing clusters, see the source of this help page. Note that the
%%%%   last, "Rest", contains all topics that are not included in any of the
%%%%   other clusters. This is to make visible forgotten topics. To prevent a
%%%%   topic appearing in "Rest", add it to cluster "ignore", which is just
%%%%   before "Rest" in the Rd source.
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{alltopics <- unique(gsub("^(([^ ]+)[^ ]*|(^[ ]))(.*)", "\\\\1", lb$info[[2]]))}
%%%%   \Sexpr[stage=build,results=hide]{alltopics <- alltopics[!grepl("^ ", alltopics)]}
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{pat <- list()}
%%%% 
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(mixAR =
%%%%     c("^mixAR-class", "^MixARGaussian-class", "^mixARgen-class",
%%%%     "^mixAR-methods", "^show_diff-methods", "^parameters-methods" )))}
%%%%   \subsection{mixAR classes}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["mixAR"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(estimation =
%%%%     c("^fit_mixAR", "^mixARreg", "^mixSARfit", "^mix_se-methods"
%%%%       )))}
%%%%   \subsection{mixAR estimation}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["estimation"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(likelihood =
%%%%     c("^est_templ", "^cond_loglik", "^make_fcond_lik-methods"  )))}
%%%%   \subsection{mixAR likelihood}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["likelihood"]])}
%%%%   }
%%%% 
%%%% 
%%%%   % "^mix_kurtosis", "^mix_variance-methods",
%%%%   % "^mix_moment-methods", "^mix_location-methods", "^mix_central_moment-methods",
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(mixARdist =
%%%%     c("^mix_pdf-methods", "^mix_moment", "^multiStep_dist-methods"   )))}
%%%%   \subsection{mixAR distributions}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["mixARdist"]])}
%%%%   }
%%%% 
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(dataandsim = c(
%%%%     "^exampleModels", "^mixARnoise_sim", "^mixAR_sim",
%%%%      "^randomArCoefficients")))}
%%%%   \subsection{Data and simulation}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["dataandsim"]])}
%%%%   }
%%%% 
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(mixFilter =
%%%%     c("^mixFilter", "^raghat1" )))}
%%%%   \subsection{mixFilter}{
%%%%     Please use whenever possible these functions, instead of programming
%%%%     the nested loops from scratch. This leads to more reliable code and
%%%%     easy maintenance. For example the functions may be replaced by
%%%%     C-code without ttouching any cother code.
%%%%     
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["mixFilter"]])}
%%%%   }
%%%%   
%%%% 
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(MixComp =
%%%%     c("^MixComp-class", "^mix_ncomp-methods", "^mix_ek-methods",
%%%%     "^mix_hatk-methods", "^\%of\%", "^inner", "^row_lengths"  )))}
%%%%   \subsection{MixComp}{
%%%%     Convenient computing with mixAR time series.
%%%%     
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["MixComp"]])}
%%%%   }
%%%% 
%%%%   
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(em =
%%%%     c("^em_[^ ]", "^tau2probhat", "^mixMstep", "^mixgenMstep",
%%%%     "^mixARemFixedPoint", "^mixAR_cond_probs", "^mixSubsolve" )))}
%%%%   \subsection{EM}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["em"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(parameters =
%%%%     c("^parameters", "^lik_params", "^predict_coef",
%%%%     "^tomarparambyComp", "^mixAR_switch", "^isStable",
%%%%     "^[.]can[^ ]*" )))}
%%%%   \subsection{Parameters}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["parameters"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(raggedCoef  =
%%%%   c("^raggedCoef", "^rag_modify", "^ragged2char" )))}
%%%%   \subsection{raggedCoef}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["raggedCoef"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(Bayesian =
%%%%     c("^bayes_mixAR", "^bx_dx", "^Choose_pk", "^err", "^label_switch",
%%%%     "^sampZpi", "^sampMuShift", "^sampMuShift", "^sampSigmaTau", "^marg_loglik")))}
%%%%   \subsection{Bayesian}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["Bayesian"]])}
%%%%   }
%%%% 
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(mixARnoise =
%%%%     c("^dist_norm", "^fnoise", "^get_edist[^ ]", "^noise[^ ]", "^stdnormmoment")))}
%%%%   \subsection{mixAR noise}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(pat[["mixARnoise"]])}
%%%%   }
%%%% 
%%%%   
%%%%   % \Sexpr[stage=build,results=hide]{pat <- c(pat, list(internal = c("^[.]can[^ ]*")))}
%%%%   % \subsection{Internal}{
%%%%   %   \Sexpr[stage=build,results=rd]{lbpaste(pat[["internal"]])}
%%%%   % }
%%%% 
%%%%   %% these are ignored, i.e. not printed anywhere
%%%%   \Sexpr[stage=build,results=hide]{pat <- c(pat, list(ignore =
%%%%   c("^mixAR-package", "^initialize-methods")))}
%%%%   
%%%%   \Sexpr[stage=build,results=hide]{rest <- alltopics[!grepl( paste0("(", paste0(unlist(pat),collapse=")|("), ")"),alltopics)]}
%%%%   \subsection{rest}{
%%%%     \Sexpr[stage=build,results=rd]{lbpaste(rest)}
%%%%   }


  % An overview of how to use the package, including the most important functions

  Package \pkg{mixAR} provides functions for modelling with mixture
  autoregressive (MixAR) models. The S4 class \code{"MixARGaussian"} can
  be used when the error distributions of the components are standard
  Gaussian. The class \code{"MixARgen"} admits arbitrary (well, within
  reason) distributions for the error components. Both classes inherit
  from the virtual class \code{"MixAR"}.

  Estimation can be done with \code{fit_mixAR}.  Currently, the EM
  algorithm is used for estimation.

  For \code{"MixARGaussian"} the M-step of the EM algorithm reduces to a
  system of linear equations. For \code{"MixARgen"} the problem is
  substantially non-linear. The implementation is fairly general but
  currently not optimised for efficiency. The specification of the error
  distributions went through several stages and may still be
  reviewed. However, backward compatibility will be kept.

}
\author{
  \Sexpr[stage=build]{pd$Author}

Maintainer: \Sexpr[stage=build]{pd$Maintainer}
}
\references{
  \insertRef{mary2013}{mixAR}

  \insertRef{boshnakov2009mar}{mixAR}

  \insertRef{Boshnakov2011on1st2nd}{mixAR}

  \insertRef{FongEtAl2007vecmar}{mixAR}

  \insertRef{ravagli2020bayesian}{mixAR}
  
  \insertRef{ravagli2020portfolio}{mixAR}

  \insertRef{shahadat2012}{mixAR}

  \insertRef{WongPhD}{mixAR}

  \insertRef{WongLi2000}{mixAR}

  \insertRef{WongLi2001a}{mixAR}

  \insertRef{WongLi2001b}{mixAR}
}
\seealso{
  fit several types of mixAR models:
  \code{\link{fit_mixAR}}, \code{\link{bayes_mixAR}},
  \code{\link{fit_mixARreg}},
  \code{\link{mixSARfit}};

  Predictive distributions and summaries:
  \code{\link{mix_pdf}}, \code{\link{mix_cdf}}, \code{\link{mix_qf}},
  \code{\link{mix_location}},
  \code{\link{mix_variance}},
  \code{\link{mix_central_moment}},
  \code{\link{mix_moment}},
  \code{\link{mix_kurtosis}},
  \code{\link{mix_ekurtosis}}

  multi-step prediction:
  \code{\link{multiStep_dist}}

}
\examples{
## object 'exampleModels' contains a number of models for examples and testing
names(exampleModels)
exampleModels$WL_ibm

## some of the models below are available in object 'exampleModels';
## the examples here show how to create them from scratch
mo_WLprob <- c(0.5439, 0.4176, 0.0385)              # model coefficients from Wong&Li
mo_WLsigma <- c(4.8227, 6.0082, 18.1716)
mo_WLar <- list(c(0.6792, 0.3208), c(1.6711, -0.6711), 1)

mo_WL <- new("MixARGaussian", prob = mo_WLprob, scale = mo_WLsigma, arcoef = mo_WLar)
                                        
mo_WL_A <- new("MixARGaussian"                         # WongLi, model A
              , prob = c(0.5, 0.5)
              , scale = c(5, 1)
              , shift = c(0, 0)
              , arcoef = list(c(0.5), c(1.1))
              )

mo_WL_B <- new("MixARGaussian"                         # WongLi, model B
              , prob = c(0.75, 0.25)
              , scale = c(5, 1)
              , shift = c(0, 0)
              , arcoef = list(c(0.5), c(1.4))
              )
                                        
mo_WL_I <- new("MixARGaussian"                          # WongLi, model I
              , prob = c(0.4, 0.3, 0.3)
              , scale = c(1, 1, 5)
              , shift = c(0, 0, -5)
              , arcoef = list(c(0.9, -0.6), c(-0.5), c(1.50, -0.74, 0.12))
              )

mo_WL_II <- new("MixARGaussian"                         # WongLi, model II
               , prob = c(0.4, 0.3, 0.3)
               , scale = c(1, 1, 5)
               , shift = c(5, 0, -5)
               , arcoef = list(c(0.9, -0.6), c(-0.7, 0), c( 0, 0.80))
               )

## MixAR models with arbitrary dist. of the components
## (user interface not finalized)

## Gaussian
mo_WLgen <- new("MixARgen", prob = mo_WLprob, scale = mo_WLsigma, arcoef = mo_WLar,
               dist = list(dist_norm))

## t_3
mo_WLt3v <- new("MixARgen", prob = mo_WLprob, scale = mo_WLsigma, arcoef = mo_WLar,
               dist = list(fdist_stdt(3, fixed = FALSE)))

## t_20, t_30, t_40 (can be used to start estimation)
mo_WLtf <- new("MixARgen", prob = mo_WLprob, scale = mo_WLsigma, arcoef = mo_WLar,
              dist = list(generator =
                              function(par)
                                  fn_stdt(par, fixed = FALSE), param = c(20, 30, 40)))

## data(ibmclose, package = "fma")  # for `ibmclose'

## The examples below are quick but some of them are marked as 'not run'
## to avoid cumulative time of more than 5s on CRAN.

## fit a MAR(2,2,1) model
a0a <- fit_mixAR(as.numeric(fma::ibmclose), c(2, 2, 1), crit = 1e-4)
## same with 2 sets of automatically generated initial values.
\donttest{
a0b <- fit_mixAR(as.numeric(fma::ibmclose), c(2, 2, 1), 2, crit = 1e-4)
}

## fix the shift parameters:
a1a <- fit_mixAR(as.numeric(fma::ibmclose), c(2, 2, 1), fix = "shift", crit = 1e-4)
## ... with 3 sets of automatically generated initial values.
\donttest{
a1b <- fit_mixAR(as.numeric(fma::ibmclose), c(2, 2, 1), 3, fix = "shift", crit = 1e-4)
}

\donttest{
## specify the model using a MixAR model object
a1c <- fit_mixAR(as.numeric(fma::ibmclose), a1a$model, init = a0a$model, fix = "shift",
           crit = 1e-4)

## fit a model like mo_WL using as initial values 2 automatically generated sets.
a2 <- fit_mixAR(as.numeric(fma::ibmclose), mo_WL, 2, fix = "shift", permute = TRUE, 
          crit = 1e-4)
}

moT_B3 <- new("MixARgen"
               , prob = c(0.3, 0.3, 0.4)
               , scale = c(2, 1, 0.5)
               , shift = c(5, -5, 0)
               , arcoef = list(c(0.5, 0.24), c(-0.9), c(1.5, -0.74, 0.12))
                                        # t4, t4, t10
               , dist = distlist("stdt", c(4,10), fixed = c(FALSE, TRUE), tr = c(1, 1, 2))
               )

moT_C1 <- new("MixARgen"
              , prob = c(0.3, 0.3, 0.4)
              , scale = c(2, 1, 0.5)
              , shift = c(5, -5, 0)
              , arcoef = list(c(0.5, 0.24), c(-0.9), c(1.5, -0.74, 0.12))
                                        # t4, t7, N(0,1)
              , dist = distlist(c("stdt", "stdt", "stdnorm"), c(4,7))
              )

## demonstrate reuse of existing models
exampleModels$WL_Bt_1
moT_C2 <- new("MixARgen"
              , model = exampleModels$WL_Bt_1
              , dist = distlist(c("stdt", "stdt", "stdnorm"), c(4,7))  # t4, t7, N(0,1)
              )
moT_C3 <- new("MixARGaussian", model = exampleModels$WL_Bt_1 )
}
\keyword{package}
\keyword{ts}

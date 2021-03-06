% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crosnma.model.R
\name{crosnma.model}
\alias{crosnma.model}
\title{Create JAGS model to synthesize cross-design evidence and cross-format data in NMA and NMR for dichotomous outcomes}
\usage{
crosnma.model(
  prt.data,
  std.data,
  trt,
  study,
  outcome,
  n,
  design,
  reference,
  trt.effect = "random",
  covariate = NULL,
  reg0.effect = "independent",
  regb.effect = "random",
  regw.effect = "random",
  split.regcoef = T,
  method.bias = NULL,
  bias = NULL,
  bias.type = NULL,
  bias.covariate = NULL,
  bias.effect = "common",
  prior = list(tau.trt = NULL, tau.reg0 = NULL, tau.regb = NULL, tau.regw = NULL,
    tau.gamma = NULL, pi.high.rct = NULL, pi.low.rct = NULL, pi.high.nrs = NULL,
    pi.low.nrs = NULL),
  run.nrs = list(var.infl = 1, mean.shift = 0, n.adapt = 500, n.iter = 10000, n.burnin
    = 4000, thin = 1, n.chains = 2)
)
}
\arguments{
\item{prt.data}{An object of class data.frame containing the individual participant dataset. Each row contains the data of a single participant.
The data frame needs to have the following columns: treatment, study identification, outcome (event and non-event), design. Additional columns might be required for certain analyses.}

\item{std.data}{An object of class data.frame containing the study-level dataset. Each row represents the information of study arm.
The data frame needs to have the following columns: treatment, study identification, outcome (number of events), sample size and design. Additional columns might be required for certain analysis.}

\item{trt}{A vector of length 2 with the names of the treatment variable (as character) in prt.data and std.data, respectively.}

\item{study}{A vector of length 2 with the names of the study variable (as character) in prt.data and std.data, respectively.}

\item{outcome}{A vector of length 2 with the names of the outcome variable (as character) in prt.data and std.data, respectively.}

\item{n}{A character of the name of the number of participants variable in std.data.}

\item{design}{A vector of length 2 with the names of the design variable (as character) in prt.data and std.data, respectively.}

\item{reference}{A character indicating the name of the reference treatment. This option must be specified otherwise the 'placebo' is set as a reference.}

\item{trt.effect}{A character defining the model for the study-specific treatment effects. Options are 'random' (default) or 'common'.}

\item{covariate}{An optional list with two elements, the first indicating the name of the covariates in prt.data and the second is for the corresponding ones in std.data, to conduct network meta regression
The covariates can be either numeric or dichotomous variables. The user can provide up to 3 covariates. The covariate needs to be provided for both prt.data and std.data, respectively.
The default option is covariate=NULL where no covariate adjustment is applied (network meta-analysis).}

\item{reg0.effect}{An optional character (when \code{covariate} is not NULL) indicating the relationship across studies for the prognostic effects expressed by the regression coefficient, ($\beta_{0,j}$), in a study \eqn{j}.
Options are 'independent' or 'random'. We recommend using 'independent' (default).}

\item{regb.effect}{An optional character (when \code{covariate} is not NULL) indicating the relationship across studies for the between-study regression coefficient. This parameter quantifies the treatment-mean covariate interaction.
Options are 'random' or 'common'. Default is 'random'.}

\item{regw.effect}{An optional character (when \code{covariate} is not NULL)  indicating the relationship across studies for the within-study regression coefficient. This parameter quantifies the treatment-covariate interaction effect at the individual level.
Options are 'random' and 'common'. Default is 'random'.}

\item{split.regcoef}{A logical value (when \code{covariate} is not NULL). If TRUE the within- and between-study coefficients will be splitted in the analysis of prt.data.
The default is TRUE. When the split.regcoef = FALSE, only a single regression coefficient will be estimated to represent both the between-studies and within-studies covariate effects.}

\item{method.bias}{An optional character for defining the method to combine randomised clinical trials (RCT) and non-randomised studies (NRS) (required when \code{design} has nrs in addition to rct).
Options are 'naive' for naive synthesize, 'prior' for using NRS to inform priors for the relative treatment effects in RCTs.
or 'adjust1' and 'adjust2' to allow a bias adjustment.}

\item{bias}{An optional vector of length 2 (required when method.bias='adjust1' or 'adjust2') indicating the name of the variable (as character) that includes the risk of bias adjustment in prt.data and std.data, respectively.
The entries of this variable should be a character with entries either low, high or unclear. These values need to be repeated for the participants that belong to the same study.}

\item{bias.type}{An optional character defining of bias on the treatment effect (required when method.bias='adjust1' or 'adjust2').
Three options are possible: 'add' for additive bias effect,'mult' for multiplicative bias effect and'both' for both an additive and a multiplicative term.}

\item{bias.covariate}{An optional vector of two characters (required when method.bias='adjust1' or 'adjust2'). It has the variable name of the variable that will be used in estimating the probability of bias.}

\item{bias.effect}{An optional character indicating the relationship for the bias coefficients across studies.
Options are 'random' or 'common' (default). It is required when method.bias='adjust1' or 'adjust2'.}

\item{prior}{An optional list to control the prior for various parameters in JAGS model. When effects are set as 'random', we can set the heterogeneity parameters for: tau.trt for the treatment effects,
tau.reg0 for the effect of prognostic covariates, tau.regb and tau.regw for within- and between-study covariate effect, respectively.
and tau.gamma for bias effect. The default of all heterogeneity parameters is 'dunif(0,2)'. Currently only the uniform distribution is supported.
When the method.bias='adjust1' or 'adjust2', the user may provide priors to control the bias probability.
For the bias probabilities, beta distributions are assumed with the following default values: RCT with low (pi.low.rct='dbeta(1,20)')/high (pi.high.rct='dbeta(3,1)') bias, NRS with low(pi.low.rct='dbeta(1,2)')/high (pi.high.rct='dbeta(30,1)') bias (pi.low.nrs, pi.high.nrs).}

\item{run.nrs}{An optional list is needed when the NRS used as a prior (method.bias='prior').
The list consists of the follwoing: (\code{var.infl}) controls the inflation of the varaince of NRS estimates and its values range between 0 (NRS do not contribute at all and the prior is vague) and 1 (the NRS evidence is used at face value, default approach).
The parameter (\code{mean.shift}) is the bias shift to be added/subtracted from the estimated NRS mean treatment effects (0 is the default). Either (\code{var.infl}) or (\code{mean.shift}) should be provided but not both.
Here you can also specify the arguments to control the MCMC chains with default value is in the parentheses: the number of adaptions n.adapt (500), number of iterations n.iter(10000), number of burn in n.burnin (4000),
number of thinning thin (1) and number of chains n.chains (2), see \code{\link{jags.model}} arguments from rjags package.}
}
\value{
\code{crosnma.model} returns an object of class \code{crosnmaModel} which is a list containing the following components:

\code{jagsmodel}  A long character string containing JAGS code that will be run in \code{\link{jags}}.

\code{data}  The data used in the JAGS code.

\code{trt.key}  A table of the treatments and its mapped integer number (used in JAGS model).

\code{trt.effect} A character defining the model for the study-specific treatment effects.

\code{method.bias}  A character for defining the method to combine randomised clinical trials (RCT) and non-randomised studies (NRS).

\code{covariate}  A list of the the names of the covariates in prt.data and std.data used in network meta-regression.

\code{split.regcoef} A logical value. If FALSE the within- and between-study regression coefficients will be considered equal.

\code{regb.effect} An optional character  indicating the model for the between-study regression coefficients across studies.

\code{regw.effect} An optional character indicating the model for the within-study regression coefficients across studies.

\code{bias.effect} An optional character indicating the model for the bias coefficients across studies.

\code{bias.type} A character indicating the effect of bias on the treatment effect; additive ('add') or multiplicative ('mult') or both ('both').
}
\description{
This function creates a JAGS model and the needed data. The JAGS code is created from the internal function \code{crosnma.code}.
}
\examples{
# An example from participant-level data and study-level data.
# data
data(prt.data)
data(std.data)
 #=========================#
  # Create a jags model  #
 #=========================#
 # We conduct a network meta-analysis assuming a random effect model.
 # The data comes from randomised-controlled trials and non-randomised studies. They will be combined naively.
 # The data has 2 different formats: individual participant data (prt.data) and study-level data (std.data).
mod <- crosnma.model(prt.data=prt.data,
                  std.data=std.data,
                  trt=c('trt','trt'),
                  study=c('study','study'),
                  outcome=c('outcome','outcome'),
                  n='n',
                  design=c('design','design'),
                  reference='A',
                  trt.effect='random',
                  covariate = NULL,
                  method.bias='naive'
                   )
 #=========================#
    # Fit jags model  #
 #=========================#
fit <- crosnma.run(model=mod,
                n.adapt = 20,
                n.iter=50,
                thin=1,
                n.chains=3)

 #=========================#
   # Display the output   #
 #=========================#
summary(fit)
plot(fit)

}
\seealso{
\code{\link{crosnma.run}}, \code{\link{jags.model}}
}

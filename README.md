# Introduction
This folder contains all the ```R``` code used to fit the proposed Mixture of Circular Regressions (MixCircReg) model

# Description of the code
In this section, we provide a brief description of the code (contained in the ```R``` script FitMixCircReg.R):

* The ```Beta_IRWLS(theta,x,weig,mu,lmd,initBeta)``` function performs the iterative reweighted least squares to estimate the component regression parameters as in Fisher and Lee (1992)
* The ```mix_circ_regression(theta,x1,x2,k,initBeta1,initBeta2,initmu,initlmd,initprop)``` function fits the MixCircReg model using the Expectation-Maximization (EM) algorithm (with the IRLWS for the regression parameters)
* The ```fit.mix_circ_regression(theta,x1,x2,k,initBeta1,initBeta2,initmu,initlmd,initprop)``` function fit multiple MixCircReg models using the function ```mix_circ_regression(...)``` to initialize the fitting algorithm
* The ```dir_con_mle(theta,x,weig,initBeta)``` function is used inside of the ```mix_circ_regression``` function to calculate the mean direction ($\mu$) and concentration paramater ($\kappa$) for each component.
* The ```g(x)``` function is the arctangent link function.

#### Arguments (inputs)

  + ```x``` is an $n\times (2q+p)$ matrix of linear and circular (sine-cosine) covariates.
  + ```x1``` is an $n\times q$ matrix of circular covariates.
  + ```x2```  is an $n\times p$ matrix of linear covariates.
  + ```theta``` is a vector of length $n$ that consists of the response variable values
  + ```mu``` is the mean direction parameter
  + ```lmd``` is the concentration parameter
  + ```weig``` is an $n-$dimensional vector of weights represented by the responsibilities of each component 
  + ```initBeta``` is the initial parameter value of the overall $\beta$ parameter for each component.
  + ```initBeta1``` is the initial parameter value of the $\beta$ parameter for the circular (sine-cosine) covariates.
  + ```initBeta2``` is the initial parameter value of the $\beta$ parameter for the linear covariates.
  + ```initmu``` is a $k-$dimensional initial parameter vector of $\mu$, the mean direction parameter
  + ```initlmd``` is a $k-$dimensional initial parameter vector of $\kappa$, the concentration parameter
  + ```initprop``` is a $k-$dimensional initial parameter vector of $\pi$, the mixing proportion

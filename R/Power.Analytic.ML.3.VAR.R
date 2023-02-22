###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for IL studies using analytic derivations
# Model with a continuous Level 1 predictor


Power.Analytic.ML.3.VAR = function(N,T.obs,
                             b00,b10,
                             sigma,rho,
                             sigma.v0,sigma.v1,rho.v,
                             mu.X,sigma.X,sigma.X.v0,rho.X, 
                             alpha,
                             side.test){
  
  tic()
  # Compute information matrix
  c_11 = ((1-rho^2)*sigma^2)^(-1)*(T.obs+rho*(T.obs-2)*(rho-2))
  c_12 = c_21 = 0
  c_22 = ((1-rho^2)*sigma^2)^(-1)*((1+rho^2)*(T.obs-1)*(sigma.X.v0^2+sigma.X^2/(1-rho.X^2)) - 2*rho*(T.obs-1)*(sigma.X.v0^2+rho.X*sigma.X^2/(1-rho.X^2)))
  Sigma.XX = cbind(c(c_11,c_21),c(c_12,c_22))
  Sigma.nu = cbind(c(sigma.v0^2,rho.v*sigma.v0*sigma.v1),c(rho.v*sigma.v0*sigma.v1,sigma.v1^2))
  Inf.beta = N*(Sigma.XX - Sigma.XX%*%solve(solve(Sigma.nu)+Sigma.XX)%*%Sigma.XX) 
  Cov.beta = solve(Inf.beta)
  StdError.beta = sqrt(diag(Cov.beta))
  
  # Compute power
  
  if (side.test == 1){ # One-side test: positive
  power.hat = 1-pnorm(qnorm(1-alpha)-(c(b00,b10)/StdError.beta)) 
  }
  
  if (side.test == 2){ # One-side test: negative
    power.hat = pnorm(-(qnorm(1-alpha)-(c(b00,b10)/StdError.beta))) 
  }
  
  if (side.test == 3){ # Two-tailed test
    power.hat = 1-pnorm(qnorm(1-alpha/2)-(c(b00,b10)/StdError.beta)) +  pnorm(qnorm(alpha/2)-(c(b00,b10)/StdError.beta))
  }

  names(StdError.beta) = names(power.hat) = c('b00','b10')
  
  toc(log = TRUE, quiet = TRUE)
  log.txt =  tic.log(format = TRUE)
  log.lst = tic.log(format = FALSE)
  tic.clearlog()
  timings.analytic.power = unlist(lapply(log.lst, function(x) x$toc - x$tic))
  
  return(c(StdError.beta=StdError.beta,
           power.hat=power.hat,
           timings.analytic.power=timings.analytic.power))
}
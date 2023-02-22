###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for IL studies using analytic derivations
# Model with a binary Level 2 predictor


Power.Analytic.ML.1 = function(N.0,N.1,T.obs,
                             b00,b01.Z,
                             sigma,rho,
                             sigma.v0,
                             alpha,
                             side.test){
  
  tic()
  # Compute the total number of participants  
  N.Total = N.0 + N.1 
  # Compute mean of Level 2 variable
  mu.Z = N.1/N.Total
  # Compute std. deviation of Level 2 variable
  sigma.Z = sqrt(mu.Z*(1-mu.Z))
  
  # Compute information matrix 
  Sigma.XX = cbind(c(N.Total,N.Total*mu.Z),c(N.Total*mu.Z,(N.Total-1)*sigma.Z^2+N.Total*mu.Z^2))
  scalar.num = T.obs+rho*(T.obs-2)*(rho-2)
  scalar.den = sigma^2*(1-rho^2)+sigma.v0^2*scalar.num
  Inf.beta = scalar.num*scalar.den^(-1)*Sigma.XX
  Cov.beta = solve(Inf.beta)
  StdError.beta = sqrt(diag(Cov.beta))
  
  # Compute power
  
  if (side.test == 1){ # One-side test: positive
  power.hat = 1-pnorm(qnorm(1-alpha)-(c(b00,b01.Z)/StdError.beta)) 
  }
  
  if (side.test == 2){ # One-side test: negative
    power.hat = pnorm(-(qnorm(1-alpha)-(c(b00,b01.Z)/StdError.beta))) 
  }
  
  if (side.test == 3){ # Two-tailed test
    power.hat = 1-pnorm(qnorm(1-alpha/2)-(c(b00,b01.Z)/StdError.beta)) +  pnorm(qnorm(alpha/2)-(c(b00,b01.Z)/StdError.beta))
  }

  names(StdError.beta) = names(power.hat) = c('b00','b01.Z')
  
  toc(log = TRUE, quiet = TRUE)
  log.txt =  tic.log(format = TRUE)
  log.lst = tic.log(format = FALSE)
  tic.clearlog()
  timings.analytic.power = unlist(lapply(log.lst, function(x) x$toc - x$tic))
  
  return(c(StdError.beta=StdError.beta,
           power.hat=power.hat,
           timings.analytic.power=timings.analytic.power))
}
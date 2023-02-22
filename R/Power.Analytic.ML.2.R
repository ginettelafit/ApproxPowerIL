###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for IL studies using analytic derivations
# Model with a continuous Level 2 predictor


Power.Analytic.ML.2 = function(N,T.obs,
                             b00,b01.W,
                             sigma,rho,
                             sigma.v0,
                             mu.W,sigma.W,isW.center,
                             alpha,
                             side.test){
  
  if (isW.center==TRUE){mu.W=0}
  
  tic()
  # Compute information matrix 
  Sigma.XX = cbind(c(N,N*mu.W),c(N*mu.W,(N-1)*sigma.W^2+N*mu.W^2))
  scalar.num = T.obs+rho*(T.obs-2)*(rho-2)
  scalar.den = sigma^2*(1-rho^2)+sigma.v0^2*scalar.num
  Inf.beta = scalar.num*scalar.den^(-1)*Sigma.XX
  Cov.beta = solve(Inf.beta)
  StdError.beta = sqrt(diag(Cov.beta))
  
  # Compute power
  
  if (side.test == 1){ # One-side test: positive
  power.hat = 1-pnorm(qnorm(1-alpha)-(c(b00,b01.W)/StdError.beta)) 
  }
  
  if (side.test == 2){ # One-side test: negative
    power.hat = pnorm(-(qnorm(1-alpha)-(c(b00,b01.W)/StdError.beta))) 
  }
  
  if (side.test == 3){ # Two-tailed test
    power.hat = 1-pnorm(qnorm(1-alpha/2)-(c(b00,b01.W)/StdError.beta)) +  pnorm(qnorm(alpha/2)-(c(b00,b01.W)/StdError.beta))
  }

  names(StdError.beta) = names(power.hat) = c('b00','b01.W')
  
  toc(log = TRUE, quiet = TRUE)
  log.txt =  tic.log(format = TRUE)
  log.lst = tic.log(format = FALSE)
  tic.clearlog()
  timings.analytic.power = unlist(lapply(log.lst, function(x) x$toc - x$tic))
  
  return(c(StdError.beta=StdError.beta,
           power.hat=power.hat,
           timings.analytic.power=timings.analytic.power))
}
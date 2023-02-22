###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for IL studies using analytic derivations
# Model with a dichotonomous Level 1 and Level 2 predictors


Power.Analytic.ML.4.VAR = function(N.0,N.1,T.obs,
                             b00,b01.Z,b10,b11.Z,
                             sigma,rho,
                             sigma.v0,sigma.v1,rho.v,
                             mu.X0,sigma.X0,mu.X1,sigma.X1,
                             sigma.v0.X0,sigma.v0.X1,
                             rho.X0,rho.X1, 
                             alpha,
                             side.test){
  
  tic()

# Compute information matrix  
c_11_0 = ((1-rho^2)*sigma^2)^(-1)*(T.obs+rho*(T.obs-2)*(rho-2))
c_12_0 = c_21_0 = 0
c_22_0 = ((1-rho^2)*sigma^2)^(-1)*((1+rho^2)*(T.obs-1)*(sigma.v0.X0^2+sigma.X0^2/(1-rho.X0^2)) - 2*rho*(T.obs-1)*(sigma.v0.X0^2+rho.X0*sigma.X0^2/(1-rho.X0^2)))
C_0 = cbind(c(c_11_0,c_21_0),c(c_12_0,c_22_0))

c_11_1 = ((1-rho^2)*sigma^2)^(-1)*(T.obs+rho*(T.obs-2)*(rho-2))
c_12_1 = c_21_1 = 0
c_22_1 = ((1-rho^2)*sigma^2)^(-1)*((1+rho^2)*(T.obs-1)*(sigma.v0.X1^2+sigma.X1^2/(1-rho.X1^2)) - 2*rho*(T.obs-1)*(sigma.v0.X1^2+rho.X1*sigma.X1^2/(1-rho.X1^2)))
C_1 = cbind(c(c_11_1,c_21_1),c(c_12_1,c_22_1))

Sigma.nu = cbind(c(sigma.v0^2,rho.v*sigma.v0*sigma.v1),c(rho.v*sigma.v0*sigma.v1,sigma.v1^2))

D_0 = C_0 %*% solve(diag(2) + Sigma.nu%*%C_0) %*% Sigma.nu %*% t(C_0)
D_1 = C_1 %*% solve(diag(2) + Sigma.nu%*%C_1) %*% Sigma.nu %*% t(C_1)

A_0 = cbind(rbind(N.0*C_0,N.0*0*C_0),rbind(N.0*0*C_0,N.0*0*C_0))
B_0 = cbind(rbind(N.0*D_0,N.0*0*D_0),rbind(N.0*0*D_0,N.0*0*D_0))

A_1 = cbind(rbind(N.1*C_1,N.1*C_1),rbind(N.1*C_1,N.1*C_1))
B_1 = cbind(rbind(N.1*D_1,N.1*D_1),rbind(N.1*D_1,N.1*D_1))

A = A_0 + A_1
B = B_0 + B_1

Inf.beta =  A - B 
Cov.beta = solve(Inf.beta)
StdError.beta = sqrt(diag(Cov.beta))
  
# Compute power
  
  if (side.test == 1){ # One-side test: positive
  power.hat = 1-pnorm(qnorm(1-alpha)-(c(b00,b10,b01.Z,b11.Z)/StdError.beta)) 
  }
  
  if (side.test == 2){ # One-side test: negative
    power.hat = pnorm(-(qnorm(1-alpha)-(c(b00,b10,b01.Z,b11.Z)/StdError.beta))) 
  }
  
  if (side.test == 3){ # Two-tailed test
    power.hat = 1-pnorm(qnorm(1-alpha/2)-(c(b00,b10,b01.Z,b11.Z)/StdError.beta)) +  pnorm(qnorm(alpha/2)-(c(b00,b10,b01.Z,b11.Z)/StdError.beta))
  }

  names(StdError.beta) = names(power.hat) = c('b00','b10','b01.Z','b11.Z')
  
  toc(log = TRUE, quiet = TRUE)
  log.txt =  tic.log(format = TRUE)
  log.lst = tic.log(format = FALSE)
  tic.clearlog()
  timings.analytic.power = unlist(lapply(log.lst, function(x) x$toc - x$tic))
  
  return(c(StdError.beta=StdError.beta,
           power.hat=power.hat,
           timings.analytic.power=timings.analytic.power))
}
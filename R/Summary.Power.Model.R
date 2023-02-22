###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################
###############################################################

# Function to compute power for multilevel models applied to IL studies as a function of the number of participants

Summary.Power.Model = function(Model,N,N.0,N.1,T.obs,  
b00,b01.Z,b01.W,b10,b11.Z,b11.W,
sigma,rho,
sigma.v0,sigma.v1,rho.v,
mu.X,sigma.X,sigma.X.v0,rho.X,
mu.X0,sigma.X0,
mu.X1,sigma.X1,
sigma.v0.X0,sigma.v0.X1,
rho.X0,rho.X1,
mu.W,sigma.W,
isX.center,isW.center,
alpha,side.test){

########################################################################################
########################################################################################
########################################################################################

if (length(alpha) == 0) {stop('The Type I error must be between 0 and 1')}
if (alpha > 1) {stop('The Type I error must be between 0 and 1')}
if (alpha == 0) {stop('The Type I error must be between 0 and 1')}
if (length(T.obs) == 0) {stop('Set the number of time points')}
if (length(sigma) == 0) {stop('Set the value of the standard error of the level-1 residual')}
if (sigma <= 0) {stop('The variance of the level-1 standard error must be positive')}

###############################################################
###############################################################
###############################################################
  
# Check parameters & and output messages
  
if (length(alpha) == 0) {stop('The Type I error must be between 0 and 1')}
if (alpha > 1) {stop('The Type I error must be between 0 and 1')}
if (alpha == 0) {stop('The Type I error must be between 0 and 1')}
  
########################################################################################
########################################################################################
########################################################################################
  
# Simulate data from Model 1: Y = b00 + b01*Z with random effects and autocorrelated errors
  
if (Model == 1){
  
N.0 = as.numeric(unlist(strsplit(N.0,",")))
N.1 = as.numeric(unlist(strsplit(N.1,",")))

if (length(N.0) == 1) {stop('The length of the vector with the number of participants in Group 0 must be larger than one')}
if (length(N.1) == 1) {stop('The length of the vector with the number of participants in Group 1 must be larger than one')}
if (abs(length(N.0) - length(N.1)) > 0) {stop('The vector with the number of participants in Group 0 must have the same length to the vector with the number of participants in Group 1')}
if (length(rho) == 0) {stop('Set the value of the autocorrelation of the Level 1 errors')}
if (abs(rho)>1) {stop('The absolute value of the autocorrelation of the Level 1 errors must be included in the interval [-1,1]')}
if (length(N.0) == 0) {stop('The number of participants in Group 0 must be positive')}
if (length(N.1) == 0) {stop('The number of participants in Group 1 must be positive')}
if (length(b00) == 0) {stop('Set the value of the fixed intercept')}
if (length(b01.Z) == 0) {stop('Set the value of the effect of the Level 2 dummy variable on the intercept')}
if (length(sigma.v0) == 0) {stop('Set the value of the standard deviation of the random intercept')}
if (sigma.v0 < 0) {stop('The standard deviation of the random intercept must be positive')}

# Compute power using the analytic approach
# Power analysis using the analytic approach    
fit.analytic = lapply(1:length(N.0), function(i)
Power.Analytic.ML.1(N.0[i],N.1[i],T.obs,
                    b00,b01.Z,
                    sigma,rho,
                    sigma.v0,
                    alpha,
                    side.test))

power.summary = list(fit.analytic=fit.analytic)
}

######################################################################################
######################################################################################
######################################################################################

# Simulate data from Model 2: Y = b00 + b01*W with random effects and autocorrelated errors
  
if (Model == 2){
    
N = as.numeric(unlist(strsplit(N,",")))

if (length(N) == 1) {stop('The length of the vector with the number of participants must be larger than one')}
if (length(rho) == 0) {stop('Set the value of the autocorrelation of the Level 1 errors')}
if (abs(rho)>1) {stop('The absolute value of the autocorrelation of the Level 1 errors must be included in the interval [-1,1]')}
if (length(N) == 0) {stop('The number of participants must be positive')}
if (length(b00) == 0) {stop('Set the value of the fixed intercept')}
if (length(b01.W) == 0) {stop('Set the value of the effect of the Level 2 continuous variable on the intercept')}
if (length(sigma.v0) == 0) {stop('Set the value of the standard deviation of the random intercept')}
if (sigma.v0 < 0) {stop('The standard deviation of the random intercept must be positive')}  

# Compute power using the analytic approach

# Power analysis using the analytic approach    
fit.analytic = lapply(1:length(N), function(i)
Power.Analytic.ML.2(N[i],T.obs,
                             b00,b01.W,
                             sigma,rho,
                             sigma.v0,
                             mu.W,sigma.W,isW.center,
                             alpha,
                             side.test))

power.summary = list(fit.analytic=fit.analytic)
}

######################################################################################
######################################################################################
######################################################################################

# Simulate data from Model 3: Y = b00 + b10*X with random effects and autocorrelated errors
# Person differences in the distribution of the Level 1 predictor 
 
if (Model == 3){
    
N = as.numeric(unlist(strsplit(N,",")))

if (length(N) == 1) {stop('The length of the vector with the number of participants must be larger than one')}
if (length(rho) == 0) {stop('Set the value of the autocorrelation of the Level 1 errors')}
if (abs(rho)>1) {stop('The absolute value of the autocorrelation of the Level 1 errors must be included in the interval [-1,1]')}
if (length(N) == 0) {stop('The number of participants must be positive')}
if (length(b00) == 0) {stop('Set the value of the fixed intercept')}
if (length(b10) == 0) {stop('Set the value of the fixed slope')}
if (length(sigma.v0) == 0) {stop('Set the value of the standard deviation of the random intercept')}
if (length(sigma.v1) == 0) {stop('Set the value of the standard deviation of the random slope')}
if (length(rho.v) == 0) {stop('Set the value of the correlation between the random intercept and the random intercept')}
if (abs(rho.v)>1) {stop('The absolute value of the correlation between the random intercept and the random intercept must be included in the interval [-1,1]')}
if (length(mu.X) == 0) {stop('Set the value of the mean of the Level 1 predictor')}
if (length(sigma.X) == 0) {stop('Set the value of the standard deviation of the level-1 predictor')}
if (sigma.v0 < 0) {stop('The standard deviation of the random intercept must be positive')}
if (sigma.v1 < 0) {stop('The standard deviation of the random slope must be positive')}

# Compute power using the analytic approach

# Level 1 predictor has a constant within-person variance  
# Power analysis using the analytic approach    
fit.analytic = lapply(1:length(N), function(i)
Power.Analytic.ML.3.VAR(N[i],T.obs,
                             b00,b10,
                             sigma,rho,
                             sigma.v0,sigma.v1,rho.v,
                             mu.X,sigma.X,sigma.X.v0,rho.X, 
                             alpha,
                             side.test))

power.summary = list(fit.analytic=fit.analytic)
}

######################################################################################
######################################################################################
######################################################################################

# Simulate data from Model 4: Y = b00 + b01*Z + b10*X + b11*X*Z with random effects and autocorrelated errors
# Person differences in the distribution of the Level 1 predictor
  
if (Model == 4){
    
N.0 = as.numeric(unlist(strsplit(N.0,",")))
N.1 = as.numeric(unlist(strsplit(N.1,",")))

if (length(N.0) == 1) {stop('The length of the vector with the number of participants in Group 0 must be larger than one')}
if (length(N.1) == 1) {stop('The length of the vector with the number of participants in Group 1 must be larger than one')}
if (abs(length(N.0) - length(N.1)) > 0) {stop('The vector with the number of participants in Group 0 must have the same length to the vector with the number of participants in Group 1')}
if (length(rho) == 0) {stop('Set the value of the autocorrelation of the Level 1 errors')}
if (abs(rho)>1) {stop('The absolute value of the autocorrelation of the Level 1 errors must be included in the interval [-1,1]')}
if (length(N.0) == 0) {stop('The number of participants in Group 0 must be positive')}
if (length(N.1) == 0) {stop('The number of participants in Group 1 must be positive')}
if (length(b00) == 0) {stop('Set the value of the fixed intercept')}
if (length(b01.Z) == 0) {stop('Set the value of the effect of the Level 2 dummy variable on the intercept')}
if (length(b10) == 0) {stop('Set the value of the fixed slope')}
if (length(b11.Z) == 0) {stop('Set the value of the effect of the Level 2 dummy variable on the level-1 slope')}
if (length(sigma.v0) == 0) {stop('Set the value of the standard deviation of the random intercept')}
if (length(sigma.v1) == 0) {stop('Set the value of the standard deviation of the random slope')}
if (length(rho.v) == 0) {stop('Set the value of the correlation between the random intercept and the random intercept')}
if (abs(rho.v)>1) {stop('The absolute value of the correlation between the random intercept and the random intercept must be included in the interval [-1,1]')}
if (length(mu.X0) == 0) {stop('Set the value of the mean of the time-varying predictor in Group 0')}
if (length(sigma.X0) == 0) {stop('Set the value of the standard deviation of the time-varying predictor in Group 0')}
if (length(mu.X1) == 0) {stop('Set the value of the mean of the Level 1 predictor in Group 1')}
if (length(sigma.X1) == 0) {stop('Set the value of the standard deviation of the Level 1 predictor in Group 1')}
if (length(sigma.v0.X0) == 0) {stop('Set the value of the standard deviation of the random intercept of the Level 1 predictor in Group 0')}
if (length(rho.X0) == 0) {stop('Set the value of the autocorrelation of the errors of the Level 1 predictor in Group 0')}
if (abs(rho.X0)>1) {stop('The absolute value of the autocorrelation of the errors of the Level 1 predictor in Group 0 must be included in the interval [-1,1]')}
if (length(sigma.v0.X1) == 0) {stop('Set the value of the standard deviation of the random intercept of the Level 1 predictor in Group 1')}
if (length(rho.X1) == 0) {stop('Set the value of the autocorrelation of the errors of the Level 1 predictor in Group 1')}
if (abs(rho.X1)>1) {stop('The absolute value of the autocorrelation of the errors of the Level 1 predictor in Group 1 must be included in the interval [-1,1]')}
if (sigma.v0 < 0) {stop('The standard deviation of the random intercept must be positive')}
if (sigma.v1 < 0) {stop('The standard deviation of the random slope must be positive')} 

# Compute power using the analytic approach

# Level 1 predictor has a constant within-person variance    
# Power analysis using the analytic approach    
fit.analytic = lapply(1:length(N.0), function(i)
Power.Analytic.ML.4.VAR(N.0[i],N.1[i],T.obs,
                             b00,b01.Z,b10,b11.Z,
                             sigma,rho,
                             sigma.v0,sigma.v1,rho.v,
                             mu.X0,sigma.X0,mu.X1,sigma.X1,
                             sigma.v0.X0,sigma.v0.X1,
                             rho.X0,rho.X1, 
                             alpha,
                             side.test))

power.summary = list(fit.analytic=fit.analytic)
}

######################################################################################
######################################################################################
######################################################################################

# Simulate data from Model 5: Y = b00 + b01*W + b10*X + b11*X*W with random effects and autocorrelated errors
# Person differences in the distribution of the Level 1 predictor
  
if (Model == 5){
    
N = as.numeric(unlist(strsplit(N,",")))

if (length(N) == 1) {stop('The length of the vector with the number of participants must be larger than one')}
if (length(rho) == 0) {stop('Set the value of the autocorrelation of the Level 1 errors')}
if (abs(rho)>1) {stop('The absolute value of the autocorrelation of the Level 1 errors must be included in the interval [-1,1]')}
if (length(N) == 0) {stop('The number of participants must be positive')}
if (length(b00) == 0) {stop('Set the value of the fixed intercept')}
if (length(b01.W) == 0) {stop('Set the value of the effect of the Level 2 continuous variable on the intercept')}
if (length(b10) == 0) {stop('Set the value of fixed slope')}
if (length(b11.W) == 0) {stop('Set the value of the effect of the Level 2 continuous variable on the slope')}
if (length(sigma.v0) == 0) {stop('Set the value of the standard deviation of the random intercept')}
if (length(sigma.v1) == 0) {stop('Set the value of the standard deviation of the random slope')}
if (length(rho.v) == 0) {stop('Set the value of the correlation between the random intercept and the random intercept')}
if (abs(rho.v)>1) {stop('The absolute value of the correlation between the random intercept and the random intercept must be included in the interval [-1,1]')}
if (length(mu.X) == 0) {stop('Set the value of the mean of the Level 1 predictor')}
if (length(sigma.X) == 0) {stop('Set the value of the standard deviation of the Level 1 predictor')}
if (length(sigma.X.v0) == 0) {stop('Set the value of the standard deviation of the random intercept of the Level 1 predictor')}
if (length(rho.X) == 0) {stop('Set the value of the autocorrelation of the errors of the Level 1 predictor')}
if (abs(rho.X)>1) {stop('The absolute value of the autocorrelation of the errors of the Level 1 predictor must be included in the interval [-1,1]')}
if (length(mu.W) == 0) {stop('Set the value of the mean of the Level 2 continuous predictor')}
if (length(sigma.W) == 0) {stop('Set the value of the standard deviation of the Level 2 continuous predictor')}
if (sigma.v0 < 0) {stop('The standard deviation of the random intercept must be positive')}
if (sigma.v1 < 0) {stop('The standard deviation of the random slope must be positive')}    

# Compute power using the analytic approach

# Level 1 predictor has a constant within-person variance   
# Power analysis using the analytic approach    
fit.analytic = lapply(1:length(N), function(i)
Power.Analytic.ML.5.VAR(N[i],T.obs,
                             b00,b01.W,b10,b11.W,
                             sigma,rho,
                             sigma.v0,sigma.v1,rho.v,
                             mu.W,sigma.W,isW.center,
                             mu.X,sigma.X,sigma.X.v0,rho.X, 
                             alpha,
                             side.test))

power.summary = list(fit.analytic=fit.analytic)
}

######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
  
return(power.summary)}

#####################################################################################



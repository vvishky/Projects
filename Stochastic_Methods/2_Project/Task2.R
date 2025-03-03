#2. Task
#A noisy signal curve is given as output of a spectroscopic analysis in the file signal.txt. A Lorentz
#model curve is to be fitted to this curve,

#a) Determine the parameters by nonlinear regression. Choose appropriate initial values. Graphically
#output the estimated signal curve and data points.

#importing data from signal.txt using File>import data set > From text(base)
signal_data = signal

n= dim(signal_data)[1]

# plotting the data points:
plot(signal_data$x,signal_data$y,xlab="Input x",ylab="Output L", main  = "Lorentz model")

# Choosing of appropriate initial values
init = list(a=2.5, x_0=255)

# nonlinear regression
model_nl = nls( y ~ a/ (1+(x-x_0)^2), data=signal_data, start=init, trace=TRUE)
summary(model_nl)

#From nonlinear regression parameters are a = 3.164e+00 and x_0 = 2.560e+02.

#plotting estimated signal curve.
y_hat = predict(model_nl)
lines(signal_data$x, y_hat, lwd=2, col="red")

#b) Examine whether the maximum signal strength (amplitude) of the estimated signal curve
#is significantly above the critical value Lc = 3.

# LS estimation for theta = (a,b)
theta_hat = as.numeric(coef(model_nl))

# estimation of error variance
sigma2_hat = deviance(model_nl) / (n-2)

# We first define a function for the evaluation of 
# j(x,theta) = nabla_theta f(x,theta)
# 
j_func = function(x, theta){
  j = c( 1/(1 + (x-theta[2])^2), (2*theta[1]*(x-theta[2])) / (1 + (x-theta[2])^2)^2)
  return(j)
}

# We use this to create the matrices J and G
J = matrix(rep(0,n*2), ncol = 2)
for(i in 1:n){
  J[i,] = j_func(signal_data$x[i], theta_hat)
}
G = solve( t(J) %*% J )

# significance level
alpha = 0.05

# we calculate the width of the asymptotic confidence band
# by means of a function
Cband = function(x){
  j = j_func(signal_data$x,theta_hat)
  # NOTE: F distribution and number of parameters = 2:
  w = sqrt(2 * sigma2_hat * qf(1-alpha,2,n-2) * t(j) %*% G %*% j)
  return(w)
}

# We calculate the width of the asymptotic prediction intervals
# by means of a function
Pint = function(x){
  j = j_func(signal_data$x,theta_hat)
  # NOTE: t distribution and + 1
  v = sqrt(sigma2_hat * (1 + t(j) %*% G %*% j)) * qt(1-alpha/2,n-2)
  return(v)
}

# We draw the graphs
par(mfrow=c(1,1)) # divide the graph output into 1 window
plot(signal_data$x,signal_data$y,xlab="Input x",ylab="Output L", main  = "Lorentz model")
y_hat = predict(model_nl)
y_true = 2.5/ (1+(signal_data$x-255)^2)
lines(signal_data$x, y_hat, lwd=2, col="green3")
lines(signal_data$x, y_true, lwd=1, lty=2)

# Adding the borders
Cwidth = rep(0,length(signal_data$x))
Pwidth = rep(0,length(signal_data$x))
for(j in 1:length(signal_data$x)){
  # width of the asymptotic confidence band 
  # and the prediction intervals at all points 
  # evaluate from x
  Cwidth[j] = Cband(signal_data$x[j]);
  Pwidth[j] = Pint(signal_data$x[j]);
}
lines(signal_data$x, y_hat - Cwidth,lty=1, col="red", lwd=2)
lines(signal_data$x, y_hat + Cwidth,lty=1, col="red", lwd=2)
lines(signal_data$x, y_hat - Pwidth,lty=2, col="purple2", lwd=2)
lines(signal_data$x, y_hat + Pwidth,lty=2, col="purple2", lwd=2)
#################################################





#3) Task
#In the text file coal data.txt you find the time intervals in days between 
#disasters in British coal mines between 1850 and 1965.

#a) Load the data into R and visualize the empirical distribution via a histogram. 
#Also draw a kernel density estimate for the underlying distribution. What do 
#you observe regarding the (shape of the) distribution?

#Loading the coal data from file
coal_data = scan("coal_data.txt")

#empirical distribution via histogram
hist (coal_data, nclass = 14, freq = FALSE)

#kernel density estimate
lines(density (coal_data),col =" red ",lwd =3)

#By observing the shape of distribution for the coal data, the distribution
#is Log normal distribution.


#b) Compute a confidence interval for the mean value of days between two 
#disasters for confidence level of 95%.

#number of data
n = 190

#Estimation of mean
mw_hat = mean(coal_data)

#Estimation of  variance 
var_hat= var(coal_data)

#Estimation of standard deviation
s_hat= sd(coal_data)

#Calculate half width of confidence interval using quantile
# of the t-distribution by qt
alpha = 0.05 #confidence level 95%

Delta_hat = s_hat/sqrt(n) * qt(1-0.5*alpha, n-1)

# lower and upper limit of the confidence interval
theta_low = mw_hat - Delta_hat
theta_upp = mw_hat + Delta_hat

confidence_interval = c(theta_low,theta_upp)

confidence_interval



#1).Task
#In an investigation of the Cd contamination of trout in a river, ten trout 
#were caught at each of two locations and their Cd content (in mg/g fresh weight)
#was determined.

#Location_A
LA = c(76.8,72.3,74.0,73.2,46.1,76.5,61.9,62.4,65.9,62.4)
#Location_B
LB = c(64.4,60.0,59.4,61.2,52.0,58.1,62.0,57.8,57.2)


#a) Draw parallel box plots. What do you observe?

#Summary and Parallel box plots

summary <- boxplot(LB,LA, horizontal = TRUE, names = c("Location_B", "Location_A"), 
                   xlab ="Cd contamination of trout in a river", ylab ="Locations")

colnames(summary$stats) <- c("Location_B","Location_A")
rownames(summary$stats)<-c("Min","First Quartile","Median","Third Quartile","Maximum")
summary$stats

#From the Box plot we can clearly observe for the Location A minimum and maximum
#values for Cd contamination are 46.1 and 76.8 respectively. The first and third
#quartiles are 62.4 and 74.0 respectively. And Median value is 69.1.

#And for the Location B minimum and maximum values for Cd contamination are 57.2 
#and 64.4 respectively. The first and third quartiles are 57.8 and 61.2 respectively. 
#Median value is 59.4. And outlier data value is 52.0.

#The distribution of data for Location A is negatively skewed because the whisker 
#and half-box are longer on the left side of the median than on the right side.
#And the center of distribution (median = 69.1) is highest as compared with Location B.

#The distribution of data for Location B is approximately symmetric because median
# is almost dividing the interquartile range in half (1.6 on the left side and 1.8
# on the right side). This location has data distribution most concentrated because 
#interquartile range is 3.4 as compared with Location A i.e 11.6. This distribution 
#include potential outliers.The interquartile range is Q3 - Q1 = 61.2-57.8 = 3.4.
#According to the definition used by the function in R software, all values lower 
#than Q1 - 1.5 x (Q3 - Q1) = 57.8 - 1.5 x 3.4 = 52.7 are outside the left whisker
#and indicated by a circle. Their is one  potential outliers in distribution B.


#b) Test for the level α = 0.05 whether the Cd contents measured at both locations 
#can be regarded as realizations of a normally distributed random variable.

#Two independent sample LA,LB.
#Test For normal distribution using Shapiro-Wilk test.

#Assumptions:- 
#H0 : Data is normally distributed ; HA : data is not normally distributed.

shapiro.test(LA)

#P-value for Location A is 0.09526.

shapiro.test(LB)
#P-value for Location B is 0.7679.

#Since, both the p-values are way above the level of significance of α/2= 0.05/2 = 0.025.
#We can not reject that both data follows normal distribution. So, data is normally distributed.


#c) Test for the level α = 0.05 whether the variances of the Cd contents are
#equal or significantly different from each other.

#Test for variances in two independent samples is done by using F-test.

#Assumptions:- 
#H0 : Variances are equal. ; HA : Variances are not equal.
var(LB)
var(LA)

var.test(LB,LA)

#P-value for both Locations is 0.01022.
#Since the p-value is below our α = 0.05 we reject the assumption that the
#variances are equal. And variances are significantly different from each other.
#(Variance of LB = 12.31444 and Variance of LA = 89.77167)

#d) Test for the level α = 0.05 whether the expected Cd content at location A 
#is significantly greater than that at location B.

#Thus, we have now identified the proper test for the task regarding the location
#data: Welch t-test.

#Assumptions:- 
#H0 : Means are equal. ; HA : Means are not equal(two sided) or Means are greater and less than(one sided).


t.test(LB,LA,var.equal = FALSE,alternative = "less")


#P-value for both Locations is 0.01433.
#Since the p-value is below our α = 0.05 we reject the assumption that the
#means are equal. And means are significantly different from each other.

#From the above test we can say that the expected Cd Content at location A  is
#significantly greater than that at location B.
#(Mean of LA=67.15 and Mean of LB=59.122)
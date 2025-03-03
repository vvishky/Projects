#1. Task
#In the data set timber.txt, the experimentally determined values of the stiffness modulus (Rigid),
#elastic modulus (Elast) and density (Dens) of air-dried wood are available for 49 samples of
#different wood species.

#a) Investigate which of the three variables has the highest linear relationship to the two other
#ones, i.e., which of the following three models is most suitable for explaining data set
#timber.txt:

#The data set:
#Importing data set by File > Import data set  > From Text (base) > timber.txt
#File name timber is visible in environment
Timber_data = timber

# number of data
n = dim(Timber_data)[1]

#Linear regression of Timber_data
#Rigid = θ0 + θ1 Elast + θ2 Dens
reg1 = lm(formula = Rigid ~ Elast+Dens, data = Timber_data);
#summary(reg1)
#Elast = θ0 + θ1 Rigid + θ2 Dens
reg2 = lm(formula = Elast ~ Rigid+Dens, data = Timber_data);
#summary(reg2)
#Dens = θ0 + θ1 Rigid + θ2 Elast
reg3 = lm(formula = Dens ~ Rigid+Elast, data = Timber_data);
#summary(reg3)

#By observing summary for all three models. I can say that the first model reg1
#has the highest linear relationship. Because the p-value of reg1 model is lowest
#among three models and also the multiple R-squared value is highest as compared
#with other two models. For reg1 model Multiple R-squared:   0.81.

#b) Using a residual analysis, check whether there are any abnormalities in the regressed model
#chosen in a).

#Plots for residual analysis
par ( mfrow = c(2, 2)) # Divide plot window in 2x2 parts
plot (reg1) # draw the four graphs

#From the first graph (Residuals vs Fitted) we can clearly observe that the red line
#varies randomly with no recognizable structure around zero and their is no strong trend in it. 
#So there is no missing of function f in the linear regression model.

#From the second graph (Normal Q-Q) we can clearly observe that the residuals
#are normally distributed.By using plot(density(reg1$residuals)) it can give one more
#proof that the residuals are normally distributed.So, the assumption for normal distribution is true.

#From the third graph(Scale - Location) we can clearly observe that their is
#no recognizable trend in red line and the points are also uniformly distributed.
#So, the assumption of homoscedasticity is true.

#From the fourth graph(Residuals vs leverage) we can clearly observe that the 46th and 41st
#data points are Outlier because standardized residual of that point is greater than 2.
#Also 40th data point is influential point because it lies behind cook's distance
#lines D = 0.5.

#Removing 40th, 41st and 46th data points
Timber_data_r = Timber_data[c(1:39,42:45,47:49),]

#Rerunning of linear regression
reg_r = lm(formula = Rigid ~ Elast+Dens, data = Timber_data_r);

#Plots for residual analysis
par ( mfrow = c(2, 2)) # Divide plot window in 2x2 parts
plot (reg_r) # draw the four graphs

#c) Investigate whether the model is significantly improved by adding an interaction term
#between the two explanatory variables.

#Adding interaction term Elas*Dens
reg_r_I = lm(formula = Rigid ~ Elast+Dens +Elast*Dens, data = Timber_data_r);

#After seeing summary of reg_r_d we can identify Elast*Dens as non significant.

#old model without interaction term
reg_r = lm(formula = Rigid ~ Elast+Dens, data = Timber_data_r)

#F-test for model comparison
#Assumption:H0 : small model is sufficient, HA : small model is not sufficient .

anova(reg_r_I,reg_r)

#From the above test we get p-value = 0.0658 which is above significance level alpha = 0.05
#which means the smaller model is sufficient. So, model without interaction term
#reg_r is significant.


#d) A 50th wood sample was measured and the following values were observed:

#Adding new data
Timber_data_new = rbind(Timber_data_r,c(2078, 237.5, 70.8))

#Model
reg_r_n = lm(formula = Rigid ~ Elast+Dens, data = Timber_data_new)

plot(reg_r_n)

#From residual analysis we can say that the new data is fitted in regressed model
#because new data is not outlier or influential point in residual vs leverage graph.
#So, new data point is in appropriate 95% range.

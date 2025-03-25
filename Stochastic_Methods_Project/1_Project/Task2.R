#2). Task
#The measured temperatures (in degree C) before and 3 hours after taking a drug can be 
#found for ten different patients.

#before
B = c(38.4,39.6,39.4,40.1,39.2,38.5,39.3,39.1,38.4,39.5)

#After
A = c(37.6,37.9,39.1,39.4,38.6,38.9,38.7,38.7,38.9,38.7)

#Difference
D = B-A

#a) Describe which sample situation (one, two, paired, multiple, ...) we have here.
#The sample situation we have here is paired samples.


#b) Generate a scatter plot of the data. What do you observe?

plot(B,A,main = " Measurement of Temperature in Celcius",xlab = " Before taking drug",
     ylab = "After taking drug",pch = 19)
#From the scatter plot we can observe that the points are spread out more. This 
#means there is no trend in data.


#c) Estimate the correlation between the temperature values before and after 
#taking the drug.

cor(B,A)

#the correlation between the temperature values before and after taking the 
#drug is 0.3485987.


#d) Test for level α = 0.05 if the correlation is significant (you can decide 
#in one- or two-sided alternative).

#Correlation significance test is done using Pearson correlation method.
#Assumption:-
#Null hypothesis H0 is equal to zero(uncorrelated), Alternative hypothesis HA is 
#not equal to zero (two sided) or greater and less than zero (one sided).

cor.test(B,A,alternative = "greater")
#P-value from the test is 0.1618.
#P-value is above given significance level α = 0.05. We can not reject the null hypothesis.Hence,
#Correlation is not significant for given data(uncorrelated).

#d)continue. To this end, research (in literature or online) which test is 
#suitable for this task and verify that its assumptions are satisfied for the data at hand.

#First step is to check whether the difference D between the data is normally distributed
#using Shapiro-Wilk test.
#Assumptions:- 
#H0 : Data is normally distributed ; HA : Data is not normally distributed.

shapiro.test(D)
#P-value is 0.3265. This p-value is above α = 0.05 we cannot reject that data follows normal
#distribution.So, difference D is normally distributed.

#As the difference D is normally distributed we can use Paired t-test.
#Assumptions:- 
#H0 : mean difference is equal to zero ; HA : mean difference is not equal to zero (two sided)
#or greater and less than zero (one sided).

t.test(B,A,alternative = "greater", paired = TRUE)
#The p-value of the test is 0.01635, which is less than the significance level α = 0.05. 
#We can then reject null hypothesis H0 and conclude that the average temperature of the patient
#before treatment is significantly different from the average temperature after treatment 
#with a p-value = 0.01635.
# Analysis of Cd Contamination in Trout, Drug Temperature Study, and Coal Mine Disasters

## Task 1: Cd Contamination in Trout

### Task Overview
This project investigates the Cd contamination of trout in a river. Ten trout were caught at each of two locations (Location A and Location B), and their Cd content (mg/g fresh weight) was measured.

### Data
- **Location A:** `76.8, 72.3, 74.0, 73.2, 46.1, 76.5, 61.9, 62.4, 65.9, 62.4`
- **Location B:** `64.4, 60.0, 59.4, 61.2, 52.0, 58.1, 62.0, 57.8, 57.2`

### Analysis

#### a) Parallel Box Plots
A box plot was created to compare the distributions of Cd contamination at both locations.

**Observations:**
- Location A has a higher median Cd contamination (69.1) than Location B (59.4).
- The distribution of Location A is negatively skewed, while Location B is approximately symmetric.
- Location B has an outlier (52.0), and its interquartile range (IQR) is smaller (3.4) compared to Location A (11.6), indicating more concentrated data in Location B.

#### b) Normality Test
A Shapiro-Wilk test was conducted to test whether the Cd contents at both locations follow a normal distribution.

**Results:**
- **Location A:** p-value = 0.09526
- **Location B:** p-value = 0.7679
- Since both p-values are greater than the significance level (0.025), we fail to reject the null hypothesis. **Both datasets follow a normal distribution.**

#### c) Variance Test
An F-test was performed to determine if the variances of Cd contents at the two locations are significantly different.

**Results:**
- **Variance of Location A:** 89.77
- **Variance of Location B:** 12.31
- **p-value:** 0.01022
- Since the p-value is less than 0.05, we reject the null hypothesis. **The variances are significantly different.**

#### d) Mean Comparison
A Welch t-test was conducted to test whether the expected Cd content at Location A is significantly greater than at Location B.

**Results:**
- **Mean of Location A:** 67.15
- **Mean of Location B:** 59.12
- **p-value:** 0.01433
- Since the p-value is less than 0.05, we reject the null hypothesis. **The expected Cd content at Location A is significantly greater than at Location B.**

## Task 2: Drug Temperature Study

### Task Overview
This study examines the temperatures of ten patients before and three hours after taking a drug. The goal is to determine whether the drug significantly affects body temperature.

### Data
- **Before Taking Drug:** `38.4, 39.6, 39.4, 40.1, 39.2, 38.5, 39.3, 39.1, 38.4, 39.5`
- **After Taking Drug:** `37.6, 37.9, 39.1, 39.4, 38.6, 38.9, 38.7, 38.7, 38.9, 38.7`
- **Difference (Before - After):** `D = B - A`

### Analysis

#### a) Sample Situation
The given data represents **paired samples** since the measurements were taken from the same individuals before and after taking the drug.

#### b) Scatter Plot
A scatter plot was generated to visualize the relationship between temperatures before and after taking the drug.

**Observations:**
- The points are spread out, indicating **no clear trend** in the data.

#### c) Correlation Estimation
The correlation between temperature values before and after taking the drug was calculated.

**Result:**
- Correlation coefficient = **0.3486**
- The correlation is weak.

#### d) Correlation Significance Test
A Pearson correlation test was performed to check if the correlation is statistically significant.

**Results:**
- **p-value:** 0.1618
- Since the p-value is greater than 0.05, we fail to reject the null hypothesis. **The correlation is not significant.**

#### e) Normality Test for Differences
A Shapiro-Wilk test was conducted to test whether the differences (D) follow a normal distribution.

**Results:**
- **p-value:** 0.3265
- Since the p-value is greater than 0.05, we fail to reject the null hypothesis. **The differences are normally distributed.**

#### f) Paired t-test
A paired t-test was conducted to determine whether there is a significant difference in mean temperatures before and after taking the drug.

**Results:**
- **p-value:** 0.01635
- Since the p-value is less than 0.05, we reject the null hypothesis. **The average temperature before taking the drug is significantly higher than after taking the drug.**

## Task 3: Coal Mine Disaster Analysis

### Task Overview
This task examines the time intervals (in days) between disasters in British coal mines from 1850 to 1965.

### Analysis

#### a) Data Visualization
A histogram and kernel density estimate were created to visualize the distribution of disaster intervals.

**Observation:**
- The data appears to follow a **log-normal distribution**.

#### b) Confidence Interval for Mean Disaster Interval
A confidence interval for the mean time between disasters was calculated using a 95% confidence level.

**Results:**
- **Lower Limit:** Computed as `theta_low`
- **Upper Limit:** Computed as `theta_upp`
- The computed confidence interval provides an estimate of the expected interval between disasters.

## Conclusion
### Task 1: Cd Contamination in Trout
- The Cd contamination at Location A is significantly higher than at Location B.
- The variances of the two locations are significantly different.
- The data for both locations follow a normal distribution.

### Task 2: Drug Temperature Study
- The correlation between temperatures before and after taking the drug is not significant.
- The differences in temperatures follow a normal distribution.
- The average temperature before taking the drug is significantly higher than after taking the drug, indicating that the drug has a significant effect.

### Task 3: Coal Mine Disaster Analysis
- The distribution of disaster intervals follows a log-normal pattern.
- A 95% confidence interval was calculated for the mean disaster interval.

## Requirements
- **R programming language**
- **Packages:** None (Base R functions used)

## Usage
Run the `Task1.R`, `Task2.R`, and `Task3.R` scripts in R to perform the analyses and visualize the results.


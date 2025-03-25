# Timber and Signal Data Analysis

## Overview
This repository contains statistical analyses performed in R for two datasets:
1. **`timber.txt`**: Experimentally determined values of stiffness modulus (`Rigid`), elastic modulus (`Elast`), and density (`Dens`) for 49 wood samples.
2. **`signal.txt`**: A noisy signal curve from spectroscopic analysis, modeled using a Lorentz function.

Each analysis investigates relationships, fits models, and evaluates their adequacy.

## Datasets

### Timber Dataset
- **File**: `timber.txt`
- **Variables**: 
  - `Rigid`: Stiffness modulus
  - `Elast`: Elastic modulus
  - `Dens`: Density
- **Sample Size**: 49 wood samples (initially)

### Signal Dataset
- **File**: `signal.txt`
- **Variables**: 
  - `x`: Input variable
  - `y`: Output signal strength
- **Purpose**: Fit a Lorentz model to the noisy signal curve.

## Tasks and Methodology

### Timber Analysis

#### Task a: Identifying the Strongest Linear Relationship
- **Objective**: Determine which variable (`Rigid`, `Elast`, or `Dens`) has the highest linear relationship with the other two.
- **Models**:
  1. `Rigid = θ0 + θ1 Elast + θ2 Dens` (`reg1`)
  2. `Elast = θ0 + θ1 Rigid + θ2 Dens` (`reg2`)
  3. `Dens = θ0 + θ1 Rigid + θ2 Elast` (`reg3`)
- **Analysis**: Linear regression models fitted using `lm()` in R.
- **Findings**: 
  - `reg1` (Rigid as response) showed the highest linear relationship.
  - Metrics: Multiple R-squared = 0.81, lowest p-value among the three models.

#### Task b: Residual Analysis of the Chosen Model
- **Objective**: Check for abnormalities in `reg1`.
- **Method**: Residual plots generated using `plot(reg1)`:
  1. Residuals vs Fitted
  2. Normal Q-Q
  3. Scale-Location
  4. Residuals vs Leverage
- **Findings**:
  - Residuals are randomly scattered with no strong trends, indicating a good fit.
  - Normal Q-Q plot and density plot confirm residuals are normally distributed.
  - Scale-Location plot supports homoscedasticity.
  - Outliers detected: 41st and 46th data points (standardized residuals > 2).
  - Influential point detected: 40th data point (beyond Cook’s distance = 0.5).
- **Action**: Removed outliers (40th, 41st, 46th rows) to create `Timber_data_r`.
- **Updated Model**: `reg_r` fitted on cleaned data, residual plots rechecked.

#### Task c: Model Improvement with Interaction Term
- **Objective**: Test if adding an interaction term (`Elast * Dens`) improves `reg_r`.
- **Models**:
  - Without interaction: `reg_r = Rigid ~ Elast + Dens`
  - With interaction: `reg_r_I = Rigid ~ Elast + Dens + Elast * Dens`
- **Analysis**: F-test (`anova(reg_r_I, reg_r)`) to compare models.
- **Findings**:
  - P-value = 0.0658 (> 0.05), indicating the interaction term is not significant.
  - Conclusion: The simpler model (`reg_r`) is sufficient.

#### Task d: Fitting a New Data Point
- **Objective**: Assess how a 50th wood sample fits into the model.
- **New Data**: `Rigid = 2078`, `Elast = 237.5`, `Dens = 70.8`
- **Method**: Added to `Timber_data_r`, new model `reg_r_n` fitted, residual plots analyzed.
- **Findings**: The new data point is not an outlier or influential (within 95% range), indicating a good fit.

---

### Signal Analysis

#### Task a: Nonlinear Regression for Lorentz Model
- **Objective**: Fit a Lorentz model (`y ~ a / (1 + (x - x_0)^2)`) to the noisy signal data and estimate parameters.
- **Initial Values**: `a = 2.5`, `x_0 = 255`
- **Method**: Nonlinear regression using `nls()` in R.
- **Findings**:
  - Estimated parameters: `a = 3.164`, `x_0 = 256.0`.
- **Visualization**: 
  - Data points plotted with `plot(signal_data$x, signal_data$y)`.
  - Estimated curve (`y_hat`) overlaid in red using `lines()`.

#### Task b: Significance of Maximum Signal Strength
- **Objective**: Test if the maximum signal strength (amplitude `a = 3.164`) is significantly above the critical value `Lc = 3`.
- **Method**:
  - Estimated parameters (`theta_hat`) extracted from `coef(model_nl)`.
  - Error variance estimated: `sigma2_hat = deviance(model_nl) / (n-2)`.
  - Jacobian matrix (`J`) and Fisher information matrix (`G`) computed for confidence bands.
  - Asymptotic confidence band width calculated using `Cband()` (F-distribution).
  - Prediction interval width calculated using `Pint()` (t-distribution).
- **Visualization**:
  - Plotted data, fitted curve (green), true curve (dashed), confidence bands (red), and prediction intervals (purple).
- **Findings**: 
  - The amplitude `a = 3.164` exceeds `Lc = 3`, and confidence bands can be used to assess significance (not explicitly quantified in the document but implied by the methodology).

## Files
- **`timber.txt`**: Dataset with 49 wood samples.
- **`signal.txt`**: Noisy signal data from spectroscopic analysis.
- **R Script**: Contains code for both analyses (embedded in the documents).

## Requirements
- **Software**: R
- **Packages**: Base R (no additional packages required).

## How to Run
1. Place `timber.txt` and `signal.txt` in your R working directory.
2. Run the R script step-by-step:
   - **Timber**: Import data (`read.table("timber.txt", header=TRUE)`), execute regression and diagnostics.
   - **Signal**: Import data (`read.table("signal.txt", header=TRUE)`), fit Lorentz model, and plot results.

## Results
- **Timber**: 
  - `reg1` (Rigid as response) is the best model (R-squared = 0.81).
  - After removing outliers (40, 41, 46), `reg_r` is robust.
  - No significant improvement with interaction term.
  - 50th sample fits well.
- **Signal**: 
  - Lorentz model parameters: `a = 3.164`, `x_0 = 256.0`.
  - Maximum signal strength exceeds `Lc = 3`, with confidence bands supporting the fit.

## Author
Analysis performed as part of a statistical modeling exercise.

---

This updated `README.md` integrates both tasks while maintaining clarity and structure. Let me know if further refinements are needed!

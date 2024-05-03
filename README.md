# Production of essential oils Time Series Project


## Overview

This project focuses on modeling and predicting an observed time series of the French economy, specifically the French Industrial Production Index (IPI). The objectives include transforming the series to achieve stationarity if necessary, identifying and estimating appropriate ARMA models, and making predictions for future values.

## Part I: The Data

- Series Description: The chosen series represents the French Industrial Production Index (IPI) corrected from seasonal variations and working days (CVS-CJO). It signifies the industrial output of France, providing insights into the overall economic activity within the industrial sector. Potential data processing may include detrending or differencing to achieve stationarity.
Stationarity Transformation: If the series is not already stationary, transformation techniques such as differencing or detrending will be applied to ensure stationarity. The choice of transformation method will be justified based on statistical tests for stationarity and visual inspection of the series.

- Graphical Representation: The chosen series will be graphically represented before and after transformation to visualize the changes in its characteristics.


## Part II: ARMA Models

- Model Selection: An ARMA(p,q) model will be selected for the corrected time series Xt based on the Akaike Information Criterion (AIC) or Bayesian Information Criterion (BIC), considering the trade-off between model complexity and goodness of fit. The chosen model will be justified based on diagnostic checks and statistical significance of model parameters.
ARIMA Model Specification: The ARIMA(p,d,q) model for the chosen series will be formulated, incorporating the differencing parameter (d) to account for non-stationarity.

# Part III: Prediction

Denote T as the length of the series, and assume Gaussian residuals.

- Confidence Region Equation: An equation for the confidence region of level α on the future values (XT+1, XT+2) will be derived, considering the Gaussian distribution of residuals and the desired confidence level.

- Hypotheses for Confidence Region: The hypotheses used to construct the confidence region will be stated, ensuring the validity of the prediction intervals.

- Graphical Representation: The confidence region for α = 95% will be graphically represented, and observations will be made regarding its width and coverage of future values.

- Open Question: The potential improvement in prediction for XT+1 using information from a stationary time series YT will be explored. Conditions under which this information can enhance prediction accuracy will be identified, along with proposed testing methodologies to validate the improvement.


Contributors
Glorieux Grégoire & Cablant Augustin - Ensae students

License
This project is licensed under the MIT License.

# **Code for "Asymmetric Models for Realized Covariances"**

Luc Bauwens, Emilija Dzuverovic, Christian Hafner

## Overview

The code in this repository is used in the empirical analyses of the paper ''Asymmetric Models for Realized Covariances". The latter were run on Matlab_R2024b, including the adoption of the MFE Toolbox of Kevin Sheppard.

The main contents of the source code are the following:
- Estimation_Tables.m: the script file leading to the main estimation results
- Forecasting_Tables.m: the script file leading to the main forecasting results
- Fun/: folders of the functions for generating the estimation and forecasting results
- Data/: folders of processed data files used as an input for generating the estimation and forecasting results
- Forecasts/: folder of processed data files used as an input for generating the forecasting results

## Instructions

The estimtion and forecasting analyses files can be run individually, in any order. Remember to specify the correct path to the local version of the repository.

## Data availability

### Solar Flare Forecasts
The data provider is the AlgoSeek company (30 Wall Street, New York, NY, 10005, USA). The data provided to us by Algoseek are the prices of the assets, observed every minute during the trading period (9:30-16), compiled from the trades that occurred in sixteen US exchanges and marketplaces. For the subsequent empirical analyses, we have constructed the time series of daily RC matrices with the corresponding decompositions (3) and (11) into positive, negative, and mixed matrices, based on a high-frequency data set for the SPDR S&P500 (SPY or Spyder) { an exchange traded fund that tracks the S&P500 index { and  ve stocks of the banking sector, i.e., Bank of America Corp. (BAC), Citigroup Inc. (C), Goldman Sachs Group Inc. (GS), JPMorgan Chase & Co. (JPM), and Wells Fargo & Co. (WFC). To avoid the measurement drawbacks due to microstructure e ects when sampling returns at very high frequencies, we compute each daily realized (semi-)covariance matrix as the sum of the outer products of the  ve-minute log-return vectors of the trading period of the day. Given the high liquidity of all the stocks, the e ect of non-synchronicity is rather negligible at the chosen frequency; the synchronization was done globally for all the stocks, using the closest prior price. The sample period is January 3, 2012 - December 31, 2021, resulting in 2517 observations.

## References

1. Mafalda (2025). MFE Toolbox - Kevin Sheppard (https://www.mathworks.com/matlabcentral/fileexchange/170381-mfe-toolbox-kevin-sheppard), MATLAB Central File Exchange.



# **Matlab Code for "Asymmetric Models for Realized Covariances"**
Luc Bauwens, Emilija Dzuverovic, Christian Hafner

## Overview
The code in this repository is used in the empirical analyses of the paper "Asymmetric Models for Realized Covariances". The latter introduces a new class of asymmetric models for RC based on the BEKK-type conditional autoregressive Wishart model (Golosnoy et al. (2012)) either by interacting the RC with the signs of the daily returns or by using the semi-covariance decomposition (Bollerslev et al. (2020)) built upon intra-daily returns. The proposed asymmetric BEKK-CAW models feature better in-sample fit and out-of-sample forecasting performance than the benchmark symmetric model. The empirical application additionally shows that the asymmetric models using the signs of the daily returns outperform the models using the signed intra-daily returns. All the empirical exercises were run on Matlab_R2024b, including the adoption of the MFE Toolbox of Kevin Sheppard.

The main contents of the code distribution are the following:
- Estimation_Tables.m: the script file leading to the main estimation results
- Forecasting_Tables.m: the script file leading to the main forecasting results
- Fun/: folders of the functions for generating the estimation and forecasting results
- Data/: folders of data files used as an input for generating the estimation and forecasting results
- Forecasts/: folder of data files used as an input for generating forecasting results

## Instructions
The estimation and forecasting analyses can be run individually, in any order. Remember to specify the correct path to the local version of the repository.

## Data availability
For the empirical analyses, we have constructed the time series of daily RC matrices with the corresponding daily and intraday returns-based decompositions by using a high-frequency data set for the SPDR S&P500 and five stocks of the banking sector, i.e., Bank of America Corp. (BAC), Citigroup Inc. (C), Goldman Sachs Group Inc. (GS), JPMorgan Chase & Co. (JPM), and Wells Fargo & Co. (WFC). The data provided to us by the AlgoSeek company (30 Wall Street, 10005 NY, USA) are the prices of the assets observed every minute during the trading period (9:30-16:00) and compiled from the trades that occurred in the major US exchanges and marketplaces. Given we are not free to distribute the AlgoSeek data, we provide the corresponding time series of computed daily realized (semi-)covariance matrices based on the 5-minute log-return vectors of the trading period of the day. The sample period is January 3, 2012 - December 31, 2021, i.e., 2517 observations.

## References
1. Mafalda (2025). MFE Toolbox - Kevin Sheppard (https://www.mathworks.com/matlabcentral/fileexchange/170381-mfe-toolbox-kevin-sheppard), MATLAB Central File Exchange.



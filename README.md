# Bayesian t-Hidden Markov Model for Financial Time Series

## Overview
This repository contains the implementation of a Bayesian Student's t-Hidden Markov Model (t-HMM) applied to a portfolio of tech assets and major indices (NVDA, MSFT, META, AAPL, GOOG, S&P 500, NASDAQ). The model is designed for regime-switching detection and quantitative risk forecasting.

## Key Features
* **Custom MCMC Sampler:** Built using `NIMBLE` in R, featuring a custom Random Walk Metropolis-Hastings sampler on a logarithmic scale for the degrees of freedom ($\nu$).
* **Label-Switching Correction:** Algorithmic correction based on market volatility (S&P 500 variance constraint) to ensure economic interpretability of the hidden states.
* **Risk Management Metrics:** Out-of-sample Bayesian forecasting (10-days ahead) with explicit computation of **Value at Risk (VaR)** and **Expected Shortfall (ES)** at a 95% confidence level.

## Tech Stack
* **Language:** R
* **Core Libraries:** `nimble`, `coda`, `quantmod`, `MASS`, `mvtnorm`, `ggplot2`

## Use Case
This framework can be used by risk management teams to stress-test portfolios under different market regimes (e.g., low-volatility bull markets vs. high-volatility crisis periods), providing more robust tail-risk estimates compared to standard Gaussian models.

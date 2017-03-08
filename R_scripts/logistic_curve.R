# logistic_curve.R
# The logistic curve model from the PMAnalyzer tool.
# Author: Daniel A. Cuevas
# Last updated: 06 March 2017

logistic_fit <- function(y0, lag.time, max.growth.rate, asymptote, time) {
    y0 + ((asymptote - y0) / (1 + exp(((4 * max.growth.rate / asymptote) * (lag.time - time)) + 2)))
}

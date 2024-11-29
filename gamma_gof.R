gamma_gof <- function(model, p) {
  y = dmd_apnoea$ahi
  mu = predict(model, newdata = dmd_apnoea, type = "response")
  weights <- rep.int(1, length(y))
  wtdmu = sum(weights * y)/sum(weights)
  dev.resids = function(y, mu, wt) {
    y1 = y + (y == 0)/(10*exp(1))
    theta <- ( y1^(1-p) - mu^(1-p) ) / (1-p)
    if (p == 2) {
      kappa <- log(y1/mu)
      theta <- (y - mu)/mu
    } else {
      kappa <- ( y1^(2-p) - mu^(2-p) ) / (2-p)
      theta <- y1 * ( y1^(1-p) - mu^(1-p) ) / (1-p)
    }
    2 * wt * (theta - kappa)
    
  }  # Gamma(link = "inverse")$dev.resids is WRONG

  nulldev = sum(dev.resids(y, wtdmu, weights))

  dev = sum(dev.resids(y, mu, weights))

  disp = dev/sum(model$prior.weights)
  p = model$rank
  AIC = -2 * sum(dgamma(ifelse(y==0, 1, y), 1/disp, scale = mu * disp, log = TRUE)) + 2 + 2*p
  # sourcode compute dev -> AIC -> logLik

  pseudo_r2 <- 1 - dev/nulldev

  result = list(
    NullDeviance = nulldev,
    ResidualDeviance = dev,
    AIC = AIC,
    pseudoR2 = pseudo_r2
  )

  return(result)

  # dispersion = sum(model$weights * model$residuals^2)/ model$df.residual
  # dispersion
}

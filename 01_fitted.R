
source("R/00_tb_data.R")
ey <- 2023
pnew <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = 2013:ey))))
pret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = 2013:ey))))

# fitted data -------------------------------------------------------------

par0  <- c(
  prop_fast = 0.1,
  beta = 2,  
  w_fast = 0.0826,
  # w_slow = 0.0006,
  sc = 0.15, 
  cfr = 0.03, 
  RR = 3 
)
lower <- c(
  prop_fast = 0.01,
  beta = 0,  
  w_fast = 0.0100,
  # w_slow = 0.0001, 
  sc = 0.00, 
  cfr = 0.01, 
  RR = 1 
)
upper <- c(
  prop_fast = 0.15,
  beta = 10, 
  w_fast = 0.2000,
  # w_slow = 0.0012, 
  sc = 0.30, 
  cfr = 0.10, 
  RR = 10
)

target_cache <- prepare_target_cache(target_df, start_year = 2013, ey = ey)

# (1) ModFit to Find Start Point
fit <- modFit(f = objective_weighted, p = par0, lower = lower, upper = upper, method = "L-BFGS-B",
              fixed_parms = fixed_parms, ey = ey, state = state, pnew = pnew, pret = pret, target_cache = target_cache)
par_best_fit <- coef(fit); par_best_fit
var0 <- summary(fit)$modVariance
sim0 <- sim_tb_model(fixed_parms = fixed_parms, parms = par_best_fit, ey = ey, state = state, pnew = pnew, pret = pret)

name <- "Incident TB"
time <- 2013:2023
eval <- sim0$e_inc_num
tt   <- target_df[target_df$target == name,]$year
tval <- target_df[target_df$target == name,]$who
tval_lower <- target_df[target_df$target == name,]$lower
tval_upper <- target_df[target_df$target == name,]$upper
ylim <- c(min(eval, tval_lower)*0.9, max(eval, tval_upper)*1.1)
plot(x = time, y = eval, xlab = "Year", ylab = name, ylim = ylim, col = "blue", type = "l")
points(x = tt, y = tval, col = "black")
lines(tt, tval_lower, col = "black", lty = 2)
lines(tt, tval_upper, col = "black", lty = 2)

name <- "Deaths from TB"
time <- 2013:2023
eval <- sim0$e_mort_num
tt   <- target_df[target_df$target == name,]$year
tval <- target_df[target_df$target == name,]$who
tval_lower <- target_df[target_df$target == name,]$lower
tval_upper <- target_df[target_df$target == name,]$upper
ylim <- c(min(eval, tval_lower)*0.9, max(eval, tval_upper)*1.1)
plot(x = time, y = eval, xlab = "Year", ylab = name, ylim = ylim, col = "blue", type = "l")
points(x = tt, y = tval, col = "black")
lines(tt, tval_lower, col = "black", lty = 2)
lines(tt, tval_upper, col = "black", lty = 2)

name <- "Incident RRTB"
time <- 2013:2023
eval <- sim0$e_inc_rr_num
tt   <- target_df[target_df$target == name,]$year
tval <- target_df[target_df$target == name,]$who
tval_lower <- target_df[target_df$target == name,]$lower
tval_upper <- target_df[target_df$target == name,]$upper
ylim <- c(min(eval, tval_lower)*0.9, max(eval, tval_upper)*1.1)
plot(x = time, y = eval, xlab = "Year", ylab = name, ylim = ylim, col = "blue", type = "l")
points(x = tt, y = tval, col = "black")
lines(tt, tval_lower, col = "black", lty = 2)
lines(tt, tval_upper, col = "black", lty = 2)

# (2) MCMC fitted
set.seed(1234) # replication
niter <- 2e4
mcmcres <- modMCMC(
  f = objective_weighted, 
  p = par_best_fit, 
  lower = lower, 
  upper = upper,
  niter = niter, 
  var0 = var0,
  burninlength = 5e3, 
  updatecov = 200,
  verbose = 2000,
  ntrydr = 2,
  # paramters for sim_tb_model
  fixed_parms = fixed_parms, 
  ey = ey, state = state, 
  pnew = pnew, 
  pret = pret, 
  target_cache = target_cache
)

save(mcmcres, target_df, file = "testresult/mcmc.RData")

# (3) Plot posterior distribution
par(mfrow = c(2,3))
# plot.new();text(0.5, 0.5, "Posterior distributions", cex = 1.2)
hist(mcmcres$pars[, "prop_fast"], main = "prop_fast", xlab = "prop_fast", breaks = 30)
hist(mcmcres$pars[, "beta"], main = "beta", xlab = "beta", breaks = 30)
# hist(mcmcres$pars[, "g"], main = "g", xlab = "g", breaks = 30)
hist(mcmcres$pars[, "w_fast"], main = "w_fast", xlab = "w_fast", breaks = 30)
# hist(mcmcres$pars[, "w_slow"], main = "w_slow", xlab = "w_slow", breaks = 30)
hist(mcmcres$pars[, "sc"], main = "sc", xlab = "sc", breaks = 30)
hist(mcmcres$pars[, "cfr"], main = "cfr", xlab = "cfr", breaks = 30)
hist(mcmcres$pars[, "RR"], main = "RR", xlab = "RR", breaks = 30)
par(mfrow = c(1,1))

# show the results --------------------------------------------------------

bind_tb <- cbind(
  apply(mcmcres$pars, 2, median),
  apply(mcmcres$pars, 2, quantile, probs = 0.025),
  apply(mcmcres$pars, 2, quantile, probs = 0.975)
)
colnames(bind_tb) <- c("Median", "Q2.5", "Q97.5")

bind_tb





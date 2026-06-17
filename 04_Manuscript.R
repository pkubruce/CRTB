
source("R/00_tb_data.R")

# load data ---------------------------------------------------------------

load("testresult/mcmc.RData")

df_burd <- as.data.table(df_burd)
inc.tg.30 <- df_burd[iso3=="CHN"&year==2015,.(e_inc_100k, e_inc_100k_lo, e_inc_100k_hi)]*(1-0.90)
inc.tg.35 <- df_burd[iso3=="CHN"&year==2015,.(e_inc_100k, e_inc_100k_lo, e_inc_100k_hi)]*(1-0.95)
dth.tg.30 <- df_burd[iso3=="CHN"&year==2015,.(e_mort_num, e_mort_num_lo, e_mort_num_hi)]*(1-0.80)
dth.tg.35 <- df_burd[iso3=="CHN"&year==2015,.(e_mort_num, e_mort_num_lo, e_mort_num_hi)]*(1-0.90)

# check fitted ------------------------------------------------------------

par(mfrow = c(1,1))

plot(mcmcres, Full = TRUE)

pairs(mcmcres, nsample = 1000)

bind_tb <- cbind(
  apply(mcmcres$pars, 2, median),
  apply(mcmcres$pars, 2, quantile, probs = 0.025),
  apply(mcmcres$pars, 2, quantile, probs = 0.975)
)
colnames(bind_tb) <- c("Median", "Q2.5", "Q97.5")

bind_tb

fwrite(bind_tb, file = "manuscript/posterior_distrib.csv", row.names = TRUE)

# Fitted Results ----------------------------------------------------------
### Baseline and predict to 2035
set.seed(1234)
nsamp <- 100
ey <- 2035
pnew <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = 2013:ey))))
pret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = 2013:ey))))
idx <- sample(1:nrow(mcmcres$pars), nsamp)
slist <- list()
for (i in 1:nsamp) {
  s <- sim_tb_model(
    fixed_parms = fixed_parms, 
    parms = mcmcres$pars[idx[i], ], 
    ey = ey, 
    state = state, 
    pnew = pnew, 
    pret = pret
  )
  s <- integerize_sim_raw_safe(s)
  s$iter <- i
  slist[[i]] <- s
}
s_df <- do.call(rbind, slist)

### The observational data
# to be data.table
s_df <- as.data.table(s_df)
df_inc <- as.data.table(df_inc)
df_burd <- as.data.table(df_burd)
df_inc_rr <- as.data.table(df_inc_rr)
# merge the observation data table
m.obs <- merge(
  df_burd[iso3 == "CHN" & year>=2013 & year<=2023, .(
    year, e_pop_num, 
    e_inc_num, e_inc_100k, 
    e_mort_num, e_mort_100k, 
    e_inc_num_lo, e_inc_num_hi, 
    e_mort_num_lo, e_mort_num_hi, 
    e_inc_100k_lo, e_inc_100k_hi, 
    e_mort_100k_lo, e_mort_100k_hi
  )], 
  df_inc_rr[,.(year, e_inc_rr_num, e_inc_rr_num_lo, e_inc_rr_num_hi)], 
  by = "year",
  all.x = TRUE
)
m.obs$e_inc_rr_100k <- m.obs$e_inc_rr_num/m.obs$e_pop_num*1e5
m.obs$e_inc_rr_100k_lo <- m.obs$e_inc_rr_num_lo/m.obs$e_pop_num*1e5
m.obs$e_inc_rr_100k_hi <- m.obs$e_inc_rr_num_hi/m.obs$e_pop_num*1e5
m.obs <- rbind(
  m.obs[,.(year, group = "inc", 
           obs.m = e_inc_num, obs.lo = e_inc_num_lo, obs.hi = e_inc_num_hi, 
           obs.rate.m = e_inc_100k, 
           obs.rate.lo = e_inc_100k_lo, 
           obs.rate.hi = e_inc_100k_hi)],
  m.obs[,.(year, group = "dth", 
           obs.m = e_mort_num, obs.lo = e_mort_num_lo, obs.hi = e_mort_num_hi,
           obs.rate.m = e_mort_100k, 
           obs.rate.lo = e_mort_100k_lo, 
           obs.rate.hi = e_mort_100k_hi)],
  m.obs[,.(year, group = "rrinc", 
           obs.m = e_inc_rr_num, obs.lo = e_inc_rr_num_lo, obs.hi = e_inc_rr_num_hi,
           obs.rate.m = e_inc_rr_100k, 
           obs.rate.lo = e_inc_rr_100k_lo, 
           obs.rate.hi = e_inc_rr_100k_hi)]
)

### The expected data
### Incidence Compare
m.inc <- s_df[time<=2023,.(year = time, pred = e_inc_num, pop = N, iter)]
m.inc$pred.rate <- m.inc$pred/m.inc$pop*1e5
m.inc$group <- "inc"
### Deaths Compare
m.dth <- s_df[time<=2023,.(year = time, pred = e_mort_num, pop = N, iter)]
m.dth$pred.rate <- m.dth$pred/m.dth$pop*1e5
m.dth$group <- "dth"
### RR Incidence Compare
m.rrinc <- s_df[time<=2023,.(year = time, pred = e_inc_rr_num, pop = N, iter)]
m.rrinc <- na.omit(m.rrinc)
m.rrinc$pred.rate <- m.rrinc$pred/m.rrinc$pop*1e5
m.rrinc$group <- "rrinc"

### Bind the data and to make plot
m.bind <- rbind(m.inc, m.dth, m.rrinc)
m.pred <- m.bind |> 
  dplyr::group_by(year, group) |> 
  dplyr::summarise(
    # Number
    pred.m = median(pred), 
    pred.lo = quantile(pred, 0.025), 
    pred.hi = quantile(pred, 0.975), 
    # Rate
    pred.rate.m = median(pred.rate), 
    pred.rate.lo = quantile(pred.rate, 0.025), 
    pred.rate.hi = quantile(pred.rate, 0.975)
  )

fwrite(m.pred, file = "tmp/model.predict.csv", row.names = FALSE)
fwrite(m.obs, file = "tmp/observation.csv", row.names = FALSE)

prop.dt <- dplyr::left_join(m.bind, m.obs, by = c("year", "group")) |> 
  na.omit() |> 
  dplyr::mutate(
    i = as.numeric(pred>obs.lo & pred<obs.hi),
    i.rate = as.numeric(pred.rate>obs.rate.lo & pred.rate<obs.rate.hi)
  ) |> 
  dplyr::group_by(group) |> 
  dplyr::summarise(
    n = n(), i = sum(i), i.rate = sum(i.rate)
  ) |> 
  dplyr::mutate(
    per.i = i/n*100,
    per.i.rate = i.rate/n*100
  )

mape.dt <- dplyr::left_join(m.pred, m.obs, by = c("year", "group")) |> 
  na.omit() |> 
  dplyr::mutate(ape = pred.m/obs.m-1, ape.rate = pred.rate.m/obs.rate.m-1) |> 
  dplyr::group_by(group) |> 
  dplyr::summarise(mape = mean(ape*100), mape.rate = mean(ape.rate*100))

write.csv(prop.dt, file = "manuscript/0_prop.csv", row.names = FALSE)
write.csv(mape.dt, file = "manuscript/0_mape.csv", row.names = FALSE)


# Intervention  -----------------------------------------------------------

rm(s_df)

baseey <- 2026 # Not the start year, but the start year - 1

init_comp_name <- c("B", "M", "S", "Efast", "Eslow", "I", "RH", "RL", "N")

# @ (0) -------------------------------------------------------------------

# Status Quo
set.seed(1234)
ey <- 2035
nsamp <- 100
idx <- sample(1:nrow(mcmcres$pars), nsamp)
slist <- list()
for (i in 1:nsamp) {
  s <- sim_tb_model(
    fixed_parms = fixed_parms, 
    parms = mcmcres$pars[idx[i], ], 
    ey = ey, 
    state = state, 
    pnew = logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = 2013:ey)))), 
    pret = logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = 2013:ey))))
  )
  s <- integerize_sim_raw_safe(s)
  s$iter <- i
  slist[[i]] <- s
}
s_df <- as.data.table(do.call(rbind, slist))

# @ (1) -------------------------------------------------------------------

# annual intervention for the recurrence high‑risk group, whereby 80% of this group would receive a vaccine that reduces their tuberculosis recurrence rate (rate ratio = 0.50, as assumed in previous studies)

schedule <- data.frame(
  start_year = c(2026),
  end_year   = c(2035),
  TPT        = c(FALSE),
  TPT_Effect = c(0.33),
  TPT_Scale  = c(0.00),
  RECH       = c(TRUE),
  RECH_Effect = c(0.50),
  RECH_Scale  = c(0.80),
  RECL       = c(FALSE),
  RECL_Effect = c(0.50),
  RECL_Scale  = c(0.00)
)
reslist <- list()
for(i in 1:nsamp){
  si <- subset(s_df, iter==i & time <=baseey); si <- as.data.frame(si)
  
  res <- sim_intervention_sequence(
    fixed_parms = fixed_parms,
    parms = mcmcres$pars[idx[i], ],
    init_state = unlist(subset(si, time == baseey)[,init_comp_name]),
    base_year = baseey,
    intervention_schedule = schedule,
    pnew = logit2p(predict(fit_new, newdata = data.frame(year = baseey:2035))),
    pret = logit2p(predict(fit_ret, newdata = data.frame(year = baseey:2035)))
  )
  
  res <- as.data.frame(rbind_fill(si, res))
  
  res$iter <- i
  
  reslist[[i]] <- res
}
intv_df_v1 <- as.data.table(do.call(rbind, reslist)) 
intv_df_v1$group <- "Intervention for the recurrence (high risk)"

# @ (2) -------------------------------------------------------------------

# vaccination for the recurrence low‑risk group, also with 80% coverage and the same assumed effect (rate ratio = 0.50)

rm(schedule); rm(reslist)

schedule <- data.frame(
  start_year = c(2026),
  end_year   = c(2035),
  TPT        = c(FALSE),
  TPT_Effect = c(0.33),
  TPT_Scale  = c(0.00),
  RECH       = c(FALSE),
  RECH_Effect = c(0.50),
  RECH_Scale  = c(0.00),
  RECL       = c(TRUE),
  RECL_Effect = c(0.50),
  RECL_Scale  = c(0.80)
)
reslist <- list()
for(i in 1:nsamp){
  si <- subset(s_df, iter==i & time <=baseey); si <- as.data.frame(si)
  
  res <- sim_intervention_sequence(
    fixed_parms = fixed_parms,
    parms = mcmcres$pars[idx[i], ],
    init_state = unlist(subset(si, time == baseey)[,init_comp_name]),
    base_year = baseey,
    intervention_schedule = schedule,
    pnew = logit2p(predict(fit_new, newdata = data.frame(year = baseey:2035))),
    pret = logit2p(predict(fit_ret, newdata = data.frame(year = baseey:2035)))
  )
  
  res <- as.data.frame(rbind_fill(si, res))
  
  res$iter <- i
  
  reslist[[i]] <- res
}
intv_df_v2 <- as.data.table(do.call(rbind, reslist)) 
intv_df_v2$group <- "Intervention for the recurrence (low risk)"


# @ (3) -------------------------------------------------------------------

# intervention every 3 years for individuals with LTBI, assuming that 80% of them would receive tuberculosis preventive treatment (TPT) using a three‑month regimen of daily rifampicin and isoniazid (3HR; odds ratio for active tuberculosis = 0.33)

rm(schedule); rm(reslist)

schedule <- data.frame(
  start_year = c(2026,  2027,  2029,  2030,  2032,  2033),
  end_year   = c(2027,  2029,  2030,  2032,  2033,  2035),
  TPT        = c(TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE),
  TPT_Effect = c(0.33,  0.33,  0.33,  0.33,  0.33,  0.33),
  TPT_Scale  = c(0.80,  0.00,  0.80,  0.00,  0.80,  0.00),
  RECH       = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  RECH_Effect = c(0.50, 0.50,  0.50,  0.50,  0.50,  0.50),
  RECH_Scale  = c(0.00, 0.00,  0.00,  0.00,  0.00,  0.00),
  RECL       = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  RECL_Effect = c(0.50, 0.50,  0.50,  0.50,  0.50,  0.50),
  RECL_Scale  = c(0.00, 0.00,  0.00,  0.00,  0.00,  0.00)
)
reslist <- list()
for(i in 1:nsamp){
  si <- subset(s_df, iter==i & time <=baseey); si <- as.data.frame(si)
  
  res <- sim_intervention_sequence(
    fixed_parms = fixed_parms,
    parms = mcmcres$pars[idx[i], ],
    init_state = unlist(subset(si, time == baseey)[,init_comp_name]),
    base_year = baseey,
    intervention_schedule = schedule,
    pnew = logit2p(predict(fit_new, newdata = data.frame(year = baseey:2035))),
    pret = logit2p(predict(fit_ret, newdata = data.frame(year = baseey:2035)))
  )
  
  res <- as.data.frame(rbind_fill(si, res))
  
  res$iter <- i
  
  reslist[[i]] <- res
}
intv_df_v3 <- as.data.table(do.call(rbind, reslist)) 
intv_df_v3$group <- "Intervention for the LTBI (every 3 years)"
# intv_df <- intv_df_v3


# @ (4) -------------------------------------------------------------------

# the combination of scenarios (1), (2) and (3)

rm(schedule); rm(reslist)

schedule <- data.frame(
  start_year = c(2026,  2027,  2029,  2030,  2032,  2033),
  end_year   = c(2027,  2029,  2030,  2032,  2033,  2035),
  TPT        = c(TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE),
  TPT_Effect = c(0.33,  0.33,  0.33,  0.33,  0.33,  0.33),
  TPT_Scale  = c(0.80,  0.00,  0.80,  0.00,  0.80,  0.00),
  RECH       = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
  RECH_Effect = c(0.50, 0.50,  0.50,  0.50,  0.50,  0.50),
  RECH_Scale  = c(0.80, 0.80,  0.80,  0.80,  0.80,  0.80),
  RECL       = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
  RECL_Effect = c(0.50, 0.50,  0.50,  0.50,  0.50,  0.50),
  RECL_Scale  = c(0.80, 0.80,  0.80,  0.80,  0.80,  0.80)
)
reslist <- list()
for(i in 1:nsamp){
  si <- subset(s_df, iter==i & time <=baseey); si <- as.data.frame(si)
  
  res <- sim_intervention_sequence(
    fixed_parms = fixed_parms,
    parms = mcmcres$pars[idx[i], ],
    init_state = unlist(subset(si, time == baseey)[,init_comp_name]),
    base_year = baseey,
    intervention_schedule = schedule,
    pnew = logit2p(predict(fit_new, newdata = data.frame(year = baseey:2035))),
    pret = logit2p(predict(fit_ret, newdata = data.frame(year = baseey:2035)))
  )
  
  res <- as.data.frame(rbind_fill(si, res))
  
  res$iter <- i
  
  reslist[[i]] <- res
}
intv_df_v4 <- as.data.table(do.call(rbind, reslist)) 
intv_df_v4$group <- "Combination"
# intv_df <- intv_df_v4

# intergrated data --------------------------------------------------------

# to calculate the rate
s_df_agg <- s_df[,.(  ## Baseline 
  # Incidence
  e_inc_num_me = quantile(e_inc_num, 0.500),
  e_inc_num_lo = quantile(e_inc_num, 0.025),
  e_inc_num_hi = quantile(e_inc_num, 0.975),
  e_inc_100k_me = quantile(e_inc_num/N*1e5, 0.500),
  e_inc_100k_lo = quantile(e_inc_num/N*1e5, 0.025),
  e_inc_100k_hi = quantile(e_inc_num/N*1e5, 0.975),
  # Mortality
  e_mort_num_me = quantile(e_mort_num, 0.500),
  e_mort_num_lo = quantile(e_mort_num, 0.025),
  e_mort_num_hi = quantile(e_mort_num, 0.975),
  e_mort_100k_me = quantile(e_mort_num/N*1e5, 0.500),
  e_mort_100k_lo = quantile(e_mort_num/N*1e5, 0.025),
  e_mort_100k_hi = quantile(e_mort_num/N*1e5, 0.975),
  # RRTB Incidence
  e_inc_rr_num_me = quantile(e_inc_rr_num, 0.500),
  e_inc_rr_num_lo = quantile(e_inc_rr_num, 0.025),
  e_inc_rr_num_hi = quantile(e_inc_rr_num, 0.975),
  e_inc_rr_100k_me = quantile(e_inc_rr_num/N*1e5, 0.500),
  e_inc_rr_100k_lo = quantile(e_inc_rr_num/N*1e5, 0.025),
  e_inc_rr_100k_hi = quantile(e_inc_rr_num/N*1e5, 0.975)
), by = "time"]

intv_df_agg <- rbind(intv_df_v1, intv_df_v2, intv_df_v3, intv_df_v4)[,.( ## Intervention 
  # Incidence
  e_inc_num_me = quantile(e_inc_num, 0.500),
  e_inc_num_lo = quantile(e_inc_num, 0.025),
  e_inc_num_hi = quantile(e_inc_num, 0.975),
  e_inc_100k_me = quantile(e_inc_num/N*1e5, 0.500),
  e_inc_100k_lo = quantile(e_inc_num/N*1e5, 0.025),
  e_inc_100k_hi = quantile(e_inc_num/N*1e5, 0.975),
  # Mortality
  e_mort_num_me = quantile(e_mort_num, 0.500),
  e_mort_num_lo = quantile(e_mort_num, 0.025),
  e_mort_num_hi = quantile(e_mort_num, 0.975),
  e_mort_100k_me = quantile(e_mort_num/N*1e5, 0.500),
  e_mort_100k_lo = quantile(e_mort_num/N*1e5, 0.025),
  e_mort_100k_hi = quantile(e_mort_num/N*1e5, 0.975),
  # RRTB Incidence
  e_inc_rr_num_me = quantile(e_inc_rr_num, 0.500),
  e_inc_rr_num_lo = quantile(e_inc_rr_num, 0.025),
  e_inc_rr_num_hi = quantile(e_inc_rr_num, 0.975),
  e_inc_rr_100k_me = quantile(e_inc_rr_num/N*1e5, 0.500),
  e_inc_rr_100k_lo = quantile(e_inc_rr_num/N*1e5, 0.025),
  e_inc_rr_100k_hi = quantile(e_inc_rr_num/N*1e5, 0.975)
), by = c("time", "group")]


# To Plot the Intervention Scenario ---------------------------------------

options(scipen = 999) # No Scientific Data

# to plot

### test data
# ind <- "e_inc_num"
# sy <- 2026
# ey <- 2035
# ylab <- "Number"
# main_name <- "TB Incident Cases"
# basedf <- s_df_agg[time>=2026]
# intvdf <- intv_df_v1_agg[time>=2026]
plot_comp <- function(
    ind = NULL, sy = 2026, ey = 2035, 
    xlab = "Year", ylab = NULL, main_name = NULL,
    basedf = NULL, intvdf = NULL,
    ylim = NULL, yseq = NULL, ylab_line = 3.5
){
  ind_me <- paste0(ind, "_me"); ind_lo <- paste0(ind, "_lo"); ind_hi <- paste0(ind, "_hi")
  if(is.null(ylim)){ylim = c(0, max(basedf[[ind_hi]], intvdf[[ind_hi]])*1.2)}
  else{ylim = ylim}
  
  plot(
    NA, xlim = c(sy, ey), ylim = ylim, xlab = "", ylab = "",
    main = main_name, xaxt = "n", yaxt = "n"
  )
  axis(side = 1, at = seq(sy, ey, by = 1), labels = seq(sy, ey, by = 1), las = 2)
  axis(side = 2, at = yseq, las = 2)
  
  mtext(xlab, side = 1, line = ylab_line, cex = 0.85)
  mtext(ylab,   side = 2, line = ylab_line, cex = 0.85)
  # Reference to color regimen of Nim (Lancet Global Health, 2019)
  base_line <- "#3F6F9F"   
  base_fill <- adjustcolor("#9DB6D3", alpha.f = 0.30)
  
  intv_line <- "#C05A6A"   
  intv_fill <- adjustcolor("#E9A1AA", alpha.f = 0.30)
  
  # Baseline uncertainty interval
  polygon(
    c(basedf[["time"]], rev(basedf[["time"]])),
    c(basedf[[ind_lo]], rev(basedf[[ind_hi]])),
    col = base_fill,
    border = NA
  )
  
  # Intervention uncertainty interval
  polygon(
    c(intvdf[["time"]], rev(intvdf[["time"]])),
    c(intvdf[[ind_lo]], rev(intvdf[[ind_hi]])),
    col = intv_fill,
    border = NA
  )
  
  # Median Value
  lines(basedf[["time"]], basedf[[ind_me]], col = base_line, lwd = 2)
  
  lines(intvdf[["time"]], intvdf[[ind_me]], col = intv_line, lwd = 2)
  
  legend(
    "topright",
    legend = c("Baseline scenario", unique(intvdf$group)),
    col = c(base_line, intv_line),
    lwd = 2,
    bty = "n"
  )
}

### test function
# plot_comp(
#   ind = "e_inc_num", sy = 2026, ey = 2035,
#   ylab = "Number", main_name = "TB Incident Cases",
#   basedf = s_df_agg[time>=2026], 
#   intvdf = intv_df_agg[time>=2026&group == "Combination"],
# )

intv_seq <- c(
  "Intervention for the recurrence (high risk)",
  "Intervention for the recurrence (low risk)",
  "Intervention for the LTBI (every 3 years)",
  "Combination"
)

pdf(file = "manuscript/1_Intervention.pdf", width = 16, height = 7)

par(
  mfrow = c(2, 4),
  mar = c(4.5, 5.3, 2.8, 0.8),   # The second number control left margin
  mgp = c(2.2, 0.6, 0),          # The second number control left line
  tcl = -0.25,
  cex.axis = 0.85,
  cex.lab = 0.9,
  cex.main = 0.95
)

for(i in 1:4){
  # TB Incidence
  plot_comp(
    ind = "e_inc_100k", sy = 2026, ey = 2035,
    ylab = "Rate (per 100,000)", main_name = "TB Incidence",
    basedf = s_df_agg[time>=2026], 
    intvdf = intv_df_agg[time>=2026&group == intv_seq[i]],
    ylim = c(0,70), yseq = seq(0, 70, 10)
  )
  error_bar(x = 2030, y = inc.tg.30$e_inc_100k, lower = inc.tg.30$e_inc_100k_lo, upper = inc.tg.30$e_inc_100k_hi, width = 0.05)
  error_bar(x = 2035, y = inc.tg.35$e_inc_100k, lower = inc.tg.35$e_inc_100k_lo, upper = inc.tg.35$e_inc_100k_hi, width = 0.05)
}

for(i in 1:4){
  plot_comp(
    ind = "e_mort_num", sy = 2026, ey = 2035,
    ylab = "Number of cases", main_name = "TB Deaths",
    basedf = s_df_agg[time>=2026], 
    intvdf = intv_df_agg[time>=2026&group == intv_seq[i]],
    ylim = c(0,36000), yseq = seq(0, 36000, 6000)
  )
  error_bar(x = 2030, y = dth.tg.30$e_mort_num, lower = dth.tg.30$e_mort_num_lo, upper = dth.tg.30$e_mort_num_hi, width = 0.05)
  error_bar(x = 2035, y = dth.tg.35$e_mort_num, lower = dth.tg.35$e_mort_num_lo, upper = dth.tg.35$e_mort_num_hi, width = 0.05)
}

dev.off()


# To compare the target ---------------------------------------------------

basedf <- s_df[time==2015, .(iter, e_inc_100k = e_inc_num/N*1e5, e_mort_num)]
setnames(basedf, old = "e_inc_100k", new = "base_e_inc_100k")
setnames(basedf, old = "e_mort_num", new = "base_e_mort_num")

intvdf <- rbind(
  s_df[time %in% c(2030, 2035), .(iter, time, group = "Baseline", e_inc_100k = e_inc_num/N*1e5, e_mort_num)],
  intv_df_v1[time %in% c(2030, 2035), .(iter, time, group, e_inc_100k = e_inc_num/N*1e5, e_mort_num)],
  intv_df_v2[time %in% c(2030, 2035), .(iter, time, group, e_inc_100k = e_inc_num/N*1e5, e_mort_num)],
  intv_df_v3[time %in% c(2030, 2035), .(iter, time, group, e_inc_100k = e_inc_num/N*1e5, e_mort_num)],
  intv_df_v4[time %in% c(2030, 2035), .(iter, time, group, e_inc_100k = e_inc_num/N*1e5, e_mort_num)]
)

m.tg <- merge(basedf, intvdf, by = c("iter"))

m.tg$inc_tg <- (1 - m.tg$e_inc_100k/m.tg$base_e_inc_100k)*100
m.tg$dth_tg <- (1 - m.tg$e_mort_num/m.tg$base_e_mort_num)*100

outm.tg <- m.tg[,.(
  # Incidence
  e_inc_100k_me = quantile(e_inc_100k, 0.500),
  e_inc_100k_lo = quantile(e_inc_100k, 0.025),
  e_inc_100k_hi = quantile(e_inc_100k, 0.975),
  # Mortality
  e_mort_num_me = quantile(e_mort_num, 0.500),
  e_mort_num_lo = quantile(e_mort_num, 0.025),
  e_mort_num_hi = quantile(e_mort_num, 0.975),
  
  # Reduction # Incidence
  inc_tg_me = quantile(inc_tg, 0.500),
  inc_tg_lo = quantile(inc_tg, 0.025),
  inc_tg_hi = quantile(inc_tg, 0.975),
  
  # Reduction # Deaths
  dth_tg_me = quantile(dth_tg, 0.500),
  dth_tg_lo = quantile(dth_tg, 0.025),
  dth_tg_hi = quantile(dth_tg, 0.975)
), by = c("time", "group")]

openxlsx::write.xlsx(outm.tg, file = "manuscript/1_Intervention.xlsx")

inc.tg.30
#    e_inc_100k e_inc_100k_lo e_inc_100k_hi
#         <num>         <num>         <num>
# 1:        6.5           5.3           7.8
inc.tg.35
#    e_inc_100k e_inc_100k_lo e_inc_100k_hi
#         <num>         <num>         <num>
# 1:       3.25          2.65           3.9
dth.tg.30
#    e_mort_num e_mort_num_lo e_mort_num_hi
#         <num>         <num>         <num>
# 1:       8400          7400          9400
dth.tg.35
#    e_mort_num e_mort_num_lo e_mort_num_hi
#         <num>         <num>         <num>
# 1:       4200          3700          4700


# To calculate Rate Ratio -------------------------------------------------

bdf <- intvdf[group == "Baseline" & time == 2035]
setnames(bdf, old = "e_inc_100k", new = "base_e_inc_100k")
setnames(bdf, old = "e_mort_num", new = "base_e_mort_num")
idf <- intvdf[group != "Baseline" & time == 2035]

merge(bdf, idf, by = c("time", "iter"))[,.(
  rateratio.inc.m = quantile(e_inc_100k/base_e_inc_100k, 0.500),
  rateratio.inc.lo = quantile(e_inc_100k/base_e_inc_100k, 0.025),
  rateratio.inc.hi = quantile(e_inc_100k/base_e_inc_100k, 0.975)
), by = c("group.y")]








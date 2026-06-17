
source("R/00_tb_data.R")

# load data ---------------------------------------------------------------

load("testresult/mcmc.RData")

# 利用一组参数 做基线 ---------------------------------------------------

baseey <- 2026

pnew <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = 2013:baseey))))
pret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = 2013:baseey))))

s_df <- sim_tb_model(
  fixed_parms = fixed_parms, 
  parms = mcmcres$pars[2000, ], 
  ey = baseey, 
  state = state, 
  pnew = pnew, 
  pret = pret
)

eys <- baseey
eye <- 2035

snoinv_df <- sim_intv(
  fixed_parms = fixed_parms, 
  parms  = mcmcres$pars[2000, ], 
  state = s_df$raw[s_df$raw[,"time"] == baseey,-1], 
  eys = eys, 
  eye = eye, 
  pnew = logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = eys:eye)))), 
  pret = logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = eys:eye)))),
  TPT = FALSE, TPT_Effect = 0.33, TPT_Scale = 0.80,
  RECH = FALSE, RECH_Effect = 1.00, RECH_Scale = 0.00,
  rtol = 1e-6, atol = 1e-6
)

# 非连续干预  --------------------------------------------------

eys <- 2026
eye <- 2028

sinv_df_before <- sim_intv(
  fixed_parms = fixed_parms, 
  parms  = mcmcres$pars[2000, ], 
  state = s_df$raw[s_df$raw[,"time"] == baseey,-1], 
  eys = eys, 
  eye = eye, 
  pnew = logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = eys:eye)))), 
  pret = logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = eys:eye)))),
  TPT = FALSE, TPT_Effect = 0.00, TPT_Scale = 0.00,
  RECH = TRUE,  RECH_Effect = 0.10, RECH_Scale = 0.90,  
  RECL = TRUE,  RECL_Effect = 0.10, RECL_Scale = 0.90,  
  rtol = 1e-6, atol = 1e-6
)

eys <- 2028
eye <- 2030

sinv_df_med <- sim_intv(
  fixed_parms = fixed_parms, 
  parms  = mcmcres$pars[2000, ], 
  state = sinv_df_before$raw[sinv_df_before$raw[,"time"] == 2028,-1], 
  eys = eys, 
  eye = eye, 
  pnew = logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = eys:eye)))), 
  pret = logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = eys:eye)))),
  # TPT = TRUE, TPT_Effect = 0.33, TPT_Scale = 0.80,
  TPT = FALSE, TPT_Effect = 0.00, TPT_Scale = 0.00,
  RECH = FALSE, RECH_Effect = 1.00, RECH_Scale = 0.00,  
  rtol = 1e-6, atol = 1e-6
)

eys <- 2030
eye <- 2035

sinv_df_after <- sim_intv(
  fixed_parms = fixed_parms, 
  parms  = mcmcres$pars[2000, ], 
  state = sinv_df_med$raw[sinv_df_med$raw[,"time"] == 2030,-1], 
  eys = eys, 
  eye = eye, 
  pnew = logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = eys:eye)))), 
  pret = logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = eys:eye)))),
  TPT = FALSE, TPT_Effect = 0.00, TPT_Scale = 0.00,
  RECH = FALSE, RECH_Effect = 1.00, RECH_Scale = 0.00,  
  rtol = 1e-6, atol = 1e-6
)

# bind data ---------------------------------------------------------------

sinv <- rbind(
  with(s_df, expr = {data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num)}),
  with(sinv_df_before, expr = {subset(data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num), year != 2026)}),
  with(sinv_df_med, expr = {subset(data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num), year != 2028)}),
  with(sinv_df_after, expr = {subset(data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num), year != 2030)})
)

snoinv <- rbind(
  with(s_df, expr = {data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num)}),
  with(snoinv_df, expr = {subset(data.frame(year, e_inc_num, e_mort_num, e_inc_rr_num), year != 2026)})
)

# show the results --------------------------------------------------------

ylim <- range(sinv$e_inc_num, snoinv$e_inc_num, finite = TRUE, na.rm = TRUE)

plot(x = snoinv$year, y = as.integer(snoinv$e_inc_num), ylim = ylim, xlab ="Year", ylab = "Number", type = "b")
points(x = sinv$year, y = sinv$e_inc_num, col = "red")
lines(x = sinv$year, y = sinv$e_inc_num, col = "red")

# 干预时间表 -------------------------------------------------------------------

schedule <- data.frame(
  start_year = c(2026,  2028, 2029,  2030),
  end_year   = c(2028,  2029, 2030,  2035),
  TPT        = c(FALSE, TRUE, FALSE, FALSE),
  TPT_Effect = c(0.00,  0.33, 0.00,  0.00),
  TPT_Scale  = c(0.00,  0.90, 0.00,  0.00),
  RECH       = c(TRUE,  TRUE, FALSE, FALSE),
  RECH_Effect = c(0.50,  0.50, 0.00,  0.00),
  RECH_Scale  = c(0.90,  0.90, 0.00,  0.00),
  RECL       = c(TRUE,  TRUE, FALSE, FALSE),
  RECL_Effect = c(0.50,  0.50, 0.00,  0.00),
  RECL_Scale  = c(0.90,  0.90, 0.00,  0.00)
)
## 只做一组
res <- sim_intervention_sequence(
  fixed_parms = fixed_parms,
  parms = mcmcres$pars[2000, ],
  init_state = s_df$raw[s_df$raw[,"time"] == baseey, -1],
  base_year = baseey,
  intervention_schedule = schedule,
  pnew = logit2p(predict(fit_new, newdata = data.frame(year = baseey:2035))),
  pret = logit2p(predict(fit_ret, newdata = data.frame(year = baseey:2035)))
)

res <- res
res_noinv <- snoinv[snoinv$year>=2027,]

ylim <- range(res$e_inc_num, res_noinv$e_inc_num, finite = TRUE, na.rm = TRUE)
plot(x = res_noinv$year, y = res_noinv$e_inc_num, ylim = ylim, xlab ="Year", ylab = "Number", type = "b")
points(x = res$year, y = res$e_inc_num, col = "red")
lines(x = res$year, y = res$e_inc_num, col = "red")

ylim <- range(res$e_mort_num, res_noinv$e_mort_num, finite = TRUE, na.rm = TRUE)
plot(x = res_noinv$year, y = res_noinv$e_mort_num, ylim = ylim, xlab ="Year", ylab = "Number", type = "b")
points(x = res$year, y = res$e_mort_num, col = "red")
lines(x = res$year, y = res$e_mort_num, col = "red")

ylim <- range(res$e_inc_rr_num, res_noinv$e_inc_rr_num, finite = TRUE, na.rm = TRUE)
plot(x = res_noinv$year, y = res_noinv$e_inc_rr_num, ylim = ylim, xlab ="Year", ylab = "Number", type = "b")
points(x = res$year, y = res$e_inc_rr_num, col = "red")
lines(x = res$year, y = res$e_inc_rr_num, col = "red")


# 画图 ----------------------------------------------------------------------

plot_intv <- function(
    rate = FALSE,
    baseline_df,
    intv_df,
    ind_name,
    main_name = ind_name,
    ey = NULL,
    xlab = "Year",
    ylab = "Number",
    col_baseline = "#134b87",
    col_intv = "red",
    alpha_traj = 0.2,
    shade_start = 2013,
    shade_end = 2023,
    legend_pos = "topright",
    ...
  ){
  
  if(rate){
    y_mat <- tapply(baseline_df[[ind_name]]/baseline_df[["N"]]*1e5, list(baseline_df$time, baseline_df$iter), identity)
    i_y_mat <- tapply(intv_df[[ind_name]]/intv_df[["N"]]*1e5, list(intv_df$time, intv_df$iter), identity)
    ylab <- "Rate"
  } else{
    i_y_mat <- tapply(intv_df[[ind_name]], list(intv_df$time, intv_df$iter), identity)
    y_mat <- tapply(baseline_df[[ind_name]], list(baseline_df$time, baseline_df$iter), identity)
  }
  
  # Prepare baseline matrix (time x iter)
  time_vec <- as.numeric(rownames(y_mat))
  median_val <- apply(y_mat, 1, median, na.rm = TRUE)
  
  # Prepare intervention matrix
  
  i_time_vec <- as.numeric(rownames(i_y_mat))
  i_median_val <- apply(i_y_mat, 1, median, na.rm = TRUE)
  
  # Determine y-axis limits (include all trajectories)
  xlim <- c(min(time_vec, i_time_vec, na.rm = TRUE), ey)
  y_lim <- range(y_mat, i_y_mat, finite = TRUE, na.rm = TRUE)
  
  # Create empty plot
  plot(
    NA,
    xlim = xlim,
    ylim = y_lim,
    xlab = xlab,
    ylab = ylab,
    main = main_name,
    xaxt = "n",
    ...
  )
  
  # Custom x-axis (every 2 years)
  at_vals <- seq(floor(xlim[1]), ceiling(xlim[2]), by = 1)
  axis(side = 1, at = at_vals, labels = at_vals, las = 2)
  
  # Baseline posterior trajectories (semi-transparent)
  matplot(
    time_vec,
    y_mat,
    type = "l",
    lty = 1,
    col = adjustcolor(col_baseline, alpha.f = alpha_traj),
    add = TRUE
  )
  # Baseline median line and points
  lines(time_vec, median_val, col = col_baseline, lwd = 2)
  points(time_vec, median_val, col = col_baseline, pch = 16, cex = 0.8)
  
  # Intervention posterior trajectories
  matplot(
    i_time_vec,
    i_y_mat,
    type = "l",
    lty = 1,
    col = adjustcolor(col_intv, alpha.f = alpha_traj),
    add = TRUE
  )
  # Intervention median line and points
  lines(i_time_vec, i_median_val, col = col_intv, lwd = 2)
  points(i_time_vec, i_median_val, col = col_intv, pch = 17, cex = 0.8)
  
  # Legend
  legend(
    legend_pos,
    legend = c("Baseline", "Intervention"),
    col = c(col_baseline, col_intv),
    pch = c(16, 17),
    lty = c(1, 1),
    lwd = 2,
    bty = "n",
    cex = 0.8
  )
}


# 基线情境 --------------------------------------------------------------------

set.seed(1234)
baseey <- 2026 # Not the start year, but the start year - 1
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

s_df <- do.call(rbind, slist)

# 最理想干预情境 -----------------------------------------------------------

schedule <- data.frame(
  start_year = c(2026,  2028, 2029,  2030),
  end_year   = c(2028,  2029, 2030,  2035),
  TPT        = c(FALSE, TRUE, FALSE, FALSE),
  TPT_Effect = c(0.00,  0.33, 0.00,  0.00),
  TPT_Scale  = c(0.00,  0.90, 0.00,  0.00),
  RECH       = c(TRUE,  TRUE, FALSE, FALSE),
  RECH_Effect = c(0.50,  0.50, 0.00,  0.00),
  RECH_Scale  = c(0.90,  0.90, 0.00,  0.00),
  RECL       = c(TRUE,  TRUE, FALSE, FALSE),
  RECL_Effect = c(0.50,  0.50, 0.00,  0.00),
  RECL_Scale  = c(0.90,  0.90, 0.00,  0.00)
)

init_comp_name <- c("B", "M", "S", "Efast", "Eslow", "I", "RH", "RL", "N")

reslist <- list()

for(i in 1:nsamp){
  si <- subset(s_df, iter==i & time <=baseey)
  
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

intv_df <- do.call(rbind, reslist)

# 跑完以后画图

par(mfrow = c(2,3))

# Number

# par(mfrow = c(1,1))

plot_intv(baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_inc_num", 
          main_name = "Total Incident TB", 
          ey = 2035)

# abline(v = 2026)

plot_intv(baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_mort_num", 
          main_name = "Deaths from TB", 
          ey = 2035)

plot_intv(baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_inc_rr_num", 
          main_name = "Incident RRTB", 
          ey = 2035)

# Rate

plot_intv(rate = TRUE,
          baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_inc_num", 
          main_name = "Total Incident TB", 
          ey = 2035)

plot_intv(rate = TRUE,
          baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_mort_num", 
          main_name = "Deaths from TB", 
          ey = 2035)

plot_intv(rate = TRUE,
          baseline_df = s_df, 
          intv_df = intv_df, 
          ind_name = "e_inc_rr_num", 
          main_name = "Incident RRTB", 
          ey = 2035)

abline(v=2026)

par(mfrow = c(1,1))























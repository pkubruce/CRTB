
source("R/00_tb_data.R")

# load data ---------------------------------------------------------------

ey <- 2023
pnew <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = 2013:ey))))
pret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = 2013:ey))))

load("testresult/mcmc.RData")


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


# simulation --------------------------------------------------------------

set.seed(1234)
nsamp <- 100
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

plot_number <- function(ind_name, target_name, main_name, 
                        col_who = "red", noti = FALSE, notidata,
                        legend_pos = "topright", y_lim = NULL){
  y_mat <- tapply(s_df[[ind_name]], list(s_df$time, s_df$iter), identity)
  time_vec <- as.numeric(rownames(y_mat))
  median_val <- apply(y_mat, 1, median, na.rm = TRUE)
  obs <- target_df[target_df$target == target_name, ]
  tt <- obs$year; val <- obs$who; val_lower <- obs$lower; val_upper <- obs$upper
  
  if(noti){noti_val <- notidata$val; noti_tt <- notidata$year} else{noti_val <- NA}
  
  # adjusted ylim
  if(is.null(y_lim)){
    y_lim <- range(y_mat, median_val, val, val_lower, val_upper, noti_val, finite = TRUE, na.rm = TRUE)
  }else{
    y_lim <- y_lim
  }
  
  plot(
    NA,
    xlim = c(2013, ey),
    ylim = y_lim,
    xlab = "Year",
    ylab = "Number",
    main = main_name,
    xaxt = "n"
  )
  
  axis(
    side = 1,
    at = seq(2013, ey + 1, by = 2),
    labels = seq(2013, ey + 1, by = 2),
    las = 2
  )
  
  # # shaded area
  # usr <- par("usr")
  # rect(
  #   xleft = 2013,
  #   ybottom = usr[3],
  #   xright = 2023,
  #   ytop = usr[4],
  #   col = adjustcolor("grey80", alpha.f = 0.5),
  #   border = NA
  # )
  
  # posterior trajectories
  matplot(
    time_vec,
    y_mat,
    type = "l",
    lty = 1,
    col = adjustcolor("#c7e0ed", alpha.f = 0.5),
    add = TRUE
  )
  
  # median 
  lines(time_vec, median_val, col = "#134b87", lwd = 2)
  points(time_vec, median_val, col = "#134b87", pch = 16, cex = 0.8)
  
  # observation value
  error_bar(x = tt, y = val, lower = val_lower, upper = val_upper, width = 0.2, col = col_who)
  
  # Notification
  if(noti){
    lines(x = noti_tt, y = noti_val, col = "black", lwd = 2)
    points(x = noti_tt, y = noti_val, col = "black", pch = 16, cex = 1)
  }
  
  # Add legend
  legend_labels <- c("CRTB Model", "WHO")
  legend_colors <- c("#134b87", col_who)
  legend_pch    <- c(16, 16)
  legend_lty    <- c(1, 1)
  
  if (noti) {
    legend_labels <- c(legend_labels, "Notification")
    legend_colors <- c(legend_colors, "black")
    legend_pch    <- c(legend_pch, 16)
    legend_lty    <- c(legend_lty, 1)
  }
  
  legend(legend_pos, legend = legend_labels, col = legend_colors,
         pch = legend_pch, lty = legend_lty, bty = "n", cex = 0.8)
}


plot_rate <- function(ind_name, target_name, main_name, 
                      col_who = "red", noti = FALSE, notidata,
                      legend_pos = "topright", y_lim = NULL){
  y_mat <- tapply(s_df[[ind_name]]/s_df[["N"]]*1e5, list(s_df$time, s_df$iter), identity)
  time_vec <- as.numeric(rownames(y_mat))
  median_val <- apply(y_mat, 1, median, na.rm = TRUE)
  obs <- target_df[target_df$target == target_name, ]
  tt <- obs$year; val <- obs$who/obs$e_pop_num*1e5; val_lower <- obs$lower/obs$e_pop_num*1e5; val_upper <- obs$upper/obs$e_pop_num*1e5
  
  if(noti){noti_val <- notidata$val; noti_tt <- notidata$year} else{noti_val <- NA}
  
  # adjusted ylim
  if(is.null(y_lim)){
    y_lim <- range(y_mat, median_val, val, val_lower, val_upper, noti_val, finite = TRUE, na.rm = TRUE)
  }else{
    y_lim <- y_lim
  }

  plot(
    NA,
    xlim = c(2013, ey),
    ylim = y_lim,
    xlab = "Year",
    ylab = "Rate (per 100,000)",
    main = main_name,
    xaxt = "n"
  )
  
  axis(
    side = 1,
    at = seq(2013, ey + 1, by = 2),
    labels = seq(2013, ey + 1, by = 2),
    las = 2
  )
  
  # # shaded area
  # usr <- par("usr")
  # rect(
  #   xleft = 2013,
  #   ybottom = usr[3],
  #   xright = 2023,
  #   ytop = usr[4],
  #   col = adjustcolor("grey80", alpha.f = 0.5),
  #   border = NA
  # )
  
  # posterior trajectories
  matplot(
    time_vec,
    y_mat,
    type = "l",
    lty = 1,
    col = adjustcolor("#c7e0ed", alpha.f = 0.5),
    add = TRUE
  )
  
  # median 
  lines(time_vec, median_val, col = "#134b87", lwd = 2)
  points(time_vec, median_val, col = "#134b87", pch = 16, cex = 0.8)
  
  # observation value
  error_bar(x = tt, y = val, lower = val_lower, upper = val_upper, width = 0.2, col = col_who)
  
  # Notification
  if(noti){
    lines(x = noti_tt, y = noti_val, col = "black", lwd = 2)
    points(x = noti_tt, y = noti_val, col = "black", pch = 16, cex = 1)
  }
  
  # Add legend
  legend_labels <- c("CRTB Model", "WHO")
  legend_colors <- c("#134b87", col_who)
  legend_pch    <- c(16, 16)
  legend_lty    <- c(1, 1)
  
  if (noti) {
    legend_labels <- c(legend_labels, "Notification")
    legend_colors <- c(legend_colors, "black")
    legend_pch    <- c(legend_pch, 16)
    legend_lty    <- c(legend_lty, 1)
  }
  
  legend(legend_pos, legend = legend_labels, col = legend_colors,
         pch = legend_pch, lty = legend_lty, bty = "n", cex = 0.8)
}

noti_subset <- subset(df_noti, iso3 == "CHN" & year>=2013 & year<=2023)
noti_subset <- noti_subset[,c("year", "c_newinc")]
names(noti_subset) <- c("year", "val")
noti_subset_rate <- merge(noti_subset, subset(df_burd, iso3=="CHN")[,c("year","e_pop_num")], by = "year")
noti_subset_rate$val <- noti_subset_rate$val/noti_subset_rate$e_pop_num*1e5

notirr_subset <- subset(df_noti, iso3 == "CHN" & year>=2013 & year<=2023)
notirr_subset <- notirr_subset[,c("year", "conf_rrmdr", "conf_rr_nfqr", "conf_rr_fqr")]
notirr_subset$val <- ifelse(notirr_subset$year<2020, notirr_subset$conf_rrmdr, notirr_subset$conf_rr_nfqr + notirr_subset$conf_rr_fqr)
notirr_subset <- notirr_subset[,c("year", "val")]
notirr_subset_rate <- merge(notirr_subset, subset(df_burd, iso3=="CHN")[,c("year","e_pop_num")], by = "year")
notirr_subset_rate$val <- notirr_subset_rate$val/notirr_subset_rate$e_pop_num*1e5


par(mfrow = c(2,3))

options(scipen = 999) # No Scientific Data

plot_number(
  ind_name = "e_inc_num",
  target_name = "Incident TB",
  main_name = "TB Incident Cases",
  noti = TRUE, notidata = noti_subset, y_lim = c(0, 1500000)
)

plot_number(
  ind_name = "e_mort_num",
  target_name = "Deaths from TB",
  main_name = "TB Deaths", y_lim = c(0, 50000)
)

plot_number(
  ind_name = "e_inc_rr_num",
  target_name = "Incident RRTB",
  main_name = "RRTB Incident Cases",
  noti = TRUE, notidata = notirr_subset, y_lim = c(0, 100000)
)

plot_rate(
  ind_name = "e_inc_num",
  target_name = "Incident TB",
  main_name = "TB Incidence",
  noti = TRUE, notidata = noti_subset_rate, y_lim = c(0, 100)
)

plot_rate(
  ind_name = "e_mort_num",
  target_name = "Deaths from TB",
  main_name = "TB Mortality", y_lim = c(0, 4)
)

plot_rate(
  ind_name = "e_inc_rr_num",
  target_name = "Incident RRTB",
  main_name = "RRTB Incidence",
  noti = TRUE, notidata = notirr_subset_rate, y_lim = c(0, 8)
)

par(mfrow = c(1,1))



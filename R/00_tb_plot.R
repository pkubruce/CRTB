plot_number <- function(
    s_df, obs, ey,
    ind_name, target_name, main_name, 
    col_who = "red", noti = FALSE, notidata,
    legend_pos = "topright"){
  y_mat <- tapply(s_df[[ind_name]], list(s_df$time, s_df$iter), identity)
  time_vec <- as.numeric(rownames(y_mat))
  median_val <- apply(y_mat, 1, median, na.rm = TRUE)
  obs <- target_df[target_df$target == target_name, ]
  tt <- obs$year; val <- obs$who; val_lower <- obs$lower; val_upper <- obs$upper
  
  if(noti){noti_val <- notidata$val; noti_tt <- notidata$year} else{noti_val <- NA}
  
  y_lim <- range(y_mat, median_val, val, val_lower, val_upper, noti_val, finite = TRUE, na.rm = TRUE)
  
  plot(
    NA,
    xlim = c(2013, ey + 1),
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
  
  # shaded area
  usr <- par("usr")
  rect(
    xleft = 2013,
    ybottom = usr[3],
    xright = 2023,
    ytop = usr[4],
    col = adjustcolor("grey80", alpha.f = 0.5),
    border = NA
  )
  
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
  legend_labels <- c("CRTB Model", "WHO Estimation")
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

plot_rate <- function(
    s_df, obs, ey,
    ind_name, target_name, main_name, 
    col_who = "red", noti = FALSE, notidata,
    legend_pos = "topright"
  ){
  y_mat <- tapply(s_df[[ind_name]]/s_df[["N"]]*1e5, list(s_df$time, s_df$iter), identity)
  time_vec <- as.numeric(rownames(y_mat))
  median_val <- apply(y_mat, 1, median, na.rm = TRUE)
  obs <- target_df[target_df$target == target_name, ]
  tt <- obs$year; val <- obs$who/obs$e_pop_num*1e5; val_lower <- obs$lower/obs$e_pop_num*1e5; val_upper <- obs$upper/obs$e_pop_num*1e5
  
  if(noti){noti_val <- notidata$val; noti_tt <- notidata$year} else{noti_val <- NA}
  
  y_lim <- range(y_mat, median_val, val, val_lower, val_upper, noti_val, finite = TRUE, na.rm = TRUE)
  
  plot(
    NA,
    xlim = c(2013, ey + 1),
    ylim = y_lim,
    xlab = "Year",
    ylab = "Rate",
    main = main_name,
    xaxt = "n"
  )
  
  axis(
    side = 1,
    at = seq(2013, ey + 1, by = 2),
    labels = seq(2013, ey + 1, by = 2),
    las = 2
  )
  
  # shaded area
  usr <- par("usr")
  rect(
    xleft = 2013,
    ybottom = usr[3],
    xright = 2023,
    ytop = usr[4],
    col = adjustcolor("grey80", alpha.f = 0.5),
    border = NA
  )
  
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
  legend_labels <- c("CRTB Model", "WHO Estimation")
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
    y_mat <- tapply(baseline_df[[ind_name]], list(baseline_df$time, baseline_df$iter), identity)
    i_y_mat <- tapply(intv_df[[ind_name]], list(intv_df$time, intv_df$iter), identity)
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


plot_comp <- function(
    ind = NULL, sy = 2026, ey = 2035, 
    xlab = "Year", ylab = NULL, main_name = NULL,
    basedf = NULL, intvdf = NULL,
    ylim = NULL, yseq = NULL, ylab_line = 3.5
){
  ind_me <- paste0(ind, "_me"); ind_lo <- paste0(ind, "_lo"); ind_hi <- paste0(ind, "_hi")
  if(is.null(ylim)){
    all_values <- c(basedf[[ind_hi]], intvdf[[ind_hi]])
    all_values <- all_values[!is.na(all_values)]
    ylim = c(0, max(all_values)*1.2)
  }else{ylim = ylim}
  
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
    legend = c("Baseline scenario", "Intervention scenario"),
    col = c(base_line, intv_line),
    lwd = 2,
    bty = "n"
  )
}

# functions for mathematical calculation ----
sigmoid <- function(x){
  x/(1-x)
}

log_sigmoid <- function(x){
  log(sigmoid(x))
}

logit2p <- function(x){
  exp(x)/(1+exp(x))
}

get_coef <- function(x, idx){
  s <- summary(x)$coefficients
  
  if(s[idx,"Pr(>|t|)"]<0.05){
    res <- s[idx,"Estimate"]
  } else{
    res <- 0
  }
  return(res)
}

safe_log <- function(x) log(pmax(x, 1e-6))

# functions for making target dataframe -----

make_target_df <- function(list_df, name_target){
  res_out <- map_dfr(1:length(list_df), function(i){
    df <- list_df[[i]]
    name_df <- names(df)
    name_bound_match <- c("e_inc_num_hi","e_inc_num_lo","e_mort_num_hi","e_mort_num_lo","e_inc_rr_num_hi","e_inc_rr_num_lo")
    name_value_match <- c("e_inc_num","e_mort_num","e_inc_rr_num")
    sum_match <- sum(name_bound_match %in% name_df)
    if(sum_match == 2){
      name_bound_select <- name_bound_match[name_bound_match %in% name_df]
      name_value_select <- name_value_match[name_value_match %in% name_df]
      df <- df |> dplyr::select("year", all_of(name_value_select), all_of(name_bound_select))
      names(df) <- c("year", "val", "hi", "lo")
      res_in <- df |> dplyr::mutate(who = val, val = log(val), se = (log(hi) - log(lo))/(2*1.96))
    } else{
      # refer to Lancet Global Health (Rebecca A Clark, 2023) 
      # lower bound = val - 1.96 se = val * 0.9
      # upper bound = val + 1.96 se = val * 1.1
      # se = val (1.0-0.9) / 1.96
      # se = val (1.1-1.0) / 1.96
      name_value_select <- name_value_match[name_value_match %in% name_df]
      df <- df |> dplyr::select("year", all_of(name_value_select))
      names(df) <- c("year", "val")
      res_in <- df |> dplyr::mutate(who = val, val = log(val), se = val * 0.1 / 1.96)
    }
    res_in <- res_in |> 
      dplyr::rename("lower" = "lo", "upper" = "hi") |> 
      dplyr::mutate(target = name_target[i]) 
    
    return(res_in)
  })
  
  return(res_out)
}

# functions for update demographic data ----
# year
# val
update_demo <- function(df){
  df$logit_val <- log_sigmoid(df$val)
  fit <- lm(logit_val ~ year, data = df)
  beta <- get_coef(fit, idx = 2)
  val0 <- as.numeric(logit2p(predict(fit, newdata = data.frame(year = 2013))))
  return(list(val0 = val0, beta = beta))
}

# functions for update RRTB data ----
update_rrtb <- function(df){
  df$logit_val <- log_sigmoid(df$val/100)
  fit <- lm(logit_val ~ year, data = df)
  return(fit)
}

error_bar <- function(x, y, lower, upper, width = 0.2, col = "black", pch = 16) {
  points(x, y, col = col, pch = pch)
  
  segments(x0 = x, y0 = lower, x1 = x, y1 = upper, col = col)
  segments(x0 = x - width, y0 = lower, x1 = x + width, y1 = lower, col = col)
  segments(x0 = x - width, y0 = upper, x1 = x + width, y1 = upper, col = col)
}


# functions for combind rows ----------------------------------------------
rbind_fill <- function(..., fillval = 0) {
  dfs <- list(...)
  if (length(dfs) == 0) return(data.frame())
  
  # 获取所有数据框的列名并集
  all_names <- unique(unlist(lapply(dfs, names)))
  
  # 对每个数据框补齐缺失列
  dfs_filled <- lapply(dfs, function(df) {
    for (nm in setdiff(all_names, names(df))) {
      df[[nm]] <- fillval
    }
    # 保证列顺序一致 按 all_names 排序 也可按原顺序，这里按 all_names 排序
    df[all_names]
  })
  
  # 合并
  do.call(rbind, dfs_filled)
}












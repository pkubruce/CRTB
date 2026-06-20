# load packages -----------------------------------------------------------

library(tidyverse)

library(deSolve)

# Simple Model ------------------------------------------------------------

tb_model <- function(t, y, parms) {
  # Compartment
  B     <- y["B"]
  M     <- y["M"]
  S     <- y["S"]
  Efast <- y["Efast"]
  Eslow <- y["Eslow"]
  I     <- y["I"]
  RH    <- y["RH"]
  RL    <- y["RL"]
  N     <- y["N"]
  
  # Parameter
  ## Demographic
  b_B        <- parms["b_B"]
  b_M        <- parms["b_M"]
  
  ## TB Epidemiology
  beta       <- parms["beta"]
  g          <- parms["g"]
  w_fast     <- parms["w_fast"]
  w_slow     <- parms["w_slow"]
  sc         <- parms["sc"]
  cfr        <- parms["cfr"]
  
  f          <- parms["f"]
  eta        <- parms["eta"]
  
  ## PostTB
  rhoH       <- parms["rhoH"]
  rhoL       <- parms["rhoL"]
  toRL       <- parms["toRL"]
  MR         <- parms["MR"]
  
  # ODE Formula
  dB <- b_B * (B) * (1-B)
  dM <- b_M * (M) * (1-M)
  
  dS     <- B*N + beta*S*I/N - M*S
  dEfast <- (1-g)*beta*S*I/N - w_fast*Efast - M*Efast 
  dEslow <-     g*beta*S*I/N - w_slow*Eslow - M*Eslow
  dI     <- w_fast*Efast + w_slow*Eslow + rhoH*RH + rhoL*RL - eta*I - sc*I - (M+cfr)*I
  dRH    <- eta*I + sc*I - rhoH*RH - toRL*RH - (M+MR)*RH
  dRL    <- toRL*RH - rhoL*RL - (M+MR)*RL
  
  dN <- dS + dEfast + dEslow +dI + dRH + dRL
  
  list(c(dB, dM, dS, dEfast, dEslow, dI, dRH, dRL, dN))
}

# Precompute target data once ---------------------------------------------

prepare_target_cache <- function(obs_data, start_year = 2013, ey) {
  years <- start_year:ey
  
  get_target <- function(target_name) {
    d <- obs_data[obs_data$target == target_name, c("year", "val", "se")]
    idx <- match(d$year, years)
    ok <- !is.na(idx) & is.finite(d$val) & is.finite(d$se) & d$se > 0
    list(idx = idx[ok], val = d$val[ok], se = d$se[ok])
  }
  
  cache <- list(
    inc  = get_target("Incident TB"),
    mort = get_target("Deaths from TB"),
    rr   = get_target("Incident RRTB")
  )
  cache$n <- length(cache$inc$val) + length(cache$mort$val) + length(cache$rr$val)
  cache
}

# Simulation --------------------------------------------------------------

sim_tb_model <- function(fixed_parms, parms, state, ey, pnew, pret,
                         rtol = 1e-6, atol = 1e-6) {
  # Parameter Initial
  prop_fast <- parms[["prop_fast"]]
  beta      <- parms[["beta"]]
  # g         <- parms[["g"]]
  w_fast    <- parms[["w_fast"]]
  # w_slow    <- parms[["w_slow"]]
  sc        <- parms[["sc"]]
  cfr       <- parms[["cfr"]]
  RR        <- parms[["RR"]]
  
  all_parms <- c(
    fixed_parms, 
    beta = beta, 
    # g = g,
    w_fast = w_fast,
    # w_slow = w_slow,
    sc =sc,
    cfr = cfr
  )
  
  Efast <- state[["E"]]*prop_fast
  Eslow <- state[["E"]] - Efast
  state <- c(state, Efast = Efast, Eslow = Eslow)
  s_order  <- c("B", "M", "S", "Efast", "Eslow", "I", "RH", "RL", "N") 
  state <- state[s_order]
  
  # run model
  times <- 2013:ey
  out <- deSolve::ode(
    y = state, times = times, 
    func = tb_model, 
    parms = all_parms, 
    method = "lsoda",
    rtol = rtol,
    atol = atol
  )
  Efast_out <- out[, "Efast"]
  Eslow_out <- out[, "Eslow"]
  M_out     <- out[, "M"]
  I_out     <- out[, "I"]
  RH_out    <- out[, "RH"]
  RL_out    <- out[, "RL"]
  
  # secondary indicator
  w_fast  <- all_parms[["w_fast"]]
  w_slow  <- all_parms[["w_slow"]]
  rhoH    <- all_parms[["rhoH"]]
  rhoL    <- all_parms[["rhoL"]]
  f       <- all_parms[["f"]]
 
  e_inc_num    <- w_fast*Efast_out + w_slow*Eslow_out + rhoH*RH_out + rhoL*RL_out
  e_recur_num  <- rhoH*RH_out + rhoL*RL_out
  r            <- e_recur_num / pmax(e_inc_num, 1e-12)
  e_inc_rr_num <- e_inc_num * ((1 - f) * pnew * ((1 - r) + r * RR) + f * pret)
  e_mort_num   <- cfr * I_out
  # e_mort_num   <- 0.04 * e_inc_num
  
  list(
    year = out[, "time"],
    e_inc_num = e_inc_num,
    e_mort_num = e_mort_num,
    e_inc_rr_num = e_inc_rr_num,
    raw = out
  )
}


# Intervention ------------------------------------------------------------

tb_intv_model <- function(
    t, y, parms, 
    TPT = FALSE, 
    TPT_Scale = NULL,
    TPT_Effect = 0.33, 
    RECH = FALSE, 
    RECH_Scale = NULL, 
    RECH_Effect = 0.50,
    RECL = FALSE, 
    RECL_Scale = NULL, 
    RECL_Effect = 0.50
  ){
  # Compartment
  B       <- y["B"]
  M       <- y["M"]
  S       <- y["S"]
  Efast   <- y["Efast"]
  Efast_T <- y["Efast_T"]
  Eslow   <- y["Eslow"]
  Eslow_T <- y["Eslow_T"]
  I       <- y["I"]
  RH      <- y["RH"]
  RH_T    <- y["RH_T"]
  RL      <- y["RL"]
  RL_T    <- y["RL_T"]
  N       <- y["N"]
  
  # Parameter
  ## Demographic
  b_B        <- parms["b_B"]
  b_M        <- parms["b_M"]
  
  ## TB Epidemiology
  beta       <- parms["beta"]
  g          <- parms["g"]
  w_fast     <- parms["w_fast"]
  w_slow     <- parms["w_slow"]
  sc         <- parms["sc"]
  cfr        <- parms["cfr"]
  
  f          <- parms["f"]
  eta        <- parms["eta"]
  
  ## PostTB
  rhoH       <- parms["rhoH"]
  rhoL       <- parms["rhoL"]
  toRL       <- parms["toRL"]
  MR         <- parms["MR"]
  
  w_fast_t <- TPT_Effect * w_fast / (1 - w_fast + TPT_Effect * w_fast)
  w_slow_t <- TPT_Effect * w_slow / (1 - w_slow + TPT_Effect * w_slow)
  rhoH_t <- RECH_Effect * rhoH / (1 - rhoH + RECH_Effect * rhoH)
  rhoL_t <- RECL_Effect * rhoL / (1 - rhoL + RECL_Effect * rhoL)
  
  if(TPT){TPT_Scale <- TPT_Scale}else{TPT_Scale <- 0}
  
  if(RECH){RECH_Scale <- RECH_Scale}else{RECH_Scale <- 0}
  
  if(RECL){RECL_Scale <- RECL_Scale}else{RECL_Scale <- 0}
  
  # ODE Formula
  dB <- b_B * (B) * (1-B)
  dM <- b_M * (M) * (1-M)
  
  dS       <- B*N + beta*S*I/N - M*S
  dEfast   <- (1-g)*beta*S*I/N - w_fast*Efast - M*Efast - TPT_Scale*Efast
  dEfast_T <- TPT_Scale*Efast - w_fast_t*Efast_T
  dEslow   <-   (g)*beta*S*I/N - w_slow*Eslow - M*Eslow - TPT_Scale*Eslow
  dEslow_T <- TPT_Scale*Eslow - w_slow_t*Eslow_T
  dI       <- w_fast*Efast + w_slow*Eslow + 
              (w_fast_t*Efast_T + w_slow_t*Eslow_T)  + 
              (rhoH*RH + rhoL*RL) + 
              (rhoH_t*RH_T + rhoL_t*RL_T) - 
              eta*I - sc*I - (M+cfr)*I
  dRH      <- eta*I + sc*I - rhoH*RH - toRL*RH - (M+MR)*RH - RECH_Scale*RH
  dRH_T    <- RECH_Scale*RH - rhoH_t*RH_T - toRL*RH_T - (M+MR)*RH_T
  dRL      <- toRL*RH - rhoL*RL - (M+MR)*RL - RECL_Scale*RL
  dRL_T    <- RECL_Scale*RL + toRL*RH_T - rhoL_t*RL_T - (M+MR)*RL_T
  
  dN <- dS + dEfast + dEslow +dI + dRH + dRL + dRH_T + dRL_T
  
  list(c(dB, dM, dS, dEfast, dEfast_T, dEslow, dEslow_T, dI, dRH, dRH_T, dRL, dRL_T, dN))
}

sim_intv <- function(fixed_parms, parms, state, eys, eye, pnew, pret,
                     TPT = FALSE, TPT_Effect = 0.33, TPT_Scale = 0,    # 1 means OR=1
                     RECH = FALSE, RECH_Effect = 0.50, RECH_Scale = 0, # 1 means OR=1
                     RECL = FALSE, RECL_Effect = 0.50, RECL_Scale = 0, # 1 means OR=1
                     rtol = 1e-6, atol = 1e-6) {
  # Parameter Initial
  prop_fast <- parms[["prop_fast"]]
  beta      <- parms[["beta"]]
  # g         <- parms[["g"]]
  w_fast    <- parms[["w_fast"]]
  # w_slow    <- parms[["w_slow"]]
  sc        <- parms[["sc"]]
  cfr       <- parms[["cfr"]]
  RR        <- parms[["RR"]]
  
  all_parms <- c(
    fixed_parms, 
    beta = beta, 
    # g = g,
    w_fast = w_fast,
    # w_slow = w_slow,
    sc =sc,
    cfr = cfr
  )
  
  name_state <- names(state)
  
  if(! "Efast_T" %in% name_state){state[["Efast_T"]] <- 0}
  if(! "Eslow_T" %in% name_state){state[["Eslow_T"]] <- 0}
  if(! "RH_T" %in% name_state){state[["RH_T"]] <- 0}
  if(! "RL_T" %in% name_state){state[["RL_T"]] <- 0}
  
  s_order  <- c("B", "M", "S", "Efast", "Efast_T", "Eslow", "Eslow_T", "I", "RH", "RH_T", "RL", "RL_T", "N") 
  state <- state[s_order]
  
  # run model
  times <- eys:eye
  out <- deSolve::ode(
    y = state, times = times, 
    func = tb_intv_model, 
    parms = all_parms, 
    method = "lsoda",
    rtol = rtol,
    atol = atol,
    TPT = TPT, 
    TPT_Scale = TPT_Scale,
    TPT_Effect = TPT_Effect, 
    RECH = RECH, 
    RECH_Scale = RECH_Scale,
    RECH_Effect = RECH_Effect,
    RECL = RECL, 
    RECL_Scale = RECL_Scale,
    RECL_Effect = RECL_Effect
  )
  
  Efast_out   <- out[, "Efast"]
  Eslow_out   <- out[, "Eslow"]
  Efast_T_out <- out[, "Efast_T"]
  Eslow_T_out <- out[, "Eslow_T"]
  M_out       <- out[, "M"]
  I_out       <- out[, "I"]
  RH_out      <- out[, "RH"]
  RH_T_out    <- out[, "RH_T"]
  RL_out      <- out[, "RL"]
  RL_T_out    <- out[, "RL_T"]
  
  # secondary indicator
  
  w_fast    <- all_parms[["w_fast"]]
  w_fast_t  <- TPT_Effect * w_fast / (1 - w_fast + TPT_Effect * w_fast)
  w_slow    <- all_parms[["w_slow"]]
  w_slow_t  <- TPT_Effect * w_slow / (1 - w_slow + TPT_Effect * w_slow)
  rhoH      <- all_parms[["rhoH"]]
  rhoH_t    <- RECH_Effect * rhoH / (1 - rhoH + RECH_Effect * rhoH)
  rhoL      <- all_parms[["rhoL"]]
  rhoL_t    <- RECL_Effect * rhoL / (1 - rhoL + RECL_Effect * rhoL)
  f         <- all_parms[["f"]]
  
  alpha     <- 0.09/100 # https://impaact4tb.org/wp-content/uploads/2020/12/2T-Ghana-Resistance-and-TPT_JH4.pdf#1#1
 
  e_inc_norr   <- (w_fast*Efast_out + w_slow*Eslow_out) + 
                  (w_fast_t*Efast_T_out*(1-alpha) + w_slow_t*Eslow_T_out*(1-alpha)) + 
                  (rhoH*RH_out +  rhoL*RL_out) + 
                  (rhoH_t*RH_T_out + rhoL_t*RL_T_out)
  e_inc_withrr <- w_fast_t*Efast_T_out*alpha + w_slow_t*Eslow_T_out*alpha
  e_inc_num    <- (w_fast*Efast_out + w_slow*Eslow_out) + 
                  (w_fast_t*Efast_T_out + w_slow_t*Eslow_T_out) + 
                  (rhoH*RH_out +  rhoL*RL_out) + 
                  (rhoH_t*RH_T_out + rhoL_t*RL_T_out)
  e_recur_num  <- (rhoH*RH_out +  rhoL*RL_out) + 
                  (rhoH_t*RH_T_out + rhoL_t*RL_T_out)
  r            <- e_recur_num / pmax(e_inc_num, 1e-12)
  e_inc_rr_num <- e_inc_norr * ((1 - f) * pnew * ((1 - r) + r * RR) + f * pret) + e_inc_withrr
  e_mort_num   <- cfr*I_out
  
  list(
    year = out[, "time"],
    e_inc_num = e_inc_num,
    e_mort_num = e_mort_num,
    e_inc_rr_num = e_inc_rr_num,
    e_inc_withrr = e_inc_withrr,
    raw = out
  )
}


# change data -------------------------------------------------------------

integerize_sim_raw_safe <- function(simlist, nonnegative = TRUE) {
  raw_df <- as.data.frame(simlist$raw)
  raw_df$e_inc_num <- simlist$e_inc_num
  raw_df$e_mort_num <- simlist$e_mort_num
  raw_df$e_inc_rr_num <- simlist$e_inc_rr_num
  
  count_cols <- intersect(
    c("S", "Efast", "Eslow", "I", "RH", "RL", "N", "e_inc_num", "e_mort_num", "e_inc_rr_num", "e_inc_withrr"),
    names(raw_df)
  )
  
  raw_df[count_cols] <- lapply(raw_df[count_cols], function(x) {
    if (nonnegative) {
      x <- pmax(x, 0)
    }
    as.integer(round(x))
  })
  
  return(raw_df)
}

# Objective function ------------------------------------------------------

objective_weighted <- function(parms, fixed_parms, state, ey, pnew, pret,
                               target_cache, penalty = 1e12) {
  bad <- rep(penalty, target_cache$n)
  
  sim_out <- tryCatch(
    suppressWarnings(sim_tb_model(parms = parms, 
                                  fixed_parms = fixed_parms, 
                                  state = state, 
                                  ey = ey, 
                                  pnew = pnew, 
                                  pret = pret)),
    error = function(e) print("Error While Simulation")
  )
  
  if (is.null(sim_out)) return(bad)
  if (any(!is.finite(sim_out$e_inc_num)) ||
      any(!is.finite(sim_out$e_mort_num)) ||
      any(!is.finite(sim_out$e_inc_rr_num))
      ) return(bad)
  
  loglik_inc <- dnorm(
    x = target_cache$inc$val,
    mean = safe_log(sim_out$e_inc_num[target_cache$inc$idx]),
    sd = target_cache$inc$se,
    log = TRUE
  )
  
  loglik_mort <- dnorm(
    x = target_cache$mort$val,
    mean = safe_log(sim_out$e_mort_num[target_cache$mort$idx]),
    sd = target_cache$mort$se,
    log = TRUE
  )
  
  loglik_rr_inc <- dnorm(
    x = target_cache$rr$val,
    mean = safe_log(sim_out$e_inc_rr_num[target_cache$rr$idx]),
    sd = target_cache$rr$se,
    log = TRUE
  )
  
  res <- -2 * c(loglik_inc, loglik_mort, loglik_rr_inc)
  if (any(!is.finite(res))) return(bad)
  res
}

# Simulation for Intervention Sequence ------------------------------------

sim_intervention_sequence <- function(
    fixed_parms,
    parms,
    init_state,
    base_year,
    intervention_schedule,
    pnew,
    pret,
    rtol = 1e-6,
    atol = 1e-6,
    end_year = NULL
) {
  
  # 确保干预时间表按开始年份排序
  schedule <- intervention_schedule[order(intervention_schedule$start_year), ]
  
  if(!is.null(end_year) & max(schedule$end_year) != end_year){
    if(max(schedule$end_year)<end_year){
      maxn <- nrow(schedule) 
      schedule_add <- schedule[maxn,]
      schedule_add$start_year <- as.numeric(schedule_add$end_year)
      schedule_add$end_year <- as.numeric(end_year)
      schedule_add$TPT <- FALSE
      schedule_add$TPT_Effect <- max(schedule$TPT_Effect)
      schedule_add$TPT_Scale <- 0
      schedule_add$RECH <- FALSE
      schedule_add$RECH_Effect <- max(schedule$RECH_Effect)
      schedule_add$RECH_Scale <- 0
      schedule_add$RECL <- FALSE
      schedule_add$RECL_Effect <- max(schedule$RECL_Effect)
      schedule_add$RECL_Scale <- 0
      schedule <- bind_rows(schedule, schedule_add)
    }
    
    if(max(schedule$end_year)>end_year){
      stop("Intervention year must be less than or equal to end year")
    }
  }
  
  # 确定需要模拟的总年份范围 从 base_year 到最后一个 end_year
  total_years <- base_year:max(schedule$end_year)
  
  pnew_vec <- pnew
  pret_vec <- pret
  
  names(pnew_vec) <- total_years
  names(pret_vec) <- total_years
  
  # 结果存储列表
  results_list <- list()
  current_state <- init_state
  current_year <- base_year
  
  # 逐段模拟
  for (i in 1:nrow(schedule)) {
    seg_start <- schedule$start_year[i]
    seg_end   <- schedule$end_year[i]
    # 开始模拟
    seg_res <- sim_intv(
      fixed_parms = fixed_parms,
      parms = parms,
      state = current_state,
      eys = seg_start,
      eye = seg_end,
      pnew = pnew_vec[as.character(seg_start:seg_end)],
      pret = pret_vec[as.character(seg_start:seg_end)],
      TPT = schedule$TPT[i],
      TPT_Effect = schedule$TPT_Effect[i],
      TPT_Scale = schedule$TPT_Scale[i],
      RECH = schedule$RECH[i],
      RECH_Effect = schedule$RECH_Effect[i],
      RECH_Scale = schedule$RECH_Scale[i],
      RECL = schedule$RECL[i],
      RECL_Effect = schedule$RECL_Effect[i],
      RECL_Scale = schedule$RECL_Scale[i],
      rtol = rtol, atol = atol
    )
    # 提取模拟结果
    seg <- cbind(
      seg_res$raw, 
      data.frame(
        year = seg_res$year, 
        e_inc_num = seg_res$e_inc_num,
        e_mort_num = seg_res$e_mort_num,
        e_inc_rr_num = seg_res$e_inc_rr_num,
        e_inc_withrr = seg_res$e_inc_withrr
      )
    )
    seg <- seg[-1, ]
    results_list <- c(results_list, list(seg))
    # 更新状态为段结束时的状态
    current_state <- seg_res$raw[nrow(seg_res$raw), -1]
  }
  # 合并所有结果
  final_df <- do.call(rbind, results_list)
  return(final_df)
}












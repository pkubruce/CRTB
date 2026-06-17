rm(list = ls())
# load default data
load("testresult/def.RData")

# load tools --------------------------------------------------------------

source("R/00_tb_data.R")

source("R/00_utils.R")

source("R/00_tb_model.R")

source("R/00_tb_plot.R")

library(FME)

library(coda)

library(tidyverse)

library(shiny)

library(rhandsontable)

cali_year <- 2023

initial_schedule <- data.frame(
  start_year  = c(2026,  2027,  2029,  2030,  2032,  2033),
  end_year    = c(2027,  2029,  2030,  2032,  2033,  2035),
  TPT         = c(TRUE,  FALSE, TRUE,  FALSE, TRUE,  FALSE),
  TPT_Effect  = c(0.33,  0.33,  0.33,  0.33,  0.33,  0.33),
  TPT_Scale   = c(0.80,  0.00,  0.80,  0.00,  0.80,  0.00),
  RECH        = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
  RECH_Effect = c(0.50,  0.50,  0.50,  0.50,  0.50,  0.50),
  RECH_Scale  = c(0.80,  0.80,  0.80,  0.80,  0.80,  0.80),
  RECL        = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
  RECL_Effect = c(0.50,  0.50,  0.50,  0.50,  0.50,  0.50),
  RECL_Scale  = c(0.80,  0.80,  0.80,  0.80,  0.80,  0.80)
)


# ui ----------------------------------------------------------------------

ui <- fluidPage(
  # Use bsCollapse for collapsible panels
  shinyBS::bsCollapse(
    id = "collapse_panels",
    multiple = TRUE,
    
    # ========== 1. Baseline Simulation Panel (existing) ==========
    shinyBS::bsCollapsePanel(
      title = "Baseline Simulation",
      value = "baseline",
      style = "info",
      wellPanel(
        fluidRow(
          column(4, fileInput("key_target_file", "Targets (.xlsx)", accept = ".xlsx")),
          column(4, fileInput("target_file", "Demo/RRTB (.xlsx)", accept = ".xlsx")),
          column(4, selectInput("pred_year", "End year", choices = seq(2025, 2050, 1), selected = 2035))
        ),
        fluidRow(
          column(2, selectInput("calibration_TB", "Calibrate TB", choices = c("NO","YES"), selected = "NO")),
          column(2, selectInput("Update_Demo", "Update Demo", choices = c("NO","YES"), selected = "NO")),
          column(2, selectInput("Update_RRTB", "Update RRTB", choices = c("NO","YES"), selected = "NO")),
          column(2, numericInput("N0", "Population 2013", value = 1367260000)),
          column(2, numericInput("prev", "Prevalence (/100k)", value = 94)),
          column(2, numericInput("noti_2013", "Notified 2013", value = 847176))
        ),
        fluidRow(
          column(2, numericInput("recur_2013", "Recurrent 2013", value = 64423)),
          column(2, numericInput("cdr", "CDR (%)", value = 88)),
          column(2, numericInput("prop_latent", "LTBI (%)", value = 18.08)),
          column(2, numericInput("rhoH", "Recurrence rate (high risk)", value = 1.87/100)),
          column(2, numericInput("rhoL", "Recurrence rate (low risk)", value = 0.47/100)),
          column(1, numericInput("succ", "Tx success (%)", value = 95.00)),
          column(1, numericInput("fail", "Tx failure (%)", value = 0.50))
        )
      )
    ),
    
    # ========== 2. Intervention Simulation Panel (new) ==========
    shinyBS::bsCollapsePanel(
      title = "Intervention Simulation",
      value = "intervention",
      style = "info",
      wellPanel(
        h5("Default Intervention Plan (Reference)"),
        tableOutput("default_schedule_table"),
        hr(),
        p("You can edit the intervention schedule in the 'Intervention' tab above."),
        helpText("Click 'Run Intervention' after editing the schedule.")
      )
    )
  ),
  
  # Run button (use actionButton instead of submitButton to avoid page refresh)
  fluidRow(
    column(12, h5("Please run simulation before intervention", align = "center")),
    column(6, actionButton("run_simulation", "Run Simulation", icon = icon("play"), class = "btn-success")),
    column(6, actionButton("run_intervention", "Run Intervention", icon = icon("chart-line"), class = "btn-primary")),
    align = "center"
  ),
  
  # Output tabs
  fluidRow(
    column(
      12,
      tabsetPanel(
        tabPanel("Posterior of parameters", plotOutput("plot_param", height = "800px")),
        tabPanel("Simulation", plotOutput("plot_predict", height = "800px")),
        tabPanel("Intervention Plan",
                 h4("Intervention Schedule (Editable)"),
                 rHandsontableOutput("intervention_hot"),
                 br(),
                 fluidRow(
                   column(6, actionButton("add_schedule_row", "Add Row", icon = icon("plus"), class = "btn-primary"))
                 ),
                 br(),
                 p("Double-click a cell to edit. Use the buttons to add/delete rows. The schedule will be used when you click 'Run Intervention'.")
        ),
        tabPanel("Intervention Effect", plotOutput("plot_intervention", height = "800px"))
    )
    )
  )
)


server <- function(input, output, session) {
  
  # ---------- 存储结果的反应式变量 ----------
  rv <- reactiveValues(
    baseline = NULL,    # 存储基线模拟的结果列表 (mcmcres, s_df, target_df 等)
    intervention = NULL # 存储干预模拟的结果 (s_df_intv)
  )
  
  # ---------- 辅助函数 ----------
  build_baseline_inputs <- function() {
    # Update Demographic Data
    if(input$Update_Demo == "YES"){
      df_birth <- openxlsx::read.xlsx(xlsxFile = input$target_file$datapath, sheet = "Natural_Birth")
      df_death <- openxlsx::read.xlsx(xlsxFile = input$target_file$datapath, sheet = "Natural_Death")
      list_birth <- update_demo(df_birth)
      list_death <- update_demo(df_death)
      demo_list <- c(B = list_birth$val0, b_B = list_birth$beta, 
                     M = list_death$val0, b_M = list_death$beta)
    } else {
      demo_list <- defdata[c("B", "b_B", "M", "b_M")]
    }
    
    # Update RRTB related Data
    if(input$Update_RRTB == "YES"){
      df_pct_new <- openxlsx::read.xlsx(xlsxFile = input$target_file$datapath, sheet = "PCT_NEW")
      df_pct_ret <- openxlsx::read.xlsx(xlsxFile = input$target_file$datapath, sheet = "PCT_RET")
      fit_new <- update_rrtb(df_pct_new)
      fit_ret <- update_rrtb(df_pct_ret)
      newyear <- 2013:as.integer(input$pred_year)
      pred_new <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = newyear))))
      pred_ret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = newyear))))
      cali_new <- pred_new[1:length(c(2013:cali_year))]
      cali_ret <- pred_ret[1:length(c(2013:cali_year))]
    } else {
      newyear <- 2013:as.integer(input$pred_year)
      pred_new <- logit2p(as.numeric(predict(fit_new, newdata = data.frame(year = newyear))))
      pred_ret <- logit2p(as.numeric(predict(fit_ret, newdata = data.frame(year = newyear))))
      cali_new <- pred_new[1:length(c(2013:cali_year))]
      cali_ret <- pred_ret[1:length(c(2013:cali_year))]
    }
    
    # Fixed parameters and initial state
    rhoH     <- input$rhoH
    rhoL     <- input$rhoL
    
    B        <- demo_list[["B"]]
    b_B      <- demo_list[["b_B"]]
    M        <- demo_list[["M"]]
    b_M      <- demo_list[["b_M"]]
    MR       <- (2.91*M - M)
    
    N       <- input$N0
    E       <- N * (input$prop_latent / 100)
    I       <- N * (input$prev/1e5)
    RH      <- (input$recur_2013/(input$cdr/100))/rhoH
    RL      <- 17174029/1426106090*input$N0 - RH
    S       <- N - E - I - RH - RL
    
    eta      <- input$noti_2013*(input$succ/100)/I
    f        <- input$fail/100
    
    fixed_parms <- c(
      # Demographic Data
      b_B = b_B,
      b_M = b_M,
      
      # TB Epidemiology
      # beta = beta,
      # sl = sl,
      g = g,
      # w_fast = w_fast,
      w_slow = w_slow,
      # sc = sc,
      # cfr = cfr,
      f = f,
      eta = eta,
      rhoH = rhoH,
      rhoL = rhoL,
      toRL = toRL,
      MR = MR
    )
    
    state <- c(
      B = B,
      M = M,
      S = S,
      E = E,
      I = I,
      RH = RH,
      RL = RL,
      N = N
    )
    
    list(
      fixed_parms = fixed_parms,
      state = state,
      pred_new = pred_new,
      pred_ret = pred_ret,
      cali_new = cali_new,
      cali_ret = cali_ret,
      newyear = newyear
    )
  }
  
  schedule_data <- reactiveVal(initial_schedule)

  # ---------- 基线模拟 ----------
  observeEvent(input$run_simulation, {
    showNotification("Running baseline simulation...", type = "message", duration = 3)
    
    inputs <- build_baseline_inputs()
    fixed_parms <- inputs$fixed_parms
    state <- inputs$state
    pred_new <- inputs$pred_new
    pred_ret <- inputs$pred_ret
    cali_new <- inputs$cali_new
    cali_ret <- inputs$cali_ret
    
    # Calibration if needed
    if (input$calibration_TB == "YES") {
      # 读取目标数据
      df_inc    <- openxlsx::read.xlsx(xlsxFile = input$key_target_file$datapath, sheet = "TotalTB_Incident")
      df_mort   <- openxlsx::read.xlsx(xlsxFile = input$key_target_file$datapath, sheet = "TotalTB_Deaths")
      df_inc_rr <- openxlsx::read.xlsx(xlsxFile = input$key_target_file$datapath, sheet = "RRTB_Incident")
      target_df <- make_target_df(list(df_inc, df_mort, df_inc_rr), name = c("Incident TB", "Deaths from TB", "Incident RRTB"))
      target_cache <- prepare_target_cache(target_df, start_year = 2013, ey = cali_year)
      
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
      
      fit <- modFit(
        f = objective_weighted, p = par0, lower = lower, upper = upper, method = "L-BFGS-B",
        fixed_parms = fixed_parms, ey = cali_year, state = state,
        pnew = cali_new, pret = cali_ret, target_cache = target_cache
      )
      par_best_fit <- coef(fit)
      var0 <- summary(fit)$modVariance
      
      set.seed(1234)
      niter <- 4e3
      mcmcres <- modMCMC(
        f = objective_weighted, p = par_best_fit, lower = lower, upper = upper,
        niter = niter, var0 = var0, burninlength = 1e3, updatecov = 200,
        verbose = 500, ntrydr = 2,
        fixed_parms = fixed_parms, ey = cali_year, state = state,
        pnew = cali_new, pret = cali_ret, target_cache = target_cache
      )
    } else {
      load("testresult/mcmc.RData")  # 应确保 mcmcres 对象存在
    }
    
    # 从 MCMC 结果中抽样进行预测
    set.seed(1234)
    nsamp <- 100
    idx <- sample(1:nrow(mcmcres$pars), nsamp)
    slist <- list()
    for (i in 1:nsamp) {
      s <- sim_tb_model(
        fixed_parms = fixed_parms,
        parms = mcmcres$pars[idx[i], ],
        ey = as.integer(input$pred_year),
        state = state,
        pnew = pred_new,
        pret = pred_ret
      )
      s <- integerize_sim_raw_safe(s)
      s$iter <- i
      slist[[i]] <- s
    }
    s_df <- do.call(rbind, slist)
    
    # 存储基线结果
    rv$baseline <- list(
      mcmcres = mcmcres,
      s_df = s_df,
      target_df = target_df,
      fixed_parms = fixed_parms,
      state = state,
      pred_new = pred_new,
      pred_ret = pred_ret,
      idx = idx,
      ey = input$pred_year
    )
    
    showNotification("Baseline simulation completed.", type = "message", duration = 3)
  })
  
  # ---------- 干预模拟 ----------
  observeEvent(input$run_intervention, {
    req(rv$baseline)  # 必须已有基线结果
    
    schedule <- schedule_data()
    req(schedule, nrow(schedule) > 0)
    
    showNotification("Running intervention simulation...", type = "message", duration = 3)
    
    fixed_parms <- rv$baseline$fixed_parms
    base_state <- rv$baseline$state
    pred_new <- rv$baseline$pred_new
    pred_ret <- rv$baseline$pred_ret
    mcmcres <- rv$baseline$mcmcres
    idx <- rv$baseline$idx
    s_df <- rv$baseline$s_df
    ey <- rv$baseline$ey
    
    set.seed(1234)
    baseey <- 2026
    nsamp <- 100
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
        pnew = pred_new[baseey:ey-2013+1],
        pret = pred_ret[baseey:ey-2013+1],
        end_year = ey
      )
      
      res <- as.data.frame(rbind_fill(si, res))
      
      res$iter <- i
      
      reslist[[i]] <- res
    }
    
    intv_df <- do.call(rbind, reslist)
    
    rv$intervention <- list(
      baseline_df = s_df,
      intv_df = intv_df
    )
    showNotification("Intervention simulation completed.", type = "message", duration = 3)
    })
  
  # ---------- 绘图输出 ----------
  output$plot_param <- renderPlot({
    req(rv$baseline$mcmcres)
    mcmcres <- rv$baseline$mcmcres
    par(mfrow = c(2,3))
    hist(mcmcres$pars[, "prop_fast"], main = "prop_fast", xlab = "prop_fast", breaks = 30)
    hist(mcmcres$pars[, "beta"], main = "beta", xlab = "beta", breaks = 30)
    hist(mcmcres$pars[, "w_fast"], main = "w_fast", xlab = "w_fast", breaks = 30)
    hist(mcmcres$pars[, "sc"], main = "sc", xlab = "sc", breaks = 30)
    hist(mcmcres$pars[, "cfr"], main = "cfr", xlab = "cfr", breaks = 30)
    hist(mcmcres$pars[, "RR"], main = "RR", xlab = "RR", breaks = 30)
    par(mfrow = c(1,1))
  }, res = 144)
  
  output$plot_predict <- renderPlot({
    req(rv$baseline$s_df, rv$baseline$target_df)
    
    req(rv$baseline$s_df, rv$baseline$target_df)
    
    s_df <- rv$baseline$s_df
    target_df <- rv$baseline$target_df
    ey <- as.integer(input$pred_year)
    
    par(mfrow = c(2,3))
    
    plot_number(
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_inc_num",
      target_name = "Incident TB",
      main_name = "Total Incident TB"
    )
    
    plot_number(
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_mort_num",
      target_name = "Deaths from TB",
      main_name = "Total TB Deaths"
    )
    
    plot_number(
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_inc_rr_num",
      target_name = "Incident RRTB",
      main_name = "Total Incident RRTB"
    )
    
    plot_rate(
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_inc_num",
      target_name = "Incident TB",
      main_name = "Total Incident TB"
    )
    
    plot_rate( 
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_mort_num",
      target_name = "Deaths from TB",
      main_name = "Total TB Deaths")
    
    plot_rate(
      s_df = s_df, obs = obs, ey = ey,
      ind_name = "e_inc_rr_num",
      target_name = "Incident RRTB",
      main_name = "Total Incident RRTB"
    )
  }, res = 144)
  
  output$default_schedule_table <- renderTable({
    initial_schedule
  }, striped = TRUE, bordered = TRUE, digits = 2)
  
  output$intervention_hot <- renderRHandsontable({
    df <- schedule_data()
    # 确保 TPT 和 REC 为逻辑值
    df$TPT <- as.logical(df$TPT)
    df$RECH <- as.logical(df$RECH)
    df$RECL <- as.logical(df$RECL)
    rhandsontable(df, stretchH = "all", rowHeaders = TRUE) %>%
      hot_col("start_year", type = "numeric", format = "0") %>%
      hot_col("end_year",   type = "numeric", format = "0") %>%
      hot_col("TPT",        type = "checkbox", checkedTemplate = TRUE, uncheckedTemplate = FALSE) %>%
      hot_col("TPT_Effect", type = "numeric", format = "0.00") %>%
      hot_col("TPT_Scale",  type = "numeric", format = "0.00") %>%
      hot_col("RECH",        type = "checkbox", checkedTemplate = TRUE, uncheckedTemplate = FALSE) %>%
      hot_col("RECH_Effect", type = "numeric", format = "0.00") %>%
      hot_col("RECH_Scale",  type = "numeric", format = "0.00") %>% 
      hot_col("RECL",        type = "checkbox", checkedTemplate = TRUE, uncheckedTemplate = FALSE) %>%
      hot_col("RECL_Effect", type = "numeric", format = "0.00") %>%
      hot_col("RECL_Scale",  type = "numeric", format = "0.00")
  })
  
  observeEvent(input$intervention_hot, {
    if (!is.null(input$intervention_hot)) {
      new_df <- hot_to_r(input$intervention_hot)
      # 确保逻辑列正确
      new_df$TPT <- as.logical(new_df$TPT)
      new_df$RECH <- as.logical(new_df$RECH)
      new_df$RECL <- as.logical(new_df$RECL)
      schedule_data(new_df)
    }
  })
  
  observeEvent(input$add_schedule_row, {
    current <- schedule_data()
    new_row <- data.frame(
      start_year  = max(current$end_year) + 1,
      end_year    = max(current$end_year) + 2,
      TPT         = FALSE,          # 逻辑值
      TPT_Effect  = 0.33,
      TPT_Scale   = 0.80,
      RECH        = FALSE,          # 逻辑值
      RECH_Effect = 0.50,
      RECH_Scale  = 0.00,
      RECL        = FALSE,          # 逻辑值
      RECL_Effect = 0.50,
      RECL_Scale  = 0.00
    )
    schedule_data(rbind(current, new_row))
  })
  
  output$plot_intervention <- renderPlot({
    req(rv$intervention)
    
    baseline_df <- rv$intervention$baseline_df
    intv_df <- rv$intervention$intv_df
    
    # to calculate the rate
    s_df_agg <- as.data.table(baseline_df)[,.(  ## Baseline 
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
    
    intv_df <- as.data.table(intv_df)
    
    intv_df_agg <- intv_df[,.( ## Intervention 
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

    par(
      mfrow = c(2,3),
      mar = c(4.5, 5.3, 2.8, 0.8),   # The second number control left margin
      mgp = c(2.2, 0.6, 0),          # The second number control left line
      tcl = -0.25,
      cex.axis = 0.85,
      cex.lab = 0.9,
      cex.main = 0.95
    )
    
    ind_num  <- c("e_inc_num", "e_mort_num", "e_inc_rr_num")
    main_num <- c("TB Incident Cases", "Deaths from TB", "RRTB Incident Cases") 
    
    ind_rate  <- c("e_inc_100k", "e_mort_100k", "e_inc_rr_100k")
    main_rate <- c("TB Incidence", "TB Mortality", "RRTB Incidence") 
    
    ey <- max(intv_df_agg$time)
    
    for(i in 1:3){
      # Number
      plot_comp(
        ind = ind_num[i], sy = 2026, ey = ey,
        ylab = "Number of cases", main_name = main_num[i],
        basedf = s_df_agg[time>=2026], 
        intvdf = intv_df_agg[time>=2026]
      )
    }
    
    for(i in 1:3){
      # Rate
      plot_comp(
        ind = ind_rate[i], sy = 2026, ey = ey,
        ylab = "Rate (per 100,000) ", main_name = main_rate[i],
        basedf = s_df_agg[time>=2026], 
        intvdf = intv_df_agg[time>=2026]
      )
    }
    
    par(mfrow = c(1,1))
  }, res = 144)
  
}

# Run the application 
shinyApp(ui = ui, server = server)

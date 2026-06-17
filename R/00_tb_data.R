rm(list = ls())

# load tools --------------------------------------------------------------

source("R/00_utils.R")

source("R/00_tb_model.R")

library(FME)

library(data.table)

library(bayesplot)

library(coda)

# load data ---------------------------------------------------------------
## Demographic Data
df_B      <- openxlsx::read.xlsx(xlsxFile = "testdata/demo_drtb_target_file.xlsx", sheet = "Natural_Birth")
df_M      <- openxlsx::read.xlsx(xlsxFile = "testdata/demo_drtb_target_file.xlsx", sheet = "Natural_Death")
## TB Data
df_inc    <- openxlsx::read.xlsx(xlsxFile = "testdata/TB_indicate_target_file.xlsx", sheet = "TotalTB_Incident")
df_mort   <- openxlsx::read.xlsx(xlsxFile = "testdata/TB_indicate_target_file.xlsx", sheet = "TotalTB_Deaths")
df_inc_rr <- openxlsx::read.xlsx(xlsxFile = "testdata/TB_indicate_target_file.xlsx", sheet = "RRTB_Incident")
## RRTB Data
df_pnew   <- openxlsx::read.xlsx(xlsxFile = "testdata/demo_drtb_target_file.xlsx", sheet = "PCT_NEW")
df_pret   <- openxlsx::read.xlsx(xlsxFile = "testdata/demo_drtb_target_file.xlsx", sheet = "PCT_RET")

## notification data
df_noti   <- read.csv("GTB2025/TB_notifications_2025-11-13.csv")
df_burd   <- read.csv("GTB2025/TB_burden_countries_2025-11-13.csv")

# set test data -----------------------------------------------------------
## Demographic Data
resB     <- update_demo(df_B)
B        <- resB$val0
b_B      <- resB$beta

resM     <- update_demo(df_M)
M        <- resM$val0
b_M      <- resM$beta

## RRTB Data

fit_new <- update_rrtb(df_pnew)
fit_ret <- update_rrtb(df_pret)

## Natural Paramter
g        <- 0.90
# sl       <- 0.03 
# w_fast   <- 0.0826
w_slow   <- 0.0006
# sc       <- 0.23
rho      <- 1.87/100 # JAMA Network Open, 2024
toRL     <- 1/2 
MR       <- (2.91*M - M)

## TB Data
N        <- 1367260000
E        <- N*0.1808

prev     <- 94/1e5 # Global Tuberculosis Report 
I        <- N*prev

noti     <- 847176 # WHO Global Tuberculosis Database
succ     <- 0.95   # WHO Global Tuberculosis Database
cdr      <- 0.88   # WHO Global Tuberculosis Database
recurnum <- 64423  # WHO Global Tuberculosis Database
RH       <- (recurnum/cdr)/rho # estimate
RL       <- (17174029/1426106090)*N - RH # Lancet Infect Dis, 2019
S        <- N - E - I - RH - RL

## Treatment 
eta      <- noti*succ/I
f        <- 0.005

## PostTB
rhoH <- rho
rhoL <- 0.47/100

# clean data --------------------------------------------------------------

target_df <- make_target_df(list(df_inc, df_mort, df_inc_rr), name = c("Incident TB", "Deaths from TB", "Incident RRTB"))
target_df <- merge(target_df, subset(df_burd,iso3=="CHN")[,c("year","e_pop_num")], by = "year")

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

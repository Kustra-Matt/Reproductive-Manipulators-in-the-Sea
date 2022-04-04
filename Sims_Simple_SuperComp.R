#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script that was used for numerical simulation across life history traits (egg size and mortality), costs, and feminization rates for 
#"On the spread of microbes that manipulate reproduction in marine invertebrates" 
#Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282. 
# Uses the simple population model. It uses “Larval_data.csv” which contains the data for larval mortality.
# The results from running this script are saved as Data_FC_5.feather” and Data_FC_A3.feather, which are used to make Figure 5 and Figure A3 by “Plotting_Sim_Heatmaps.R”.
# 1. Setting up ------------------------------------------------
# * 1.a Loading up required libraries -------------------------------------
library(tidyverse)
library(doParallel)
library(feather)
# * 1.b Clear working space -------------------------------------
rm(list = ls())
# 2. Functions ------------------------------------------------

# * 2.a Number of zygote function -----------------------------------------
# Males is male adult density (Males per m^2)
# sigma is egg size (cross sectional area of the egg; mm^2)
# v is sperm speed (mm/s)
# Fe is fertilization efficiency
# tb is time for polyspermy block (s)
# sr is sex ratio (proportion of population that is male)
# S is sperm density per unit male density (sperm/uL/individual/m^2)
# E is egg density per unit female density (egg/uL/individual/m^2)
# tau is half life of sperm (s)
# c is cost due to infection.
# nif is density of infected females (females per m^2)
# Fertilization dynamics function
# psi is settlment constant (uL/m^2)
# Fertilization dynamics function
number_zygotes <- function(Males, sigma, v, Fe, tb, sr, S, E, tau, nif, c, psi = 1) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
  S0 <- Males * S
  Females <- (Males / sr) - Males
  # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
  E0 <- E * (Females - nif) + E * (1 - c) * nif
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizaing sperm
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b)))
  return(prop_mono * E0 * psi)
}
# * 2.b Population recursion equations ------------------------------------
# Population dynamic functions for feminization
# Nz is density of settled zygotes
# Nif is density of in infected female
# Nuf is density of uninfected females
# Nim is density of infected males
# Num is density of uninfected males
# K is adult carrying capacity (individuals/m^2)
# d is mortality rate of changee
# Ma is base adult mortality
# ml is larval mortality
# R is feminzation rate
# Feminization models
# (nuf/Females) + (1-c)*(nuf/Females)
# infected females
Dens_F_i_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, R, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzif <- Nz * (R * (((1 - c) * Nif) / (Nif * (1 - c) + Nuf)))
  return((Nzif * (1 - Ml) + Nif) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# uninfected females
Dens_F_u_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzuf <- Nz * (Nuf / (Nif * (1 - c) + Nuf)) * 0.5
  return((Nzuf * (1 - Ml) + Nuf) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# infected males
Dens_M_i_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, R, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzim <- Nz * ((1 - R) * ((Nif * (1 - c)) / (Nif * (1 - c) + Nuf)))
  return((Nzim * (1 - Ml) + Nim) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# uninfected females
Dens_M_u_t1f <- function(Nz, Nif, Nuf, Nim, Num, K, d, Ma, Ml, c) {
  Nt <- Nif + Nuf + Nim + Num
  Nzmf <- Nz * (Nuf / (Nif * (1 - c) + Nuf)) * 0.5
  return((Nzmf * (1 - Ml) + Num) * (1 - Ma / (1 + exp(-d * (
    Nt + Nz * (1 - Ml) - K / 2
  )))))
}
# * 2.c Vectorized equal function -----------------------------------------
# vectorized equal
# tol is how close doubles can be together
# taken from https://stackoverflow.com/questions/35097815/vectorized-equality-testing
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}
# 3.Run simulations for Figure 5. -----------------------------------------
# simulate function
# initialize things
# register number of cores for parallel computing
registerDoParallel(23)
# setting working directory
setwd("/hb/home/...")
# do a smaller subset of eggs size. these egg sizes were generally infected when infection was possible in a less fine scale analysis
eggs <- seq(450, 800, 50)
eggs <- (pi * (eggs / 1000)^2) / 4
# read in mortality data set
mortality <- read.csv("Data/Larval_data.csv")
# do a smaller subset of mortalities,these larval days were generally infected when infection was possible in a less fine scale analysis
mortality <- mortality[c(3, 7, 11, 15, 20, 32), ]

# make list of feminization values to iterate through
femsr <- seq(0.501, 1, 0.001)
# make list of cost values to iterate through
costs <- seq(0, 0.6, 0.001)

fcsims <-
  foreach(c = costs, .combine = rbind) %:%
  foreach(e = eggs, .combine = rbind) %:%
  foreach(
    mli = 1:length(mortality$Mortality),
    .combine = rbind
  ) %:%
  foreach(r = femsr, .combine = rbind) %dopar% {
    # set flag to false, this lets us know whether population sizes of equilibriated
    flags <- c(FALSE, FALSE, FALSE)
    # keep track of any numerical errors
    er <- "None"
    # intilize starting denisties
    nt0 <- c(.01, .01, .99, .99)
    # Larval day
    mld <- mortality[mli, 1]
    # actual mortality rate
    ml <- mortality[mli, 2]
    # calculate number of eggs based on egg size. Number taken from same overal reproductive value as H. tub.
    numegg <- 1.132059 / e
    # set generation time to time
    time <- 1
    # Keep going if flags have not stabilized
    while (!all(flags)) {
      nz <-
        number_zygotes(
          Males = nt0[2] + nt0[4],
          sigma = e,
          v = 0.125,
          Fe = 0.09444,
          tb = 1,
          sr = nt0[2] + nt0[4] / sum(nt0),
          S = 700,
          E = numegg,
          tau = 5400,
          nif = nt0[1],
          c = c
        )
      nt1 <-
        c(
          Dens_F_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c
          ),
          Dens_M_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c
          ),
          Dens_F_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          ),
          Dens_M_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          )
        )
      # test if there is any change
      flags <- is_equal_tol(nt0, nt1, tol = .Machine$double.eps^0.5)
      if (any(is.na(flags))) {
        er <- "NAs"
        break
      }
      # some simulations result in cycles, this breaks the cycle
      if (time > 100000) {
        er <- "Exceeded Time"
        break
      }
      nt0 <- nt1
      time <- time + 1
    }
    dens <- sum(nt0)
    sr <- (nt0[1] + nt0[3]) / dens
    ir <- (nt0[1] + nt0[2]) / dens
    result <-
      data.frame(
        Egg_Size = e,
        Larval_d = mld,
        Larval_m = ml,
        Egg_number = numegg,
        Feminization = r,
        Time = time,
        Density = dens,
        Infection = ir,
        Sex_ratio = sr,
        Error = er,
        nif0 = nt0[1],
        nim0 = nt0[2],
        nuf0 = nt0[3],
        num0 = nt0[4],
        cost = c
      )
    return(result)
  }
write_feather(fcsims, "Data_FC_5.feather")

# 4.Run simulations for Figure A3. -----------------------------------------
# simulate function
# initialize things
# register number of cores for parallel computing
registerDoParallel(23)
# do a smaller subset of eggs size. these egg sizes were generally infected when infection was possible in a less fine scale analysis
# mortalityFC<-mortality[3:52,]

eggsA3 <- seq(40, 1000, 10)
eggsA3 <- (pi * (eggsA3 / 1000)^2) / 4
# read in mortality data set
mortalityA3 <- read.csv("Data/Larval_data.csv")
# do a subset of mortalities that capture a large range of mortalities.

mortalityA3 <- mortalityA3[3:52, ]

# make list of feminization values to iterate through
femsrA3 <- c(0.65, 0.7, 0.85, 0.95)
# make list of cost values to iterate through
costsA3 <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

fcsimsA3 <-
  foreach(c = costsA3, .combine = rbind) %:%
  foreach(e = eggsA3, .combine = rbind) %:%
  foreach(
    mli = 1:length(mortalityA3$Mortality),
    .combine = rbind
  ) %:%
  foreach(r = femsrA3, .combine = rbind) %dopar% {
    # set flag to false, this lets us know whether population sizes of equilibriated
    flags <- c(FALSE, FALSE, FALSE)
    # keep track of any numerical errors
    er <- "None"
    # intilize starting denisties
    nt0 <- c(.01, .01, .99, .99)
    # Larval day
    mld <- mortalityA3[mli, 1]
    # actual mortality rate
    ml <- mortalityA3[mli, 2]
    # calculate number of eggs based on egg size. Number taken from same overal reproductive value as H. tub.
    numegg <- 1.132059 / e
    # set generation time to time
    time <- 1
    # Keep going if flags have not stabilized
    while (!all(flags)) {
      nz <-
        number_zygotes(
          Males = nt0[2] + nt0[4],
          sigma = e,
          v = 0.125,
          Fe = 0.09444,
          tb = 1,
          sr = nt0[2] + nt0[4] / sum(nt0),
          S = 700,
          E = numegg,
          tau = 5400,
          nif = nt0[1],
          c = c
        )
      nt1 <-
        c(
          Dens_F_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c
          ),
          Dens_M_i_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            R = r,
            c = c
          ),
          Dens_F_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          ),
          Dens_M_u_t1f(
            Nz = nz,
            Nif = nt0[1],
            Nuf = nt0[3],
            Nim = nt0[2],
            Num = nt0[4],
            K = 100,
            d = 0.1,
            Ma = 0.99,
            Ml = ml,
            c = c
          )
        )
      # test if there is any change
      flags <- is_equal_tol(nt0, nt1, tol = .Machine$double.eps^0.5)
      if (any(is.na(flags))) {
        er <- "NAs"
        break
      }
      # some simulations result in cycles, this breaks the cycle
      if (time > 100000) {
        er <- "Exceeded Time"
        break
      }
      nt0 <- nt1
      time <- time + 1
    }
    dens <- sum(nt0)
    sr <- (nt0[1] + nt0[3]) / dens
    ir <- (nt0[1] + nt0[2]) / dens
    result <-
      data.frame(
        Egg_Size = e,
        Larval_d = mld,
        Larval_m = ml,
        Egg_number = numegg,
        Feminization = r,
        Time = time,
        Density = dens,
        Infection = ir,
        Sex_ratio = sr,
        Error = er,
        nif0 = nt0[1],
        nim0 = nt0[2],
        nuf0 = nt0[3],
        num0 = nt0[4],
        cost = c
      )
    return(result)
  }

write_feather(fcsimsA3, "Data_FC_A3.feather")

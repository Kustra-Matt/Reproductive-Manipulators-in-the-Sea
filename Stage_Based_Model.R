#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script file that runs and plots the population model with different size classes for
#"On the spread of microbes that manipulate reproduction in marine invertebrates" 
#Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282. 
# Currently coded to simulate and plot Figure 3, Figure A7, Figure A8, and Figure A9
# but users can run/plot the model under different parameter combinations.
# Creates/uses “Sensitivity_Stage_Based.csv” to make Figure A8.
# 1. Setting up -----------------------------------------------------------

# * 1.a Load up libraries -------------------------------------------------
library(viridis)
library(patchwork)
library(scales)
library(tidyverse)
library(deSolve)
library(svglite)
library(doParallel)

# * 1.b Setting my default plotting theme -----------------------
mytheme <-
  theme_linedraw() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(size = 12, color = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    axis.text = element_text(color = "black")
  ) # makes all the text larger and bold
theme_set(mytheme)


# 2. Functions ------------------------------------------------------------
# function to calculate number of zygotes
number_zygotes <-
  function(f1i, #density of infected females of 20-39mm size class
           f1u,#density of uninfected females of 20-39mm size class
           m1i,#density of infected males of 20-39mm size class
           m1u,#density of uninfected males of 20-39mm size class
           f2i,#density of infected females of 40-59mm size class
           f2u,#density of uninfected females of 40-59mm size class
           m2i,#density of infected males of 40-59mm size class
           m2u,#density of uninfected males of 40-59mm size class
           f3i,#density of infected females of 60+mm size class
           f3u,#density of uninfected females of 60+mm size class
           m3i,#density of infected males of 60+mm size class
           m3u,#density of uninfected males of 60+mm size class
           sigma,#egg size (cross sectional area of the egg; mm^2)
           v,#sperm speed (mm/s)
           Fe,#fertilization efficiency
           tb,#time for polyspermy block (s)
           S,#sperm density per unit male density (sperm/uL/individual/m^2)
           E,# E is egg density per unit female density (egg/uL/individual/m^2)
           tau,#half life of sperm (s)
           c,#cost due to infection.
           ded1, #reduction in fecundity for 20-39mm relative to 60+mm
           ded2,#reduction in fecundity for 40-59mm relative to 60+mm
           psi=1# psi is settlment constant (uL/m^2)
           ) {
    # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
    S0 <- S * ((m3i + m3u) + (m2i + m2u) * (ded2) + (m1i + m1u) * (ded1))
    # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
    E0 <-
      E * ((f3i * (1 - c) + f3u) + (f2i * (1 - c) + f2u) * (ded2) + (f1i * (1 -
        c) + f1u) * (ded1))
    beta0 <- sigma * v
    tau <- tau
    # x=Average number of potential fertilizaing sperm
    x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
    # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
    b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
    prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
      (1 - exp(-1 * b)))
    return(prop_mono * E0*psi)
  }
# birth -> death -> transition
# function to run simulations
pop_next <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Nz <-
      number_zygotes(
        f1i,
        f1u,
        m1i,
        m1i,
        f2i,
        f2u,
        m2i,
        m2u,
        f3i,
        f3u,
        m3i,
        m3u,
        sigma,
        v,
        Fe,
        tb,
        S,
        E,
        tau,
        c,
        ded1,
        ded2
      )
    # total sperm
    st <- (m3i + m3u) + (m2i + m2u) * (ded2) + (m1i + m1u) * (ded1)
    # total Egg
    et <-
      ((f3i * (1 - c) + f3u) + (f2i * (1 - c) + f2u) * (ded2) + (f1i * (1 - c) +
        f1u) * (ded1))
    # infected egg
    ei <- (f3i * (1 - c)) + (f2i * (1 - c) * ded2) + (f1i * (1 - c) * ded1)
    # uninfected egg
    eu <- f3u + f2u * ded2 + f1u * ded1
    # infected sperm
    si <- m3i + m2i * ded2 + m1i * ded1
    # uninfected sperm
    su <- m3u + m2u * ded2 + m1u * ded1
    Nzim <- (Nz * (ei / et) * 0.5) * (1 - Ml) * (1 - Mk)
    Nzum <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -
      Ml)
    Nzif <- (Nz * (ei / et) * 0.5) * (1 - Ml)
    Nzuf <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -
      Ml)
    Nzt <- Nzim + Nzum + Nzif + Nzuf
    # total adult density
    Nta <- f1i + f1u + m1i + m1u + f2i + f2u + m2i + m2u + f3i + f3u + m3i +
      m3u
    # death
    # mortality
    mort <- (1 - Ma / (1 + exp(-d * (Nta + Nzt - K / 2))))
    f1i_plus <- (Nzif + f1i) * mort
    f1u_plus <- (Nzuf + f1u) * mort
    m1i_plus <- (Nzim + m1i) * mort
    m1u_plus <- (Nzum + m1u) * mort
    f2i_plus <- f2i * mort
    f2u_plus <- f2u * mort
    m2i_plus <- m2i * mort
    m2u_plus <- m2u * mort
    f3i_plus <- f3i * mort
    f3u_plus <- f3u * mort
    m3i_plus <- m3i * mort
    m3u_plus <- m3u * mort
    # transition
    f1i_plus1 <- f1i_plus * (1 - tif)
    f1u_plus1 <- f1u_plus * (1 - t)
    m1i_plus1 <- m1i_plus * (1 - t)
    m1u_plus1 <- m1u_plus * (1 - t)
    f2i_plus1 <- f1i_plus * tif + f2i_plus * (1 - tif)
    f2u_plus1 <- f1u_plus * t + f2u_plus * (1 - t)
    m2i_plus1 <- m1i_plus * t + m2i_plus * (1 - t)
    m2u_plus1 <- m1u_plus * t + m2u_plus * (1 - t)
    f3i_plus1 <- f2i_plus * tif + f3i_plus
    f3u_plus1 <- f2u_plus * t + f3u_plus
    m3i_plus1 <- m2i_plus * t + m3i_plus
    m3u_plus1 <- m2u_plus * t + m3u_plus
    list(
      c(
        f1i_plus1,
        f1u_plus1,
        m1i_plus1,
        m1u_plus1,
        f2i_plus1,
        f2u_plus1,
        m2i_plus1,
        m2u_plus1,
        f3i_plus1,
        f3u_plus1,
        m3i_plus1,
        m3u_plus1
      )
    )
  })
}

# vectorized equal
# tol is how close doubles can be together
# taken from https://stackoverflow.com/questions/35097815/vectorized-equality-testing
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}
# 3. Running simulations --------------------------------------------------


# starting states
# f/m indicates sex (m= male, f= female)
# 1 indicates size class where 1 is smalles 3 is largest
# i/u indicates infection status (i = infected, u = uninfected)
states_t <-
  c(
    f1i = 0.01,
    f1u = 0.99,
    m1i = 0.01,
    m1u = 0.99,
    f2i = 0,
    f2u = 0,
    m2i = 0,
    m2u = 0,
    f3i = 0,
    f3u = 0,
    m3i = 0,
    m3u = 0
  )


# *3.a.1 HE no CI, no MK ----------------------------------------------------------
# paramaters
# equation 2 Levitan 1991 Biol. Bull. 181:261-268
# using this to estimate deduction
# Ded1 60mm vs 30
Ded1he <- (exp(3.56) * 30 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
# Ded2 60mm vs 50
Ded2he <- (exp(3.56) * 50 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
Heparms <-
  c(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0,
    Mk = 0
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HE <-
  ode(
    y = states_t,
    times = 0:60000,
    func = pop_next,
    parms = Heparms,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )
# plot densities of each subpop
Hep <- HE %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
Hep
# get sex ratio at each stage
HESR <- HE %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )
# make sex ratio plot
HESRP <- ggplot(HESR, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HESRP
# calculate infection status
HEI <- HE %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")
# plot infection
HEIP <- ggplot(HEI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HEIP
# *3.a.2 HE CI ----------------------------------------------------------
# paramaters
HeparmsCI <-
  c(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0.7,
    Mk = 0
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HECI <-
  ode(
    y = states_t,
    times = 0:60000,
    func = pop_next,
    parms = HeparmsCI,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )
HeCIP <- HECI %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HeCIP
# get sex ratio
HESRCI <- HECI %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HESRCIP <- ggplot(HESRCI, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HESRCIP
# infection status
HEICI <- HECI %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HEICIP <- ggplot(HEICI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HEICIP
# *3.a.3 HE no CI but MK ----------------------------------------------------------
# paramaters
HeparmsMk <-
  c(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0,
    Mk = 0.7
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HEMk <-
  ode(
    y = states_t,
    times = 0:60000,
    func = pop_next,
    parms = HeparmsMk,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )

HepMk <- HEMk %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HepMk
# calculate sex ratio
HESRMK <- HEMk %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HESRPMK <- ggplot(HESRMK, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HESRPMK
# get infection status
HEIMK <- HEMk %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HEIMKP <- ggplot(HEIMK, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HEIMKP

# *3.a.3 HE CI AND MK ----------------------------------------------------------
# paramaters
HeparmsMkCI <-
  c(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0.7,
    Mk = 0.7
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HEMkCI <-
  ode(
    y = states_t,
    times = 0:60000,
    func = pop_next,
    parms = HeparmsMkCI,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )

HepMkCI <- HEMkCI %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HepMkCI
HESRMKCI <- HEMkCI %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HESRPMKCI <- ggplot(HESRMKCI, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HESRPMKCI

HEIMKCI <- HEMkCI %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HEIMKPCI <-
  ggplot(HEIMKCI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HEIMKPCI

# 4. Putting dataframes together -----------------------------------------------
# Adding new column to keep track of what model data belong to
HEIMK$Model <- "Enhanced growth + MK"
HESRMK$Model <- "Enhanced growth + MK"
HEIMKCI$Model <- "Enhanced growth + MK + CI"
HESRMKCI$Model <- "Enhanced growth + MK + CI"
HEICI$Model <- "Enhanced growth + CI"
HESRCI$Model <- "Enhanced growth + CI"
HEI$Model <- "Enhanced growth"
HESR$Model <- "Enhanced growth"
HEMk$Model <- "Enhanced growth + MK"
HEMkCI$Model <- "Enhanced growth + MK + CI"
HECI$Model <- "Enhanced growth + CI"
HE$Model <- "Enhanced growth"
# putting the dataframes together for sex ratio
SRall <- rbind(HESRMK, HESR, HESRMKCI)
# ordering the way I want to
SRall$Model <-
  factor(
    SRall$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )
# putting dataframes together for infection %
Iall <- rbind(HEIMK, HEI, HEIMKCI)
Iall$Model <-
  factor(
    Iall$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )
# putting dataframes together for density
Dall <- rbind(HEMk, HE, HEMkCI)
Dall$Model <-
  factor(
    Dall$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )

# 5.a Making Figure 3 -----------------------------------------------------
# putting infection rate plot together
Iallp <- ggplot(Iall, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Time") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#C783AA", "#666A9E", "#2F4E2F"),
    name = "Size Class"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  scale_y_continuous(label = percent) +
  facet_grid(. ~ Model)
# putting sex ratio plots together
SRallp <- ggplot(SRall, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Time") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#C783AA", "#666A9E", "#2F4E2F"),
    name = "Size Class"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  facet_grid(. ~ Model) +
  ylim(0.2, 0.8)
# Figure 3
Allp <-
  ((SRallp / Iallp) &
    theme(legend.position = "right")) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
Allp

# ggsave("Plots/Stage_based_4_29.pdf",Allp,height = 140,width=190,units="mm")
# ggsave("Plots/Stage_based_4_29.svg",Allp,height = 140,width=190,units="mm")


# 5.b Making Figure A7 ----------------------------------------------------
# Putting Density plots together
DallP <-
  Dall %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Time", color = "Subpopulations") +
  scale_x_continuous(
    label =
      comma, limits = c(20, 60000)
  ) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  facet_grid(. ~ Model) +
  ylim(0, 0.35)
DallP
# ggsave("Plots/Stage_based_Dens_4_28_21.pdf",DallP,height = 140,width=190,units="mm")



# 6. Sensitivity analysis for Figure A8 ------------------------------------------------
# can skip simulation as it takes a few hours to run by reading in the result file.
# starting states
states_t <-
  c(
    f1i = 0.01,
    f1u = 0.99,
    m1i = 0.01,
    m1u = 0.99,
    f2i = 0,
    f2u = 0,
    m2i = 0,
    m2u = 0,
    f3i = 0,
    f3u = 0,
    m3i = 0,
    m3u = 0
  )

# paramaters
# equation 2 Levitan 1991 Biol. Bull. 181:261-268
# using this to estimate deduction
# Deduction in gamate production 60mm vs 30
Ded1he <- (exp(3.56) * 30 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
# Deduction in gamate production 60mm vs 30
Ded2he <- (exp(3.56) * 50 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
# vector of male killing values to iterate through
Mks <- seq(0, 0.99, 0.01)
# vector of enhanced growth rates to iterate through.
Tifs <- seq(0.11, 0.35, 0.01)
# registering cores to perform calculations
registerDoParallel(10)
# simulatee until equilibrium across all male killing rates and many enhanced growth rates.
m <- foreach(mk = Mks, .combine = rbind) %:%
  foreach(tf = Tifs, .combine = rbind) %dopar% {
    states_t <-
      c(
        f1i = 0.01,
        f1u = 0.99,
        m1i = 0.01,
        m1u = 0.99,
        f2i = 0,
        f2u = 0,
        m2i = 0,
        m2u = 0,
        f3i = 0,
        f3u = 0,
        m3i = 0,
        m3u = 0
      )
    Heparms <-
      c(
        sigma = 0.203,
        v = 0.125,
        Fe = 0.09444,
        tb = 1,
        S = 700,
        E = 0.05,
        tau = 5400,
        c = 0,
        ded1 = Ded1he,
        ded2 = Ded2he,
        # fertilization parms
        tif = tf,
        t = 0.1,
        K = 100,
        d = 0.1,
        Ma = 0.99,
        Ml = 0.66,
        CI = 0,
        Mk = mk
      )
    # actual simulation
    HES <-
      ode(
        y = states_t,
        times = 0:20000,
        func = pop_next,
        parms = Heparms,
        method = "iteration"
      ) %>% data.frame()
    # iterating at increments of 20,000 generations
    # counter keeps track of how many 20,000 increments we've done.
    counter <- 1
    # stop when there is no change
    while (!all(is_equal_tol(HES[HES$time == 19999, ][2:13], HES[HES$time ==
      20000, ][2:13]))) {
      # keep going until max time
      states_temp <- as.numeric(HES[HES$time == 20000, ][2:13])
      names(states_temp) <- names(HES[HES$time == 20000, ][2:13])
      HES <-
        ode(
          y = states_temp,
          times = 0:20000,
          func = pop_next,
          parms = Heparms,
          method = "iteration"
        ) %>% data.frame()
      counter <- counter + 1
    }
    Temp <-
      rbind(HES[HES$time == 20000, ][2:13]) %>% mutate(
        Generation = 20000 * counter,
        Mk = mk,
        tif = tf
      )
    return(Temp)
  }
# save output
# write.csv(m,"Data/Sensitivity_Stage_Based.csv")
m <- read.csv("Data/Sensitivity_Stage_Based.csv")
glimpse(m)
# equilibrium sexratio at different stages
srm <-
  m %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  ) %>%
  group_by(Mk, tif, Generation) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )
# equilbirum infection rate
inf <-
  m %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  ) %>%
  group_by(Mk, tif, Generation) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

# 7. Plotting Figure A8 -------------------------------------------------
mytheme <-
  theme_classic() + theme(
    legend.position = "right",
    # this puts legend on the bottom
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(size = 12, color = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) # makes all the text larger and bold
theme_set(mytheme)
# labels for stages
srlabels <- c("20 - 39 mm", "40 - 59 mm", "60+ mm")
names(srlabels) <- c("sexratio1", "sexratio2", "sexratio3")
# labels for infectiioin
Ilabels <- c("20 - 39 mm", "40 - 59 mm", "60+ mm")
names(Ilabels) <- c("I1", "I2", "I3")
# divide enhanced growth rate by 0.1 to convert it into relative enhanced growth rate.
srm$tif2 <- srm$tif / 0.1

# make the actual plot.
Srp <- ggplot(srm, aes(x = Mk, y = tif2, fill = `Sex Ratio`)) +
  geom_tile() +
  facet_wrap(vars(Stage), labeller = labeller(Stage = srlabels)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_distiller(
    palette = "PuOr",
    na.value = "grey",
    direction = -1,
    limits = c(0, 1),
    labels = percent
  ) +
  ylab("Relative Enhanced Growth") +
  xlab("Male Killling Rate") +
  theme(aspect.ratio = 1)

Srp
# ggsave("SI_egmk.png",Srp,height = 150,width=183,units="mm",path="Plots")


# 8.a.1 HT no CI, no MK ----------------------------------------------------------
# paramaters
# equation 2 Levitan 1991 Biol. Bull. 181:261-268
# using this to estimate deduction
# Ded1 60mm vs 30
Ded1he <- (exp(3.56) * 30 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
# Ded2 60mm vs 50
Ded2he <- (exp(3.56) * 50 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
Htparms <-
  c(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0,
    Mk = 0
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HT <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = Htparms,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )
# plot densities of each subpop
Htp <- HT %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
Htp
# get sex ratio at each stage
HTSR <- HT %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )
# make sex ratio plot
HTSRP <- ggplot(HTSR, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HTSRP
# calculate infection status
HTI <- HT %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")
# plot infection
HTIP <- ggplot(HTI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HTIP
# *8.a.2 HT CI ----------------------------------------------------------

# paramaters
HtparmsCI <-
  c(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0.7,
    Mk = 0
  )
# pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HTCI <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HtparmsCI,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )
HtCIP <- HTCI %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HtCIP
# get sex ratio
HTSRCI <- HTCI %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HTSRCIP <- ggplot(HTSRCI, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HTSRCIP
# infection status
HTICI <- HTCI %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HTICIP <- ggplot(HTICI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HTICIP
# *8.a.3 HT no CI but MK ----------------------------------------------------------
# paramaters

HtparmsMk <-
  c(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0,
    Mk = 0.7
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HTMk <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HtparmsMk,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )

HtpMk <- HTMk %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HtpMk
# calculate sex ratio
HTSRMK <- HTMk %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HTSRPMK <- ggplot(HTSRMK, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HTSRPMK
# get infection status
HTIMK <- HTMk %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HTIMKP <- ggplot(HTIMK, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HTIMKP

# *8.a.3 HT CI AND MK ----------------------------------------------------------
# paramaters
HtparmsMkCI <-
  c(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400,
    c = 0,
    ded1 = Ded1he,
    ded2 = Ded2he,
    # fertilization parms
    tif = 0.25,
    t = 0.1,
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0.7,
    Mk = 0.7
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HTMkCI <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HtparmsMkCI,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )

HtpMkCI <- HTMkCI %>%
  separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
  ggplot(aes(
    x = time,
    y = Density,
    color = Class,
    linetype = Infection
  )) +
  geom_line() +
  # change labels of our axis
  labs(y = "Density", x = "Generation", color = "Subpopulations") +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(label = comma) +
  scale_linetype(labels = c("Infected", "Uninfected")) +
  scale_color_manual(
    labels = c(
      "Female 20 - 39 mm",
      "Female 40 - 59 mm",
      "Female 60+ mm",
      "Male 20 - 39 mm",
      "Male 40 - 59 mm",
      "Male 60+ mm"
    ),
    values = c(
      "#564157",
      "#796FBD",
      "#13D1E2",
      "#FF0015",
      "#B2D75E",
      "#ACD9C3"
    )
  )
HtpMkCI
HTSRMKCI <- HTMkCI %>%
  group_by(time) %>%
  summarise(
    sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) /
      (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population ==
        "m1i"] + Density[Population == "m1u"]),
    sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) /
      (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population ==
        "m2i"] + Density[Population == "m2u"]),
    sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) /
      (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population ==
        "m3i"] + Density[Population == "m3u"])
  ) %>%
  pivot_longer(c(sexratio1, sexratio2, sexratio3),
    names_to = "Stage",
    values_to = "Sex Ratio"
  )

HTSRPMKCI <- ggplot(HTSRMKCI, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  ylim(0.2, 0.8)
HTSRPMKCI

HTIMKCI <- HTMkCI %>%
  group_by(time) %>%
  summarise(
    I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population ==
      "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population ==
      "m1u"]),
    I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population ==
      "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population ==
      "m2u"]),
    I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population ==
      "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population ==
      "m3u"])
  ) %>%
  pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")

HTIMKPCI <-
  ggplot(HTIMKCI, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Generation") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#13D1E2", "#796FBD", "#FF0015"),
    name = "Size Class"
  ) +
  theme(legend.position = "right", aspect.ratio = 1) +
  scale_x_continuous(
    label =
      comma
  ) +
  scale_y_continuous(label = percent)
HTIMKPCI

# 9. Putting dataframes together -----------------------------------------------
# Adding new column to keep track of what model data belong to
HTIMK$Model <- "Enhanced growth + MK"
HTSRMK$Model <- "Enhanced growth + MK"
HTIMKCI$Model <- "Enhanced growth + MK + CI"
HTSRMKCI$Model <- "Enhanced growth + MK + CI"
HTICI$Model <- "Enhanced growth + CI"
HTSRCI$Model <- "Enhanced growth + CI"
HTI$Model <- "Enhanced growth"
HTSR$Model <- "Enhanced growth"
HTMk$Model <- "Enhanced growth + MK"
HTMkCI$Model <- "Enhanced growth + MK + CI"
HTCI$Model <- "Enhanced growth + CI"
HT$Model <- "Enhanced growth"
# putting the dataframes together for sex ratio
SRallT <- rbind(HTSRMK, HTSR, HTSRMKCI)
# ordering the way I want to
SRallT$Model <-
  factor(
    SRallT$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )
# putting dataframes together for infection %
IallT <- rbind(HTIMK, HTI, HTIMKCI)
IallT$Model <-
  factor(
    IallT$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )
# putting dataframes together for density
DallT <- rbind(HTMk, HT, HTMkCI)
DallT$Model <-
  factor(
    DallT$Model,
    levels = c(
      "Enhanced growth",
      "Enhanced growth + MK",
      "Enhanced growth + MK + CI"
    )
  )

# 10.a Making Figure A9 -----------------------------------------------------
# putting infection rate plot together
IallpT <- ggplot(IallT, aes(x = time, y = `Infection Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Infection %", x = "Time") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#C783AA", "#666A9E", "#2F4E2F"),
    name = "Size Class"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  scale_y_continuous(label = percent) +
  facet_grid(. ~ Model)
# putting sex ratio plots together
SRallpT <- ggplot(SRallT, aes(x = time, y = `Sex Ratio`, color = Stage)) +
  geom_line() +
  # change labels of our axis
  labs(y = "Sex Ratio", x = "Time") +
  scale_color_manual(
    labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"),
    values = c("#C783AA", "#666A9E", "#2F4E2F"),
    name = "Size Class"
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  facet_grid(. ~ Model) +
  ylim(0.2, 0.8)
# Figure 3
AllpT <-
  ((SRallpT / IallpT) &
    theme(legend.position = "right")) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
AllpT

# ggsave("Stage_based_11_17_HT_A9.pdf",AllpT,height = 140,width=190,units="mm")
# ggsave("Stage_based_11_17_HT_A9..svg",Allp,height = 140,width=190,units="mm")

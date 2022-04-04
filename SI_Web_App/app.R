#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# Supplemental web app for "On the spread of microbes that manipulate reproduction in marine invertebrates" Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282
# 1. Setting up -----------------------------------------------------------
rm(list = ls())

#* 1.a Loading up libraries ----------------------------------------------
library(shiny)
library(plotly)
library(scales)
library(tidyverse)
library(patchwork)
library(viridis)
library(feather)
library(shinythemes)
library(deSolve)

#* 1.b Setting up plotting themes ----------------------------------------


mytheme <-
  theme_linedraw() + theme(
    legend.position = "top",
    # this puts legend on the bottom
    axis.title = (element_text(face = "bold")),
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(
      size = 18, face = "bold", color =
        "black"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    )
  ) # makes all the text larger and bold
theme_set(mytheme)

# 2. Functions ---------------------------------------------------------
# sigma=egg size (cros section of the area mm^2).
# v = sperm velocity(mm/sec^-1);
# Fe=Fertilization efficiency (5% as a starting point do we have an estimate?)
# S=Sprem density per unit male density (uL-1)
# E=Egg densiity pere unit female density (uL-1)
# tau=sperm half life (substitued with t, time(s), eggs are exposed to serm it t<tau)...apparently dependent on sperm density...
# tb=time for a block to polyspermy to be established. 0 assumes perfect blocking
# Fertilization dynamics function

# * 2.a Styan Fert functions ----------------------------------------------

# * 2.b Sex ratio Fert functions ------------------------------------------
prob_monof <- function(Dens, sigma, v, Fe, tb, sr, S, E, tau) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
  Males <- Dens * sr
  S0 <- Males * S
  Females <- Dens * (1 - sr)
  # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
  E0 <- E * Females
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizaing sperm
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b)))
  return(prop_mono)
}

# * 2.c Number of zygotes produced ----------------------------------------
number_zygotesf <- function(Dens, sigma, v, Fe, tb, sr, S, E, tau) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
  Males <- Dens * sr
  S0 <- Males * S
  Females <- Dens * (1 - sr)
  # E0<-S0*E*(Females/Males)#Egg density should be proportional to #females
  E0 <- E * Females
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizaing sperm
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b)))
  return(prop_mono * E0)
}

# plots probability of fertilization
plot_functionf <- function(sigma, v, Fe, tb, S, E, tau) {
  Data <- data.frame(Dens = 0)
  p <-
    ggplot(data = Data, mapping = aes(Dens = Dens)) +
    xlab(bquote(bold('Density (individuals per m'^2*")"))) +
    ylab("Fertilization (%)") +
    labs(color = "Functions") +
    scale_x_continuous(
      trans = "log10",
      limits = c(0.00001, 100),
      breaks = c(0.0001, 0.01, 1, 100),
      labels = c("0.0001", "0.01", "1", "100")
    ) +
    scale_y_continuous(labels = label_percent())
  p + stat_function(
    fun = prob_monof,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.5,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.5")
  ) + stat_function(
    fun = prob_monof,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.7,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.7")
  ) + stat_function(
    fun = prob_monof,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.9,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.9")
  ) + stat_function(
    fun = prob_monof,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.3,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.3")
  ) + stat_function(
    fun = prob_monof,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.1,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.1")
  ) + scale_color_manual(
    name = "Sex ratio",
    breaks = c("0.1", "0.3", "0.5", "0.7", "0.9"),
    values = viridis(5)
  ) + guides(color = guide_legend(override.aes = list(size = 2))) #+ggtitle(paste0("Egg Diameter Size = ", (sqrt(sigma*4/pi)*1000)," mm;", "# Eggs = ", E))
}

# * 2.d plotting function for Zygote density ------------------------------
plot_function2f <- function(sigma, v, Fe, tb, S, E, tau) {
  Data <- data.frame(Dens = 0)
  p <-
    ggplot(data = Data, mapping = aes(Dens = Dens)) +
    xlab(bquote(bold('Density (individuals per m'^2*")"))) +
    ylab((bquote(bold("Zygote density (zygotes per"~ mu*L~")")))) +
    labs(color = "Functions") +
    scale_x_continuous(
      trans = "log10",
      limits = c(0.00001, 100),
      breaks = c(0.0001, 0.01, 1, 100),
      labels = c("0.0001", "0.01", "1", "100")
    )
  p + stat_function(
    fun = number_zygotesf,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.5,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.5")
  ) + stat_function(
    fun = number_zygotesf,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.7,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.7")
  ) + stat_function(
    fun = number_zygotesf,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.3,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.3")
  ) + stat_function(
    fun = number_zygotesf,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.9,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.9")
  ) + stat_function(
    fun = number_zygotesf,
    args = (list(
      sigma = sigma,
      v = v,
      Fe = Fe,
      tb = tb,
      sr = 0.1,
      S = S,
      E = E,
      tau = tau
    )),
    mapping = aes(color = "0.1")
  ) + scale_color_manual(
    name = "Sex ratio",
    breaks = c("0.1", "0.3", "0.5", "0.7", "0.9"),
    values = viridis(5)
  ) + ggtitle(paste0("Egg Diameter Size = ", (sqrt(sigma * 4 / pi) * 1000), " mm;", " # Eggs = ", E)) +
    theme(legend.position = "none")
}

#  * 2.e Larval death rate ------------------------------------------------


# * 2.f Population Simulation functions -----------------------------------
number_zygotes <- function(Males, sigma, v, Fe, tb, sr, S, E, tau, c, nif) {
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
  return(prop_mono * E0)
}
# Population dynamic functions
# Nz is number of zygotes
# Nif is number of in infected female
# Nuf is number of uninfected females
# Nim is number of infected males
# Num is number of uninfected males
# K is adult carrying capacity (individuals/m^2)
# d is mortality rate of changee
# Ma is base adult mortality
# ml is larval mortality
# R is feminzation rate
# Mk is male killing mortality
# Model asks what kind of model
pop_next <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Nta <- fi + fu + mi + mu
    sr <- (mi + mu) / Nta
    Nz <-
      number_zygotes(
        Males = (mi + mu),
        sigma = sigma,
        v = v,
        Fe = Fe,
        tb = tb,
        sr = sr,
        S = S,
        E = E,
        tau = tau,
        c = c,
        nif = fi
      )*psi
    # total sperm
    st <- mi + mu
    # total Egg
    et <- (fi * (1 - c)) + fu
    # infected egg
    ei <- (fi * (1 - c))
    # uninfected egg
    eu <- fu
    # infected sperm
    si <- mi
    # uninfected sperm
    su <- mu
    Nzim <- (Nz * (ei / et) * (1 - r)) * (1 - Ml) * (1 - Mk)
    Nzum <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -
      Ml)
    Nzif <- (Nz * (ei / et) * r) * (1 - Ml)
    Nzuf <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -
      Ml)
    Nzt <- Nzim + Nzum + Nzif + Nzuf
    # total adult density
    # death
    # mortality
    mort <- (1 - Ma / (1 + exp(-d * (Nta + Nzt - K / 2))))
    fi_plus1 <- (Nzif + fi) * mort
    fu_plus1 <- (Nzuf + fu) * mort
    mi_plus1 <- (Nzim + mi) * mort
    mu_plus1 <- (Nzum + mu) * mort
    list(c(fi_plus1, fu_plus1, mi_plus1, mu_plus1))
  })
}
# original simulation


# * 2.g Stage based simulation --------------------------------------------
number_zygotesS <-
  function(f1i,
           f1u,
           m1i,
           m1u,
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
           ded2) {
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
    return(prop_mono * E0)
  }
# birth -> death -> transition
# function to run simulations
pop_nextS <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    Nz <-
      number_zygotesS(
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
      )*psi
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
states_tS <-
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
Ded1he <- (exp(3.56) * 30 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
# Ded2 60mm vs 50
Ded2he <- (exp(3.56) * 50 - exp(6.59)) / (exp(3.56) * 60 - exp(6.59))
HeparmsS <-
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
    Mk = 0.7,
    psi=1
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
Data_1S <-
  ode(
    y = states_tS,
    times = 0:30000,
    func = pop_nextS,
    parms = HeparmsS,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(
    c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u),
    names_to = "Population",
    values_to = "Density"
  )
# * 2. h loading up data --------------------------------------------------
states_t <- c(
  fi = 0.01,
  fu = 0.99,
  mi = 0.01,
  mu = 0.99
)
# *2.a.1 HE no CI, no MK, Fem ----------------------------------------------------------
# paramaters
Heparmsf <-
  c(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400,
    c = 0,
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0,
    Mk = 0,
    r = 0.7,
    psi=1
  )

Data_1 <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = Heparmsf,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  )

# * adult mortality -------------------------------------------------------
death <- function(Ma, d, Nt, k) {
  return((1 - Ma / (1 + exp(-d * (
    Nt - k / 2
  )))))
}

#  * larval mortality -----------------------------------------------------

deathl <- function(Mz, t) {
  return((1 - exp(-Mz * t)))
}
# figure 4 A
Datat <- read_feather("data/Data_3_3_21_01_99_all.feather")
Datat <- Datat %>%
  mutate(Egg_Size = (((Egg_Size * 4) / pi)^0.5) * 1000) %>%
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
Data2t <- Datat %>% complete(Egg_Size, Feminization, Larval_d)
species <- read.csv("data/species2.csv")
# Figure 4 B
DatatB <- read_feather("data/Data_4b_8_22_21.feather")
DatatB <- DatatB %>%
  mutate(Egg_Size = (((Egg_Size * 4) / pi)^0.5) * 1000) %>%
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
Data2tB <- DatatB %>% complete(Egg_Size, Feminization, Larval_d)
# species data
species <- read.csv("data/species2.csv")


eggs <- as.numeric(unique(Data2t$Egg_Size))
days <- as.numeric(unique(Data2t$Larval_d))
which.min(abs(-eggs))
species$EggSize
matchclose <- function(x, vec) {
  vec[which.min(abs(x - vec))]
}
species$EggSize2 <- sapply(species$EggSize, FUN = matchclose, vec = eggs)
species$Days <- sapply(species$PlanktonicTime, FUN = matchclose, vec = days)

foo <- function(data) {
  rowname <- rep(0, length(row.names(data)))
  for (i in 1:length(rowname)) {
    rowname[i] <-
      paste(colnames(data)[c(1, 3, 4, 8, 9)], ":", data[i, c(1, 3, 4, 8, 9)], collapse = "\n")
  }
  return(rowname)
}
species$labls <- foo(species)

# Define UI for application that draws a histogram

# 3. User interaction section ---------------------------------------------


ui <- fluidPage(
  theme = shinytheme("superhero"), class = "home",
  # this code should make it look allright on a phone.
  HTML('<meta name="viewport" content="width=1024">'),
  navbarPage(
    "", # need to make this to say we are having multiple tabs

    # * 3.a Homepage ----------------------------------------------------------
    tabPanel(
      "Homepage",
      fluidRow(
        titlePanel(h1(
          "Supporting information for:",
          align = "center"
        ))
      ), fluidRow(h1("'On the spread of microbes that manipulate reproduction in marine invertebrates'", align = "center")),
      fluidRow(h3("Matthew C. Kustra & Tyler J. Carrier", align = "center")),
      fluidRow(
        column(
          12,
          h2("Abstract:"),
          br(),
          p(
            "Bacterial symbionts are functionally integral to animal reproduction and development, some of which have evolved additional mechanisms to override these host programs. One habitat that is increasingly recognized to contain phylogenetically related lineages of reproductive manipulators is the ocean. The reproduction of marine invertebrates often occurs by free-spawning instead of by the physical contact of copulation in terrestrial systems. We developed an integrated model to understand whether and when microbes that manipulate host reproduction by cytoplasmic incompatibility, feminization, and male killing spread within populations of free-spawning marine invertebrates. Our model support three primary findings. First, sex ratio distortion leads to suboptimal fertilization and zygote production in planktotrophs (feeding larvae), but enhance these processes in lecithotrophs (nonfeeding larvae). Second, feminization and a combination of male killing plus enhanced growth are effective at spreading reproductive manipulators while also inducing a female-biased sex ratio. Third, the majority of free-spawning marine invertebrates could be infected across a range of life-history combinations, with infections harming species with smaller eggs and longer pelagic durations while benefiting species with larger eggs and shorter pelagic durations. Together, this supports the general premise that microbes may manipulate the reproduction of free-spawning marine invertebrates (e.g., by inducing changes in developmental life-history) and that these types of manipulations overlap considerably with terrestrial systems."
          )
        )
      ),
      fluidRow(
        column(
          12, br(),
          h2("Description:"),
          br(),
          p(
            "This is a supplemental web application for the paper: 'On the spread of microbes that manipulate reproduction in marine invertebrates' (https://www.journals.uchicago.edu/doi/10.1086/720282). At each tab, one can make supplemental figures using the model from this paper. One can get to different pages of the web application by using the tabs above. Description of the figures are given on each page, and short descriptions of each tab are given below."
          )
        )
      ),
      fluidRow(
        column(
          4,
          br(),
          h3("Fertilization dynamics (Eqs. 4 & 5)"),
          br(),
          p(
            "Here one can generate graphs of how species reproductive biology (e.g. egg size) interact with density and sex ratio to affect fertilization dynamics. Equations 4 and 5 in the paper used to generate these figures were adopted from the Styan (1998) polyspermy model."
          ),
          br(),
          h3("Larval mortality (Eq. 6)"),
          br(),
          p(
            "Here one can change the instantaneous mortality rate to see how the length of larval duration determines larval mortality used in the population models. Equation 6 of the paper is used to generate the plot and is adapted from Rumrill (1990)."
          )
        ),
        column(
          4,
          br(),
          h3("Adult mortality (Eq. 7)"),
          br(),
          p(
            "Here one can see how different parameters in Equation 7 of the paper influence the density dependence of adult mortality."
          ),
          br(),
          h3("Simulations: no population structure"),
          br(),
          p(
            "Here one can simulate infection dynamics of a population with no size structure under different parameter values (Equations 8 - 11)."
          )
        ),
        column(
          4, br(),
          h3("Simulations: population structure"),
          br(),
          p(
            "Here one can simulate infection dynamics of a size structured population under different parameter values (Equations 12 - 20)."
          ), br(),
          h3("Marine invertebrates (Fig. 4)"),
          br(),
          p(
            "Here one can see which known marine invertebrate species are predicted to be infected based on Figure 4 in the main paper."
          )
        )
      )
    ),

    # * 3.d Fertilzation Dynamics Sex ratio ----------------------------------
    tabPanel(
      "Fertilization dyanmics (Eqs. 4 & 5)", # name of third tab
      titlePanel("Simulate fertilzation dyanmics here"), # gives us the title of the webpage
      # Sidebarlayout
      sidebarLayout( # saying to use a sidebarLayout. T
        sidebarPanel(
          h3("Description:"), p("Here you can explore how sex ratio (proportion males), density, and various reproductive traits (e.g. egg size) affect fertilization dynamics in broadcast spawners. Enter in values above and watch how these parameters affect the probability of fertilization and density of zygotes production. These graphs are produced using Equations 4 and 5 from the main paper based on the polyspermy model from Styan (1998). Initial graphs produced use species-specific values for Heliocidaris tuberculata and are the same as the top row of Figure 1 from the main paper."), br(), # usually where you put controls that user interacts with
          numericInput("Sigmaf",
            "Egg size (mm^2)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 2, # max of slider input
            value = 0.017
          ),
          numericInput("vf",
            "Sperm velocity (mm per sec)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 2, # max of slider input
            value = 0.14
          ),
          # Fe=0.05,tb=1,S=700,E=0.05,tau=5400,generations=5000
          numericInput("Fef",
            "Fertilization efficiency constant", # what you are displaying
            min = 0, # miniumum of slider input
            max = 0.5, # max of slider input
            value = 0.09444
          ),
          numericInput("tbf",
            "Time to block polysmery (s)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 800, # max of slider input
            value = 1
          ),
          numericInput("Sf", # n
            "Sperm/uL/male/m^2", # what you are displaying
            min = 0, # miniumum of slider input
            max = 10000, # max of slider input
            value = 700
          ), # starting value
          numericInput("Ef",
            "Egg/uL/female/m^2", # what you are displaying
            min = 0, # miniumum of slider input
            max = 100, # max of slider input
            value = 66.59171429
          ),
          numericInput("tauf",
            "Sperm half life (s)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 10000, # max of slider input
            value = 5400
          )
        ),
        mainPanel( # usually where the plot or main thing is
          plotOutput("fertf", height = 800) # where
        )
      )
    ),


    # * 3.e1 Larval Mortality -------------------------------------------------
    tabPanel(
      "Larval mortality (Eq. 6)", # name of the tab
      titlePanel("Larval mortality"), # gives us the title of the webpage

      # Sidebar with a slider input for number of bins
      sidebarLayout( # saying to use a sidebarLayout. Two parts: sidebarPanel and mainPanel
        sidebarPanel(
          h3("Description:"), p("Here you can explore how the pelagic larval duration relates to larval mortality by varying the instantaneous mortality rate parameter (Rumrill 1990). Starting graph shows the instantaneous mortality rate used in the paper."), br(), # usually where you put controls that user interacts with
          numericInput("Mzl", # name of the input that you access in server wiht input$xxx
            "Instantenous mortality rate (0 to 1)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 1, # max of slider input
            value = 0.226,
            step = 0.05
          )
        ),
        mainPanel( # usually where the plot or main thing is
          plotOutput("al_func") # where you specify the name of the object you are plotting
        )
      )
    ),


    # * 3.e Adult mortality --------------------------------------------------
    tabPanel(
      "Adult mortality (Eq. 7)", # name of the tab
      titlePanel("Adult mortality"), # gives us the title of the webpage

      # Sidebar with a slider input for number of bins
      sidebarLayout( # saying to use a sidebarLayout. Two parts: sidebarPanel and mainPanel
        sidebarPanel(
          h3("Description:"), p("Here you can explore how the d parameter in the model affects density-dependent mortality in the model. Starting graph shows the d value used for numerical simulations in the model."), br(), # usually where you put controls that user interacts with
          numericInput("Ma", # name of the input that you access in server wiht input$xxx
            "Max adult mortality (0 to 1)", # what you are displaying
            min = 0, # miniumum of slider input
            max = 1, # max of slider input
            value = 0.99,
            step = 0.05
          ), # starting value
          numericInput("d", # name of the input that you access in server wiht input$xxx
            "Shape parameter (0 to 1):", # what you are displaying
            min = 0, # miniumum of slider input
            max = 1, # max of slider input
            value = 0.1,
            step = 0.05
          ),
          numericInput("k", # name of the input that you access in server wiht input$xxx
            "Carrying capacity (individuals per m^2):", # what you are displaying
            min = 1, # miniumum of slider input
            max = 200,
            step = 10, # max of slider input
            value = 100
          )
        ),
        mainPanel( # usually where the plot or main thing is
          plotOutput("am_func") # where you specify the name of the object you are plotting
        )
      )
    ),

    # * 3.f Simulations ------------------------------------------------------


    tabPanel(
      "Simulations: no population structure", # name of third tab
      titlePanel("Simulate infection dyanmics here"), # gives us the title of the webpage
      # Sidebarlayout
      sidebarLayout( # saying to use a sidebarLayout. T
        sidebarPanel(
          h3("Description:"), p("Here you can simulate the model to see how infection may spread in a marine invertebrate species without population structure under different parameter combinations (press 'Run Simulation' to make graphs). (A) the density of different sub-populations, (B) the proportion of males in the population, and (C) the percentage of infected individuals.  The simulations look at three potential mechanisms: cytoplasmic incompatibility, feminization, and male-killing. First, use the checkboxes to display what parameters you want to change, then change parameters to your choosing, and finally click run simulation. Please note that some subpopulations may overlap; hence they do not appear on the figure (i.e. uninfected females and uninfected males). Simulations can take a few minutes to run on this website. We recommend downloading the R code from dryad to run the model when the number of generations exceeds 20000. Initial values and graphs generated are for Heliocidaris erythrogramma with feminization only and are in the text: (A) Figure A6 center bottom panel, (B) Figure 2 center top panel, (C) Figure 2 center bottom panel."), br(), # usually where you put controls that user interacts with
          checkboxGroupInput("parms",
            h3("What parameters do you want to change"),
            choices = list("Starting densities" = "dens", "Fertilization" = "ferts", "Population" = "pops", "Infection" = "inf"), inline = TRUE
          ),
          conditionalPanel(
            condition = "input.parms.includes('dens')",
            numericInput("Nif0", # n
              "Starting density of infected females (individuals per m^2)",
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.01, step = 0.05
            ), # starting value
            numericInput("Nuf0",
              "Starting density of uninfected females (individuals per m^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.99, step = 0.05
            ),
            numericInput("Nim0",
              "Starting density of infected males (individuals per m^2)",
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.01, step = 0.05
            ),
            numericInput("Num0", # n
              "Starting density of uninfected males (individuals per m^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.99, step = 0.05
            )
          ), # starting value
          conditionalPanel(
            condition = "input.parms.includes('ferts')",
            numericInput("Sigma",
              "Egg size (mm^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 2, # max of slider input
              value = 0.203, step = 0.05
            ),
            numericInput("v",
              "Sperm velocity (mm per sec)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 2, # max of slider input
              value = 0.125, step = 0.05
            ),
            # Fe=0.05,tb=1,S=700,E=0.05,tau=5400,generations=5000
            numericInput("Fe",
              "Fertilization efficiency constant", # what you are displaying
              min = 0, # miniumum of slider input
              max = 0.5, # max of slider input
              value = 0.09444, step = 0.01
            ),
            numericInput("tb",
              "Time to block polysmery (s)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 800, # max of slider input
              value = 1, step = 10
            ),
            numericInput("S", # n
              "Sperm/uL/male/m^2", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1000, # max of slider input
              value = 700, step = 20
            ), # starting value
            numericInput("E",
              "Egg/uL/female/m^2", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.05, step = 1
            ),
            numericInput("tau",
              "Sperm half life (s)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 10000, # max of slider input
              value = 5400, step = 100
            )
          ),
          conditionalPanel(
            condition = "input.parms.includes('pops')",
            numericInput("ks",
              "Carrying capacity (individuals per m^2)", # what you are displaying
              min = 10, # miniumum of slider input
              max = 200, # max of slider input
              value = 100, step = 10
            ),
            numericInput("ds",
              "Shape parameter for adult mortality", # what you are displaying
              min = 0, # miniumum of slider input
              max = 10, # max of slider input
              value = 0.1, step = 0.05
            ),
            numericInput("Ma", # n
              "Max adult mortality",
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.99, step = 0.05
            ), # starting value
            numericInput("Ml",
              "Larval mortality", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.66, step = 0.05
            ),
            numericInput("psi",
                         "Settlement constant (uL per m^2)", # what you are displaying
                         min = 0, # miniumum of slider input
                         max = 10, # max of slider input
                         value = 1, step = 0.05
            )
          ),
          conditionalPanel(
            condition = "input.parms.includes('inf')",
            numericInput("R",
              "Rate of Feminization", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.7, step = 0.05
            ),
            numericInput("Mk", # n
              "Male killing rate", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0, step = 0.05
            ), # starting value
            numericInput("CI", # n
              "Cytoplasmic incompatibility mortality rate", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0, step = 0.05
            ),
            numericInput("cost", # n
              "Cost of infection as proportion of fecundity reduction", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0, step = 0.05
            )
          ),
          numericInput("gen",
            "Number of generations to run simulation:",
            min = 0, # miniumum of slider input
            max = 100000, # max of slider input
            value = 20000, step = 1000
          ),
          actionButton("go", "Run Simulation")
        ),
        mainPanel( # usually where the plot or main thing is
          plotOutput("sim", height = 800) # where
        )
      )
    ),


    # * 3.g Simulations stage based -------------------------------------------


    tabPanel(
      "Simulations: population structure", # name of third tab
      titlePanel("Simulate infection dyanmics here"), # gives us the title of the webpage
      # Sidebarlayout
      sidebarLayout( # saying to use a sidebarLayout. T
        sidebarPanel(
          h3("Description:"), p("Here you can simulate the model to see how infection may spread in a marine invertebrate species under different parameter combinations with a size-based population structure (press 'Run Simulation' to make graphs). (A) the density of different sub-populations, (B) the proportion of males at different size classes, and (C) the percentage of infected individuals at different size classes. The simulations look at three potential mechanisms: enhanced growth rate in infected females, cytoplasmic incompatibility, and male-killing. First, use the checkboxes to display what parameters you want to change, then change parameters to your choosing, and finally click run simulation. Simulations can take a few minutes to run on this website. We recommend downloading the R code from dryad to run the model when the number of generations exceeds 20000. Initial values and graphs generated are for Heliocidaris erythrogramma with enhanced growth, male killing, and cytoplasmic incompatibility and are in the text: (A) Figure A7 right panel, (B) Figure 3 top right panel, (C) Figure 3 bottom right panel."), br(), # usually where you put controls that user interacts with
          checkboxGroupInput("parmsS",
            h3("What parameters do you want to change"),
            choices = list("Starting densities" = "densS", "Fertilization" = "fertsS", "Population" = "popsS", "Infection" = "infS"), inline = TRUE
          ),
          conditionalPanel(
            condition = "input.parmsS.includes('densS')",
            numericInput("Nif0S", # n
              "Starting density of infected females (individuals per m^2)",
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.01, step = 0.05
            ), # starting value
            numericInput("Nuf0S",
              "Starting density of uninfected females (individuals per m^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.99, step = 0.05
            ),
            numericInput("Nim0S",
              "Starting density of infected males (individuals per m^2)",
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.01, step = 0.05
            ),
            numericInput("Num0S", # n
              "Starting density of uninfected males (individuals per m^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.99, step = 0.05
            )
          ), # starting value
          conditionalPanel(
            condition = "input.parmsS.includes('fertsS')",
            numericInput("SigmaS",
              "Egg size (mm^2)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 2, # max of slider input
              value = 0.203, step = 0.05
            ),
            numericInput("vS",
              "Sperm velocity (mm per s)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 2, # max of slider input
              value = 0.125, step = 0.05
            ),
            # Fe=0.05,tb=1,S=700,E=0.05,tau=5400,generations=5000
            numericInput("FeS",
              "Fertilization efficiency constant", # what you are displaying
              min = 0, # miniumum of slider input
              max = 0.5, # max of slider input
              value = 0.09444, step = 0.01
            ),
            numericInput("tbS",
              "Time to block polysmery (s)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 800, # max of slider input
              value = 1, step = 10
            ),
            numericInput("SS", # n
              "Sperm/uL/male/m^2", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1000, # max of slider input
              value = 700, step = 20
            ), # starting value
            numericInput("ES",
              "Egg/uL/female/m^2", # what you are displaying
              min = 0, # miniumum of slider input
              max = 100, # max of slider input
              value = 0.05, step = 1
            ),
            numericInput("tauS",
              "Sperm half life (s)", # what you are displaying
              min = 0, # miniumum of slider input
              max = 10000, # max of slider input
              value = 5400, step = 100
            )
          ),
          conditionalPanel(
            condition = "input.parmsS.includes('popsS')",
            numericInput("ksS",
              "Carrying capacity (individuals per m^2)", # what you are displaying
              min = 10, # miniumum of slider input
              max = 200, # max of slider input
              value = 100, step = 10
            ),
            numericInput("dsS",
              "Shape parameter for adult mortality", # what you are displaying
              min = 0, # miniumum of slider input
              max = 10, # max of slider input
              value = 0.1, step = 0.05
            ),
            numericInput("MaS", # n
              "Max adult mortality",
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.99, step = 0.05
            ), # starting value
            numericInput("MlS",
              "Larval mortality", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.66, step = 0.05
            ),
            numericInput("psiS",
                         "Settlement constant (uL per m^2)", # what you are displaying
                         min = 0, # miniumum of slider input
                         max = 10, # max of slider input
                         value = 1, step = 0.05
            ),
            numericInput("ded1S",
              "Deduction in fecundity for smallest class ", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.237, step = 0.05
            ),
            numericInput("ded2S",
              "Deduction in fecundity for middle class", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.746, step = 0.05
            )
          ),
          conditionalPanel(
            condition = "input.parmsS.includes('infS')",
            numericInput("GIS",
              "Growth rate of infected females", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.25, step = 0.05
            ),
            numericInput("GuS",
              "Growth rate of other individuals", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.1, step = 0.05
            ),
            numericInput("MkS", # n
              "Male killing rate", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.7, step = 0.05
            ), # starting value
            numericInput("CIS", # n
              "Cytoplasmic incompatibility mortality rate", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0.7, step = 0.05
            ),
            numericInput("costS", # n
              "Cost of infection as proportion of fecundity reduction", # what you are displaying
              min = 0, # miniumum of slider input
              max = 1, # max of slider input
              value = 0, step = 0.05
            )
          ),
          numericInput("genS",
            "Number of generations to run simulation:",
            min = 0, # miniumum of slider input
            max = 100000, # max of slider input
            value = 20000, step = 1000
          ),
          actionButton("goS", "Run Simulation")
        ),
        mainPanel( # usually where the plot or main thing is
          plotOutput("simS", height = 800) # where
        )
      )
    ),

    # * 3.h Marine Invert species ---------------------------------------------


    tabPanel(
      "Marine invertebrates (Fig. 4)", # name of third tab
      titlePanel("Species plotted on Figure 4"),
      sidebarLayout( # saying to use a sidebarLayout. Two parts: sidebarPanel and mainPanel
        sidebarPanel(
          h3("Description:"), conditionalPanel(condition = "input.data=='fem'", p("Pelagic larval duration and egg size influence infection spread. Here you can explore where known marine invertebrates fall on these two life-history axes (Figure 4 left). Hover over data points to get more information on those species. Shown is a heatmap of how long it took for the model to reach equilibrium, the proportion of infected individuals at equilibrium, or population density at equilibrium. For these graphs, we assumed no cost to infection, feminization rate = 0.7, and other model parameters were held constant (Tables 1 and 2).")), conditionalPanel(condition = "input.data=='mk'", p("Pelagic larval duration and egg size influence infection spread. Here you can explore where known marine invertebrates fall on these two life-history axes (Figure 4 right). Hover over data points to get more information on those species. Shown is a heatmap of either how long it took for the model to reach equilibrium, the proportion of infected individuals at equilibrium, or the population density at equilibrium. For these graphs, we assumed no cost to infection, male killing = 0.7 and cost = -1.25 representing a benefit, and other parameters were held constant (Tables 1 and 2).")),
          br(), # usually where you put controls that user interacts with
          radioButtons("data", "Which simulation results?", choices = c("Feminization" = "fem", "Male Killing + Enhanced Growth" = "mk"), selected = "fem"),
          radioButtons(
            "graph", # name of input
            "What to show?", # the header
            choices = list( # where you put the choices
              "Generations to equilibrium" = "gens", # name to show = variable to use in server
              "Infection rate at equilibrium" = "inf", # name to show  = variable to use in server
              "Density at equilibrium" = "dens"
            ),
            selected = "inf" # this is where you put the default
          )
        ),
        mainPanel( # usually where the plot or main thing is
          plotlyOutput("speciesf", height = 800) # where you specify the name of the object you are plotting
        )
      )
    )
  )
)

# 4. Server Side ----------------------------------------------------------
# Define server logic required to draw a histogram
server <- function(input, output) { # server is where the actual plotting/data crunching happens

  # * 4.b Fertilization sex ratio ------------------------------------------
  output$fertf <- renderPlot({
    a <- plot_functionf(sigma = input$Sigmaf, v = input$vf, Fe = input$Fef, tb = input$tbf, S = input$Sf, E = input$Ef, tau = input$tauf) + ggtitle("")
    b <- plot_function2f(sigma = input$Sigmaf, v = input$vf, Fe = input$Fef, tb = input$tbf, S = input$Sf, E = input$Ef, tau = input$tauf) + ggtitle("")
    c <- a + b + plot_layout(guides = "collect")
    return(c)
  })
  # * 4.c Juv mortality ---------------------------------------------------
  output$al_func <- renderPlot({ # use renderPlot to actually make the plot
    ggplot(data = data.frame(t = 0), mapping = aes(t = t)) +
      xlab("Pelagic larval duration (days)") +
      ylab("Mortality") +
      theme(text = element_text(size = 15)) +
      labs(color = "Functions") +
      stat_function(fun = deathl, args = (list(Mz = input$Mzl))) +
      xlim(0, 200)
  })
  # * 4.c Adult mortality ---------------------------------------------------
  output$am_func <- renderPlot({ # use renderPlot to actually make the plot
    ggplot(data = data.frame(Nt = 0), mapping = aes(Nt = Nt)) +
      xlab(bquote(bold('Density (individuals per m'^2*")"))) +
      ylab("Probability of Survival") +
      theme(text = element_text(size = 15)) +
      labs(color = "Functions") +
      stat_function(fun = death, args = (list(Ma = input$Ma, d = input$d, k = input$k))) +
      xlim(0, input$k + 10)
  })

  # * 4.d Simulation part ---------------------------------------------------
  sim_reactive <- eventReactive(input$go, {
    states_t <- c(fi = input$Nif0, fu = input$Nuf0, mi = input$Nim0, mu = input$Num0)
    
    # *3.a.1 HE no CI, no MK, Fem ----------------------------------------------------------
    # paramaters
    parmst <- c(
      sigma = input$Sigma, v = input$v, Fe = input$Fe, tb = input$tb, S = input$S, E = input$E, tau = input$tau, c = input$cost, # fertilization parms
      K = input$ks, d = input$ds, Ma = input$Ma, Ml = input$Ml, CI = input$CI, Mk = input$Mk, r = input$R,psi=input$psi
    )
    maxtime <- input$gen
    Data_1 <- ode(y = states_t, times = 0:maxtime, func = pop_next, parms = parmst, method = "iteration") %>%
      data.frame() %>%
      pivot_longer(c(fi, fu, mi, mu), names_to = "Population", values_to = "Density")
  })

  output$sim <- renderPlot({
    if (input$go == 0) {
      Data_1_Sr <- Data_1 %>%
        group_by(time) %>%
        summarise(SexRatio = (Density[Population == "mu"] + Density[Population == "mi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population == "mu"] + Density[Population == "mi"]))
      Data_1_I <- Data_1 %>%
        group_by(time) %>%
        summarise(InfRatio = (Density[Population == "mi"] + Density[Population == "fi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population == "mu"] + Density[Population == "mi"]))
      p1 <- ggplot(Data_1, aes(x = time, y = Density, color = Population)) +
        geom_line() +
        theme(legend.position = "top") +
        ylab(bquote(bold('Density (individuals per m'^2*")"))) +
        scale_x_continuous(label = comma) +
        scale_color_manual(labels = c("Infected Females", "Infected Males", "Uninfected Females", "Uninfected Males"), breaks = c("fi", "mi", "fu", "mu"), values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      sr1 <- ggplot(Data_1_Sr, aes(x = time, y = SexRatio)) +
        geom_line(color = "#6D5B97") +
        theme(legend.position = "right") +
        ylab("Sex ratio") +
        scale_x_continuous(label = comma)

      I1 <- ggplot(Data_1_I, aes(x = time, y = InfRatio)) +
        geom_line(color = "#6D5B97") +
        theme(legend.position = "right") +
        ylab("Infected") +
        scale_x_continuous(label = comma) +
        scale_y_continuous(label = percent, limits = c(0, 1))
      ((p1 / sr1 / I1) & xlab("Time")) + plot_annotation(tag_levels = "A")
    } else {
      Data_1 <- sim_reactive()
      Data_1_Sr <- Data_1 %>%
        group_by(time) %>%
        summarise(SexRatio = (Density[Population == "mu"] + Density[Population == "mi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population == "mu"] + Density[Population == "mi"]))
      Data_1_I <- Data_1 %>%
        group_by(time) %>%
        summarise(InfRatio = (Density[Population == "mi"] + Density[Population == "fi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population == "mu"] + Density[Population == "mi"]))
      p1 <- ggplot(Data_1, aes(x = time, y = Density, color = Population)) +
        geom_line() +
        theme(legend.position = "top") +
        ylab(bquote(bold('Density (individuals per m'^2*")"))) +
        scale_x_continuous(label = comma) +
        scale_color_manual(labels = c("Infected Females", "Infected Males", "Uninfected Females", "Uninfected Males"), breaks = c("fi", "mi", "fu", "mu"), values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      sr1 <- ggplot(Data_1_Sr, aes(x = time, y = SexRatio)) +
        geom_line(color = "#6D5B97") +
        theme(legend.position = "right") +
        ylab("Sex Ratio") +
        scale_x_continuous(label = comma)
      I1 <- ggplot(Data_1_I, aes(x = time, y = InfRatio)) +
        geom_line(color = "#6D5B97") +
        theme(legend.position = "right") +
        ylab("Infected") +
        scale_y_continuous(label = percent, limits = c(0, 1)) +
        scale_x_continuous(label = comma)
      ((p1 / sr1 / I1) & xlab("Time")) + plot_annotation(tag_levels = "A")
    }
  }) %>% bindCache(input$go, sim_reactive())

  # * 4.e Stage based simulation --------------------------------------------
  sim_reactiveS <- eventReactive(input$goS, {
    states_tS <- c(f1i = input$Nif0S, f1u = input$Nuf0S, m1i = input$Nim0S, m1u = input$Num0S, f2i = 0, f2u = 0, m2i = 0, m2u = 0, f3i = 0, f3u = 0, m3i = 0, m3u = 0)
    # *3.a.1 HE no CI, no MK, Fem ----------------------------------------------------------
    # paramaters
    parmstS <- c(
      sigma = input$SigmaS, v = input$vS, Fe = input$FeS, tb = input$tbS, S = input$SS, E = input$ES, tau = input$tauS, c = input$costS, # fertilization parms
      K = input$ksS, d = input$dsS, Ma = input$MaS, Ml = input$MlS, CI = input$CIS, Mk = input$MkS, t = input$GuS, tif = input$GIS, ded1 = input$ded1S, ded2 = input$ded2S,psi=input$psiS
    )
    maxtimeS <- input$genS
    Data_1S <- ode(y = states_tS, times = 0:maxtimeS, func = pop_nextS, parms = parmstS, method = "iteration") %>%
      data.frame() %>%
      pivot_longer(c(f1i, f1u, m1i, m1u, f2i, f2u, m2i, m2u, f3i, f3u, m3i, m3u), names_to = "Population", values_to = "Density")
  })

  output$simS <- renderPlot({
    if (input$goS == 0) {
      Data_1_SrS <- Data_1S %>%
        group_by(time) %>%
        summarise(sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) / (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population == "m1u"]), sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) / (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population == "m2u"]), sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) / (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population == "m3u"])) %>%
        pivot_longer(c(sexratio1, sexratio2, sexratio3), names_to = "Stage", values_to = "Sex Ratio")
      Data_1_IS <- Data_1S %>%
        group_by(time) %>%
        summarise(I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population == "m1u"]), I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population == "m2u"]), I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population == "m3u"])) %>%
        pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")
      p1S <- Data_1S %>%
        separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
        ggplot(aes(x = time, y = Density, color = Class, linetype = Infection)) +
        geom_line() +
        # change labels of our axis
        labs(y = bquote(bold('Density (individuals per m'^2*")")), x = "Generation", color = "Subpopulations") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        scale_linetype(labels = c("Infected", "Uninfected")) +
        scale_color_manual(labels = c("Female 20 - 39 mm", "Female 40 - 59 mm", "Female 60+ mm", "Male 20 - 39 mm", "Male 40 - 59 mm", "Male 60+ mm"), values = c("#564157", "#796FBD", "#13D1E2", "#FF0015", "#B2D75E", "#ACD9C3")) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      sr1S <- ggplot(Data_1_SrS, aes(x = time, y = `Sex Ratio`, color = Stage)) +
        geom_line() +
        # change labels of our axis
        labs(y = "Sex Ratio", x = "Generation") +
        scale_color_manual(labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"), values = c("#C783AA", "#666A9E", "#2F4E2F"), name = "Size Class") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        ylim(0, 1) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      I1S <- ggplot(Data_1_IS, aes(x = time, y = `Infection Ratio`, color = Stage)) +
        geom_line() +
        # change labels of our axis
        labs(y = "Infection %", x = "Generation") +
        scale_color_manual(labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"), values = c("#C783AA", "#666A9E", "#2F4E2F"), name = "Size Class") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        scale_y_continuous(label = percent) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      ((p1S / sr1S / I1S) & xlab("Time")) + plot_annotation(tag_levels = "A")
    } else {
      Data_1S <- sim_reactiveS()
      Data_1_SrS <- Data_1S %>%
        group_by(time) %>%
        summarise(sexratio1 = (Density[Population == "m1i"] + Density[Population == "m1u"]) / (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population == "m1u"]), sexratio2 = (Density[Population == "m2i"] + Density[Population == "m2u"]) / (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population == "m2u"]), sexratio3 = (Density[Population == "m3i"] + Density[Population == "m3u"]) / (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population == "m3u"])) %>%
        pivot_longer(c(sexratio1, sexratio2, sexratio3), names_to = "Stage", values_to = "Sex Ratio")
      Data_1_IS <- Data_1S %>%
        group_by(time) %>%
        summarise(I1 = (Density[Population == "m1i"] + Density[Population == "f1i"]) / (Density[Population == "f1i"] + Density[Population == "f1u"] + Density[Population == "m1i"] + Density[Population == "m1u"]), I2 = (Density[Population == "m2i"] + Density[Population == "f2i"]) / (Density[Population == "f2i"] + Density[Population == "f2u"] + Density[Population == "m2i"] + Density[Population == "m2u"]), I3 = (Density[Population == "m3i"] + Density[Population == "f3i"]) / (Density[Population == "f3i"] + Density[Population == "f3u"] + Density[Population == "m3i"] + Density[Population == "m3u"])) %>%
        pivot_longer(c(I1, I2, I3), names_to = "Stage", values_to = "Infection Ratio")
      p1S <- Data_1S %>%
        separate(Population, c("Class", "Infection"), "(?<=\\G..)") %>%
        ggplot(aes(x = time, y = Density, color = Class, linetype = Infection)) +
        geom_line() +
        # change labels of our axis
        labs(y = bquote(bold('Density (individuals per m'^2*")")), x = "Generation", color = "Subpopulations") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        scale_linetype(labels = c("Infected", "Uninfected")) +
        scale_color_manual(labels = c("Female 20 - 39 mm", "Female 40 - 59 mm", "Female 60+ mm", "Male 20 - 39 mm", "Male 40 - 59 mm", "Male 60+ mm"), values = c("#564157", "#796FBD", "#13D1E2", "#FF0015", "#B2D75E", "#ACD9C3")) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      sr1S <- ggplot(Data_1_SrS, aes(x = time, y = `Sex Ratio`, color = Stage)) +
        geom_line() +
        # change labels of our axis
        labs(y = "Sex Ratio", x = "Generation") +
        scale_color_manual(labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"), values = c("#C783AA", "#666A9E", "#2F4E2F"), name = "Size Class") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        ylim(0, 1) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      I1S <- ggplot(Data_1_IS, aes(x = time, y = `Infection Ratio`, color = Stage)) +
        geom_line() +
        # change labels of our axis
        labs(y = "Infection %", x = "Generation") +
        scale_color_manual(labels = c("20 - 39 mm", "40 - 59 mm", "60+ mm"), values = c("#C783AA", "#666A9E", "#2F4E2F"), name = "Size Class") +
        theme(legend.position = "right") +
        scale_x_continuous(label = comma) +
        scale_y_continuous(label = percent) +
        guides(color = guide_legend(override.aes = list(size = 2)))
      ((p1S / sr1S / I1S) & xlab("Time")) + plot_annotation(tag_levels = "A")
    }
  }) %>% bindCache(input$goS, sim_reactiveS())
  # * 4. e Species part -----------------------------------------------------
  output$speciesf <- renderPlotly({
    if (input$data == "mk") {
      if (input$graph == "gens") {
        ggplotly(ggplot(Data2tB[Data2tB$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Time)) +
          scale_fill_viridis(option = "magma", direction = -1, labels = comma, na.value = "darkgrey") +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Years to \n equilibrium") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      } else if (input$graph == "inf") {
        ggplotly(ggplot(Data2tB[Data2tB$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Infection)) +
          scale_fill_distiller(palette = "Blues", na.value = "grey", direction = 1, limits = c(0, 1)) +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Proportion \n   infected") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg Diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      } else {
        ggplotly(ggplot(Data2tB[Data2tB$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Density)) +
          scale_fill_distiller(palette = "Greens", na.value = "grey", direction = 1, breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), limits = c(0, 30)) +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Density") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg Diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      }
    } else {
      if (input$graph == "gens") {
        ggplotly(ggplot(Data2t[Data2t$Feminization == 0.7 & Data2t$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Time)) +
          scale_fill_viridis(option = "magma", direction = -1, labels = comma, na.value = "darkgrey") +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Years to \n equilibrium") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      } else if (input$graph == "inf") {
        ggplotly(ggplot(Data2t[Data2t$Feminization == 0.7 & Data2t$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Infection)) +
          scale_fill_distiller(palette = "Blues", na.value = "grey", direction = 1, limits = c(0, 1)) +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Proportion \n   infected") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      } else {
        ggplotly(ggplot(Data2t[Data2t$Feminization == 0.7 & Data2t$Larval_d > 1, ], aes(x = Egg_Size, y = Larval_d)) +
          geom_tile(aes(fill = Density)) +
          scale_fill_distiller(palette = "Greens", na.value = "grey", direction = 1, breaks = c(0, 10, 20, 30), labels = c(0, 10, 20, 30), limits = c(0, 30)) +
          labs(x = "", y = "Pelagic larval duration (days)", fill = "Density (individuals per m^2)") +
          theme(legend.position = "right") +
          geom_point(data = species, aes(x = EggSize2, y = Days, shape = DevelopmentalMode, text = labls), size = 1.5, fill = "grey", color = "black") +
          xlab("Egg diameter (um)") +
          scale_x_continuous(expand = c(0.0, 0.0)) +
          scale_y_continuous(expand = c(0.0, 0.0)) +
          scale_shape_manual(name = "", values = c(22, 21), labels = c("Lecithotrophic", "Planktotrophic")) +
          guides(shape = guide_legend(override.aes = list(size = 3))) +
          theme(axis.line = element_line(size = 1)), tooltip = "text") %>%
          add_annotations(
            text = "Developmental \n mode", xref = "paper", yref = "paper",
            x = 1.00, xanchor = "left",
            y = 0.4, yanchor = "bottom", # Same y as legend below
            legendtitle = TRUE, showarrow = FALSE
          ) %>%
          layout(legend = list(y = 0.4, yanchor = "top"))
      }
    }
  }) %>% bindCache(input$graph, input$data)
}




# 5. Runing the shiny app -------------------------------------------------


# Run the application
shinyApp(ui = ui, server = server)

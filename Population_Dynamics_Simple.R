#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script file that runs and plots the simple population model for 
#"On the spread of microbes that manipulate reproduction in marine invertebrates" #Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282. 
#Currently coded to simulate and plot Figure 2 and Figure A6
#but users can run/plot the model under different parameter combinations. 
# 1. Setting up -----------------------------------------------------------

# * 1.a Load up libraries -------------------------------------------------
library(viridis)
library(patchwork)
library(scales)
library(tidyverse)
library(deSolve)
library(svglite)

# * 1.b Setting my default plotting theme -----------------------
mytheme <-
  theme_linedraw() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
    axis.title = (element_text(face = "bold")),
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(
      size = 12, face = "bold", color =
        "black"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 2
    )
  ) # makes all the text larger and bold
theme_set(mytheme)


# 2. Functions ------------------------------------------------------------
# function to calculate number of zygotes
#slightly differs from "fertilization.R" code because for population models we need to know how many infected female there are.
#Males is male adult density (Males per m^2)
#sigma is egg size (cross sectional area of the egg; mm^2)
#v is sperm speed (mm/s)
#Fe is fertilization efficiency
#tb is time for polyspermy block (s)
#sr is sex ratio (proportion of population that is male)
#S is sperm density per unit male density (sperm/uL/individual/m^2)
#E is egg density per unit female density (egg/uL/individual/m^2)
#tau is half life of sperm (s)
#c is cost due to infection.
#nif is density of infected females (females per m^2)
number_zygotes <- function(Males, sigma, v, Fe, tb, sr, S, E, tau, c, nif) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)
  S0 <- Males * S
  Females <- (Males / sr) - Males
  E0 <- E * (Females - nif) + E * (1 - c) * nif
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizing sperm
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b)))
  return(prop_mono * E0)
}
# birth -> death -> transition
# function to run simulations
#needs to be in this format to work with deSolve package
pop_next <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #Calculate total adult density
    Nta <- fi + fu + mi + mu
    #calculate sex ratio (proportion males)
    sr <- (mi + mu) / Nta
    #calculate number of zygotes produced
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
      ) * psi #settlement coefficient
    #slopes (S and E) are not included since they cancel out
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
    #density of infected male zygotes (eq 11)
    Nzim <- (Nz * (ei / et) * (1 - r)) * (1 - Ml) * (1 - Mk)
    #density of uninfected male zygotes (eq 9)
    Nzum <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -Ml)
    #density of infected female zygotes (eq 10)
    Nzif <- (Nz * (ei / et) * r) * (1 - Ml)
    #density of uninfected female zygotes (eq 9)
    Nzuf <- ((Nz * (eu / et) * 0.5) * (si * (1 - CI) + su) / st) * (1 -Ml)
    #total zygote density
    Nzt <- Nzim + Nzum + Nzif + Nzuf
    # total adult density
    # death
    # mortality (eq 7)
    mort <- (1 - Ma / (1 + exp(-d * (Nta + Nzt - K / 2))))
    #calculate the next generation
    fi_plus1 <- (Nzif + fi) * mort
    fu_plus1 <- (Nzuf + fu) * mort
    mi_plus1 <- (Nzim + mi) * mort
    mu_plus1 <- (Nzum + mu) * mort
    list(c(fi_plus1, fu_plus1, mi_plus1, mu_plus1))
  })
}

# 3. Running simulations --------------------------------------------------


# starting states
# fi = infected female density
# fu = uninfected female density
# mi = infected male density
# mu= uninfected male density
states_t <- c(
  fi = 0.01,
  fu = 0.99,
  mi = 0.01,
  mu = 0.99
)
# HE= Heliocidaris erythrogramma,
# *3.a.1 HE no CI, no MK, Fem ----------------------------------------------------------
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
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HEf <-
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
  ) %>%
  mutate(Model = "Feminization")
# *3.a.2 HE  CI, no MK, No Fem ----------------------------------------------------------
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
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0.7,
    Mk = 0,
    r = 0.5,
    psi=1
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HECI <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HeparmsCI,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  ) %>%
  mutate(Model = "Cytoplasmic Incompatability")

# *3.a.1 HE  No CI, MK, No Fem ----------------------------------------------------------
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
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.66,
    CI = 0,
    Mk = 0.7,
    r = 0.5,
    psi=1
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HEMK <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HeparmsMk,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  ) %>%
  mutate(Model = "Male Killing")
# *3.b HE putting it together ----------------------------------------------------------
#merging data frames
HE_all <- rbind(HEf, HECI, HEMK)
#calculating sex ratio
Data_HE_Naive_Sr <- HE_all %>%
  group_by(Model, time) %>%
  summarise(SexRatio = (Density[Population == "mu"] + Density[Population ==
    "mi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population ==
    "mu"] + Density[Population == "mi"]))
#calculated infection rate
Data_HE_Naive_I <- HE_all %>%
  group_by(Model, time) %>%
  summarise(InfRatio = (Density[Population == "mi"] + Density[Population ==
    "fi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population ==
    "mu"] + Density[Population == "mi"]))

# HT is Heliocidaris tuberculata
# 4.a.1 HT no CI, no MK, Fem ----------------------------------------------------------
# paramaters

Htparmsf <-
  c(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400,
    c = 0,
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0,
    Mk = 0,
    r = 0.7,
    psi=1
  ) # pop dyanmics equation
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HTf <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = Htparmsf,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  ) %>%
  mutate(Model = "Feminization")
# 4.a.2 HT  CI, no MK, No Fem ----------------------------------------------------------
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
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0.7,
    Mk = 0,
    r = 0.5,
    psi=1
  ) # pop dyanmics equation
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
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  ) %>%
  mutate(Model = "Cytoplasmic Incompatability")

# 4.b.1 HE  No CI, MK, No Fem ----------------------------------------------------------
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
    # fertilization parms
    K = 100,
    d = 0.1,
    Ma = 0.99,
    Ml = 0.998863725,
    CI = 0,
    Mk = 0.7,
    r = 0.5,
    psi=1
  )
# this gives a matrix of ODE solver output...Used method iteration for discrete time model
HTMK <-
  ode(
    y = states_t,
    times = 0:30000,
    func = pop_next,
    parms = HtparmsMk,
    method = "iteration"
  ) %>%
  data.frame() %>%
  pivot_longer(c(fi, fu, mi, mu),
    names_to = "Population", values_to =
      "Density"
  ) %>%
  mutate(Model = "Male Killing")
# *4.b HT putting it together ----------------------------------------------------------
HT_all <- rbind(HTf, HTCI, HTMK)
#calculating sex ratio
Data_HT_Naive_Sr <- HT_all %>%
  group_by(Model, time) %>%
  summarise(SexRatio = (Density[Population == "mu"] + Density[Population ==
    "mi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population ==
    "mu"] + Density[Population == "mi"]))
#calculating infection rate
Data_HT_Naive_I <- HT_all %>%
  group_by(Model, time) %>%
  summarise(InfRatio = (Density[Population == "mi"] + Density[Population ==
    "fi"]) / (Density[Population == "fi"] + Density[Population == "fu"] + Density[Population ==
    "mu"] + Density[Population == "mi"]))

# 5 putting it all together to make Figure 2 ----------------------------------------------
# *5.a putting data together ----------------------------------------------
#adding in a species column so we can later facet
Data_HE_Naive_Sr$Species <- "H. erythrogramma"
Data_HT_Naive_Sr$Species <- "H. tuberculata"
Data_HE_Naive_I$Species <- "H. erythrogramma"
Data_HT_Naive_I$Species <- "H. tuberculata"

Sr_all <- rbind(Data_HE_Naive_Sr, Data_HT_Naive_Sr)
I_all <- rbind(Data_HE_Naive_I, Data_HT_Naive_I)

# * 5.b Ploting -----------------------------------------------------------
# Plot of sex ratio
All_Naive_sr <-
  ggplot(Sr_all, aes(x = time, y = SexRatio, color = Species)) +
  geom_line(
    size =
      1
  ) +
  facet_grid(. ~ Model) +
  ylab("Sex Ratio") +
  xlab("") +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  ylim(0.2, 0.8) +
  xlab("Time") +
  scale_color_manual(values = c("#6D5B97", "#A8332A"))
All_Naive_sr
# Plot of infection rate
All_Naive_I <-
  ggplot(I_all, aes(x = time, y = InfRatio, color = Species)) +
  geom_line(
    size =
      1
  ) +
  theme(legend.position = "right") +
  facet_grid(. ~ Model) +
  ylab("Infection %") +
  xlab("Time") +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  scale_x_continuous(label = comma) +
  scale_y_continuous(label = percent, limits = c(0, 1)) +
  scale_color_manual(
    values =
      c("#6D5B97", "#A8332A")
  )
All_Naive_I
# Putting it all together
Figure2all <-
  ((All_Naive_sr / All_Naive_I) &
    theme(legend.position = "right")) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A")
Figure2all

# ggsave("Figure_2_all.svg",Figure2all,height = 150,width=200,units="mm",path="Plots")




# 6. Density plots (SI Figure 6) --------------------------------------------------------

# * 6a putting data together  ----------------------------------------------
HE_all$Species <- "H. erythrogramma"
HT_all$Species <- "H. tuberculata"


# * 6b Plotting -----------------------------------------------------------
HEP <- ggplot(HE_all, aes(color = Population, x = time, y = Density)) +
  geom_line() +
  facet_grid(. ~ Model) +
  scale_x_continuous(label = comma) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  xlab("Time") +
  scale_color_manual(
    labels = c(
      "Infected Females",
      "Infected Males",
      "Uninfected Females",
      "Uninfected Males"
    ),
    breaks = c("fi", "mi", "fu", "mu"),
    values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))

HTP <- ggplot(HT_all, aes(color = Population, x = time, y = Density)) +
  geom_line() +
  facet_grid(. ~ Model) +
  scale_x_continuous(
    label = comma, limits =
      c(0, 10000)
  ) +
  theme(
    legend.position = "right",
    aspect.ratio = 1,
    strip.background = element_blank(),
    strip.text = element_text(color = "black", size = 9, face = "italic")
  ) +
  xlab("Time") +
  scale_color_manual(
    labels = c(
      "Infected Females",
      "Infected Males",
      "Uninfected Females",
      "Uninfected Males"
    ),
    breaks = c("fi", "mi", "fu", "mu"),
    values = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
  ) +
  guides(color = guide_legend(override.aes = list(size = 2)))

SIFigure6 <-
  ((HTP / HEP) &
    theme(legend.position = "right")) + plot_layout(guides = "collect")
SIFigure6

# ggsave("SIFigure6.pdf",SIFigure3,height = 150,width=200,units="mm",path="Plots")

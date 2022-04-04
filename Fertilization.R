#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script file that runs and plots the fertilization dynamics for "On the spread of microbes that manipulate reproduction in marine invertebrates" Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282
#Runs analyses and makes plots for Figure 1, Figure A4, and Figure A5.
#Makes/uses “Fertilization_opt_f.csv” and “Fertilization_opt_z.csv”.
# 1. Setting up -----------------------------------------------------------
# * 1.a Loading up libraries ----------------------------------------------
library(viridis)
library(patchwork)
library(scales)
library(svglite)
library(foreach)
library(doParallel)
library(tidyverse)
library(ggrepel)
library(optimx)
# * 1.b Setting my default plotting theme ----------------------
mytheme <-
  theme_linedraw() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
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
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold.italic")
  ) # makes all the text larger and bold
theme_set(mytheme)

# 2. Functions of fertilization dynamics ----------------------------------
# * 2.a Probabilty of successful fertilization  ---------------------------
#Dens is total adult density (individuals per m^2)
#sigma is egg size (cross sectional area of the egg; mm^2)
#v is sperm speed (mm/s)
#Fe is fertilization efficiency
#tb is time for polyspermy block (s)
#sr is sex ratio (proportion of population that is male)
#S is sperm density per unit male density (sperm/uL/individual/m^2)
#E is egg density per unit female density (egg/uL/individual/m^2)
#tau is half life of sperm (s)
prob_mono <- function(Dens, sigma, v, Fe, tb, sr, S, E, tau) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v) (eq 2 in text)
  Males <- Dens * sr
  S0 <- Males * S
  Females <- Dens * (1 - sr)
  E0 <- E * Females
  beta0 <- sigma * v 
  tau <- tau
  # x=Average number of potential fertilizaing sperm (eq 1 in text)
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted (eq 3 in text)
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b))) #(eq 4 in text)
  return(prop_mono)
}


# * 2.b Number of zygotes produced ----------------------------------------
#paramaters are the same as above we just multiply result by egg density to get zygote density.
number_zygotes <- function(Dens, sigma, v, Fe, tb, sr, S, E, tau) {
  # beta0=collision rate constant; estimated by egg size (sigma) and velocity(v)(eq 2 in text)
  Males <- Dens * sr
  S0 <- Males * S
  Females <- Dens * (1 - sr)
  E0 <- E * Females
  beta0 <- sigma * v
  tau <- tau
  # x=Average number of potential fertilizaing sperm (eq 1 in text)
  x <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tau)))
  # b=mean number of extr fertilizing sperm that will contact an egg in a time period tb in a population of eggs already contacted (eq 3 in text)
  b <- Fe * (S0 / E0) * (1 - exp(-1 * (beta0 * E0 * tb)))
  prop_mono <- 1 - exp(-1 * x) - ((1 - exp(-1 * x) - x * exp(-1 * x)) *
    (1 - exp(-1 * b))) #(eq 4 in text)
  return(prop_mono * E0)
}


# 3. Plotting functions ---------------------------------------------------
#functions using the equations above to make plots of fertilization dynamics
# * 3.a Plotting function for % Fertilization success ---------------------
plot_function <- function(sigma, v, Fe, tb, S, E, tau) {
  Data <- data.frame(Dens = 0)
  p <-
    ggplot(data = Data, mapping = aes(Dens = Dens)) +
    xlab("Density") +
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
    fun = prob_mono,
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
    fun = prob_mono,
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
    fun = prob_mono,
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
    fun = prob_mono,
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
    fun = prob_mono,
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
  )
}

# * 3.b Plotting function for Zygote density ------------------------------
plot_function2 <- function(sigma, v, Fe, tb, S, E, tau) {
  Data <- data.frame(Dens = 0)
  p <-
    ggplot(data = Data, mapping = aes(Dens = Dens)) +
    xlab("Density") +
    ylab("Zygote Density") +
    labs(color = "Functions") +
    scale_x_continuous(
      trans = "log10",
      limits = c(0.00001, 100),
      breaks = c(0.0001, 0.01, 1, 100),
      labels = c("0.0001", "0.01", "1", "100")
    )
  p + stat_function(
    fun = number_zygotes,
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
    fun = number_zygotes,
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
    fun = number_zygotes,
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
    fun = number_zygotes,
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
    fun = number_zygotes,
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
  )
}

# 4. Making species specific plots ----------------------------------------

# * 4.a Heliocidaris erythrogramma plots ----------------------------------
#Fertiliztion percentage
ae <-
  plot_function(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400
  ) + ggtitle("") + theme(aspect.ratio = 1)
#zygote density
be <-
  plot_function2(
    sigma = 0.203,
    v = 0.125,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 0.05,
    tau = 5400
  ) + ggtitle("") + theme(aspect.ratio = 1)

# * 4.b Heliocidaris tuberculata plots ----------------------------------
#Fertilization percentage
at <-
  plot_function(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400
  ) + ggtitle("") + theme(aspect.ratio = 1)
#Zygote density
bt <-
  plot_function2(
    sigma = 0.017,
    v = 0.14,
    Fe = 0.09444,
    tb = 1,
    S = 700,
    E = 66.59171429,
    tau = 5400
  ) + ggtitle("") + theme(aspect.ratio = 1)
# * 4.c Putting the plots together to make figure 1 ----------------------------------
Figure1 <- (at + bt) / (ae + be) + plot_layout(guides = "collect")
Figure1
# 5. Saving figure 1 ----------------------------------
# ggsave("Figure_1_tdens.png",Figure1,height = 140,width=183,units="mm",path="Plots/")
# 6. Describe the fertilization pattern more broadly with a heatmap FigureA4 --------------------------------------
# make vector of sex ratios to calculate values for
srs <- seq(0.005, 1, 0.005)
# make vector of densities to calculate values for
dens <- seq(0.005, 100, 0.005)
# registering cores to use for analyssi
registerDoParallel(6)
# caluclate density of zygotes and fertilizatiion percentage across range of values for HTub
dataSRDHT <- foreach(s = srs, .combine = rbind) %dopar% {
  htod <-
    number_zygotes(
      Dens = dens,
      sigma = 0.017,
      v = 0.14,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 66.59171429,
      tau = 5400,
      sr = s
    )
  htof <-
    prob_mono(
      Dens = dens,
      sigma = 0.017,
      v = 0.14,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 66.59171429,
      tau = 5400,
      sr = s
    )
  return(
    data.frame(
      Species = "Heliocidaris tuberculata",
      Dens = dens,
      Fert = htof,
      Zygote = htod,
      sexratio = s
    )
  )
}
# caluclate density of zygotes and fertilizatiion percentage across range of values for HEry
dataSRDHE <- foreach(s = srs, .combine = rbind) %dopar% {
  heod <-
    number_zygotes(
      Dens = dens,
      sigma = 0.203,
      v = 0.125,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 0.05,
      tau = 5400,
      sr = s
    )
  heof <-
    prob_mono(
      Dens = dens,
      sigma = 0.203,
      v = 0.125,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 0.05,
      tau = 5400,
      sr = s
    )
  return(
    data.frame(
      Species = "Heliocidaris erythrogramma",
      Dens = dens,
      Fert = heof,
      Zygote = heod,
      sexratio = s
    )
  )
}
#need to use a different theme for heatmaps
mythemehp <-
  theme_classic() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
    # this makes the axis titles in bold,
    axis.line = element_line(
      color = "black", size =
        0
    ),
    # Makes the axis line black and  thicker
    text = element_text(
      size = 12, color =
        "black"
    ),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(
      color = "black", fill = "white", size =
        2
    ),
    strip.text = element_text(color = "black"),
    aspect.ratio = 1
  ) # makes all the text larger and bold

# He Fertilization percentage
FHE <-
  ggplot(dataSRDHE[dataSRDHE$sexratio > 0.1 &
    dataSRDHE$sexratio < 0.9 &
    dataSRDHE$Dens < 5, ], aes(x = sexratio, y = Dens, fill = Fert)) +
  geom_tile() +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_fill_viridis(
    option = "magma",
    na.value = "grey",
    labels = percent,
    limits = c(0, 1),
    direction = -1
  ) +
  labs(x = "Sex ratio", y = "Adult Densty", fill = "Fertilization (%)") +
  mythemehp
# He Zygote
ZHE <-
  ggplot(dataSRDHE[dataSRDHE$sexratio > 0.1 &
    dataSRDHE$sexratio < 0.9 &
    dataSRDHE$Dens < 5, ], aes(x = sexratio, y = Dens, fill = Zygote)) +
  geom_tile() +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_fill_viridis(
    option = "viridis",
    na.value = "grey",
    direction = -1
  ) +
  labs(x = "Sex ratio", y = "Adult Densty", fill = "Zygote Density") +
  mythemehp
# HT Fertilization percentage
FHT <-
  ggplot(dataSRDHT[dataSRDHT$sexratio > 0.1 &
    dataSRDHT$sexratio < 0.9, ], aes(x = sexratio, y = Dens, fill = Fert)) +
  geom_tile() +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_fill_viridis(
    option = "magma",
    na.value = "grey",
    labels = percent,
    limits = c(0, 1),
    direction = -1
  ) +
  labs(x = "Sex ratio", y = "Adult Densty", fill = "Fertilization (%)") +
  mythemehp
# HT Zygote
ZHT <-
  ggplot(dataSRDHT[dataSRDHT$sexratio > 0.1 &
    dataSRDHT$sexratio < 0.9, ], aes(x = sexratio, y = Dens, fill = Zygote)) +
  geom_tile() +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_fill_viridis(
    option = "viridis",
    na.value = "grey",
    direction = -1
  ) +
  labs(x = "Sex ratio", y = "Adult Densty", fill = "Zygote Density") +
  mythemehp

# putting it all together
FigureSI4 <-
  ((FHT + ZHT) / (FHE + ZHE) &
    theme(legend.position = "right")) + plot_layout(guides = "collect")
FigureSI4
# ggsave("FigureSI4.png",FigureSI4,height =183,width=183,units="mm",dpi=720)
#please note order of legends and subtitles were added manually.

# 7. Supplemental analysis finding maximum  possible zygote density and fertilization percentage Figure A5 --------------------------------

# * 7.a Calculating optimal sex ratio -------------------------------------
# can skip this bit since it takes a long time to run. The output is saved in the data folder so it can be read instead of reran.
# list of sex ratios to calculate
srs <- seq(0.001, 1, 0.001)
# list of densities to calculate
dens <- seq(0.0001, 100, 0.0001)
# register number of cores used for analysis
registerDoParallel(10)
# Calculate maximum zygote production at different sex ratios
maxdata <- foreach(s = srs, .combine = rbind) %dopar% {
  htod <-
    number_zygotes(
      Dens = dens,
      sigma = 0.017,
      v = 0.14,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 66.59171429,
      tau = 5400,
      sr = s
    )
  imt <- which.max(htod)
  if (length(imt) < 1) {
    tz <- NA
    td <- NA
  } else {
    tz <- htod[imt]
    td <- dens[imt]
  }
  heod <-
    number_zygotes(
      Dens = dens,
      sigma = 0.203,
      v = 0.125,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 0.05,
      tau = 5400,
      sr = s
    )
  ime <- which.max(heod)
  if (length(ime) < 1) {
    ez <- NA
    ed <- NA
  } else {
    ez <- heod[ime]
    ed <- dens[ime]
  }
  return(data.frame(
    Species = c("Heliocidaris tuberculata", "Heliocidaris erythrogramma"),
    Dens = c(td, ed),
    Zygote = c(tz, ez),
    sexratio = s
  ))
}

# Calculate maximum fertilization percentage at different sex ratios
maxdatasr <- foreach(s = srs, .combine = rbind) %dopar% {
  htod <-
    prob_mono(
      Dens = dens,
      sigma = 0.017,
      v = 0.14,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 66.59171429,
      tau = 5400,
      sr = s
    )
  imt <- which.max(htod)
  if (length(imt) < 1) {
    tz <- NA
    td <- NA
  } else {
    tz <- htod[imt]
    td <- dens[imt]
  }
  heod <-
    prob_mono(
      Dens = dens,
      sigma = 0.203,
      v = 0.125,
      Fe = 0.09444,
      tb = 1,
      S = 700,
      E = 0.05,
      tau = 5400,
      sr = s
    )
  ime <- which.max(heod)
  if (length(ime) < 1) {
    ez <- NA
    ed <- NA
  } else {
    ez <- heod[ime]
    ed <- dens[ime]
  }
  return(data.frame(
    Species = c("Heliocidaris tuberculata", "Heliocidaris erythrogramma"),
    Dens = c(td, ed),
    Fert = c(tz, ez),
    sexratio = s
  ))
}

# write.csv(maxdatasr,"Data/Fertilization_opt_f.csv")
# write.csv(maxdata,"Data/Fertilization_opt_z.csv")

# * 6.b Reading in and processing results ---------------------------------
# Maximum zyogte density production at diff sex ratio
maxdata <- read.csv("Data/Fertilization_opt_z.csv") %>%
  mutate(Fert = ifelse(
    Species == "Heliocidaris erythrogramma",
    Zygote / (0.05 * (Dens * (1 - sexratio))),
    Zygote / (66.59171429 * (Dens * (1 - sexratio)))
  ))
#factoring the species column and reordering it.
maxdata$Species <-
  factor(maxdata$Species,
    levels = c("Heliocidaris tuberculata", "Heliocidaris erythrogramma")
  )
# Maximum fertilization percentage at diff sex ratio
maxdatasr <- read.csv("Data/Fertilization_opt_f.csv") %>%
  mutate(Zygote = ifelse(
    Species == "Heliocidaris erythrogramma",
    Fert * (0.05 * (Dens * (1 - sexratio))),
    Fert * (66.59171429 * (Dens * (1 - sexratio)))
  ))
# change the orderiing of factors.
maxdatasr$Species <-
  factor(
    maxdatasr$Species,
    levels = c("Heliocidaris tuberculata", "Heliocidaris erythrogramma")
  )

# 7. Plotting -------------------------------------------------------------
# HE maximum zygote density production at diff sex ratios
zyge <-
  ggplot(maxdata[maxdata$Species == "Heliocidaris erythrogramma", ], aes(
    x =
      sexratio, y = Zygote
  )) +
  geom_line() +
  ylab("Zygote density") +
  xlab("Sex ratio") +
  theme(aspect.ratio = 1)
zyge
# HE maximum fertilizaton percentage at diff sex ratios
fertesr <-
  ggplot(maxdatasr[maxdatasr$Species == "Heliocidaris erythrogramma", ], aes(
    x =
      sexratio, y = Fert
  )) +
  geom_line() +
  ylab("Fertilization %") +
  xlab("Sex ratio") +
  theme(aspect.ratio = 1) +
  scale_y_continuous(labels = label_percent())
fertesr
# HT maximum zygote density production at diff sex ratios
zygt <-
  ggplot(maxdata[maxdata$Species == "Heliocidaris tuberculata", ], aes(
    x =
      sexratio, y = Zygote
  )) +
  geom_point(size = 0.0001) +
  ylab("Zygote density") +
  xlab("Sex ratio") +
  theme(aspect.ratio = 1)
zygt
# HT maximum fertilizaton percentage at diff sex ratios
ferttsr <-
  ggplot(maxdatasr[maxdatasr$Species == "Heliocidaris tuberculata", ], aes(
    x =
      sexratio, y = Fert
  )) +
  geom_line() +
  ylab("Fertilization %") +
  xlab("Sex ratio") +
  theme(aspect.ratio = 1) +
  scale_y_continuous(labels = label_percent())
ferttsr
# putting it all together
SIFigure5 <- (ferttsr + zygt) / (fertesr + zyge)
SIFigure5
# ggsave("SIFigure_5.pdf",SIFigure5,height =183,width=183,units="mm",path="Plots/")

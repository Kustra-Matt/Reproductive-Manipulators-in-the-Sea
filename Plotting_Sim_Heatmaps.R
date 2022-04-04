#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
#R script file that plots the parameter simulation results for the simple population model (Figures 4, 5, and A3) for "On the spread of microbes that manipulate reproduction in marine invertebrates" Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282. 
#Uses “Data_F_4A.feather” and “Data_EGMK_4B.feather”,  to make Figure 4A. 
#Uses “Data_F_4A_U.feather”, and  “Data_EGMK_4B_U.feather” to run density analyses which generates “Fem_D.feather”  and “MKEG_D.feather” and Figure 4B. 
#Uses “Data_FC_5.feather” to make Figure 5. Uses “Data_FC_A3.feather” to make Figure A3.
# 1. Setting up -----------------------------------------------------------

# * 1.a Loading up libraries ----------------------------------------------
library(tidyverse)
library(feather)
library(viridis)
library(patchwork)
library(scales)
# * 1.b clearing workspace ------------------------------------------------
rm(list = ls())
# * 1.c setting working directory ------------------------------------------------
# setwd("~/Desktop/my_papers/Projects/SeaUrchinInfection/Dryad")
# * 1.d Setting my default plotting theme ----------------------
mytheme <-
  theme_classic() + theme(
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
    strip.background = element_rect(
      color = "black", fill = "white", size =
        2
    ),
    strip.text = element_text(color = "black")
  ) # makes all the text larger and bold
theme_set(mytheme)
# * 1.e Vectorized equal function -----------------------------------------
# This is needed to more accurately apply logical == to numerical values in a fast vectorized tidyverse way.
# vectorized equal
# tol is how close doubles can be together (stands for tolerance)
# taken from https://stackoverflow.com/questions/35097815/vectorized-equality-testing
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}

# 2. Reading and processing data ------------------------------------------------
# * 2.a reading in data as a feather (file format used for reading in large files faster) ------------------------------------------------
# Figure 4a feminization
DatatA <- read_feather("Data/Data_F_4A.feather")
# Figure 4b enhanced growth + male killing
DatatB <- read_feather("Data/Data_EGMK_4B.feather")
# * 2.b processing data ------------------------------------------------
# Data for figure 4A
DatatA <- DatatA %>% # Converting egg size to diameter
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  # filtering out simulation runs that resulted in negative values due to chaos
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
# Add those parameter combinations back in, but as NA
Data2tA <- DatatA %>% complete(Egg_Size, Feminization, Larval_d)

# Data for figure 4b enhanced growth + male killing
DatatB <- DatatB %>% # Converting egg size to diameter
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  # filtering out simulation runs that resulted in negative values due to chaos
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
# Add those parameter combinations back in, but as NA
Data2tB <- DatatB %>% complete(Egg_Size, cost, Larval_d)
# 3. Making Figure 4a ------------------------------------------------
# Making the feminization plot for infection
p4ai <-
  ggplot(Data2tA[Data2tA$Larval_d >= 1, ], aes(x = Egg_Size, y = Larval_d)) +
  geom_tile(aes(fill = Infection)) +
  scale_fill_distiller(
    palette = "Blues",
    na.value = "grey",
    direction = 1,
    limits = c(0, 1),
    labels = percent
  ) +
  labs(y = "Pelagic larval duration (days)", fill = "Percent \n   infected") +
  theme(legend.position = "right") +
  xlab(bquote("Egg Diameter (" * mu * "m)")) +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80)
  ) +
  theme(aspect.ratio = 1) +
  ggtitle("Feminization") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 12))
# making figure 4b (male killing + enhanced growth infection)
p4bi <-
  ggplot(Data2tB[Data2tB$Larval_d >= 1, ], aes(x = Egg_Size, y = Larval_d)) +
  geom_tile(aes(fill = Infection), show.legend = F) +
  scale_fill_distiller(
    palette = "Blues",
    na.value = "grey",
    direction = 1,
    limits = c(0, 1)
  ) +
  labs(y = "Pelagic larval duration (days)", fill = "Proportion \n   infected") +
  theme(legend.position = "right") +
  xlab(bquote("Egg Diameter (" * mu * "m)")) +
  scale_x_continuous(expand = c(0.0, 0.0)) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80)
  ) +
  theme(aspect.ratio = 1) +
  ggtitle("Male killing + enhanced growth") +
  theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 12))

# putting it together to make panel A of figure 4
Figure4a <-
  ((p4ai + p4bi) &
    theme(legend.position = "right")) + plot_layout(guides = "collect")

Figure4a

# 4. Making Figure 5. Feminization/Cost -----------------------------------
# 4.a Reading in/processing data ------------------------------------------------
# reading in feather file
Data5 <- read_feather("Data/Data_FC_5.feather")
# converting egg size and filtering out runs that resulted in weird numerical values due to very high growth rates.
Data5 <- Data5 %>%
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  filter(
    Density > 0 |
      is_equal_tol(Density, 0),
    Sex_ratio > 0 |
      is_equal_tol(Sex_ratio, 0),
    Infection > 0 |
      is_equal_tol(Infection, 0),
    Infection < 1 |
      is_equal_tol(Infection, 1),
    Sex_ratio < 1 | is_equal_tol(Sex_ratio, 1)
  )
# * 4.b Summarize data into categories of infection -----------------------------------------------------
# See what combinations of cost and feminization result in complete or mostly infected (greater than 99.9%)
sum_data_inf5 <- Data5 %>%
  group_by(cost, Feminization) %>%
  summarize(
    Infectedall = ifelse(any(is_equal_tol(Infection, 1)), "Yes", "No"),
    InfectedMost = ifelse(any(Infection >= 0.999), "Yes", "No")
  )

# 4c. Plots ----------------------------------------------------------------
Fig5 <- ggplot(sum_data_inf5, aes(x = Feminization, y = cost)) +
  geom_tile(aes(fill = InfectedMost)) +
  labs(x = "Feminization rate", y = "Fitness cost (%)", fill = ">99.9% infection possible?") +
  scale_fill_manual(values = c("#deebf7", "#08519c")) +
  theme(legend.position = "right") +
  scale_x_continuous(
    expand = c(0, 0), limits =
      c(0.501, 1)
  ) +
  scale_y_continuous(expand = c(0, 0), labels = percent_format(accuracy = 1)) +
  theme(aspect.ratio = 1)
Fig5

# ggsave("Figure_5_cost.png",Fig5,height = 183,width=183,units="mm",path="")

# 5. Feminization/cost supplement Figure A3 ----------------------------------
Datafc <- read_feather("Data/Data_FC_A3.feather")
# * 5.a processing data ------------------------------------------------
Datafc <- Datafc %>% # Converting egg size to diameter
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  # filtering out simulation runs that resulted in negative values due to chaos
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
# Add those parameter combinations back in, but as NA
Data2fc <- Datafc %>% complete(Egg_Size, Feminization, Larval_d, cost)

# * 5.b Plot -------------------------------------------------------------
#make cost a factor so we can facet it
Data2fc$cost <- as.factor(Data2fc$cost)
piFC <-
  ggplot(Data2fc[Data2fc$Larval_d > 1 &
    Data2fc$cost %in% c(0., 0.1, 0.2, 0.3, 0.4, 0.5), ], aes(
    x = Egg_Size, y =
      Larval_d
  )) +
  geom_tile(aes(fill = Infection)) + #+facet_grid(.~Feminization)+
  scale_fill_distiller(
    palette = "Blues",
    na.value = "grey",
    direction = 1,
    limits = c(0, 1),
    labels = percent
  ) +
  labs(y = "Pelagic larval duration (days)", fill = "Proportion \n   infected") +
  theme(legend.position = "right") +
  xlab(bquote(bold("Egg Diameter (" * mu * "m)"))) +
  scale_x_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_y_continuous(expand = c(0.0, 0.0)) +
  facet_grid(cost ~ Feminization)
piFC


# 6. Making Figure 4B (comparing density from uninfected populations to infected populations) --------------------------------------------------------------
#  Reading and processing data ------------------------------------------------
# * 6.a reading in data as a feather ------------------------------------------------
# Figure 4b feminization uninfected data 
DatatAu <- read_feather("Data/Data_F_4A_U.feather")
# Figure 4b enhanced growth + male killing
DatatBu <- read_feather("Data/Data_EGMK_4B_U.feather")

# * 6.b processing data ------------------------------------------------
# Data for figure 4B Feminization
DatatAu <- DatatAu %>% # Converting egg size to diameter
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  # filtering out simulation runs that resulted in negative values due to chaos
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
# Add those parameter combinations back in, but as NA
Data2tAu <- DatatAu %>% complete(Egg_Size, Feminization, Larval_d)

# Data for figure 4B enhanced growth + male killing
DatatBu <- DatatBu %>% # Converting egg size to diameter
  mutate(Egg_Size = round((((
    Egg_Size * 4
  ) / pi)^0.5) * 1000)) %>%
  # filtering out simulation runs that resulted in negative values due to chaos
  filter(
    Density >= 0,
    Sex_ratio >= 0,
    Infection >= 0,
    Infection <= 1,
    Sex_ratio <= 1
  )
# Add those parameter combinations back in, but as NA
Data2tBu <- DatatBu %>% complete(Egg_Size, cost, Larval_d)

# Merging uninfected data with infected data for feminzation to compare densities
#<100001 indicates that it didn't reach max timelimit due to cycles. Also focusing on combinations that resulted in infections
Data3tAu <- Data2tA %>%
  select(Egg_Size, Larval_d, Time, Density, Infection) %>%
  left_join(Data2tAu,
    by = c("Egg_Size", "Larval_d"),
    suffix = c("I", "U")
  ) %>%
  mutate(
    DiffDensity = DensityI - DensityU,
    l2DiffDensity = log2(DensityI / DensityU)
  ) %>%
  mutate(
    DiffDensity = ifelse(InfectionI > 0.25 &
      (TimeI < 100001 & TimeU < 100001), DiffDensity, NA),
    l2DiffDensity = ifelse(
      InfectionI > 0.25 &
        (TimeI < 100001 & TimeU < 100001),
      l2DiffDensity,
      NA
    )
  )

# Merging uninfected data with infected data for male killing + enhanced growth to compare densities
#<100001 indicates that it didn't reach max timelimit due to cycles. Also focusing on combinations that resulted in infections
Data3tBu <- Data2tB %>%
  select(Egg_Size, Larval_d, Time, Density, Infection) %>%
  left_join(Data2tBu,
    by = c("Egg_Size", "Larval_d"),
    suffix = c("I", "U")
  ) %>%
  mutate(
    DiffDensity = DensityI - DensityU,
    l2DiffDensity = log2(DensityI / DensityU)
  ) %>%
  mutate(
    DiffDensity = ifelse(InfectionI > 0.25 &
      (TimeI < 100001 & TimeU < 100001), DiffDensity, NA),
    l2DiffDensity = ifelse(
      InfectionI > 0.25 &
        (TimeI < 100001 & TimeU < 100001),
      l2DiffDensity,
      NA
    )
  )

#write_feather(Data3tAu, "Data/Fem_D.feather")
#write_feather(Data3tBu, "Data/MKEG_D.feather")
# * 6.c Figure 4b  Benefit/cost of infection -----------------------------------
# getting limits for outliers for feminization and rounding to nearest whole number
limA <-
  round(max(
    IQR(Data3tAu$l2DiffDensity, na.rm = T) + quantile(Data3tAu$l2DiffDensity,
      na.rm =
        T, c(0.25, 0.75)
    )
  ), 0)
# getting limits for outliers for male killing and rounding to nearest whole number
limB <-
  round(max(
    IQR(Data3tBu$l2DiffDensity, na.rm = T) + quantile(Data3tBu$l2DiffDensity,
      na.rm =
        T, c(0.25, 0.75)
    )
  ), 0)

p4au <-
  ggplot(Data3tAu[Data3tAu$Larval_d >= 1, ], aes(x = Egg_Size, y = Larval_d)) +
  geom_tile(aes(fill = l2DiffDensity)) + #+facet_grid(.~Feminization)+
  scale_fill_distiller(
    palette = "PRGn",
    na.value = "grey",
    direction = 1,
    limits = c(-1, 1) * max(limA, limB),
    oob = squish,
    breaks = c(-3, -1.5, 0, 1.5, 3),
    labels = c("<-3", "-1.5", "0", "1.5", ">3")
  ) +
  labs(y = "Pelagic larval duration (days)", fill = bquote("Log"[2] ~ "FC"[Density])) +
  theme(legend.position = "right") +
  xlab(bquote("Egg Diameter (" * mu * "m)")) +
  scale_x_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80)
  ) +
  theme(aspect.ratio = 1) +
  ggtitle("Feminization") +
  theme(plot.title = element_text(
    hjust = 0.5, face =
      "italic", size = 12
  ))

p4bu <-
  ggplot(Data3tBu[Data3tBu$Larval_d >= 1, ], aes(x = Egg_Size, y = Larval_d)) +
  geom_tile(aes(fill = l2DiffDensity)) + #+facet_grid(.~Feminization)+
  #+facet_grid(.~Feminization)+
  scale_fill_distiller(
    palette = "PRGn",
    na.value = "grey",
    direction = 1,
    limits = c(-1, 1) * max(limA, limB),
    oob = squish,
    breaks = c(-3, -1.5, 0, 1.5, 3),
    labels = c("<-3", "-1.5", "0", "1.5", ">3")
  ) +
  labs(y = "Pelagic larval duration (days)", fill = bquote("Log"[2] ~ "FC"[Density])) +
  theme(legend.position = "right") +
  xlab(bquote("Egg Diameter (" * mu * "m)")) +
  scale_x_continuous(
    expand =
      c(0.0, 0.0)
  ) +
  scale_y_continuous(
    expand = c(0.0, 0.0),
    limits = c(0, 80),
    breaks = c(0, 20, 40, 60, 80)
  ) +
  theme(aspect.ratio = 1) +
  ggtitle("Male killing + enhanced growth") +
  theme(plot.title = element_text(
    hjust = 0.5, face =
      "italic", size = 12
  ))

# 7. Putting Figure 4A and Figure 4B together --------------------------------------------------
#infection heat maps of feminization and male killing + enhanced growth
Figure4 <-
  ((p4ai + p4bi) &
    theme(legend.position = "right")) + plot_layout(guides = "collect")
#Change in density heat maps for feminization and male killing + enhanced growth
Figure4u <-
  ((p4au + p4bu) &
    theme(legend.position = "right")) + plot_layout(guides = "collect")

(Figure4 / Figure4u) + plot_annotation(tag_levels = c("A"))
#Note that we edited the tag levels manually for the figure in the paper as well as reformated the legends.
# ggsave("Figure_4_All.pdf",height = 183,width=183,units="mm")
# ggsave("Figure_4_All.svg",height = 183,width=183,units="mm")

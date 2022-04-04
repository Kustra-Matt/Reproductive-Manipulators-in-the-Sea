#Code by Matthew Kustra. For questions: mkustra@ucsc.edu
# R script file that performs the species analysis described in the paper:
#"On the spread of microbes that manipulate reproduction in marine invertebrates" 
#Kustra & Carrier 2022 https://www.journals.uchicago.edu/doi/10.1086/720282. 
# Uses results from “Data_F_4A.feather”, “Data_EGMK_4B.feather”,  “Fem_D.feather” and  “MKEG_D.feather” to see whether known marine invertebrate species (“Species_data.csv”) are likely to be infected,
# whether or not infection is associated with phylum or developmental mode,
# and if the change in density is influenced by phylum or developmental mode.
# This script makes Figure A10 and Figure A11.
# 1.a Loading up libraries -------------------------------------------------
library(MASS)
library(tidyverse)
library(reshape2)
library(car)
library(scales)
library(statmod)
library(feather)
# * 1.b Plotting theme ----------------------------------------------------
mytheme <-
  theme_classic() + theme(
    legend.position = "bottom",
    # this puts legend on the bottom
    axis.line = element_line(
      color = "black", size =
        0
    ),
    axis.text = element_text(color = "black"),
    # Makes the axis line black and  thicker
    text = element_text(size = 12, color = "black"),
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_blank(),
    strip.text = element_text(color = "black", face = "bold")
  ) # makes all the text larger and bold
theme_set(mytheme)

# * 1.c Making function needed to match----------------------------------------------------
# gets the value that is closest to a given value
matchclose <- function(x, vec) {
  vec[which.min(abs(x - vec))]
}
# vectorized equal
# tol is how close doubles can be together
# taken from https://stackoverflow.com/questions/35097815/vectorized-equality-testing
is_equal_tol <- function(x, y, tol = .Machine$double.eps) {
  abs(x - y) < tol
}
# 2 Load up and process data ---------------------------------------------------------
# Figure4a feminization results
DatatA <-
  read_feather("Data/Data_F_4A.feather") %>% # Converting egg size to diameter
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
  ) %>%
  # Add those parameter combinations back in, but as NA
  complete(Egg_Size, Feminization, Larval_d) %>%
  mutate(across(c("nif0", "nim0", "nuf0", "num0"), ~
  replace(., is_equal_tol(., 0), 0))) %>%
  mutate(
    Sex_ratio = (nim0 + num0) / (nim0 + num0 + nif0 + nuf0),
    Infection = (nif0 + nim0) / (nim0 + num0 + nif0 + nuf0)
  )

# Figure 4b enhanced growth + male killing results
DatatB <-
  read_feather("Data/Data_EGMK_4B.feather") %>% # Converting egg size to diameter
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
  ) %>%
  # Add those parameter combinations back in, but as NA
  complete(Egg_Size, Feminization, Larval_d) %>%
  mutate(across(c("nif0", "nim0", "nuf0", "num0"), ~
  replace(., is_equal_tol(., 0), 0))) %>%
  mutate(
    Sex_ratio = (nim0 + num0) / (nim0 + num0 + nif0 + nuf0),
    Infection = (nif0 + nim0) / (nim0 + num0 + nif0 + nuf0)
  )
# Species Data
SpData <- read.csv("Data/Species_data.csv")
# 3. Merging species data -------------------------------------------------------
# create vector of unique egg sizes
eggs <- as.numeric(unique(DatatA$Egg_Size))
# create vector of unique larval days
days <- as.numeric(unique(DatatA$Larval_d))
# create new column with eggsize now made to be closest egg size simulated at
SpData$EggSize2 <-
  as.numeric(sapply(SpData$EggSize, FUN = matchclose, vec = eggs))
# create new column with eggsize now made to be closest planktonic duration
SpData$Days <- sapply(SpData$PlanktonicTime, FUN = matchclose, vec = days)

# *3.A Merging feminization sims (Fig.4A) ---------------------------------
# creating summary of infection for Feminization sims (4A)
sum_data_infA <- DatatA %>%
  group_by(Egg_Size, Larval_d) %>%
  summarize(InfectedMost = ifelse(any(Infection >= 0.999), "Yes", "No"))
# merge data frames together
togetherA <-
  left_join(sum_data_infA,
    SpData,
    by = c("Egg_Size" = "EggSize2", "Larval_d" = "Days")
  ) %>%
  drop_na(-InfectedMost) %>%
  mutate(InfectedMost = ifelse(is.na(InfectedMost), "No", InfectedMost)) %>% # filter out the Phyla we want
  filter(Phylum %in% c("Annelida", "Mollusca", "Echinodermata"))
# summarize the data
Sum_dataA <- togetherA %>%
  group_by(Phylum, DevelopmentalMode) %>%
  summarize(
    Num_infected = length(InfectedMost[InfectedMost == "Yes"]),
    Num_uninfected = length(InfectedMost[InfectedMost == "No"]),
    Total = length(InfectedMost),
    Prop = Num_infected / Total
  )

# *3.A Merging Malekilling + growth sims (Fig.4) ---------------------------------
# creating summary of infection for Feminization sims (4)
sum_data_infB <- DatatB %>%
  group_by(Egg_Size, Larval_d) %>%
  summarize(InfectedMost = ifelse(any(Infection >= 0.999), "Yes", "No"))
# merge data frames together
togetherB <-
  left_join(sum_data_infB,
    SpData,
    by = c("Egg_Size" = "EggSize2", "Larval_d" = "Days")
  ) %>%
  drop_na(-InfectedMost) %>%
  mutate(InfectedMost = ifelse(is.na(InfectedMost), "No", InfectedMost)) %>% # filter out the Phyla we want
  filter(Phylum %in% c("Annelida", "Mollusca", "Echinodermata"))
# summarize the data
Sum_dataB <- togetherB %>%
  group_by(Phylum, DevelopmentalMode) %>%
  summarize(
    Num_infected = length(InfectedMost[InfectedMost == "Yes"]),
    Num_uninfected = length(InfectedMost[InfectedMost == "No"]),
    Total = length(InfectedMost),
    Prop = Num_infected / Total
  )


# 4.a Analysis + Statistics for Feminization -----------------------------------------------
# make dummy column of 0 or 1
togetherA$InfectedMostDum <- ifelse(togetherA$InfectedMost == "Yes", 1, 0)
# fit binomial glm
modA <-
  glm(
    data = togetherA,
    InfectedMostDum ~ Phylum * DevelopmentalMode,
    family = binomial
  )
# check residuals of the model
qresA <- qresid(modA)
qqnorm(qresA, las = 1)
qqline(qresA)
# looks good
summary(modA)
# Want to do type 2 test
AmodA <- Anova(modA, test = "LR", type = "II")
AmodA

# 4.B Analysis + Statistics for Male killing + growth -----------------------------------------------
# make dummy column of 0 or 1
togetherB$InfectedMostDum <- ifelse(togetherB$InfectedMost == "Yes", 1, 0)
# fit binomial glm
modB <-
  glm(
    data = togetherB,
    InfectedMostDum ~ Phylum * DevelopmentalMode,
    family = binomial
  )
# check residuals of the model
qresB <- qresid(modB)
qqnorm(qresB, las = 1)
qqline(qresB)
summary(modB)
# Want to do type 2 test
AmodB <- Anova(modB, test = "LR", type = "II")
AmodB
# 5.a Plotting Feminization results (Figure A10A) ---------------------------
Plot_dataA <- togetherA %>%
  group_by(Phylum, DevelopmentalMode, InfectedMost) %>%
  summarize(Count = n()) %>%
  mutate(County = Count / sum(Count), Label = paste0("N = ", Count))

FigA10A <-
  ggplot(
    Plot_dataA,
    aes(x = DevelopmentalMode, y = Count, fill = InfectedMost)
  ) +
  geom_bar(
    stat = "identity",
    color = "black",
    position = "fill"
  ) +
  geom_text(aes(label = Label, y = County), position = position_stack(
    vjust =
      0.5
  )) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  facet_grid(. ~
  Phylum) +
  scale_fill_manual(values = c("#deebf7", "#08519c"), name = "Infection?") +
  scale_x_discrete(expand = c(0, 0)) +
  ylab("Percent") +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  ) +
  xlab("Developmental Mode")
FigA10A

# 5.b Plotting Male Killing + growth results (Figure A10B) ---------------------------
Plot_dataB <- togetherB %>%
  group_by(Phylum, DevelopmentalMode, InfectedMost) %>%
  summarize(Count = n()) %>%
  mutate(County = Count / sum(Count), Label = paste0("N = ", Count))

FigA10B <-
  ggplot(
    Plot_dataB,
    aes(x = DevelopmentalMode, y = Count, fill = InfectedMost)
  ) +
  geom_bar(
    stat = "identity",
    color = "black",
    position = "fill"
  ) +
  geom_text(aes(label = Label, y = County), position = position_stack(
    vjust =
      0.5
  )) +
  scale_y_continuous(labels = percent, expand = c(0, 0)) +
  facet_grid(. ~
  Phylum) +
  scale_fill_manual(values = c("#deebf7", "#08519c"), name = "Infection?") +
  scale_x_discrete(expand = c(0, 0)) +
  ylab("Percent") +
  theme(
    panel.spacing = unit(1, "lines"),
    legend.position = "right"
  ) +
  xlab("Developmental Mode")
FigA10B

FigA10A/ FigA10B + plot_annotation(tag_levels = "A")

# 6.Density diff analysis ---------------------------------------------------
Fem_D <- read_feather("Data/Fem_D.feather") %>%
  dplyr::select(Egg_Size, Larval_d, DiffDensity, l2DiffDensity)
MKEG_D <- read_feather("Data/MKEG_D.feather") %>%
  dplyr::select(Egg_Size, Larval_d, DiffDensity, l2DiffDensity)

# *6.a Merging species data -------------------------------------------------------
# create vector of unique egg sizes
eggs <- as.numeric(unique(DatatA$Egg_Size))
# create vector of unique larval days
days <- as.numeric(unique(DatatA$Larval_d))
# create new column with eggsize now made to be closest egg size simulated at
SpData$EggSize2 <-
  as.numeric(sapply(SpData$EggSize, FUN = matchclose, vec = eggs))
# create new column with eggsize now made to be closest planktonic duration
SpData$Days <- sapply(SpData$PlanktonicTime, FUN = matchclose, vec = days)

# *6.b Merging feminization sims (Fig.4B) ---------------------------------
# creating summary of infection for Feminization sims (4B)

# merge data frames together and filter out extreme outliers
togetherFD <-
  left_join(Fem_D, SpData, by = c("Egg_Size" = "EggSize2", "Larval_d" = "Days")) %>% # filter out the Phyla we want
  filter(
    Phylum %in% c("Annelida", "Mollusca", "Echinodermata"),
    abs(l2DiffDensity) < 20
  ) %>%
  mutate(Phylum = droplevels(factor(Phylum)))

#  *6.c plotting Feminization results (Figure A11A) -----------------------


a11A <-
  ggplot(
    togetherFD,
    aes(y = l2DiffDensity, x = DevelopmentalMode, fill = DevelopmentalMode)
  ) +
  geom_jitter(alpha = 0.3, aes(color = DevelopmentalMode)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~Phylum) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold")
  ) +
  scale_x_discrete(labels = c("Lecithotrophy", "Planktotrophy")) +
  xlab("Developmental Mode") +
  ylab(bquote(bold("Log"[2] ~ "FC"[Density]))) +
  scale_color_manual(values = c("#6D5B97", "#A8332A")) +
  scale_fill_manual(values = c("#6D5B97", "#A8332A"))

# *6.d Merging Malekilling + enhanced growth sims (Fig.4B) ---------------------------------
# creating summary of infection for Feminization sims (4B)

# merge data frames together
togetherMKEG_D <-
  left_join(MKEG_D, SpData, by = c("Egg_Size" = "EggSize2", "Larval_d" = "Days")) %>% # filter out the Phyla we want
  filter(
    Phylum %in% c("Annelida", "Mollusca", "Echinodermata"),
    abs(l2DiffDensity) < 20
  ) %>%
  mutate(Phylum = droplevels(factor(Phylum)))
#  *6.e Plotting Malekilling + enhanced growth sims (Figure A11B) ---------
a11B <-
  ggplot(
    togetherMKEG_D,
    aes(y = l2DiffDensity, x = DevelopmentalMode, fill = DevelopmentalMode)
  ) +
  geom_jitter(alpha = 0.3, aes(color = DevelopmentalMode)) +
  geom_boxplot(alpha = 0.2, outlier.shape = NA) +
  facet_wrap(~Phylum) +
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold")
  ) +
  scale_x_discrete(labels = c("Lecithotrophy", "Planktotrophy")) +
  xlab("Developmental Mode") +
  ylab(bquote(bold("Log"[2] ~ "FC"[Density]))) +
  scale_color_manual(values = c("#6D5B97", "#A8332A")) +
  scale_fill_manual(values = c("#6D5B97", "#A8332A"))


a11A / a11B + plot_annotation(tag_levels = "A")

ggsave("FigA11.pdf",
  height = 183,
  width = 183,
  units = "mm"
)


# 7.a Analysis + Statistics for Feminization -----------------------------------------------
# permutation difference test function
# g1 is the vector of one group, g2 is vector of group 2, and n is the number of reps
perm_diff_means <- function(g1, g2, n) {
  # calculate real difference
  real_diff <- mean(g1) - mean(g2)
  # concatenate the two vectors together
  all <- c(g1, g2)
  # create a vector of the possible indicies
  all_indices <- seq(1, length(all), 1)
  # calculate the length of one vector, allows us to use this function for data with unequal group sample size
  g1_length <- length(g1)
  # make a vector to store the fake test statistics (differences in the mean between two fake groups)
  sim_diffs <- numeric(length = n)
  for (i in 1:n) {
    # does the loop for as many simulations as specified
    # Get the indicies to make into fake group 1
    fake_g1_indices <- sample(all_indices, g1_length)
    # make fake group 1
    fake_g1 <- all[fake_g1_indices]
    # make fake group 2 by getting all other data that wasn't assigned into group 1
    fake_g2 <- all[setdiff(all_indices, fake_g1_indices)]
    # calculate the difference in the average between these two groups
    sim_diffs[i] <- mean(fake_g1) - mean(fake_g2)
  }
  # make histogram of fake data
  hist(sim_diffs)
  # plot the real data
  abline(v = real_diff, col = "red")
  # calculate pvalue by looking at whether or not observed data was as extreme
  # this assumes lower.tail, but to make upper tail subtract this value from 1.
  pvalue <- length(sim_diffs[sim_diffs <= real_diff]) / n
  # ways of estimating 95% CI on simulated p values from Higgens Introduction to Modern Nonparametric Statistics
  lower <- max(0, pvalue - (2 * sqrt(pvalue * (1 - pvalue) / n)))
  upper <- min(1, pvalue + (2 * sqrt(pvalue * (1 - pvalue) / n)))
  return(
    list(
      "sim_diffs" = sim_diffs,
      "Real Diff" = real_diff,
      "p-value" = pvalue,
      "95% CI of pvalue estimation" = paste("[", lower, ",", upper, "]",
        sep =
          ""
      )
    )
  )
}
# filter out for only echinoderms
togetherFD_E <- togetherFD[togetherFD$Phylum == "Echinodermata", ]

# look at distribution
ggplot(togetherFD_E, aes(x = l2DiffDensity)) +
  geom_histogram() +
  facet_grid(~DevelopmentalMode)
# definitely not normally distributed, so need to perform a permutation test
#Set random seed so it can be replicated
set.seed(12)
FD_mod <-
  perm_diff_means(
    togetherFD_E$l2DiffDensity[togetherFD_E$DevelopmentalMode ==
      "Planktotrophic"],
    togetherFD_E$l2DiffDensity[togetherFD_E$DevelopmentalMode == "Lecithotrophic"],
    10000
  )
# 7.b Analysis + Statistics for Malekilling + enhanced growth -----------------------------------------------
# filter for only echinoderms
togetherMKEG_D_E <-
  togetherMKEG_D[togetherMKEG_D$Phylum == "Echinodermata", ]
# look at histogram
ggplot(togetherMKEG_D_E, aes(x = l2DiffDensity)) +
  geom_histogram() +
  facet_grid(~DevelopmentalMode)
# definintely not normally distriubted so doing a permutation test
set.seed(12)
MKEGD_mod <-
  perm_diff_means(
    togetherMKEG_D_E$l2DiffDensity[togetherMKEG_D_E$DevelopmentalMode ==
      "Planktotrophic"],
    togetherMKEG_D_E$l2DiffDensity[togetherMKEG_D_E$DevelopmentalMode == "Lecithotrophic"],
    10000
  )

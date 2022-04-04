# Reproductive Manipulators in the Sea
 Code associated with Kustra & Carrier. 2022. On the Spread of microbes that manipulate reproduction in marine invertebrates. American Naturalist. <https://www.journals.uchicago.edu/doi/10.1086/720282>.

  For data files, either email me or go to dryad repository: <https://doi.org/10.7291/D15Q4X>.

**[Code:]{.underline}**

-   **"Fertilization.R":** R script file that runs and plots the
    fertilization dynamics. Runs analyses and makes plots for Figure 1,
    Figure A4, and Figure A5. Makes/uses "**Fertilization_opt_f.csv"**
    and **"Fertilization_opt_z.csv".**

-   **"Population_Dynamics_Simple.R":** R script file that runs and
    plots the simple population model. Currently coded to simulate and
    plot Figure 2 and Figure A6, but users can run/plot the model under
    different parameter combinations.

-   **"Stage_Based_Model.R":** R script file that runs and plots the
    population model with different size classes. Currently coded to
    simulate and plot Figure 3, Figure A7, Figure A8, and Figure A9 but
    users can run/plot the model under different parameter combinations.
    Creates/uses "**Sensitivity_Stage_Based.csv"** to make Figure A8.

-   **"Sims_Simple_Local.R**": R script that was used for numerical
    simulation across life history traits (egg size and mortality) using
    the simple population model. It uses **"Larval_data.csv**" which
    contains the data for larval mortality. The results from running
    this script are saved as "**Data_F\_4A.feather**",
    **"Data_EGMK_4B.feather",** "**Data_F\_4A_U.feather**", and
    **"Data_EGMK_4B_U.feather"**, which are used to make Figure 4 by
    "**Plotting_Sim_Heatmaps.R**".

-   **"Sims_Simple_SuperComp.R**": R script that was used for numerical
    simulation across life history traits (egg size and mortality),
    costs, and feminization rates using the simple population model. It
    uses **"Larval_data.csv**" which contains the data for larval
    mortality. The results from running this script are saved as
    **Data_FC_5.feather**" and **Data_FC_A3.feather**, which are used to
    make Figure 5 and Figure A3 by "**Plotting_Sim_Heatmaps.R**".

-   "**Plotting_Sim_Heatmaps.R**": R script file that plots the
    parameter simulation results for the simple population model. Uses
    **"Data_F\_4A.feather**" and **"Data_EGMK_4B.feather",** to make
    Figure 4A. Uses "**Data_F\_4A_U.feather**", and
    **"Data_EGMK_4B_U.feather"** to run density analyses which generates
    **"Fem_D.feather"** and **"MKEG_D.feather"** and Figure 4B. Uses
    **"Data_FC_5.feather**" to make Figure 5. Uses
    **"Data_FC_A3.feather**" to make Figure A3.

-   "**Species_Analysis.R**": R script file that performs the species
    analysis described in the paper. Uses results from
    **"Data_F\_4A.feather**", **"Data_EGMK_4B.feather**",
    **"Fem_D.feather**" and **"MKEG_D.feather**" to see whether known
    marine invertebrate species (**"Species_data.csv**") are likely to
    be infected, whether or not infection is associated with phylum or
    developmental mode, and if the change in density is influenced by
    phylum or developmental mode. This script also makes Figure A10 and
    Figure A11.

**[Data:]{.underline}** Folder with all the different data files used to
make plots/analyze data

-   "**Fertilization_opt_z.csv**": Data file of maximum zygote
    production at each sex ratio for both species and the density at
    that sex ratio that maximizes zygote production. The R code file
    "**Fertilization.R**" creates this file but since computation time
    is long the R code can load this csv to make Figure A5. Here are the
    descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Species*:** Species name that data is for: either *Heliocidaris
    erythrogramma* or *Heliocidaris tuberculata.*

-   ***Dens*:** Adult density that maximizes zygote production at a
    given sex ratio.

-   ***Zygote*:** Maximum zygote production at a given sex ratio.

-   ***sexratio*:** Proportion of males

```{=html}
<!-- -->
```
-   **"Fertilization_opt_f.csv**": Data file of maximum fertilization
    percentage at each sex ratio for both species and the density at
    that sex ratio that maximizes fertilization percentage. The R code
    file "**Fertilization.R**" creates this file but since computation
    time is long the R code can load this csv to make Figure A5. Here
    are the descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Species*:** Species name that data is for: either *Heliocidaris
    erythrogramma* or *Heliocidaris tuberculata.*

-   ***Dens*:** Adult density that maximizes fertilization percentage at
    a given sex ratio.

-   ***Fert*:** Maximum fertilization percentage

-   ***sexratio*:** Proportion of males

```{=html}
<!-- -->
```
-   **"Sensitivity_Stage_Based.csv**": Data file of equilibrium
    densities of different subpopulations for stage-based model
    simulations across different male killing rates and infected female
    transition (growth) rates. The R code file "**Stage_Based_Model.R**"
    creates this file but since computation time is long the R code can
    load this csv to make Figure A8. Here are the descriptions of the
    column variables:

```{=html}
<!-- -->
```
-   ***f1i:*** Equilibrium density of 20-39mm infected females

-   ***f1u:*** Equilibrium density of 20-39mm uninfected females

-   ***m1i:*** Equilibrium density of 20-39mm infected males

-   ***m1u:*** Equilibrium density of 20-39mm uninfected males

-   ***f2i:*** Equilibrium density of 40-59mm infected females

-   ***f2u* :** Equilibrium density of 40-59mm uninfected females

-   ***m2i:*** Equilibrium density of 40-59mm infected males

-   ***m2u:*** Equilibrium density of 40-59mm uninfected males

-   ***f3i:*** Equilibrium density of 60+ mm infected females

-   ***f3u:*** Equilibrium density of 60+ mm uninfected females

-   ***m3i:*** Equilibrium density of 60+ mm infected males

-   ***m3u:*** Equilibrium density of 60+ mm uninfected males

-   ***Generation*:** Approximate generation where equilibrium was
    reached

-   ***Mk*:** Male killing rate

-   ***tif*:** Infected female transition (growth) rate

```{=html}
<!-- -->
```
-   **"Larval_data.csv**": Data file of larval duration and larval
    mortality. The R code file "**Simple_Sims.R**" uses this file. Here
    are the descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Days*:** Days of pelagic larval duration.

-   ***Mortality*:** Expected larval mortality

```{=html}
<!-- -->
```
-   **"Data_F\_4A.feather**": Data file of the results from simulations
    of the simple population model across different egg sizes, and
    larval duration. The R code file "**Plotting_Sim_Heatmaps.R**" uses
    this data file to make figure 4A. The R code file
    "**Species_Analysis.R"** uses this to make figure A10. This data
    file was saved as a feather file to allow faster loading into R. The
    file was generated using "**Sims_Simple_Local.R".** Here are the
    descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium total population density

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection

```{=html}
<!-- -->
```
-   **"Data_F\_4A_U.feather**": Data file of the results from
    simulations of the simple population model across different egg
    sizes, and larval duration when there is no infection. The R code
    file "**Plotting_Sim_Heatmaps.R**" uses this data file to make
    figure 4B. This data file was saved as a feather file to allow
    faster loading into R. The file was generated using
    "**Sims_Simple_Local.R".** Here are the descriptions of the column
    variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium total population density

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection

```{=html}
<!-- -->
```
-   **"Fem_D.feather**": Data file of the difference in density between
    uninfected an infected from simulations of the simple population
    model across different egg sizes, and larval duration. The R code
    file "**Plotting_Sim_Heatmaps.R**" uses this data file to make
    figure 4B. The R code file "**Species_Analysis.R"** uses this to
    make figure A11. This data file was saved as a feather file to allow
    faster loading into R. The file was generated using
    "**Plotting_Sim_Heatmaps.R**". Here are the descriptions of the
    column variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***TimeI*:** Time it took for equilibrium to be reached or an error
    for infected runs

-   ***DensityI*:** Equilibrium total population density of infected
    runs

-   ***InfectionI*:** Proportion of infected individuals at end of
    simulation for infected runs

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error
    for uninfected runs

-   ***Density*:** Equilibrium total population density of uninfected
    runs

-   ***Infection*:** Proportion of infected individuals at end of
    simulation of uninfected runs

-   ***Sex_ratio*:** Sex ratio at end of simulation of uninfected runs

-   ***Error*:** What kind of errors were produced of uninfected runs

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females of uninfected
    runs

-   ***nim0*:** Equilibrium density of infected males of uninfected runs

-   ***nuf0*:** Equilibrium density of uninfected females of uninfected
    runs

-   ***num0*:** Equilibrium density of uninfected males of uninfected
    runs

-   ***cost*:** Cost of infection of uninfected runs

-   ***DiffDensity*:** Difference in density between infected and
    uninfected populations

-   ***l2DiffDensity*:** Log~2~ fold change in density between infected
    and uninfected populations

```{=html}
<!-- -->
```
-   **"Data_EGMK_4B.feather**": Data file of the results from
    simulations of the simple population model across different egg
    sizes, and larval duration. The R code file
    "**Plotting_Sim_Heatmaps.R**" uses this data file to make figure 4A.
    The R code file "**Species_Analysis.R"** uses this to make figure
    A10. This data file was saved as a feather file to allow faster
    loading into R. The file was generated using
    "**Sims_Simple_Local.R".** Here are the descriptions of the column
    variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium total population density

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection (values are negative to create a
    benefit of infection)

```{=html}
<!-- -->
```
-   **"Data_EGMK_4B_U.feather**": Data file of the results from
    simulations of the simple population model across different egg
    sizes, and larval duration. The R code file
    "**Plotting_Sim_Heatmaps.R**" uses this data file to make figure 4B.
    The R code file "**Species_Analysis.R"** uses this to make figure
    A11. This data file was saved as a feather file to allow faster
    loading into R. The file was generated using
    "**Sims_Simple_Local.R".** Here are the descriptions of the column
    variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium total population density

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection (values are negative to create a
    benefit of infection)

```{=html}
<!-- -->
```
-   **"MKEG_D.feather**": Data file of the difference in density between
    uninfected an infected from simulations of the simple population
    model across different egg sizes, and larval duration. The R code
    file "**Plotting_Sim_Heatmaps.R**" uses this data file to make
    figure 4B. The R code file "**Species_Analysis.R"** uses this to
    make figure A11B. This data file was saved as a feather file to
    allow faster loading into R. The file was generated using
    "**Plotting_Sim_Heatmaps.R**" Here are the descriptions of the
    column variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***TimeI*:** Time it took for equilibrium to be reached or an error
    for infected runs

-   ***DensityI*:** Equilibrium total population density of infected
    runs

-   ***InfectionI*:** Proportion of infected individuals at end of
    simulation for infected runs

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error
    for uninfected runs

-   ***Density*:** Equilibrium total population density for uninfected
    populations

-   ***Infection*:** Proportion of infected individuals at end of
    simulation for uninfected runs

-   ***Sex_ratio*:** Sex ratio at end of simulation of uninfected runs

-   ***Error*:** What kind of errors were produced of uninfected runs

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females of uninfected
    runs

-   ***nim0*:** Equilibrium density of infected males of uninfected runs

-   ***nuf0*:** Equilibrium density of uninfected females of uninfected
    runs

-   ***num0*:** Equilibrium density of uninfected males of uninfected
    runs

-   ***cost*:** Cost of infection of uninfected runs

-   ***DiffDensity*:** Difference in density between infected and
    uninfected populations

-   ***l2DiffDensity*:** Log~2~ fold change in density between infected
    and uninfected populations

```{=html}
<!-- -->
```
-   **"Data_FC_5.feather**": Data file of the results from simulations
    of the simple population model across different egg sizes, larval
    duration, feminization rates, and costs. The R code file
    "**Plotting_Sim_Heatmaps.R**" uses this data file to make Figure 5.
    This data file was saved as a feather file to allow faster loading
    into R. The file was generated using "**Sims_Simple_SuperComp.R".**
    Here are the descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium density of 40-59mm infected males

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection

```{=html}
<!-- -->
```
-   **"Data_FC_A3.feather**": Data file of the results from simulations
    of the simple population model across different egg sizes, larval
    duration, feminization rates, and costs. The R code file
    "**Plotting_Sim_Heatmaps.R**" uses this data file to make figure A3.
    This data file was saved as a feather file to allow faster loading
    into R. The file was generated using "**Sims_Simple_SuperComp.R".**
    Here are the descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Egg_Size*:** Egg size parameter

-   ***Larval_d*:** Larval duration

-   ***Larval_m*:** Larval mortality

-   ***Egg number*:** Egg density produced by a female

-   ***Feminization:*** Feminization rate

-   ***Time*:** Time it took for equilibrium to be reached or an error

-   ***Density*:** Equilibrium density of 40-59mm infected males

-   ***Infection*:** Proportion of infected individuals at end of
    simulation

-   ***Sex_ratio*:** Sex ratio at end of simulation

-   ***Error*:** What kind of errors were produced

    -   *"NAs"*: Na's were produced due to too high growth

    -   *"Exceeded Time*": simulation went past limit of 100,000
        generations, which was generally indicative of a limit cycle.

    -   ***"**None"*: means no errors occurred.

-   ***nif0*:** Equilibrium density of infected females

-   ***nim0*:** Equilibrium density of infected males

-   ***nuf0*:** Equilibrium density of uninfected females

-   ***num0*:** Equilibrium density of uninfected males

-   ***cost*:** Cost of infection

-   ***init_inf_dens*:** gives original ratio of infected to uninfected
    densities.

```{=html}
<!-- -->
```
-   **"Species_data.csv**": Data file of species data for marine
    invertebrates (Álvarez-Noriega et al. 2020; Marshall et al. 2012).
    The R code file "**Species_Analysis.R**" uses this file to make
    Figure A9. Here are the descriptions of the column variables:

```{=html}
<!-- -->
```
-   ***Phylum*:** Phylum of the species

-   ***Class*:** Class of the species

-   ***Genus:*** Genus of the species

-   ***Species*:** Species name of the species

-   ***DevelopmentalMode*:** Larval developmental mode. Either
    *Lecithotrophic or Planktotrophic*

-   ***EggSize*:** Egg diameter (um) of the species

-   ***PlanktonicTime*:** The number of days larvae are planktonic

-   ***Longitude*:** Longitude of species distribution

-   ***Latitutde*:** Latitude of species distribution

-   ***Reference*:** Reference for this information.

**[SI_Web_App:]{.underline}** Folder that contains all the files and
code needed to generate the supplemental web app.

-   "**app.R**": The R code file that creates the supplemental web app.
    Code is not needed to run the web app but can be run locally if
    needed.

-   "**data**": Folder that contains data files that are used to run the
    web app.

    -   "**Data_3\_3_21_01_99_all.feather":** Data file used for
        plotting species on feminization analysis (Figure 4A) in the
        panel: "Marine invertebrates Fig. 4"

    -   "**Data_4b_8\_22_21.feather":** Data file used for plotting
        species on Male Killing + Enhanced growth analysis (Figure 4B)
        in the panel: "Marine invertebrates Fig. 4"

    -   **"species2.csv":** Data file with species data used for
        plotting species in the panel: "Marine invertebrates Fig. 4"

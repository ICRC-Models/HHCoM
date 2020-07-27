Generic CEA code
=======================

Administered by Marina Antillon, PhD.
Health Systems Research Unit  
Swiss Tropical and Public Health (Swiss TPH) Institute  
(An affiliated institute of the University of Basel)  
Basel, Switzerland

COPYRIGHT 2020, Swiss TPH

---

# Project Objective 

To provide tools to pursue cost-effectiveness analyses in R.

---

# Brief description

## Analysis tools

The analysis is broadly defined in four parts, each of these parts is executed in various R code files:  
I. Projecting clinical outcomes under alternative strategies with existing and upcoming treatments.  
II. Costing: for clinical activities as well as screening and prevention (vector control) activities under alternative strategies.  
III. Cost-effectiveness analysis of alternative strategies.  
IV. Value-of-information analysis, tipping point (threshold analysis)   

The tools for the economic model are coded in R. I would recommend tracking the project via git, though as of this writing (June 2020) I still find it challenging to get used to that. I would recommend that any final code posted in the internet include a package library by using Packrat, but again, I don't follow my own advice in this respect.

While I would highly recommend using the code within RStudio environment (in part because of it's features to manage the project with .RProj, git, and packrat) this is not strictly necessary and the benefits of git and packrat are available from a classic R interface and a shell command line. I personally use and love RStudio.

## File structure and functions

There is a master file that runs all the code.
The master file is separated into sections and those sections are separated into subsections:
- A) Set up for the section of the analysis
- B) Load output from previous sections of the analysis (sometimes NA)
- C) Run code
- D) Save output
- E) Summary tables (sometimes NA)

## Purpose of each section of the master file
__Section 0:__ Set the working directory, make the directories where the output will be stored, set the seed of the analysis, and make labels of the interventions
__Section 1:__ Set up of parameters for analysis without uncertainty 
__Section 2:__ Do CEA without uncertainty
__Section 3:__ Set up analysis with uncertainty 
__Section 4:__ Simulate outcomes with uncertainty
__Section 5:__ CEAs with uncertainty
__Section 6:__ VOI/EVPPI
__Section 7:__ Threshold analysis

 ## Outputs
 __Section 0:__ Set the working directory, make the directories where the output will be stored, set the seed of the analysis, and make labels of the interventions
     - output_data, output_fig_tbls: directory paths to store outputs
     - myseed: seed of the analysis (for reproducibility)
     - lbl: list of labels (for now, only the labels of the interventions)
     - keep: the list of objects in the r environment that I want to have all the time (this is to clean out other things with the rm command)
__Section 1:__ Set up of parameters for analysis without uncertainty - 
     - fixedpars: list of parameters for the analysis without uncertainty
__Section 2:__ Do CEA without uncertainty
     - prob: probabilities for the probability tree
     - ylls, ylds, dalys: measures of health effect
     - cost_meds, cost_treat, costs: costs - both the breakdown and the total
     - sum.df: summaries of intermediate outcomes
     - ICER.df: decision analysis
 __Section 3:__ Set up analysis with uncertainty 
     - iter: number of iterations
     - par: parameters for analysis with uncertainty
 __Section 4:__ Simulate outcomes with uncertainty
     - prob: probabilities for the probability tree
     - ylls, ylds, dalys: measures of health effect
     - cost_meds, cost_treat, costs: costs - both the breakdown and the total
     - sum.df: summaries of intermediate outcomes
 __Section 5:__ CEAs with uncertainty
     - thresholds
     - nmb, nmb_means: net monetary benefits to make ceac's and ceaf's (respectively)
     - ceac, ceac_pre, ceac_df, ceac_pie: arrays for cost-effectiveness acceptability CURVES
     - ceaf_prob: arrays for cost-effectiveness acceptability FRONTIERS
 __Section 6:__ VOI/EVPPI
     - partial.evpi, partial.evpi_df: array and data frame to make EVPPI graph
 __Section 7:__ Threshold analysis
(nothing)

---

# Installation to-do list (all free)
- [R](https://www.r-project.org) (required)
- [Packrat](https://rstudio.github.io/packrat/) (package installed within R; required)
- [RStudio](https://www.rstudio.com/) (highly recommended though not required; it integrates well with the other management tools used in the project. The free version is more than enough.)
- Latex engine (not required, though some documentation created automatically won't run). See [here](https://support.rstudio.com/hc/en-us/articles/200532257-Customizing-LaTeX-Options) for the recommended engines.

# Helpful Tutorials

I would recommend the [RStudio Essential webinars](https://www.rstudio.com/resources/webinars/rstudio-essentials-webinar-series-part-1/) if you are completely new to RStudio.

# Documentation with Rmarkdown, Bookdown, and Latex

Some of the documentation for the project is produced automatically with the code in the project (or it will be, once I get to do some more things). This is done via Rmarkdown (using bookdown) and Latex. See the help links later in this document for more information, but suffice it to say here that you may have to install an Latex engine (preferably pdfLaTeX and XeLaTeX). 

The clinical and economic model is a probability tree model, diagrammed in using Latex. This model takes the output of the dynamic model developed and operated by the Warwick team and projects the clinical outcomes and the accompanying costs of treating patients. In addition, prevention campaigns (active screening and vector control) are modeled and integrated into the analysis (to be diagrammed).

- Rmarkdown [website](https://rmarkdown.rstudio.com) and [documentation](https://bookdown.org/yihui/rmarkdown/).
- Bookdown [website](https://bookdown.org) and [documentation](https://bookdown.org/yihui/bookdown/).
- Latex [website](https://www.latex-project.org) and [documentation](https://en.wikibooks.org/wiki/LaTeX/Basics).
- Latex in RStudio [website](https://support.rstudio.com/hc/en-us/articles/200532257-Customizing-LaTeX-Options).

Much of the documentation of this project (in particular the documentation that will result from a query to the parameter database) will be done via rmarkdown (using the bookdown package, to be precise) and rendered in PDF via latex, and often, rendered as an html as well. Rmarkdown and bookdown are available on the private library of the project, so you have nothing to do by hand to set up Rmarkdown

**However, you have to install a latex engine separately.** Since we are using RStudio I would go for pdfLaTeX and XeLaTeX; see more information [here](https://support.rstudio.com/hc/en-us/articles/200532257-Customizing-LaTeX-Options).

**Keep in mind that if you do not want to install latex you should be able to run all of the code in this project, just some of the documentation will not compile.** (Note on 17 July 2020: This should not be an issue at the moment).

---

# Troubleshooting

**These are related to packrat, which is currently disabled for this project until further notice.**

1. Once packrat is set for a project, it usually always checks the package library against the same CRAN mirror that was indicated when the library was created. To figure out that what repo are we using, use the R command in `getOption("repos")`. I once had an issue with a package that seemed not to be installing from the ETH-Zurich repository, so I had to change the mirror that the project tapped into. To do this: go to ./packrat/packrat.lock and change the line (probably the 4th) that says 'Repos: CRAN='. Here are all the [CRAN mirrors](https://cran.r-project.org/mirrors.html) and [their status](https://cran.r-project.org/mirmon_report.html). Copy paste the address of the CRAN mirror to the packrat.lock file. 
% note to self: this issue talks about the crudeness of this solution https://github.com/rstudio/packrat/pull/429
2. Adding a package to the Packrat library - just `install.packages()` it. Note that it won't install for other projects
3. Packrat just installs for everything, even half-broken code within the project folder, so be careful what you set as the target for packrat (the hatmepp_cea directory was too broad for packrat)

---

# TODO (for Marina to do in a future iteration of this code)
1. Create a directory tree in Latex: http://tug.ctan.org/macros/generic/dirtree/dirtree.pdf
2. Finish the two sections in this document related to steps to sync with bitbucket repo (under construction).
3. Re-apply pack rat for the project.
4. Have some steps that are error checking or sanity checks.
5. Change code that uses reshape2 or plyr to use tidyr and dplyr.
6. Set some outputs up to streamline into an rmarkdown file.
7. Make into a tutorial.

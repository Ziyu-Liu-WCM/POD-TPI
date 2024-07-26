## PoD-TPI: Probability-of-Decision Toxicity Probability Interval Design to Accelerate Phase I Trials

## Introduction

We have proposed the PoD-TPI design to accelerate phase I trials and have developed a statistical methodology to calculate the posterior distribution of a dose assignment decision in the presence of pending toxicity outcomes. The posterior distribution directly reflects the confidence of all possible decisions, and the optimal decision is computed under a decision-theoretic framework. The PoD-TPI design is built upon the mTPI-2 design. Nevertheless, the proposed strategy can be applied to other model-free or model-assisted designs such as 3+3, BOIN, keyboard or i3+3.


## User-provided inputs

In dev


## Output
Dose decision will be decided from four possible choices(Suspend, De-escalation, Stay, Escalation). The recommended dose assignment for the current patient will also be given.


## Algortihm of POD-TPI
![](Algorithm.png)



## Package Requirement
Run the following code in R:

```
install_missing_packages <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

packages_to_check <- c("shiny", "dfcrm", "BOIN", "truncdist", 
                        "MCMCpack", "poisbinom", "zoo")

install_missing_packages(packages_to_check)
```

## Run POD-TPI

```
library(shiny)
runGitHub("POD-TPI", "Ziyu-Liu-WCM")
```
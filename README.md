## PoD-TPI: Probability-of-Decision Toxicity Probability Interval Design to Accelerate Phase I Trials

### Introduction

We have proposed the PoD-TPI design to accelerate phase I trials and have developed a statistical methodology to calculate the posterior distribution of a dose assignment decision in the presence of pending toxicity outcomes. The posterior distribution directly reflects the confidence of all possible decisions, and the optimal decision is computed under a decision-theoretic framework. The PoD-TPI design is built upon the mTPI-2 design. Nevertheless, the proposed strategy can be applied to other model-free or model-assisted designs such as 3+3, BOIN, keyboard or i3+3.


### Example Dataset 
We will use the dataset in the Fig.1(a) of the PoD-TPI paper. The input dataset 'df' should include at least three columns: 'dose', 'event_time' and 'event'. Suppose that you are assigning a dose level to the nth patient of the trial, the 'dose' column should contain the dose level history of all past (n-1) patients and the 'event_time' column should contain the time to event(DLT here) observed for all past (n-1) patients, and the 'event' column should contain whether event is observed on each patient.

There are several restrictions here. First, The input dataset should be arranged according to the order of patient recruitment. Patients recruited later need to be placed later in the dataset. Second, the 'event_time' variable must be less than or equal to the DLT assessment window. If DLT is not observed on a patient within the assessment window, the event_time should be equal to assessment window. Third, for patients with pending outcomes, their 'event' variable should be NA and their 'event_time' variables are the pending follow-up times.

```
df <- data.frame(dose = c(1, 2, 2, 2, 2, 2),
                 event_time = c(28, 28, 9, 26, 15, 8),
                 event = c(0, 0, 1, 1, NA, NA))
```


### Other Parameters
**current_dose:** current dose assignment for the incoming patient.

**W:** DLT assessment window.

**MaxDose:** max dose level of the trial.

**p_T:** target DLT rate.

**epsilon:** 2 dimensional vector, the bounds of the equivalence interval.

**niter:** number of iterations for MCMC sampling, greater the better, > 1000 recommended.

**type_dose_decision:** Can be "mTPI" or "i3+3", the base complete data dose decision


**type_suspension:** Three types of suspension rules.

= 0: Never suspends

= 1: Suspend if r_d > (n_d + m_d + r_d) / 2 

= 2: Suspend based on probability thresholds

= 3: Suspend if Pr(decision) < 1, i.e., "look-ahead" design

Suppose n_d is the number of DLT patients observed and m_d is the number of non-DLT patients observed and r_d is the number of pending patients

(1) if no observed outcome at current dose (n_d + m_d = 0), always stay.
    when suspension type is 2, if r_d >= 3 and n_d + m_d = 0, also suspend
    
(2) if m_d = 0, no non-DLTs, for suspension rules 0 or 1, if decision is E, change to S.
    for suspension rule 2, if decision is E, suspend

**type_safety:** Three types of safety rules.

= 0: No safety rule

= 1: Only early termination

= 2: Both safety rules


**q1:** pi_E for Trial Suspension, must be within [0.33, 1], > 0.8 recommended. A larger pi_E represents more conservative dose escalations.

**q2:** pi_D for Trial Suspension, must be within [0, 0.5], <= 0.25 recommended. A smaller pi_D represents more conservative stays.


### Run POD-TPI
```
POD_TPI_df(df, current_dose = 2, W = 28, MaxDose = 4,
           p_T = 0.3, epsilon = c(0.05, 0.05), 
           niter = 1000,
           type_dose_decision = "mTPI",
           type_p_prior = "uniform",
           type_t_model = "pwuniform",
           type_suspension = 2,
           type_safety = 2,
           q1 = 1, q2 = 0.15)

## [1] "Dose decision: De-escalation. Dose assignment for this patient: 1"
```

### Output
Dose decision will be decided from four possible choices(Suspend, De-escalation, Stay, Escalation). The recommended dose assignment for the current patient will also be given.


### Algortihm of POD-TPI
![](Algorithm.png)



### Package Requirement
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

### POD-TPI R Shiny App

```
library(shiny)
runGitHub("POD-TPI", "Ziyu-Liu-WCM")
```
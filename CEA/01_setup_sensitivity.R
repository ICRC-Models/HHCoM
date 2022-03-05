# Setup Sensitivity
# Update: July 16, 2021
# MSahu

## SENSITIVITY ANALYSIS SETUP (turn off) -------------------------------------------------------------

# Primary settings

vs_scalar = T # include the scalar to get from ART + VS to ART
vmmc_cost = F # do not include VMMC costs


# Comparison scenario - 2 or 2a? (in 05_icer_bia Script)

alt_scen_2 = T
alt_scen_2a = F


# Discount rate (in 05_icer_bia Script)

sn_dr_main = T
sn_dr_lo = F
sn_dr_hi = F 

bia = F

# Percentage on cb-ART (for Scenario 2 -- lines 425-447 of 03_costs Script)

sn_pct_cbART_main = T
sn_pct_cbART_lo = F
sn_pct_cbART_hi = F

# Maximum Allowable Cost

max_allowable_cost = "OFF"

# ---------------------------------------------------------------------------------------------------


# Set up DF with list of sensitivitity analyses

lengthSN = 4 # Number of one-way sensitivity analyses

snDF <- data.frame(sn_type = rep("Cost", lengthSN),
                   sn_name = c("Home testing", "Hospitalization", "Community ART", "Clinic ART"),
                   sn_no = c(1:lengthSN))

# Bounds

bound <- c("lower", "upper")
snDF <- snDF %>% crossing(bound) %>% mutate(bound_abb = ifelse(bound == "lower", "lb", "ub")) %>% 
  arrange(sn_no, bound)

# Labels

snDF <- snDF %>% 
  mutate(name_abb = paste0( "sn" , sn_no, ".", bound_abb),
         sn_name_full = paste(sn_type, "for", sn_name))

# Set on or off!

snDF$ON_or_OFF = F


# ---------------------------------------------------------------------------------------------------


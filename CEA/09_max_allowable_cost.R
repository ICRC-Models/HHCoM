# Maximum Allowable Cost
# March 4, 2021
# msahu

# Depencies: 03_costs (lines 90+) ; 01_setup sensitivity

# Guess and check..

# ----------------------------------------------------------------------------------


max_allowable_cost = "ON"
MAX_COST = 1616

source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

# Using lower bound

max_allowable_cost = "ON"
MAX_COST = 1604

snDF[7,"ON_or_OFF"] = T
source(paste0(cea_path, "03_costs.R"))  # MUST BE CONNECTED TO VPN, or will get error
source(paste0(cea_path, "05_icer_bia.R"))

# Turn off

max_allowable_cost = "OFF"
snDF[7,"ON_or_OFF"] = F
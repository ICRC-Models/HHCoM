# ICER Table
# March 3, 2021

ICER_TABLE = "ON"

# ===============================================================================

print("Scenario 2 RESULTS")

# Scenario 2, Discounting 0% [costs and health gains]

sn_dr_lo = T
sn_dr_main = F
alt_scen_2 = T
alt_scen_2a = F

source(paste0(cea_path, "05_icer_bia.R"))

# Scenario 2, Discounting 3% [ICERs]

sn_dr_lo = F
sn_dr_main = T

source(paste0(cea_path, "05_icer_bia.R"))

# ===============================================================================

print("Scenario 2a RESULTS")

# Scenario 2a, Discounting 0% [costs and health gains]

sn_dr_lo = T
sn_dr_main = F
alt_scen_2 = F
alt_scen_2a = T

source(paste0(cea_path, "05_icer_bia.R"))

# Scenario 2a, Discounting 3% [ICERs]

sn_dr_lo = F
sn_dr_main = T

source(paste0(cea_path, "05_icer_bia.R"))

# Clean up

sn_dr_main = T
sn_dr_lo = F
alt_scen_2 = T
alt_scen_2a = F

ICER_TABLE = "OFF"
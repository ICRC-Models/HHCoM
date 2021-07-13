# Calculate population sizes for on ART - each arm

#------------------------------------------------------------------------------------

# Setup params

version <- "DoArtOutputs/"
names(version) <- "main"

discount_rate <- 0
names(discount_rate) <- "dr0"

yrs <- c(1,6, 11) # year = 2020, 2025, 2030 


#-----------------------------------------------------------------------------------

# WITH SCALAR ADJUSTMENT

vs_scalar = "on"

## 

scenarios <- "Scenario1"
source(paste0(cea_path, "05_costs.R"), 345:422)


fpop_soc = (Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_soc   
mpop_soc = (Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art)) / MpctVS_soc 

fpop_cb_art = (Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_cbArt 
mpop_cb_art =  (Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art)) / MpctVS_cbArt

mpop_soc[yrs, 1]
fpop_soc[yrs, 1] 
mpop_cb_art[yrs, 1]  
fpop_cb_art[yrs, 1]  

rm(mpop_soc, fpop_soc, fpop_cb_art, mpop_cb_art)

##

scenarios <- "Scenario2"
source(paste0(cea_path, "05_costs.R"), 345:422)

fpop_soc = (Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_soc   
mpop_soc = (Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art)) / MpctVS_soc 

fpop_cb_art = (Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art)) / FpctVS_cbArt 
mpop_cb_art =  (Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art)) / MpctVS_cbArt

mpop_soc[yrs, 1]
fpop_soc[yrs, 1] 
mpop_cb_art[yrs, 1]  
fpop_cb_art[yrs, 1]  

rm(mpop_soc, fpop_soc, fpop_cb_art, mpop_cb_art)

# -----------------------------------------------------------------------------------------

# WITHOUT SCALAR ADJUSTMENT

vs_scalar = "off"

##

scenarios <- "Scenario1"
source(paste0(cea_path, "05_costs.R"), 345:422)


fpop_soc = Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)   
mpop_soc = Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art) 

fpop_cb_art = Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art) 
mpop_cb_art =  Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art) 

mpop_soc[yrs, 1]
fpop_soc[yrs, 1] 
mpop_cb_art[yrs, 1]  
fpop_cb_art[yrs, 1]  

rm(mpop_soc, fpop_soc, fpop_cb_art, mpop_cb_art)

##

scenarios <- "Scenario2"
source(paste0(cea_path, "05_costs.R"), 345:422)

fpop_soc = Fscalar_soc * (Fpop_prev_art + 0.5 * Fpop_init_art)   
mpop_soc = Mscalar_soc * (Mpop_prev_art + 0.5 * Mpop_init_art) 

fpop_cb_art = Fscalar_cbArt * (Fpop_prev_art + 0.5 * Fpop_init_art) 
mpop_cb_art =  Mscalar_cbArt * (Mpop_prev_art + 0.5 * Mpop_init_art) 

mpop_soc[yrs, 1]
fpop_soc[yrs, 1] 
mpop_cb_art[yrs, 1]  
fpop_cb_art[yrs, 1]  

rm(mpop_soc, fpop_soc, fpop_cb_art, mpop_cb_art)
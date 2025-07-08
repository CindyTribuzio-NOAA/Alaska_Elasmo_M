

# setup ----
libs <- c("tidyverse", "FishLife", 'janitor')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function


# fishlife
edge_names = c( FishBase_and_Morphometrics$tree$tip.label,
                FishBase_and_Morphometrics$tree$node.label[-1] ) # Removing root

which_g = match( "Squalus suckleyi", edge_names )
Table2023 = cbind( 
  Mean = FishBase_and_Morphometrics$beta_gv[which_g,],
  SE = sqrt(diag(FishBase_and_Morphometrics$Cov_gvv[which_g,,]))
)

SD_fishlife <- as.data.frame(Table2023) %>% 
  rownames_to_column(var = 'parameter') %>% 
  clean_names() %>% 
  mutate(param_val = exp(mean))

which_g = match( "Lamna ditropis", edge_names )
Table2023 = cbind( 
  Mean = FishBase_and_Morphometrics$beta_gv[which_g,],
  SE = sqrt(diag(FishBase_and_Morphometrics$Cov_gvv[which_g,,]))
)

SS_fishlife <- as.data.frame(Table2023) %>% 
  rownames_to_column(var = 'parameter') %>% 
  clean_names() %>% 
  mutate(param_val = exp(mean))

which_g = match( "Somniosus pacificus", edge_names )
Table2023 = cbind( 
  Mean = FishBase_and_Morphometrics$beta_gv[which_g,],
  SE = sqrt(diag(FishBase_and_Morphometrics$Cov_gvv[which_g,,]))
)
PSS_fishlife <- as.data.frame(Table2023) %>% 
  rownames_to_column(var = 'parameter') %>% 
  clean_names() %>% 
  mutate(param_val = exp(mean))

# trying to find correct sci name for skate species
which_g = match( "Bathyraja", edge_names )
Table2023 = cbind( 
  Mean = FishBase_and_Morphometrics$beta_gv[which_g,],
  SE = sqrt(diag(FishBase_and_Morphometrics$Cov_gvv[which_g,,]))
)
PSS_fishlife <- as.data.frame(Table2023) %>% 
  rownames_to_column(var = 'parameter') %>% 
  clean_names() %>% 
  mutate(param_val = exp(mean))

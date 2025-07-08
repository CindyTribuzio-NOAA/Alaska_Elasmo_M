# M analysis
# cindy.tribuzio@noaa.gov
# June 2025

# 1) pull in life history parameters from google sheets
# 2) append tmax estimates obtained from norpac and racebase specimen data
# 3) use selected empirical estimators for M
# 4) format output, save

# setup ----
libs <- c("tidyverse", "googlesheets4", "writexl", "FishLife", 'janitor')
if(length(libs[which(libs %in% rownames(installed.packages()) == FALSE )]) > 0) {
  install.packages(libs[which(libs %in% rownames(installed.packages()) == FALSE)])}
lapply(libs, library, character.only = TRUE)
'%nin%'<-Negate('%in%') #this is a handy function

out_path <- paste0("results/skates")
dir.create(out_path)
fig_path <- paste0(out_path, '/figures')
dir.create(fig_path)

source('R/M_fxns.R')

# data ----

# see elasmo_M_updates google sheet in SAFE>Natural_Mortality google folder
lh <- read_sheet("https://docs.google.com/spreadsheets/d/1cj5mtrbQ4w2GCJTK-mRY-rKKBnhm2e1yedlw9KMfXuo/edit?gid=0#gid=0",
           sheet = 'life_history_parameters') %>% 
  filter(species %nin% c('spiny_dogfish', 'blue_shark', 'pacific_sleeper_shark', 'salmon_shark'))
lh$species <- sub('_', ' ', lh$species)

nrow(lh)
lh
lh$use

# tmax estimates based on combined fishery and survey age data. see tmax.R
tmax <- read_csv(paste0(getwd(), '/results/tmax_summary.csv')) %>% 
  filter(Species %nin% c('spiny dogfish', 'blue shark', 'pacific sleeper shark', 'salmon shark')) %>% 
  select(common_name = `Species`, area = Region, 
       max_age = `Maximum age observed`, 
       mean_top5 = `Mean top 5 ages`)

# test to see if the spreadsheet has any errors identifying which region the
# data represent. Should be zero rows, and if it's not GOA:WC should be 0 and
# multi_region should be 1
lh %>% 
  mutate(multi_region_test = NWP+GOA+BS+AI+BC+WC+multi_region) %>% 
  filter(multi_region_test > 1)

# bind tmax estimates from AFSC fishery and survey data to user-defined life
# history parameter data
lhfull <- tmax %>% 
  select(species = common_name, data_area = area, 
         # afsc_allages_q99 = q99, # decided as a group not to use this one
         `AFSC max age` = max_age, `AFSC mean top 5` = mean_top5) %>% 
  pivot_longer(cols = -c(species, data_area), names_to = 'method_or_source', values_to = 'lh_param_estimate') %>% 
  mutate(lh_param_type = 'maxage_yr',
         use = 1,
         NWP = 0,
         GOA = ifelse(data_area == 'GOA', 1, 0),
         BS = ifelse(data_area == 'BS', 1, 0),
         AI = ifelse(data_area == 'AI', 1, 0),
         BC = 0,
         WC = ifelse(data_area == 'WOC', 1, 0),
         multi_region = 0) %>% 
  bind_rows(lh %>% 
              select(-notes)) %>% 
  rename(lh_param = lh_param_type,
         estimate = lh_param_estimate)

unique(lhfull$lh_param)

# omit any data that aren't going to get used only keep parameter inputs that
# are relevant to the four selected methods (see M_fxns.R)
input_data <- lhfull %>%
  filter(use == 1) %>% 
#  filter(lh_param %in% c('maxage_yr',
#                         'gsi',
#                         'dry_weight_g',
#                         'temp_C',
#                         'vbgf_k_cm-1',
#                         'vbgf_linf_cm')) %>% 
  select(-use) 

#input_data %>% filter(BC == 1 | WC == 1)

# create new column for input_area, get rid of 0 values, create a version number
# for each species, input_area, and lh_param
input_data <- input_data %>% 
  pivot_longer(cols = c('NWP', 'GOA', 'BS', 'AI', 'BC', 'WC', 'multi_region'),
               names_to = 'input_area', values_to = 'use') %>% # print(n=Inf)
  filter(use == 1) %>% 
  select(-use) %>% 
  tidyr::complete(species, input_area, lh_param) %>% 
  group_by(species, input_area, lh_param) %>% 
  mutate(version = row_number()) %>% 
  ungroup()


input_data %>% print(n=Inf)

# Run analysis ----

#setup for fishlife
# uses scientific name as input and doesn't have region specificity
# get scientific names for things
stocks <- read_sheet("https://docs.google.com/spreadsheets/d/1cj5mtrbQ4w2GCJTK-mRY-rKKBnhm2e1yedlw9KMfXuo/edit?gid=0#gid=0",
                     sheet = 'stocks_of_interest') %>% 
  filter(common_name %nin% c('spiny_dogfish', 'blue_shark', 'pacific_sleeper_shark', 'salmon_shark')) %>% 
  rename(species = common_name)
stocks$species <- sub('_', ' ', stocks$species)

#setup fishbase info
edge_names = c( FishBase_and_Morphometrics$tree$tip.label,
                FishBase_and_Morphometrics$tree$node.label[-1] )

# calculate four M estimators for each input species and area combo

species_ls <- unique(input_data$species)
area_ls <- unique(input_data$input_area)

# save inputs so they can be version controlled as we refine our spreadsheet
write_csv(input_data, paste0(getwd(), '/data/M_analysis_input_data_skates.csv'))

for(i in 1:length(species_ls)) {
  for(j in 1:length(area_ls)) {
    
    # i = 1; j = 1; k = 1
    df <- input_data %>% 
      filter(species == species_ls[i] & input_area == area_ls[j]) 
    
    # fishlife estimator ----
    fl_out <- data.frame(method = 'fishlife',
                         version = 1,
                         Mest = calcM_fishlife(i = i),
                         method_or_source = 'Thorson')
    
    # max age estimator ----
    tmpdf <- df %>% filter(lh_param == 'maxage_yr')
    
    # occasionally there are more than one version of the area and species
    # specific life history parameters. loop through these
    amax_out <- data.frame(method = 'amax',
                           version = unique(tmpdf$version),
                           Mest = NA,
                           method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      amax_out[k, 'Mest'] <- calcM_amax(input_amax = v_tmpdf$estimate)
    }
    
    # gsi estimator ----
    tmpdf <- df %>% filter(lh_param == 'gsi')
    
    gsi_out <- data.frame(method = 'gsi',
                          version = unique(tmpdf$version),
                          Mest = NA,
                          method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      gsi_out[k, 'Mest'] <- calcM_gsi(input_gsi = v_tmpdf$estimate)
    }
    
    # temp/dry mass estimator ----
    #tmpdf <- df %>% filter(lh_param %in% c('dry_weight_g', 'temp_C'))
    
    #masstemp_out <- data.frame(method = 'temp',
    #                           version = unique(tmpdf$version),
    #                           Mest = NA,
    #                           method_or_source = unique(tmpdf$method_or_source))
    
    #for(k in 1:length(unique(tmpdf$version))) {
    #  v_tmpdf <- tmpdf %>% filter(version == k)
    #  masstemp_out[k, 'Mest'] <- calcM_mass_temp(input_dry_mass = v_tmpdf$estimate[v_tmpdf$lh_param == 'dry_weight_g'],
    #                                             input_temp = v_tmpdf$estimate[v_tmpdf$lh_param == 'temp_C'])
    #}
    
    # then linf/k estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('vbgf_k_cm-1', 'vbgf_linf_cm'))
    
    then_out <- data.frame(method = 'then',
                          version = unique(tmpdf$version),
                          Mest = NA,
                          method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      then_out[k, 'Mest'] <- calcM_lvb(input_linf = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_linf_cm'],
                                      input_k = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_k_cm-1'])
    }
    
    # hamel_k estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('vbgf_k_cm-1'))
    
    hamelk_out <- data.frame(method = 'hamel_k',
                           version = unique(tmpdf$version),
                           Mest = NA,
                           method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      hamelk_out[k, 'Mest'] <- calcM_hamelk(input_k = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_k_cm-1'])
    }
    
    # Jensen_k1 estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('vbgf_k_cm-1'))
    
    jensenk1_out <- data.frame(method = 'jensen_k1',
                             version = unique(tmpdf$version),
                             Mest = NA,
                             method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      jensenk1_out[k, 'Mest'] <- calcM_jensenk1(input_k = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_k_cm-1'])
    }
    
    # Jensen_k2 estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('vbgf_k_cm-1'))
    
    jensenk2_out <- data.frame(method = 'jensen_k2',
                               version = unique(tmpdf$version),
                               Mest = NA,
                               method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      jensenk2_out[k, 'Mest'] <- calcM_jensenk2(input_k = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_k_cm-1'])
    }
    
    # Frisk_k estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('vbgf_k_cm-1'))
    
    friskk_out <- data.frame(method = 'frisk_k',
                               version = unique(tmpdf$version),
                               Mest = NA,
                               method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      friskk_out[k, 'Mest'] <- calcM_friskk(input_k = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_k_cm-1'])
    }

    # Frisk_tmat estimator ----
    tmpdf <- df %>% filter(lh_param %in% c('agemat_yr'))
    
    frisktmat_out <- data.frame(method = 'frisk_tmat',
                             version = unique(tmpdf$version),
                             Mest = NA,
                             method_or_source = unique(tmpdf$method_or_source))
    
    for(k in 1:length(unique(tmpdf$version))) {
      v_tmpdf <- tmpdf %>% filter(version == k)
      frisktmat_out[k, 'Mest'] <- calcM_frisktmat(input_tmat = v_tmpdf$estimate[v_tmpdf$lh_param == 'agemat_yr'])
    } 
    # Hisano_t0 estimator ----
    #tmpdf <- df %>% filter(lh_param %in% c('agemat_yr', 'vbgf_t0_yr'))
    
    #hisanotmat_out <- data.frame(method = 'hisano_tmat',
    #                            version = unique(tmpdf$version),
    #                            Mest = NA,
    #                            method_or_source = unique(tmpdf$method_or_source))
    
    #for(k in 1:length(unique(tmpdf$version))) {
    #  v_tmpdf <- tmpdf %>% filter(version == k)
    #  hisanotmat_out[k, 'Mest'] <- calcM_hisanotmat(input_tmat = v_tmpdf$estimate[v_tmpdf$lh_param == 'agemat_yr'],
    #                                              input_t0 = v_tmpdf$estimate[v_tmpdf$lh_param == 'vbgf_t0_yr'])
    #} 
    
       # Bind and output data
    out <- bind_rows(fl_out, amax_out, gsi_out, then_out, hamelk_out, jensenk1_out, jensenk2_out,
                     friskk_out, frisktmat_out) %>% #, hisanotmat_out taken out for now
      mutate(species = species_ls[i],
             input_area = area_ls[j]) %>% 
      select(species, input_area, M_method = method, M_estimate = Mest, 
             version, method_or_source)
    
    if(i == 1 & j == 1) {
      fullout <- out 
    } else {
      fullout <- bind_rows(fullout, out)
    }
  }
}

# M results ----
nrow(fullout)
fullout %>% select(-method_or_source)
write_csv(fullout, paste0(out_path, "/l_M_estimates_skates.csv"))
l_fullout <- fullout

fullout <- l_fullout %>% 
  select(-method_or_source) %>% 
  filter(!is.na(M_estimate)) %>% 
  pivot_wider(id_cols = c(species, version, M_method),
              names_from = input_area,
              values_from = M_estimate) %>% 
  select(species, version, M_method, GOA, BS, BC, WC, multi_region) %>% #note: no data from region = AI
  arrange(species, M_method, version) 

fullout #%>% View()

write_csv(as.data.frame(fullout), paste0(out_path, "/M_estimates_skates.csv"))

# join inputs to M outputs
l_fullout %>% as_tibble()
unique(l_fullout$M_method)
unique(input_data$lh_param)

summ <- l_fullout %>% 
  left_join(input_data, by = c("species", "input_area", "version", "method_or_source")) %>% 
  filter(!is.na(M_estimate)) %>% 
  mutate(keep = ifelse(M_method == 'amax' & lh_param == 'maxage_yr', 1, 
                       ifelse(M_method %in% c('hamel_k', 'jensen_k1', 'jensen_k2', 'frisk_k') 
                              & lh_param == 'vbgf_k_cm-1', 1, 
                              ifelse(M_method == 'gsi' & lh_param == 'gsi', 1, 
                                     ifelse(M_method == 'frisk_tmat' & lh_param == 'agemat_yr', 1,
                                            ifelse(M_method == 'then' & lh_param %in% c('vbgf_k_cm-1', 'vbgf_linf_cm'), 1, 
                                                   ifelse(M_method == 'fishlife', 1, 0))))))) %>% 
  filter(keep == 1)


#summ <- input_data %>% 
#  mutate(M_method = case_when(lh_param %in% c('maxage_yr') ~ 'amax',
#                              lh_param %in% c('gsi') ~ 'gsi',
#                              lh_param %in% c('vbgf_linf_cm') ~ 'lvb', #'vbgf_k_cm-1', 
#                              lh_param %in% c('agemat_yr') ~ 'amat'),)

########
##############
###########
#add fishlife to this
summ <- summ %>% 
  filter(M_method %in% c('gsi', 'amax')) %>% 
  mutate(lh_param = ifelse(lh_param == 'gsi', 'GSI', 'Max age (y)')) %>% 
  rename(input_names = lh_param, input_values = estimate) %>%
  mutate(input_values = ifelse(input_names == 'GSI',
                               formatC(round(input_values, 4), format = 'f', digits = 4),
                               ifelse(grepl('mean top five', method_or_source),
                                      formatC(round(input_values, 1), format = 'f', digits = 1),
                                      formatC(round(input_values, 0), format = 'f', digits = 0)))) %>% 
  bind_rows(summ %>% 
              filter(M_method %in% c('then')) %>% 
              mutate(lh_param = ifelse(lh_param == 'vbgf_k_cm-1', 'k', 'linf')) %>% 
              pivot_wider(id_cols = c('species', 'input_area', 'data_area', 'version', 'M_method', 'method_or_source', 'M_estimate'),
                          names_from = lh_param, values_from = estimate) %>% 
              mutate(input_names = 'VBGF Linf (cm) / k (cm-1)',
                     input_values = ifelse(is.na(k), NA,
                                           paste0(formatC(round(linf, 1), digits = 1, format = 'f'), 
                                                  " / ", 
                                                  formatC(round(k, 3), format = 'f', digits = 3)))) %>% 
              select(-k, -linf)) %>% 
  bind_rows(summ %>% 
              filter(M_method %in% c('frisk_tmat')) %>% 
              #mutate(lh_param = 'amat') %>% 
              rename(input_names = lh_param, input_values = estimate) %>%
              mutate(input_names = 'Age at Maturity (yr)',
                     input_values = ifelse(is.na(input_values), NA,
                                           paste0(formatC(round(input_values, 1), format = 'f', digits = 1))))) %>%
  bind_rows(summ %>% 
              filter(M_method %in% c('fishlife')) %>% 
              #mutate(lh_param = 'amat') %>% 
              rename(input_names = lh_param, input_values = estimate) %>%
              mutate(input_names = 'species name',
                     data_area = input_area,
                     input_values = ifelse(is.na(input_values), NA,
                                           paste0(formatC(round(input_values, 1), format = 'f', digits = 1))))) %>%
  bind_rows(summ %>% 
              filter(M_method %in% c('hamel_k', 'jensen_k1', 'jensen_k2', 'frisk_k', 'zhou_k')) %>% 
              rename(input_names = lh_param, input_values = estimate) %>%
              mutate(input_names = 'VBGF k (cm-1)',
                     input_values = ifelse(is.na(input_values), NA,
                                           paste0(formatC(round(input_values, 3), format = 'f', digits = 3))))) %>% 
  mutate(M_estimate = formatC(round(M_estimate, 3), format = 'f', digits = 3),
         M_method = paste0(M_method, '.v', version),
         area = ifelse(input_area == 'multi_region', data_area, input_area)) %>% 
  select(species, area, M_method, data_input = input_names, 
         data_input_values = input_values, M_estimate,
         references = method_or_source) 

summ




#summ <- summ %>% 
#  mutate(species = factor(species, 
#                          labels = c('dusky rockfish', 'harlequin rockfish',
#                                     'rebs rockfish','rougheye rockfish','blackspotted rockfish',
#                                     'redbanded rockfish', 'redstripe rockfish', 
#                                     'sharpchin rockfish', 'shortraker rockfish','silvergray rockfish',
#                                     'yelloweye rockfish', 'shortspine thornyhead'),
#                          levels = c('dusky rockfish', 'harlequin rockfish',
#                                     'rebs rockfish','rougheye rockfish','blackspotted rockfish',
#                                     'redbanded rockfish', 'redstripe rockfish', 
#                                     'sharpchin rockfish', 'shortraker rockfish','silvergray rockfish',
#                                     'yelloweye rockfish', 'shortspine thornyhead'),
#                          ordered = TRUE)) %>% 
#  arrange(species, area, data_input, version)
#summ
summ %>% write_csv(paste0(out_path, '/formatted_M_results_skates.csv'))

# write results to separate excel sheets using writexl and purrr ----
summls <- summ %>% 
  dplyr::group_split(species)

# these become the names of each sheet in excel
names(summls) <- summls %>% 
  purrr::map(~pull(., species)) %>% 
  map(~as.character(.)) %>% 
  map(~unique(.))

summls %>% 
  writexl::write_xlsx(path = paste0(out_path, '/M_species_tables_skates.xlsx'))

reference_lkup <- summ %>% 
  distinct(references) %>% 
  mutate(reference_number = row_number()) 
# 
reference_lkup %>% 
  mutate(txt = paste0(reference_number, ') ', references))

summ %>% 
  left_join(reference_lkup) %>% 
  select(Species = species, Region = area, `Parameter(s)` = data_input,
         `Parameter values(s)` = data_input_values, `M estimate` = M_estimate, Reference = references) %>% 
  write_xlsx(path = paste0(out_path, '/M_species_tables_full_skates.xlsx'))



# Figures ----

l_fullout <- l_fullout %>% 
  select(-method_or_source) %>% 
  mutate(M_method2 = paste0(M_method, '.v', version)) %>% 
  mutate(input_area = factor(input_area,
                             labels = c('GOA', 'BS', 'AI', 'BC', 'WC', 'multi_region'),
                             levels = c('GOA', 'BS', 'AI', 'BC', 'WC', 'multi_region'),
                             ordered = TRUE))

mybarplot <- function(df = plot_data,
                      title = 'your plot title (e.g., species name)') {
  df %>%
    filter(!is.na(M_estimate)) %>% 
    ggplot(aes(x = M_method2, y = M_estimate, fill = M_method, col = M_method)) +
    geom_bar(stat = 'identity') +
    facet_grid(species~input_area) +
    labs(x = NULL, y = 'Natural mortality (M)',
         title = title) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
}

theme_set(theme_bw(base_size = 12))

unq_key <- unique(l_fullout$species)

for(i in 1:length(unq_key)){
  plot_data <- l_fullout %>% 
    filter(species == unq_key[i])
  
  mybarplot(title = unq_key[i])
  ggsave(paste0(fig_path, '/', unq_key[i], '_M_results.png'),
         dpi = 300, units = 'in', width = 8, height = 6)
}


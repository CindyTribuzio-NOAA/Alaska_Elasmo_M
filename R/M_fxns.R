# Formulas for the natural mortality estimators selected for this study.
# See further rationale for the selection of these estimators at
# https://docs.google.com/spreadsheets/d/1cj5mtrbQ4w2GCJTK-mRY-rKKBnhm2e1yedlw9KMfXuo/edit?gid=0#gid=0

# Last updated June 2025

# References ----
# Gunderson, D. R. 1997. Trade-off between reproductive effort and adult
# survival in oviparous and viviparous fishes. Canadian Journal of Fisheries and
# Aquatic Sciences, 54: 990 –998.
# Hamel, O. S. 2015. A method for calculating a meta-analytical prior for the
# natural mortality rate using multiple life history correlates. – ICES Journal
# of Marine Science, 72: 62–69.
# Hamel, O. S. and J. M. Cope. 2022. Development and considerations for application
# of a longevity-based prior for the natural mortality rate. Fisheries Research 256: 106477
# Then, A. Y., Hoenig, J. M., Hall, N. G., and Hewitt, D. A. 2015. Evaluating the
# predictive performance of empirical estimators of natural mortality rate using
# information on over 200 fish species. – ICES Journal of Marine Science, 72:
# 82–92.
# Jensen 1996 (pull code from NMT)
# Zhang and Megrey 2006
# Frisk et al. 2001
# Hisano et al. 2011
# Natural Mortality Tool (NMT) Jason Cope https://github.com/shcaba/Natural-Mortality-Tool

# M estimator functions ----

# Hamel 2022 update to Then et al 2015. Use log-log transformation on the linear
# regression due to underlying heteroscedasticity in the data.
# M_vals[3]=Hamel_amax in Cope's NMT
calcM_amax <- function(input_amax # years
) {
  if(is.na(input_amax)) {
    out <- NA
    } else { 
  out <- 5.4 / input_amax
  }
  return(out)
}


# Gunderson (1997) used linear regression in real space found a coefficent of
# 1.79. Hamel (2015) updated this coefficient to 1.8717, which is used here. This
# updated value is based on a log-log trasformation on the regression
# coefficient, which was most appropriate given heteroscedasticity in the
# underlying data. Since M is positively correlated with GSI, the log-log
# transformation forces the relationship through the origin in real space. The
# regression was forced to have a slope of 1 in log space, which makes them
# linear in real space.
calcM_gsi <- function(input_gsi) {
  if(is.na(input_gsi)) {
    out <- NA
  } else { 
  out <- 1.817 * input_gsi # GSI = weight oavry / somatic weight. GnD_GSI_M in Cope's NMT
  }
  return(out)
}

###### Omitted, not applicable to elasmos
# McCoy and Gillooly 2008 and reformulation in Hamel 2015. 
#calcM_mass_temp <-function(input_dry_mass, # grams
#                           input_temp # Celsius
#) {
  # see equations 1 and 2 in McCoy and Gillooly 2008  and equation 17 in Hamel
  # 2015 for definitions. The one difference is that the value of -7540 is
  # expressed as -eV/K in equation 17 in Hamel 2015, where K=Boltzmann's
  # constant=8.62e-5, and
  # eV=0.65.
#  if(is.na(input_dry_mass) | is.na(input_temp)) {
#    out <- NA
#  } else { 
#  out <- 3.2 * ((input_dry_mass / 4)^ -0.27) * exp(-7540 * ((1 / (273 + input_temp))- (1 / 293.15)))
#  }
#  return(out)
#}

# Then et al 2015 revision of the Pauly (1980) estimator that omits
# temperature (see Table 1 in Then et al 2015, M=ak^(b)*Linf^(c)). Then et al
# 2015 estimated the coefficients to be equal to
# Mest=4.111*K^(0.73)*Linf^(-0.333) with length in cm using non-linear least
# squares regression methods. Recommendation to use this estimator IF RELIABLE
# ESTIMATES OF A_MAX WERE NOT AVAILABLE. However, Hamel and Cope (2023) argue that 
# Then et al. methods are incorrect because of log-space transformations.
calcM_lvb <-function(input_linf, # cm
                 input_k # cm^-1
) {
  if(is.na(input_linf) | is.na(input_k)) {
    out <- NA
  } else { 
  out <- 4.11 * input_k ^ 0.73 * input_linf ^ -0.33
  }
  return(out)
}

# Zhang and Megrey (2006). This estimator is a modification of the Alverson and
# Carney 1975 method and is reformulated to include beta and t0 parameters.
# Emailed Jason C 2022-03-14 about whether or not its ok to use length-based
# VBGF k and t0 parameters instead of weight-based VBGF parameters. Jane is
# concerned about using this estimator because of the Then et al 2015 recommendation
# to not use the Alverson-Carney estimator because it has no additional
# advantage over the estimators that rely solely on Amax. Perhaps this
# recommendation is not relevant since Zhang and Megrey (2006) revised Alverson
# and Carney 1975. Also in email correspondence with Jason Cope 2022-03-14, he
# said that this estimator routinely performs poorly

# M_ZM_AC<-function(k,Amax,t0) {
#   Mvals.out<-c(NA,NA)
#  Mvals.out[1]<- (3*k)/(exp(k*((0.302 * Amax)-t0))-1) # pelagic - not relevant for our study
#   Mvals.out[2]<- (3*k)/(exp(k*((0.44 * Amax)-t0))-1) # demersal
#   return(Mvals.out) 
# }

 ######
 #add in new functions below here

 
# Hamel_k (Hamel 2015). This is the log-transofrmed approach to Jensen_k2,
# thus an alternative to the 1.6 Jensen value. Hamel_M_k and Jensen_k1, k2 are
# are basically sensitivities of the same method.
 
 calcM_hamelk<-function(input_k) {
   if(is.na(input_k)) {
     out <- NA
   } else {   
   out <- input_k*1.753
   }
   return(out)
 }
 
# Jensen_k1, Jensen_k2 (Jensen 1996). Relies on the M/k ratio = 1.5. Derived value.
# The 1.6 is the regression coefficient from the linear regression of M to k data
calcM_jensenk1<-function(input_k) {
  if(is.na(input_k)) {
    out <- NA
  } else {   
    out <- input_k*1.5
  }
  return(out)
}
calcM_jensenk2<-function(input_k) {
  if(is.na(input_k)) {
    out <- NA
  } else {   
    out <- input_k*1.6
  }
  return(out)
}

# Frisk_k (2001). Demonstrated that elasmobranchs as a group
# violate the M/k ratio assumptions, however, the data in Fig 3
# suggests that Rajidae do not. Skates may have a "faster"
# life history? The k model is done in log space, so fits 
# the Hamel arguments without concern. Note that the original 
# equation was a linear regression, but Zhou et al. put it in 
# log space, with no change in the constants.

calcM_friskk <- function(input_k) {
  if(is.na(input_k)) {
    out <- NA
  } else {   
    out <- 0.436*(input_k^0.42)
  }
  return(out)
}

# Fristk_tmat. Avoids the M/k assumptions, assumes 
 # M=1/(0.44tmat+1.87)
# FEMALE ONLY
calcM_frisktmat <- function(input_tmat) {
  if(is.na(input_tmat)) {
    out <- NA
  } else {   
    out <- 1/(0.44*input_tmat + 1.87)
  }
  return(out)
}

# Hisano_tmat Hisano et al. 2011. Argued that Jensen ignored t0
# which is OK for teleost, but not for elasmos, which are not
# length = at 'birth'. Supercedes Jensen_tmat equation.
# this one is problematic due to the input data.mismatches between
# sources that supply information and/or t0 being kinda meaningless
# for sharks

calcM_hisanotmat <- function(input_tmat,
                             input_t0) {
  if(is.na(input_tmat) | is.na(input_tmat)) {
    out <- NA
  } else {   
    out <- 1.65/(input_tmat - input_t0)
  }
  return(out)
}

 
 

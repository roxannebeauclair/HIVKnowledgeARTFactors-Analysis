# =================
# Importing dataset
# =================

# Author: Roxanne Beauclair
# Description: This scripts contains functions used for the HERStory 
# HIV knowledge and ART factors analysis


# ================================================
# Functions for loading and detaching dependencies
# ================================================
InstallLoad <- function(package1, ...) {
  # This function will load, and install if necessary, libraries needed for script
  
  # Convert arguments to vector
  pkgs <- c(package1, ...)
  
  # Start loop to determine if each package is installed
  for(p in pkgs) {
    
    # If package is installed locally, load
    if(p %in% rownames(installed.packages()))
      do.call("library", list(p))
    
    # If package is not installed locally, download, then load
    else {
      install.packages(p)
      do.call("library", list(p))
    }
  }
}


# Remove libraries
RemoveLibraries <- function(package1, ...) {
  # This function will remove all libraries that have been loaded
  # Packages must be entered as character strings
  
  # Convert argument into a vector
  pkgs <- c(package1, ...)
  
  #Start loop
  for(p in pkgs) {
    item <- paste("package", p, sep = ":")
    
    #Check to see that package is loaded before detaching
    while(item %in% search()) {
      detach(item, character.only = TRUE)
    }
  }
}

# ==========================================
# Functions for producing summary statistics
# In particular table format
# ==========================================

# This function must be used with the tidyverse and
# svryr packages

# NOTES ON THE PROGRAMMING FUNCTIONS USED BELOW
# enquo() looks at the argument to see what I typed, and returns that
# value as a quosure. A quosure is a special type of formula
# The "!!" unquotes the quosure so it can be evaluated
# quo_name() creates a string out of the quosure


# ==============================================
# Univarite summary statistics for Cascade Paper
# ==============================================

# This function returns a table for a single categorical variable
# The table contains the frequencies (Freq) and proportions (Prop)
# in each category
recap_uni_catvar_svy <- function(svy, var) {
  
  df <- svy
  group_var <- enquo(var)
  group_var_str <- quo_name(group_var)
  
  c_levels <- unique(df$variables[[group_var_str]])
  
  tab_n <- df %>%
    group_by(!! group_var, .drop = FALSE) %>%
    summarise(Freq = unweighted(n()))
  
  tab_prop <- map_df(c_levels, function(cat) {
    df %>%
      summarise(Prop = survey_mean(!! group_var == cat,
                                   na.rm = T,
                                   vartype = "ci",
                                   proportion = TRUE, 
                                   prop_method = "beta")) %>%
      mutate(!! group_var_str := cat) # := is a helper for creating new names
    
  })
  
  tab_prop <- tab_prop %>%
    mutate_at(vars(starts_with("Prop")), 
              list(~format(round(. * 100, 1), nsmall = 1))) %>%
    unite(CI, Prop_low, Prop_upp, sep = " — ") 
  
  tab <- tab_n %>%
    left_join(tab_prop, by = group_var_str) %>%
    data.frame(row.names = 1) %>%
    as.matrix()
  
  
  return(tab)
  
}

# This function returns a table for a single numeric variable
# The table contains the medians (Freq), IQR (Prop), and 95% CI
# for the median. 
# The table has the column names "n", "Prop" and "CI", so that 
# it can easily be appended to univariate summaries of
# categorical variables

recap_uni_numvar_svy <- function(svy, var) {
  
  df <- svy
  numvar <- enquo(var)
  
  tab <- df %>%
    summarise(med = survey_quantile(!! numvar,
                                    na.rm = T,
                                    quantiles = c(0.50, 0.25, 0.75),
                                    vartype = "ci")) %>%
    mutate_at(vars(starts_with("med")), list(~round(., 1))) %>%
    select(-med_q25_low, -med_q25_upp, -med_q75_low, -med_q75_upp) %>%
    unite(Prop, med_q25, med_q75, sep = " — ") %>%
    unite(CI, med_q50_low, med_q50_upp, sep = " — ") %>%
    rename(Freq = med_q50) %>%
    as.matrix()
  
  row.names(tab) <- "Med and IQR"
  
  return(tab)
  
}


# This function returns a table that summarizes variables overall
# and by another variable, using survey data. It also performs
# The chi-square test uses the svychisq() function from the 
# survey package. It computes first and second-order Rao-Scott 
# corrections to the Pearson chisquared test (defaults "F")
# This function returns frequencies (Freq) and proportions (Prop)
# and 95% intervals 
# The percent argument specifies if you want the percents to be as 
# a fraction of the column variable ("col") or row variable ("row")

recap_bi_catvar_svy <- function(svy, var, byvar, percent = "col") {
  
  df <- svy
  rowvar <- enquo(var)
  colvar <- enquo(byvar)
  rowvar_str <- quo_name(rowvar)
  colvar_str <- quo_name(colvar)
  formula <- as.formula(paste("~", rowvar_str, "+", colvar_str))
  
  overall_tab <- df %>%
    group_by(!! rowvar, .drop = FALSE) %>%
    summarise(Freq = unweighted(n()),
              Prop = survey_mean(na.rm = T,
                                 vartype = "ci")) %>%
    mutate_at(vars(starts_with("Prop")), 
              list(~format(round(. * 100, 1), nsmall = 1))) %>%
    unite(CI, Prop_low, Prop_upp, sep = " — ") %>%
    rownames_to_column() # Need to do this so can join p-value table (wrong number of rows)
                         # for bind_cols
  
  if(percent == "col") {
    
    by_tab <- df %>%
      group_by(!! colvar, !! rowvar, .drop = FALSE) %>%
      summarise(Freq = unweighted(n()),
                Prop = survey_mean(na.rm = T,
                                   vartype = "ci")) %>%
      mutate_at(vars(starts_with("Prop")), 
                list(~format(round(. * 100, 1), nsmall = 1))) %>%
      unite(CI, Prop_low, Prop_upp, sep = " — ") %>%
      gather(key = "stat", value = "value", Freq, Prop, CI) %>%
      unite(var_stat, !! colvar, stat, sep = "_") %>%
      spread(key = var_stat, value = value) %>%
      select(!! rowvar, No_Freq, No_Prop, No_CI, Yes_Freq, Yes_Prop, Yes_CI) %>%
      rownames_to_column() 
    
  } else {
    
    by_tab <- df %>%
      group_by(!! rowvar, !! colvar, .drop = FALSE) %>%
      summarise(Freq = unweighted(n()),
                Prop = survey_mean(na.rm = T,
                                   vartype = "ci")) %>%
      mutate_at(vars(starts_with("Prop")), 
                list(~format(round(. * 100, 1), nsmall = 1))) %>%
      unite(CI, Prop_low, Prop_upp, sep = " — ") %>%
      gather(key = "stat", value = "value", Freq, Prop, CI) %>%
      unite(var_stat, !! colvar, stat, sep = "_") %>%
      spread(key = var_stat, value = value) %>%
      select(!! rowvar, No_Freq, No_Prop, No_CI, Yes_Freq, Yes_Prop, Yes_CI) %>%
      rownames_to_column()
  }
  
  
  p <- svychisq(formula, 
                df, 
                na.rm = T, 
                statistic = "adjWald") %>%
    tidy() %>%
    select(p.value) %>%
    mutate(p.value = format(round(p.value, 4), nsmall = 4)) %>%
    rownames_to_column()
  
  tab <- overall_tab %>%
    left_join(by_tab, by = c("rowname", rowvar_str)) %>%
    left_join(p, by = "rowname") %>%
    select(-rowname) %>%
    data.frame(row.names = 1, fix.empty.names = F) %>%
    as.matrix()
  
  return(tab)
  
}


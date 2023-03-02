

# This code was written by Alexander Y. Mitrophanov, PhD. Finalized in February 2023.

# This script generates confidence interval (CI) data.

# !!!!!!!!!! RUN THIS CODE before running f_main_plots.R


library(readxl)


# ------------------------- CODE PARAMS ----------------------------------

# It is assumed that the project scripts, including this one, will be run from the current folder ("curr").
# It is also assumed that the current folder has a subfolder DATA for all the (input and output) data files.
# !!!!! Users: change directory names and locations as needed.
curr <- getwd()
DATA_FOLDER <- paste(curr, "/DATA/", sep = "")

num_reps <- 10000 # The number of noise-perturbed data sets (per model).
set.seed(2)

st_link <- "probit"  #standard link function for this project


# -------------------------- FUNCTIONS ------------------------------------


pl <- function(x) { # A convenient output function.
  print(x, quote = FALSE)
}


glm_boot <- function(dset_predict, dset_labeled, link_func, reg_form) {
  # This function fits the statistical model to one data set and predicts on the other.
  # dset_labeled is for training, dset_predict is for prediction
  # The data sets are assumed to be NOT transformed to type *factor* yet (we do that in this function).
  
  formula_string <- paste("classification ~ ", reg_form, sep = "")
  
  dset_labeled$classification <- factor(dset_labeled$classification)
  fit <- glm(as.formula(formula_string), control = list(maxit = 500), data = dset_labeled, family = binomial(link = link_func))
  PIF_list <- predict(fit, dset_predict, type = "response") # here, not all variants may be labeled
  
  return(PIF_list)
  
}



PIF_noise <- function(data, num_reps, link_func, reg_form) {
  # This function adds Gaussian noise to the variable values in the training set, refits the model,
  # and generates "perturbed" PIFs.
  
  ssdd <- 0.05 # SD for the Gaussian noise (mean is 0 by default)
  
  output_PIFs <- list() # one list element is one PIF set for a perturbed data set (num_reps of them)
  
  data_pathog <- subset(data, classification == 1)
  data_neut <- subset(data,  classification == 0)
  # no need to factor classification here, we do that in glm_boot
  
  data_labeled_sample_1 <- rbind(data_pathog,data_neut)
  NN <- nrow(data_labeled_sample_1)
  
  for (I in 1:num_reps) {
    
    data_labeled_sample <- data_labeled_sample_1
    data_labeled_sample$DMSO <- data_labeled_sample$DMSO + rnorm(NN, sd = ssdd)
    data_labeled_sample$Cis <- data_labeled_sample$Cis + rnorm(NN, sd = ssdd)
    data_labeled_sample$Ola <- data_labeled_sample$Ola + rnorm(NN, sd = ssdd)
    
    output_PIFs[[I]] <- glm_boot(data,data_labeled_sample, link_func, reg_form)
    
  }
  
  return(output_PIFs) 
  
}



quant_intervals <- function(boot_output) { # This is, in fact, "noise output" here, from PIF_noise.
  # This function calculates and returns CI information.
  
  # these upper and lower numbers correspond to 95% CI
  lower_quant <- .025
  upper_quant <- .975
  
  num_variants <- length(boot_output[[1]])
  num_replicates <- length(boot_output)
  
  upper <- 1:num_variants
  lower <- 1:num_variants
  
  for (I in 1:num_variants){
    
    samples <- c()
    for (J in 1:num_replicates) {
      samples[J] <- boot_output[[J]][I]
    }
    
    upper[I] <- quantile(samples, probs = upper_quant)
    lower[I] <- quantile(samples, probs = lower_quant)
    
  }
  
  return(data.frame(lo = lower,up = upper))
  
}





# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# --------------------- MAIN BODY OF THE CODE ---------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------



# -----------------------------------------------------------------------------
# ------------- data loading and processing -----------------------------------
# -----------------------------------------------------------------------------

file_in <- paste(DATA_FOLDER, "FS_599BRCA2_final.xlsx", sep = "")
datafr <- read_excel(file_in)

datafr_rep1 <- datafr[c(1,2,3,5,7)]
datafr_rep2 <- datafr[c(1,2,4,6,8)]


# ----------------- basic stats ----------------------


pl(" ")
pl("Basic stats")
pl("Total variants:")
print(nrow(datafr))
pl("Functional variants:")
print(nrow(subset(datafr, CLASS == "Functional")))
pl("Non-functional variants:")
print(nrow(subset(datafr, CLASS == "Non-functional")))
pl("Total labeled variants:")
print( nrow(subset(datafr, CLASS == "Functional")) + nrow(subset(datafr, CLASS == "Non-functional")) )
pl(" ")


# ---------- averaging over 2 replicates --------------


datafr[,3:5] <- (datafr_rep1[,3:5] + datafr_rep2[,3:5])/2
datafr <- datafr[1:5]



# --------- change names for compatibility------

datafr[datafr == "Non-functional"] <- "Pathogenic"
datafr[datafr == "Functional"] <- "Neutral"
colnames(datafr) <- c("variant", "classification", "DMSO", "Cis", "Ola")


# ----------- Replacing P and N with 1 and 0, respectively -- for classification

classification <- c()
for (I in 1:nrow(datafr)) {
  if (datafr[I,2] == "Pathogenic") (classification[I] <- 1)
  else if (datafr[I,2] == "Neutral") (classification[I] <- 0)
  else (classification[I] <- 0.5)
}

datafr_classif <- data.frame(datafr$variant, classification, datafr$DMSO, datafr$Cis, datafr$Ola)
colnames(datafr_classif) <- c("variant", "classification", "DMSO", "Cis", "Ola")



# ------------------------ CI generation and saving ------------------------------------------


reg_formula <- "DMSO"
noise_PIF_out_file <- "f_noise_PIFs_DMSO.RData"
noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
boot_out <- quant_intervals(noise_out)
file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
save(boot_out, file = file_out)

reg_formula <- "DMSO + Cis + Ola"
noise_PIF_out_file <- "f_noise_PIFs_DMSO_Cis_Ola.RData"
noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
boot_out <- quant_intervals(noise_out)
file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
save(boot_out, file = file_out)

reg_formula <- "Cis"
noise_PIF_out_file <- "f_noise_PIFs_Cis.RData"
noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
boot_out <- quant_intervals(noise_out)
file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
save(boot_out, file = file_out)

reg_formula <- "Ola"
noise_PIF_out_file <- "f_noise_PIFs_Ola.RData"
noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
boot_out <- quant_intervals(noise_out)
file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
save(boot_out, file = file_out)


# # -------------------- additional combinations (2 variables) ---------------------------
# 
# reg_formula <- "DMSO + Cis"
# noise_PIF_out_file <- "f_noise_PIFs_DMSO_Cis.RData"
# noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
# boot_out <- quant_intervals(noise_out)
# file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
# save(boot_out, file = file_out)
# 
# reg_formula <- "DMSO + Ola"
# noise_PIF_out_file <- "f_noise_PIFs_DMSO_Ola.RData"
# noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
# boot_out <- quant_intervals(noise_out)
# file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
# save(boot_out, file = file_out)
# 
# reg_formula <- "Cis + Ola"
# noise_PIF_out_file <- "f_noise_PIFs_Cis_Ola.RData"
# noise_out <- PIF_noise(datafr_classif, num_reps, st_link, reg_formula)
# boot_out <- quant_intervals(noise_out)
# file_out <- paste(DATA_FOLDER, noise_PIF_out_file, sep = "")
# save(boot_out, file = file_out)



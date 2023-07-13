

# Code developed by Alexander Y. Mitrophanov, PhD. Finalized in June 2023.

# This is the main script that calculates and plots PIFs, and does cross-validation.

# !!!!! Before running it, run f_noise.R, which generates and saves CI (confidence interval) information



library(dplyr)
library(readxl)
library(writexl)


# -------------------------------- code parameters ----------------------------

# It is assumed that the project scripts, including this one, will be run from *this* folder ("curr").
# It is also assumed  that this folder has a subfolder DATA for all the (input and output) data files 
# and subfolder FIGURES for all code-generated figures.
# !!! Users: change directory names and locations as needed.
curr <- getwd()
DATA_FOLDER <- paste(curr, "/DATA/", sep = "")
# !!!!!!!!! It is assumed that the main input-data file has the name "FS_NNNSGE_BRCA2.xlsx"
main_input_data_file <- "FS_NNNSGE_BRCA2.xlsx"

FIG_FOLDER <- paste(curr, "/FIGURES/", sep = "")
# This script generates figures as PDF files, and they are saved in this FIGURES folder.


pathog_thresh <- 0.99 # greater than that is Non-functional
neutral_thresh <- 0.05 # less than or equal to that is Functional
threshs <- c(pathog_thresh,neutral_thresh)
st_link <- "probit" # main type of classifier we are using here
set.seed(2) # randomization is needed for K fold selection


# noise PIF file names (hard-coded in the f_noise.R script)
boot_PIF_file_all <- "f_noise_PIFs_DMSO_Cis_Ola.RData"
boot_PIF_file_DMSO <- "f_noise_PIFs_DMSO.RData"
boot_PIF_file_Cis <- "f_noise_PIFs_Cis.RData"
boot_PIF_file_Ola <- "f_noise_PIFs_Ola.RData"



# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# --------------------------------- FUNCTIONS ---------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


pl <- function(x) { # A convenient text-output function.
  print(x, quote = FALSE)
}


set_graphics_params <- function() {
  
  par(cex.main = 1.5, cex.lab = 1.3, cex.axis = 1.2, font.lab = 2)
  par(font = 2)
  
  return(NULL)
}



LOOCV_glm <- function(datafr_classif,thresh_vec, link_func, reg_form, output_flag, out_flag1) {
  # This function performs leave-one-out cross-validation.
  # 
  # reg_form is the regression-formula string in terms of variables
  # output_flag is for extra output of misclassified variant names
  # # output_flag1 is for additional text output
  
  formula_string <- paste("classification ~ ", reg_form, sep = "")

  pathog_thresh <- thresh_vec[1]
  neutral_thresh <- thresh_vec[2]

  correct_positives <- 0
  correct_negatives <- 0
  correct_predictions <- 0
  dset <- subset(datafr_classif, classification == 0 | classification == 1)
  dset$classification <- factor(dset$classification)

  num_variants <- nrow(dset)
  num_pathog <- nrow(subset(dset, classification == 1))
  num_neut   <- nrow(subset(dset, classification == 0))

  for (I in 1:num_variants) {

    test_set <- dset[I,]
    training_set <- dset[-I,]

    fit <- glm(as.formula(formula_string), control = list(maxit = 500), data = training_set, family = binomial(link = link_func))
    ppp <- predict(fit, test_set, type = "response")


    if (ppp > pathog_thresh & test_set$classification == 1) {
      correct_positives <- correct_positives + 1
      correct_predictions <- correct_predictions + 1
    }

    if (ppp <= pathog_thresh & test_set$classification == 1 & output_flag == TRUE) {
      print(test_set[[1]])
      print("Misclassified non-functional")
      pl(" ")
    }

    if (ppp <= neutral_thresh & test_set$classification == 0) {
      correct_negatives <- correct_negatives + 1
      correct_predictions <- correct_predictions + 1
    }
    if (ppp > neutral_thresh & test_set$classification == 0 & output_flag == TRUE) {
      print(test_set[[1]])
      print("Misclassified functional")
      pl(" ")
    }

  }

  sensitivity <- correct_positives/num_pathog*100
  specificity <- correct_negatives/num_neut*100
  predict_accuracy <- correct_predictions/(num_pathog + num_neut)*100
  
  if (out_flag1 == TRUE) {
    
    pl(' ----------------------------------- ')
    pl("LOOCV output")
    pl(" ")
    pl(paste("Regression formula:", reg_form))
    pl(" ")
    pl("Sn, Sp, Acc (%):")
    print(sensitivity)
    print(specificity)
    print(predict_accuracy)
    pl(" ")
    # pl(' ----------------------------------- ')
    
  }
  
  return( c(sensitivity,specificity,predict_accuracy))
  
}



fold_CV_dataset_gen <- function(train_set,num_folds) {
  # This function generates num_folds data sets for num_folds-fold CV. These data sets are then used in
  # the function that actually does the CV and returns its accuracy (averaged over the folds).
  #
  # It is assumed that the variants in the data set are either functional or non-functional -- i.e., labeled --
  # obtained using train_set <- subset(data_set, classification == 0 | classification == 1)
  
  N_train <- nrow(train_set)
  rnd_order <- sample(1:N_train)
  
  dsets <- list()
  fold_size <- floor(N_train/num_folds)
  
  fold_start <- 0
  fold_end <- 0
  
  for (I in 1:(num_folds - 1) ) { # fold number
    
    fold_start <- (I - 1)*fold_size + 1
    fold_end <- fold_start + fold_size - 1
    
    dsets[[I]] <- train_set[rnd_order[fold_start:fold_end],]
    
  }
  dsets[[num_folds]] <- train_set[ rnd_order[(fold_end + 1):N_train], ]
  
  # pl(" ~~~~~~~~~~~~~~~ ")
  # pl("K folds generation stats:")
  # pl("K is:")
  # print(num_folds)
  # pl("Standard fold size:")
  # print(fold_size)
  # pl("Last fold size:")
  # print(N_train - fold_end)
  # pl(" ~~~~~~~~~~~~~~~~ ")
  
  return(dsets)
  
}


K_fold_CV_glm <- function(data_set,num_folds,thresh_vec, link_func, reg_form) {
  # This function performs K-fold cross-validation.
  
  formula_string <- paste("classification ~ ", reg_form, sep = "")
  
  data_set <- subset(data_set, classification == 0 | classification == 1)
  data_set$classification <- factor(data_set$classification)
  
  fold_data <- fold_CV_dataset_gen(data_set,num_folds)
  fold_size_standard <- nrow(fold_data[[1]])
  fold_size_last <- nrow(fold_data[[num_folds]])
  fold_sizes <- rep_len(fold_size_standard, num_folds - 1)
  fold_sizes <- c(fold_sizes, fold_size_last)
  
  
  accuracies <- rep_len(0, num_folds)
  sensitivities <- rep_len(0, num_folds)
  specificities <- rep_len(0, num_folds)
  
  
  for (I in 1:num_folds) { # all possible test sets
    
    test_set <- fold_data[[I]]
    train_set <- c()
    ind <- 1:num_folds
    ind <- ind[-I]
    
    for (J in ind) {
      train_set <- rbind(train_set,fold_data[[J]])
    }
    
    fit <- glm(as.formula(formula_string), control = list(maxit = 500), data = train_set, family = binomial(link = link_func))
    PIFs <- predict(fit, test_set, type = "response")
    
    correct_predictions <- 0
    correct_positives <- 0
    correct_negatives <- 0
    
    for (J in 1:fold_sizes[I]) { # J goes through all variants in the test set
      
      if ( PIFs[J] > thresh_vec[1] & test_set$classification[J] == 1 ) {
        correct_predictions <- correct_predictions + 1
        correct_positives <- correct_positives + 1
      }
      else if ( PIFs[J] <= thresh_vec[2] & test_set$classification[J] == 0 ) {
        correct_predictions <- correct_predictions + 1
        correct_negatives <- correct_negatives + 1
      }
      
    }
    
    num_N <- nrow(subset(test_set, classification == 0))
    num_P <- nrow(subset(test_set, classification == 1))
    
    if (num_N > 0) { specificities[I] <- correct_negatives/num_N }
    else { specificities[I] <- NA }
    
    if (num_P > 0) { sensitivities[I] <- correct_positives/num_P }
    else { sensitivities[I] <- NA }
    
    accuracies[I] <- correct_predictions/fold_sizes[I]
    
    
  }
  
  Sn <- mean(sensitivities, na.rm = TRUE)*100
  Sp <- mean(specificities, na.rm = TRUE)*100
  Acc <- mean(accuracies)*100
  
  pl(' ----------------------------------- ')
  pl("K-fold CV output")
  pl(" ")
  pl(paste("Regression formula:", reg_form))
  pl(" ")
  pl("K is:")
  print(num_folds)
  pl("Sn, Sp, Acc (%):")
  print(Sn)
  print(Sp)
  print(Acc)
  pl(" ")
  # pl(' ----------------------------------- ')
  
  
  return(c(Sn,Sp,Acc))
  
}



acc_on_data_glm_with_CIs <- function(full_dset, thr, link_func, reg_form, boot_PIF_file_with_folder_name, fig_file_with_folder_name) {
  # This function does (probit) binomial regression on a data set, calculates the (fitting) accuracy, and 
  # makes a sorted-PIF plot.
  
  formula_string <- paste("classification ~ ", reg_form, sep = "")
  pathog_thresh <- thr[1]
  neutral_thresh <- thr[2]
  
  dset <- subset(full_dset, classification == 0 | classification == 1)
  dset$classification <- factor(dset$classification)
  fit <- glm(as.formula(formula_string), control = list(maxit = 500), data = dset, family = binomial(link = link_func))
  fit_null <- glm(classification ~ 1, control = list(maxit = 500), data = dset, family = binomial(link = link_func))
  
  
  PIFs <- predict(fit, dset, type = "response")
  # this is needed for accuracy calculation
  num_v <- nrow(dset)
  counter <- 0
  for (I in 1:num_v) {
    if ( (dset[[I,2]] == 0 & PIFs[I] <= neutral_thresh)
         | (dset[[I,2]] == 1 & PIFs[I] > pathog_thresh)) {counter <- counter + 1}
  }
  
  
  # -----------------
  # "plotting insert"
  load(boot_PIF_file_with_folder_name)
  
  # this is needed for plotting PIFs
  PIFs <- predict(fit, full_dset, type = "response")
  data_set <- cbind(full_dset,PIFs)
  data_set_out <- data_set # this will be returned by the function as part of its output
  data_set_out <- cbind(data_set_out, boot_out)
  
  data_set <- cbind(data_set, boot_out)
  data_set <- arrange(data_set, PIFs)
  
  # dev.new(width = 5.5, height = 2.8, unit = "in")
  pdf(fig_file_with_folder_name, width = 8, height = 4.5) # inches
  set_graphics_params()
  CI_col_unlabeled <- "gray60"
  plot(data_set$PIFs, cex = 1.5, 
       xlab = expression( paste( bolditalic("BRCA2 "), bold ("variants in order of increasing PIF" ) )  ), 
       ylab = "Probability of impact on function (PIF)" ) #, main = reg_form)
  abline(h = thr[1], lty = 2)
  abline(h = thr[2], lty = 2)
  
  large_CI_counter <- 0
  large_CI_thresh <- .01
  large_CI_counter_labeled <- 0
  
  
  CI_thickness <- 1
  
  for (I in 1:nrow(data_set)) {
    
    if (data_set[I,2] == 0) {points(I,data_set$PIFs[I], col = "blue", cex = 1.5)}
    else if (data_set[I,2] == 1) { points(I,data_set$PIFs[I], col = "red", cex = 1.5) }
    
    lines(  c(I,I), c(data_set$lo[I], data_set$up[I]), col = CI_col_unlabeled, lwd = CI_thickness )
    
    # if (reg_form == "DMSO + Cis + Ola") {
    #   if (data_set$PIFs[I] <= thr[1] & data_set$PIFs[I] > thr[2]) {
    #     text(I - 65, data_set$PIFs[I], labels = data_set[I,1], cex = .7, col = "black")
    #   }
    # }
    
    if ( data_set$up[I] - data_set$lo[I] > large_CI_thresh  )  
    { large_CI_counter <- large_CI_counter + 1}
    
    if ( (data_set$up[I] - data_set$lo[I] > large_CI_thresh) & !(data_set[[I,2]] == .5) )  
      { large_CI_counter_labeled <- large_CI_counter_labeled + 1}

  }
  
  # plotting PIF points on top of CIs, "for beauty" 
  points(data_set$PIFs, cex = 1.5)
  for (I in 1:nrow(data_set)) {
    
    if ( data_set[I,2] == .5 & data_set$PIFs[I] > neutral_thresh & data_set$PIFs[I] <= pathog_thresh ) {  
      points(I,data_set$PIFs[I], col = "cyan", cex = 1.5)
    }
    
  }
  for (I in 1:nrow(data_set)) {
    
    if (data_set[I,2] == 0) {
      points(I,data_set$PIFs[I], col = "blue", cex = 1.5)
      lines(  c(I,I), c(data_set$lo[I], data_set$up[I]), col = "blue", lwd = CI_thickness )
    }
    
    if (data_set[I,2] == 1) { 
      points(I,data_set$PIFs[I], col = "red", cex = 1.5)
      lines(  c(I,I), c(data_set$lo[I], data_set$up[I]), col = "red", lwd = CI_thickness )
    }
    
  }
  
  
  legend(-5,.95, legend=c("Non-functional PIF", "Functional PIF", "Experimental PIF", "Experimental intermediate PIF",
                        "Non-functional PIF CI","Functional PIF CI","Experimental PIF CI"), bg = "white",
         col=c("red", "blue", "black", "cyan", "red", "blue", CI_col_unlabeled), 
         lty=c(0,0,0,0,1,1,1), pch = c(1,1,1,1,NA,NA,NA), cex=.8)
  
  
  dev.off()
  
  # end of "plotting insert"
  # -----------
  
  
  acc <- counter/num_v
  out_list <- list()
  out_list[[1]] <- acc*100
  
  colnames(data_set_out)[colnames(data_set_out) == "lo"] <- "CI_lower"
  colnames(data_set_out)[colnames(data_set_out) == "up"] <- "CI_upper"
  colnames(data_set_out)[colnames(data_set_out) == "DMSO"] <- "Avg. DMSO_FS"
  colnames(data_set_out)[colnames(data_set_out) == "Cis"] <- "Avg. Cis_FS"
  colnames(data_set_out)[colnames(data_set_out) == "Ola"] <- "Avg. Ola_FS"
  colnames(data_set_out)[colnames(data_set_out) == "PIFs"] <- "PIF"
  data_set_out[[2]] <- as.character(data_set_out[[2]])
  for (I in 1:nrow(data_set_out)) {
    if (data_set_out[I,2] == "1") { data_set_out[I,2] <- "Non-functional" }
    else if (data_set_out[I,2] == "0") { data_set_out[I,2] <- "Functional" }
    else { data_set_out[I,2] <- "Experimental"  }
  }
  
  out_list[[2]] <- data_set_out
  out_list[[3]] <- fit
  out_list[[4]] <- fit_null
  
  pl(" ")
  pl(" ")
  #pl(" ------------------------------------- ")
  pl(paste("Regression formula: ", reg_form))
  pl(" ")
  pl("Fitting results")
  pl("Accuracy on the training set (%):")
  print(acc*100)
  pl(" ")
  pl("Number of variants with above-threshold CIs:")
  print(large_CI_counter)
  pl("...fracion of the data set (%):")
  print(large_CI_counter/nrow(full_dset)*100)
  pl("Number of labeled variants with above-threshold CIs:")
  print(large_CI_counter_labeled) 
  pl("...% of all labeled variants:")
  print(large_CI_counter_labeled/num_v*100)
  #pl(" ------------------------------------- ")
  
  return(out_list)
  
}





# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# --------------------- MAIN BODY OF THE CODE ---------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


plot(1,1)
text(1,1.1,"empty plot by design") 


# -----------------------------------------------------------------------------
# ------------- data loading and processing -----------------------------------
# -----------------------------------------------------------------------------


file_in <- paste(DATA_FOLDER, main_input_data_file, sep = "")
datafr <- read_excel(file_in)
datafr_main_out <- datafr # this is the output data frame, where we combine the input data and PIFs/CIs (see end of script)

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


pl(" ")
pl("Means:")
print(mean(datafr$DMSO))
print(mean(datafr$Cis))
print(mean(datafr$Ola))
pl(" ")


# ----------- Replacing P and N with 1 and 0, respectively -- for classification

classification <- c()
for (I in 1:nrow(datafr)) {
  if (datafr[I,2] == "Pathogenic") (classification[I] <- 1)
  else if (datafr[I,2] == "Neutral") (classification[I] <- 0)
  else (classification[I] <- 0.5)
}

datafr_classif <- data.frame(datafr$variant, classification, datafr$DMSO, datafr$Cis, datafr$Ola)
colnames(datafr_classif) <- c("variant", "classification", "DMSO", "Cis", "Ola")




# -----------------------------------------------------------------------------
# ------------- plotting PIFs for 3 vars (HAT, Cis, Ola) ----------------------
# -----------------------------------------------------------------------------

pl(" ")
pl(" ")

pl(" ------------------------- Accuracy on the full data set (labeled variants) for combined DMSO + Cis + Ola --------------------- ")

file_in <- paste(DATA_FOLDER, boot_PIF_file_all, sep = "")
fig_file <- paste(FIG_FOLDER, "f_PIFs_DMSO_Cis_Ola.pdf", sep = "")
mm <- acc_on_data_glm_with_CIs(datafr_classif, threshs, st_link, "DMSO + Cis + Ola", file_in, fig_file)
out_DMSO_Cis_Ola <- mm[[2]]
# out_file <- paste(DATA_FOLDER, "f_PIFs_DMSO_Cis_Ola_only.xlsx", sep = "")
# write_xlsx(out, path = out_file)



# -----------------------------------------------------------------------------
# ------------ plotting PIFs for DMSO -----------------------------------------
# -----------------------------------------------------------------------------


pl(" ")

pl(" ------------------------- Accuracy on the full data set (labeled variants) for DMSO ------------------------------- ")

file_in <- paste(DATA_FOLDER, boot_PIF_file_DMSO, sep = "")
fig_file <- paste(FIG_FOLDER, "f_PIFs_DMSO.pdf", sep = "")
mm <- acc_on_data_glm_with_CIs(datafr_classif, threshs, st_link, "DMSO", file_in, fig_file)
out_DMSO <- mm[[2]]



# -----------------------------------------------------------------------------
# ------------ plotting PIFs for Cis and Ola ----------------------------------
# -----------------------------------------------------------------------------

pl(" ")

pl(" ------------------------- Accuracy on the full data set (labeled variants) for Cis --------------------- ")

file_in <- paste(DATA_FOLDER, boot_PIF_file_Cis, sep = "")
fig_file <- paste(FIG_FOLDER, "f_PIFs_Cis.pdf", sep = "")
mm <- acc_on_data_glm_with_CIs(datafr_classif, threshs, st_link, "Cis", file_in, fig_file)
out_Cis <- mm[[2]]



pl(" ")

pl(" ------------------------- Accuracy on full data set (labeled variants) for Ola --------------------- ")

file_in <- paste(DATA_FOLDER, boot_PIF_file_Ola, sep = "")
fig_file <- paste(FIG_FOLDER, "f_PIFs_Ola.pdf", sep = "")
mm <- acc_on_data_glm_with_CIs(datafr_classif, threshs, st_link, "Ola", file_in, fig_file)
out_Ola <- mm[[2]]



# -----------------------------------------------------------------------------
# ---------- cross-validation for 3 vars -------------------------------------
# -----------------------------------------------------------------------------

K1 <- 5 
K2 <- 10 
K_max <- 50 # 50 labeled variants in the full data set



pl(" ")
pl(" ")
pl('----------------------------------------------------------------------------')
pl('--------------- cross-validation for combined DMSO + Cis + Ola --------------')
pl('----------------------------------------------------------------------------')
# both LOOCV and K-fold
pl(" ")

out <- LOOCV_glm(datafr_classif,threshs, st_link, "DMSO + Cis + Ola", FALSE, TRUE)
# verification
out <- K_fold_CV_glm(datafr_classif, K_max, threshs, st_link, "DMSO + Cis + Ola")

out1 <- K_fold_CV_glm(datafr_classif, K1, threshs, st_link, "DMSO + Cis + Ola")
out2 <- K_fold_CV_glm(datafr_classif, K2, threshs, st_link, "DMSO + Cis + Ola")


# -----------------------------------------------------------------------------
# ---------- cross-validation for DMSO ---------------------------------------
# -----------------------------------------------------------------------------
# both LOOCV and K-fold

pl(" ")
pl('----------------------------------------------------------------------------')
pl('--------------------- cross-validation for DMSO ----------------------------')
pl('----------------------------------------------------------------------------')
pl(" ")

out <- LOOCV_glm(datafr_classif,threshs, st_link, "DMSO", FALSE, TRUE)
# verification
out <- K_fold_CV_glm(datafr_classif, K_max, threshs, st_link, "DMSO")

out1 <- K_fold_CV_glm(datafr_classif, K1, threshs, st_link, "DMSO")
out2 <- K_fold_CV_glm(datafr_classif, K2, threshs, st_link, "DMSO")




# -----------------------------------------------------------------------------
# -------------------------- saving data --------------------------------------
# -----------------------------------------------------------------------------

colnames(out_DMSO_Cis_Ola)[6:8] <- c("PIF_DMSO_Cis_Ola","CI_lower_DMSO_Cis_Ola","CI_upper_DMSO_Cis_Ola")
colnames(out_DMSO)[6:8] <- c("PIF_DMSO","CI_lower_DMSO","CI_upper_DMSO")
colnames(out_Cis)[6:8] <- c("PIF_Cis","CI_lower_Cis","CI_upper_Cis")
colnames(out_Ola)[6:8] <- c("PIF_Ola","CI_lower_Ola","CI_upper_Ola")

our_class <- out_DMSO_Cis_Ola[2]
colnames(our_class)[1] <- "our_classification"
for (I in 1:nrow(out_DMSO_Cis_Ola)) {
  if (out_DMSO_Cis_Ola[I,6] > pathog_thresh) { our_class[I,1] <- "Non-functional" }
  else if ( out_DMSO_Cis_Ola[I,6] <= neutral_thresh ) { our_class[I,1] <- "Functional" }
  else ( our_class[I,1] <- "Intermediate" )
}

datafr_main_out <- cbind(datafr_main_out, out_DMSO_Cis_Ola[6:8], out_DMSO[6:8], out_Cis[6:8], out_Ola[6:8], our_class)
out_file <- paste(DATA_FOLDER, "FS_NNNSGE_BRCA2_withPIFs.xlsx", sep = "")


# saving project data
write_xlsx(datafr_main_out, path = out_file)













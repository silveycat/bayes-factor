# Code for running Case 1 simulation
# Reported in paper "Bayes factors for mixed-effects models"


# load libraries

# dplyr for data manipulation
library(dplyr)
# lme4 for mixed-effects models
library(lme4)
# ggplot2 for plotting
library(ggplot2)
# MASS for generating multivariate random distributions (for subject effects)
library(MASS)
# boot for inv.logit
library(boot)
# plotrix for standard error
library(plotrix)

# function to load lab libraries from Github
loadFunctionsGithub <-function(urlFolder, urlRaw, listFunctions){
  if (!require(httr)) {
    stop("httr not installed")
  } 
  else {
    print('----Downloading. Please wait----')
  };
  httr::GET(urlFolder)-> req
  stop_for_status(req)
  filelist <- unlist(lapply(content(req)$tree, "[", "path"), use.names = F)
  urlFunctions <- grep("docs/tools/", filelist, value = TRUE, fixed = TRUE)
  gsub("docs/tools/", "", urlFunctions) -> functions
  if (length(listFunctions) == 0){ #load all
    for (i in 1:length(functions)){
      httr::GET(paste0(urlRaw, functions[i]), ssl.verifypeer = FALSE)-> temp
      content(temp)->p3
      eval(parse(text = p3), envir = .GlobalEnv)
    } 
  } else {
    functions[functions %in% listFunctions]-> functionsIlike
    for (i in 1:length(functionsIlike)){
      httr::GET(paste0(urlRaw, functionsIlike[i]), ssl.verifypeer = FALSE)-> temp
      content(temp)->p3
      eval(parse(text = p3), envir = .GlobalEnv)
    }
  };
  print('----Download completed----')
}

# call function to load relevant libraries
urlFolder <- 'https://api.github.com/repos/n400peanuts/languagelearninglab/git/trees/master?recursive=1'
urlRaw <- 'https://raw.githubusercontent.com/n400peanuts/languagelearninglab/master/tools/'

# list functions to be loaded
listFunctions <- c("Bf.R")

loadFunctionsGithub(urlFolder, urlRaw, listFunctions)

# function for getting mode
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# function for categorising BFs
bf_cat <- function(bf) {
  # categorise one-tailed BFs
  if (bf > 100) {
    # extreme evidence for H1
    cat = 9
  } else if (bf > 30) {
    # very strong evidence for H1
    cat = 8
  } else if (bf > 10) {
    # strong evidence for H1
    cat = 7
  } else if (bf > 3) {
    # moderate evidence for H1
    cat = 6
  } else if (bf > 1/3) {
    # no evidence to speak of
    cat = 5
  } else if (bf > 1/10) {
    # moderate evidence for H0
    cat = 4
  } else if (bf > 1/30) {
    # strong evidence for H0
    cat = 3
  } else if (bf > 1/100) {
    # very strong evidence for H0
    cat = 2
  } else {
    # extreme evidence for H0
    cat = 1
  }
  return(cat)
}


# ** PARAMETERS FOR SIMULATION **

# number of subjects
n_subj = 40
# number of observations per subject - must be divisible by 2
n_obs = 20

# random effect SDs
subj_tau = c(0.4,0.9)
# correlation between subject intercepts and slopes
subj_corr = .2


# function to generate multilevel data with binary DV

generate_bin <- function(n_subj, n_obs, alpha, beta, subj_corr, subj_tau) {
  # make data frame
  data <- data.frame(matrix(0, ncol=0, nrow = n_subj * n_obs))
  # make subject vector
  data$subject <- as.factor(rep(seq(1:n_subj), each = n_obs))
  # make condition vector
  data$cond <- as.factor(rep(c(0,1), each = n_obs/2))
  # centre this
  data$c_cond <- as.numeric(data$cond) - mean(as.numeric(data$cond))
  # for subject effects
  # we specify the correlation between intercept and slope and put it in a matrix
  corr_matrix <- matrix(c(1,subj_corr,subj_corr,1), nrow = 2)
  # Variance covariance matrix for subject effects
  # multiplies the subject effect sds (in matrix form) by the correlation matrix
  # and then again by the subject effect sds
  # so we end up with the sds squared (on the diagonal) = variance, and covariance on the off-diagonal
  subj_v_cov <- diag(subj_tau,2,2) %*% corr_matrix  %*% diag(subj_tau,2,2)
  # Create the correlated subject effects using mvrnorm to sample from multivariate normal distribution
  # means of subject intercepts and slopes are 0
  u <- mvrnorm(n = n_subj, c(0,0), subj_v_cov)
  # check the correlation
  # print(cor(u))
  # check the SDs
  # print(sd(u[,1]))
  # print(sd(u[,2]))
  # generate mus for each subject and y for each observation
  data <- data %>%
    mutate(
      # We first calculate the linear predictor
      eta = alpha + u[data$subject,1] +
        data$c_cond * (beta + u[data$subject,2]),
      # then transform by inverse logit to a probability (for this subject/condition)
      mu = inv.logit(eta),
      # and generate a 0 or 1 for this row based on the probability for this subject/condition
      y = rbinom(nrow(data),1,mu))
  return(data)
}


# range of baseline values (performance in lower of 2 cells) to test
baseline_range <- seq(0, 3, 0.1)
# range of upper values (performance in higher of 2 cells) to test
values_range <- seq(0, 3, 0.1)

# set up data frame for collecting results

results_collect <- data.frame(matrix(0, ncol=14, nrow=0))
names(results_collect) <- c("N_subj","N_obs", "cell_1", "cell_2", "alpha", 
                            "beta",
                            "run", "grand_mean", "model_effect", 
                            "model_se",
                            "singular", 
                            "conv_error", "BF", "BF_cat")

# main simulation loop
for (b in baseline_range) {
  print(sprintf('current cell 1: %s', b))
  # define values range for other cell - starting from baseline and increasing to 3
  values_range <- seq(0, 3, 0.1)
  for (v in values_range) {
    print(sprintf('current cell 2: %s', v))
    for (run in 1:20) {
      # sometimes in the case of extreme evidence for the null, the BF is so small it comes out NaN
      # to avoid this, we initially set the BF to NaN and generate the data within a while loop
      BF = NaN
      # indicators for singularity/convergence error
      is_sing = "N"
      conv_error = "N"
      # calculate alpha and beta
      alpha = (b + v)/2
      beta = (v - b)
      # while loop to catch NaN BF
      while (is.nan(BF)) {
        # generate data
        data <- generate_bin(n_subj, n_obs, alpha, beta, subj_corr, subj_tau)
        # analyse using mixed-effects model
        model1 <- glmer(y ~ c_cond + (1 + c_cond|subject), family=binomial, data, control=glmerControl(optimizer = "bobyqa"))
        output <- summary(model1)
        # our H1 estimate here is the intercept
        mean_intercept <- output$coefficients["(Intercept)", "Estimate"]
        # we also save the effect estimate and SE as our summary of the data
        mean_effect <- output$coefficients["c_cond", "Estimate"]
        se_effect <- output$coefficients["c_cond", "Std. Error"]
        # calculate Bayes factor
        BF <- Bf(se_effect, mean_effect, likelihood = "normal", modeloftheory = "normal", modeoftheory = 0, 
                       scaleoftheory = mean_intercept, tail = 1, method = "new")
      }
      # categorise BF
      BF_cat <- bf_cat(BF)
      # check for singularity and convergence warnings and log them
      if (isSingular(model1)) {
        is_sing = "Y"
      }
      conv <- output$optinfo$conv$lme4$messages
      # collect if not duplicate of singularity warning
      if (!is.null(conv) && conv != "boundary (singular) fit: see ?isSingular") {
        conv_error = "Y"
      }
      # add all this to the data
      results_collect <- rbind(results_collect, data.frame(N_subj = n_subj, N_obs = n_obs, cell_1 = b,
                                                           cell_2 = v, alpha = alpha,
                                                           beta = beta,
                                                           run = run, grand_mean = mean_intercept,
                                                           model_effect = mean_effect, model_se = se_effect,
                                                           singular = is_sing,
                                                           conv_error = conv_error,
                                                           BF = BF,
                                                           BF_cat = BF_cat))
    }
  }
}

# write to file
write.csv(results_collect, "sim_1_results_40.csv", row.names = FALSE)

# plot by performance in cell1 (x-axis) and effect size (y-axis)

# calculate difference between cell2 and cell1 as a variable
results_collect$effect <- round(as.numeric(results_collect$cell_2) - as.numeric(results_collect$cell_1), 1)

# group by cell1 and effect size
grouped_data <- group_by(results_collect, cell_1, effect)

# get modal Bayes factor and make an ordered factor
summarised_data <- summarise(grouped_data, mode = getmode(BF_cat))
summarised_data$mode <- ordered(summarised_data$mode)

# colour brewer values for heatmap
color_values = c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b")

# plot
heatmap <- ggplot(summarised_data, aes(x=cell_1, y=effect, fill=mode))
heatmap + geom_tile() + scale_fill_manual("Evidence category", limits=c(1,2,3,4,5,6,7,8,9), 
                                          values = color_values,
                                          labels=c("Extreme H0", "V strong H0", "Strong H0", "Mod H0", 
                                                   "No evidence", 
                                                   "Mod H1", "Strong H1", "V strong H1", "Extreme H1")) +
  scale_x_continuous("\nPerformance in cell 1", limits=c(-0.1,3.1)) + 
  scale_y_continuous("Effect (difference between cell2 and cell1) \n", limits = c(-2.1,2.1)) +
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("Case 1 (N = 40)")
  

# save this
ggsave("Case_1_40.png")




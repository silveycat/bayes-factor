# Code for running Case 2 simulation
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
  # first check if this is one we skipped
  if (typeof(bf) == "character") {
    cat = NA
  } else {
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
    } else if (bf <= 1/100) {
      # extreme evidence for H0
      cat = 1
    }
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

# function to generate multilevel data with 2 interacting predictors and binary DV

generate_bin <- function(n_subj, n_obs, alpha, beta1, beta2, beta3, subj_corr, subj_tau) {
  # make data frame
  data <- data.frame(matrix(0, ncol=0, nrow = n_subj * n_obs))
  # make subject vector
  data$subject <- as.factor(rep(seq(1:n_subj), each = n_obs))
  # make condition 1 vector - within subjects
  data$cond1 <- as.factor(rep(c(0,1), each = n_obs/2))
  # centre this
  data$c_cond1 <- as.numeric(data$cond1) - mean(as.numeric(data$cond1))
  # make condition 2 vector - between subjects
  data$cond2 <- as.factor(rep(c(0,1), each = n_obs))
  # centre this
  data$c_cond2 <- as.numeric(data$cond2) - mean(as.numeric(data$cond2))
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
        data$c_cond1 * (beta1 + u[data$subject,2]) + data$c_cond2 * beta2 +
        (data$c_cond1 * data$c_cond2) * beta3,
      # then transform by inverse logit to a probability (for this subject/condition)
      mu = inv.logit(eta),
      # and generate a 0 or 1 for this row based on the probability for this subject/condition
      y = rbinom(nrow(data),1,mu))
  return(data)
}


# set values of cell1 and cell2
# this is set up so you can test multiple values within one simulation loop
# but be warned this makes the simulation extremely long!

cell1_values = c(0.5)
cell2_values = c(0.5)

# for 3rd and 4th cell, same range of values as for Case 1, i.e. from p 0.5 to around 0.95
cell3_values = seq(0, 3, 0.1)
cell4_values = seq(0, 3, 0.1)

# main simulation loop
for (a in cell1_values) {
  cell1 = a
  print(sprintf('current cell 1: %s', cell1))
  for (b in cell2_values) {
    cell2 = b
    print(sprintf('current cell 2: %s', cell2))
    # data frame for collecting results
    results_collect <- data.frame(matrix(0, ncol=21, nrow=0))
    names(results_collect) <- c("N_subj","N_obs", "cell_1", "cell_2", "cell_3", "cell_4",
                                "alpha", "beta1", "beta2", "beta3",
                                "run", "model_grand_mean", "model_beta1",
                                "model_beta3", "model_beta3_se",
                                "singular", 
                                "conv_error",
                                "BF_GM", "BF_GM_cat", "BF_ME", "BF_ME_cat")
    for (c in cell3_values) {
      cell3 = c
      print(sprintf('current cell 3: %s', cell3))
      for (x in cell4_values) {
        cell4 = x
        print(sprintf('current cell 4: %s', cell4))
        for (run in 1:20) {
          # set BFs initially to NaN for while loop
          BF_GM = NaN
          BF_ME = NaN
          is_sing = "N"
          conv_error = "N"
          # calculate true values of parameters
          alpha = (cell1 + cell2 + cell3 + cell4)/4
          beta1 = (cell3 + cell4)/2 - (cell1 + cell2)/2
          beta2 = (cell2 + cell4)/2 - (cell1 + cell3)/2
          beta3 = (cell4 - cell3) - (cell2 - cell1)
          # while loop to catch NaN BFs
          while (is.nan(BF_GM) || is.nan(BF_ME)) {
            data <- generate_bin(n_subj, n_obs, alpha, beta1, beta2, beta3, subj_corr, subj_tau)
            model1 <- glmer(y ~ c_cond1 * c_cond2 + (1 + c_cond1|subject), family=binomial, data, control=glmerControl(optimizer = "bobyqa"))
            output <- summary(model1)
            # get interaction estimate and SE
            mean_effect <- output$coefficients["c_cond1:c_cond2", "Estimate"]
            se_effect <- output$coefficients["c_cond1:c_cond2", "Std. Error"]
            # get BF for interaction. estimate here is 1) 2 * intercept
            grand_mean <- output$coefficients["(Intercept)", "Estimate"]
            BF_GM <- Bf(se_effect, mean_effect, likelihood = "normal", modeloftheory = "normal", modeoftheory = 0, 
                        scaleoftheory = 2 * grand_mean, tail = 1, method = "new")
            # 2) the main effect of the within-subjects factor (beta1)
            # skip this where beta1 is negative
            main_effect <- output$coefficients["c_cond1", "Estimate"]
            if (beta1 >= 0) {
              BF_ME <- Bf(se_effect, mean_effect, likelihood = "normal", modeloftheory = "normal", modeoftheory = 0, 
                          scaleoftheory = main_effect, tail = 1, method = "new")
            } else {
              BF_ME <- "skip"
            }
          }
          # categorise BFs
          BF_GM_cat <- bf_cat(BF_GM)
          BF_ME_cat <- bf_cat(BF_ME)
          # check for singularity and convergence warnings
          if (isSingular(model1)) {
            is_sing = "Y"
          }
          conv <- output$optinfo$conv$lme4$messages
          # collect if not duplicate of singularity warning
          if (!is.null(conv) && conv != "boundary (singular) fit: see ?isSingular") {
            conv_error = "Y"
          }
          # add all this to the data
          # we use stringsAsFactors = false to avoid running into problems with 'skip' rows
          current_row <- data.frame(N_subj = n_subj, N_obs = n_obs, cell_1 = cell1,
                                    cell_2 = cell2, cell_3 = cell3, cell_4 = cell4,
                                    alpha = alpha,
                                    beta1 = beta1, beta2 = beta2, beta3 = beta3,
                                    run = run, model_grand_mean = grand_mean,
                                    model_beta1 = main_effect,
                                    model_beta3 = mean_effect, 
                                    model_beta3_se = se_effect,
                                    singular = is_sing,
                                    conv_error = conv_error,
                                    BF_GM = BF_GM,
                                    BF_GM_cat = BF_GM_cat,
                                    BF_ME = BF_ME,
                                    BF_ME_cat = BF_ME_cat, stringsAsFactors = FALSE)
          results_collect <- rbind(results_collect, current_row)
        }
      }
    }
    # write to file
    filename <- sprintf('sim_2_me_results_%s_%s_%s.csv', n_subj, cell1, cell2)
    write.csv(results_collect, filename, row.names = FALSE)
  }
}


# plot - puts results from both estimates on one 2-facet graph

# calculate interaction effect as a variable
results_collect$effect <- round((as.numeric(results_collect$cell_4) - as.numeric(results_collect$cell_3)) - 
                                  (as.numeric(results_collect$cell_2) - as.numeric(results_collect$cell_1)), 1)

# group data
grouped_data <- group_by(results_collect, cell_3, effect)

# summarise ME and GM Bayes factor category
summarised_data <- summarise(grouped_data, mode = getmode(BF_ME_cat))
summarised_data_GM <- summarise(grouped_data, mode = getmode(BF_GM_cat))

# combine these
summarised_data$type <- "ME"
summarised_data_GM$type <- "GM"
summarised_data <- rbind(summarised_data, summarised_data_GM)

# mode is ordered factor

summarised_data$mode <- ordered(summarised_data$mode)

# facet labels

facet.labs <- c("Main effect", "2 * intercept")
names(facet.labs) <- c("ME", "GM")

# colour brewer values
color_values = c("#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b")


heatmap <- ggplot(summarised_data, aes(x=cell_3, y=effect, fill=mode, group=type))
heatmap + geom_tile() + scale_fill_manual("Evidence category", limits=c(1,2,3,4,5,6,7,8,9), 
                                          values = color_values,
                                          labels=c("Extreme H0", "V strong H0", "Strong H0", "Mod H0", 
                                                   "No evidence", 
                                                   "Mod H1", "Strong H1", "V strong H1", "Extreme H1")) +
  scale_x_continuous("\nLV post-test performance", limits=c(-0.1,3.1), 
                     breaks = c(0,1,2,3), 
                     labels = c("0 (.5)", "1 (.73)", "2 (.88)", "3 (.95)")) + 
  scale_y_continuous("Interaction effect \n", limits = c(-2.1,2.1)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(~ type, labeller = labeller(type = facet.labs)) + 
  ggtitle("Performance in pre-test at chance")

# save this
ggsave("Case_2_chance_40.png")

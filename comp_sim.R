# Simulation to compare brms bayesfactor results with Dienes calculator method
# Reported in paper "Bayes factors for mixed-effects models"

# load libraries ** check which of these are actually needed
library(dplyr)
library(RColorBrewer)
library(lme4)
library(brms)
library(ggplot2)
library(bayesplot)
library(BayesFactor)
library(MASS)
library(boot)

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


# define parameters for the data simulation

# number of subjects
n_subj = 40
# number of observations per subject - must be divisible by 2
n_obs = 20

# SDs of subject intercepts and slopes
subj_tau = c(0.4,0.9)
# correlation between subject intercepts and slopes
subj_corr = .2


# function for generating multilevel binary data with one within-subjects predictor
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

# first, one toy example to check timing

cell1 = 0
cell2 = 0.5

# calculate true intercept and condition effect
alpha = (cell1 + cell2)/2
beta = (cell2 - cell1)

# generate data
data <- generate_bin(n_subj, n_obs, alpha, beta, subj_corr, subj_tau)

# lmer/Dienes analysis function

run_dienes <- function() {
  start.time <- Sys.time()
  model1 <- glmer(y ~ c_cond + (1 + c_cond|subject), family=binomial, data, control=glmerControl(optimizer = "bobyqa"))
  output <- summary(model1)
  # our BF estimate here is the grand mean
  mean_intercept <- output$coefficients["(Intercept)", "Estimate"]
  mean_effect <- output$coefficients["c_cond", "Estimate"]
  se_effect <- output$coefficients["c_cond", "Std. Error"]
  onetail_BF <- Bf(se_effect, mean_effect, likelihood = "normal", modeloftheory = "normal", modeoftheory = 0, 
                   scaleoftheory = mean_intercept, tail = 1, method = "new")
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  result <- list('time' = time.taken, 'BF' = onetail_BF)
  return(result)
}

# run Dienes version
run_dienes()

# brms/bridge sampling analysis
# get intercept first

model1 <- glmer(y ~ c_cond + (1 + c_cond|subject), family=binomial, data, control=glmerControl(optimizer = "bobyqa"))
output <- summary(model1)
# our BF estimate here is the intercept
mean_intercept <- output$coefficients["(Intercept)", "Estimate"]

# brms/bridge sampling function
run_mc <- function() {
  start.time <- Sys.time()
  h1_string <- sprintf("normal(0, %f)", mean_intercept)
  h1_prior <- set_prior(h1_string, class="b", lb = 0)
  model1_full <- brm(y ~ c_cond + (1 + c_cond|subject), family = bernoulli, data = data, 
                     prior = h1_prior, warmup = 2000, iter = 20000, sample_prior = TRUE,
                     save_pars = save_pars(all = TRUE))
  
  model1_null <- brm(y ~ (1 + c_cond|subject), family = bernoulli, data = data, 
                     warmup = 2000, iter = 20000,
                     save_pars = save_pars(all = TRUE))
  brms_bf_example <- bayes_factor(model1_full, model1_null)$bf
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  result <- list('time' = time.taken, 'BF' = brms_bf_example)
  return(result)
}

# run brms version
run_mc()

# define ranges to test in simulation

cell1_range <- seq(0, 3, 1)
cell2_range <- seq(0, 3, 1)

# 16 cases, more manageable. but x10 runs for each case = still 320 brms models (null and full)...
# so in the morning we're not even halfway through (3rd run in 8th case of 16)
# end of 2nd day we're at 13th case of 16, run 3
# morning of 3rd day we're on run 4 of last case

# number of runs - try 10 to start
runs = 10


# collect results
# here, Dienes BF is 1, brms BF is 2
results_collect <- data.frame(matrix(0, ncol=14, nrow=0))
names(results_collect) <- c("N_subj","N_obs", "cell_1", "cell_2", "alpha", 
                            "beta",
                            "run", "grand_mean", "model_effect", 
                            "model_se", "BF_1", "BF_2", "BF_1_cat", "BF_2_cat")


# simulation loop

for (b in cell1_range) {
  print(sprintf('current cell 1: %s', b))
  for (v in cell2_range) {
    print(sprintf('current cell 2: %s', v))
    for (run in 1:runs) {
      # calculate alpha and beta
      alpha = (b + v)/2
      beta = (v - b)
      # initially set ZD BF to NaN
      ZD_BF = NaN
      # inside while loop, generate data and calculate BFs (to avoid NaN BFs)
      while (is.nan(ZD_BF)) {
        # generate data
        data <- generate_bin(n_subj, n_obs, alpha, beta, subj_corr, subj_tau)
        # lmer model
        model1 <- glmer(y ~ c_cond + (1 + c_cond|subject), family=binomial, data, control=glmerControl(optimizer = "bobyqa"))
        output <- summary(model1)
        # our BF scale estimate here is the grand mean
        # we take absolute value in order to account for the few cases where the estimate is negative
        mean_intercept <- abs(output$coefficients["(Intercept)", "Estimate"])
        mean_effect <- output$coefficients["c_cond", "Estimate"]
        se_effect <- output$coefficients["c_cond", "Std. Error"]
        # calculate and record Dienes BF
        ZD_BF <- Bf(se_effect, mean_effect, likelihood = "normal", modeloftheory = "normal", modeoftheory = 0, 
                    scaleoftheory = mean_intercept, tail = 1, method = "new")
        # define prior for brms analysis - half-normal with intercept from lmer as SD
        h1_string <- sprintf("normal(0, %f)", mean_intercept)
        h1_prior <- set_prior(h1_string, class="b", lb = 0)
        # run brms full model
        print(sprintf('running full model: cell1 %s, cell2 %s, run %d', b, v, run))
        model1_full <- brm(y ~ c_cond + (1 + c_cond|subject), family = bernoulli, data = data, 
                           prior = h1_prior, warmup = 2000, iter = 20000,
                           save_all_pars = TRUE)
        # run brms null model
        print(sprintf('running null model: cell1 %s, cell2 %s, run %d', b, v, run))
        model1_null <- brm(y ~ (1 + c_cond|subject), family = bernoulli, data = data, 
                           warmup = 2000, iter = 20000,
                           save_all_pars = TRUE)
        # run bridge sampling and record brms BF
        print(sprintf('bridge sampling: cell1 %s, cell2 %s, run %d', b, v, run))
        brms_BF <- bayes_factor(model1_full, model1_null)$bf
      }
      # categorise BFs
      ZD_BF_cat <- bf_cat(ZD_BF)
      brms_BF_cat <- bf_cat(brms_BF)
      # add results to data frame
      results_collect <- rbind(results_collect, data.frame(N_subj = n_subj, N_obs = n_obs, cell_1 = b,
                                                           cell_2 = v, alpha = alpha,
                                                           beta = beta,
                                                           run = run, grand_mean = mean_intercept,
                                                           model_effect = mean_effect, model_se = se_effect,
                                                           BF_1 = ZD_BF,
                                                           BF_2 = brms_BF,
                                                           BF_1_cat = ZD_BF_cat,
                                                           BF_2_cat = brms_BF_cat))
    }
  }
}

# write to file
write.csv(results_collect, "brms_sim_results_40.csv", row.names = FALSE)

# create table of results

results_collect <- mutate(results_collect, match = ifelse(BF_1_cat == BF_2_cat, TRUE, FALSE))

# group by cell1 and cell2
grouped_data <- group_by(results_collect, cell_1, cell_2)

# summarise
summarised_data <- summarise(grouped_data, matches = sum(match), 
                             modal_d = getmode(BF_1_cat), 
                             modal_b = getmode(BF_2_cat))

print(summarised_data)

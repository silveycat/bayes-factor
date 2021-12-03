This repository contains three R scripts:

1) `case_1_simulation.R` - simulation code generating Bayes factors for a range of veridical effect sizes in a design with one within-subjects predictor
2) `case_2_simulation.R` - simulation code generating Bayes factors for a range of veridical effect sizes in a design with an interaction between one within-subjects and one between-subjects predictor
3) `comp_sim.R` - simulation code generating Bayes factors via the Dienes method and the brms model comparison method for a range of veridical effect sizes in a design with one within-subjects predictor

The repository also includes, in the folder 'case_study_data', the data from the case study reported in the paper, and 
in the folder 'output', pre-generated .csv files containing the results from the simulations reported in the paper.

## Running the scripts
Please note that since comp_sim involves running 320 brms models, it takes a long time to run (several days on a new MacBook Air).

## Citation
> Silvey, C., Dienes, Z., & Wonnacott, E. (2021, December 3). Bayes factors for mixed-effects models. Retrieved from psyarxiv.com/m4hju
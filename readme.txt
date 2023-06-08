# THIS CODE HAS THE FOLLOWING SECTIONS
# 1. case 1: varying sample size
# # a. Non-optimal GMM
# # b. Non-optimal GMM with known intercept
# # c. Optimal GMM
# # d. Optimal GMM with known intercept
# # e. Pseudo-likelihood
#---------------------------
# 2. case 2: varying rho
# # a. Non-optimal GMM with known intercept
# # b. Optimal GMM with known intercept
# # c. Pseudo-likelihood
#---------------------------
# 3. case 3: model misspecification
# # a. Non-optimal GMM
# # b. Non-optimal GMM with known intercept
# # c. Optimal GMM
# # d. Optimal GMM with known intercept
# # e. Pseudo-likelihood
#---------------------------
# 4. tables and figures
#---------------------------


# 1. case 1: varying sample size
source("./sample-size.R")

# 2. case 2: varying rho
source("./varying-correlation.R")

# 3. case 3: model misspecification
source("./model-misspecification.R")

# 4. tables and figures
source("./table-figure.R")
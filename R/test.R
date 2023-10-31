library(devtools)

libs <- c("usethis", "devtoos", "gitcreds")
if (!requireNamespace(libs, quietly = TRUE)){
  install.packages(libs)
}

library(gitcreds)
library(usethis)

# Running git
usethis::use_git()
# Configurating my user id
use_git_config(user.name = "AraujoTheus", user.email = "matheusaraujolima@live.com")
# Creating a git project
usethis::use_github()

# Getting cred
gitcreds::gitcreds_set()


# Creating R files for my functions
use_r("nodes_config")

load_all()

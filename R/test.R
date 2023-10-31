library(devtools)

libs <- c("usethis", "devtoos", "gitcreds")
if (!requireNamespace(libs, quietly = TRUE)){
  install.packages(libs)
}


library(usethis)
usethis::use_git()

use_git_config(user.name = "AraujoTheus", user.email = "matheusaraujolima@live.com")

library(usethis)
create_github_token()

library(gitcreds)
gitcreds::gitcreds_set()

use_git()
usethis::use_github()
use_git()

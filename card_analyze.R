library(tidyverse)
library(googlesheets)

##aws s3 sync  s3://timpawsanalysis/171012_vre/kraken/ ~/Data/kraken/
workdir="~/Data/kraken"
plotdir="~/Dropbox/Data/Nanopore/171116_vre"
dir.create(plotdir, recursive=T)

gs_auth(token = "~/googlesheets_token.rds")

fullsheet=gs_url("https://docs.google.com/spreadsheets/d/1sOa1AP7K9mwNjgPX4qeBbRlqxAa6-gUBazj_QLvcAcE/edit?usp=sharing")
dataloc=gs_read(fullsheet, ws="VRE_171011")

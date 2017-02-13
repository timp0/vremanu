##As sudo
##source("https://raw.githubusercontent.com/fbreitwieser/pavian/master/inst/shinyapp/install-pavian.R")

datadir="/atium/Data/NGS/Aligned/161110_VREkraken/for_pavian"

options(pavian.server_access = TRUE)
options(pavian.server_dir = "/atium/Data/NGS/Aligned/161110_VREkraken/for_pavian")
options(pavian.start_data_input_with = "VREStudy")  # Tab name to start with - could also be "Example data" or "Upload files"



#library(rsconnect)
#rsconnect::deployApp('path/to/your/app')


pavian::runApp(host="128.220.138.134", port=5000, server_access=TRUE,
               server_dir=datadir)



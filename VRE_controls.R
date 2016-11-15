library(reshape)
library(ggplot2)

#read file into R here
controlsgg = read.csv("C:/Users/Stephanie/OneDrive/Documents/TimpLab/VREpaper/Illumina_controls_gg.csv")
controlsout <- split(controlsgg,controlsgg$Control)

spikein = c("0","10^2","10^4","10^6")

#reordered to match how we spiked into controls
KPNcontrol = c(controlsout$KPN$Percent[4],controlsout$KPN$Percent[6],controlsout$KPN$Percent[7],controlsout$KPN$Percent[5])
ENFMcontrol = c(controlsout$ENFM$Percent[4],controlsout$ENFM$Percent[7],controlsout$ENFM$Percent[5],controlsout$ENFM$Percent[6])
ENFScontrol = c(controlsout$ENFS$Percent[4], controlsout$ENFS$Percent[5],controlsout$ENFS$Percent[6],controlsout$ENFS$Percent[7])
combinedcontrols = data.frame(spikein,KPNcontrol,ENFMcontrol,ENFScontrol)

meltcombinedcontrols = melt(combinedcontrols)

controlsspikein = meltcombinedcontrols$spikein
controlsvalue = meltcombinedcontrols$value
controlsvariable = meltcombinedcontrols$variable

controlgraph = data.frame(controlsspikein,controlsvalue,controlsvariable)

ggplot(meltcombinedcontrols,aes(x=controlsspikein,y=controlsvalue,group=controlsvariable)) + geom <- line(aes(color = controlsvariable)) + geom <- point(aes(color=controlsvariable))


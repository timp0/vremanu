## script adapted from Yunfan's mito_tree.R file

library(ape)

clustdir="/home/shao4/VRE_KPC/"
plotdir="/home/shao4/VRE_KPC/"
tree=read.tree(paste0(clustdir, "VRE7KPC.dnd"))

pdf(file=paste0(plotdir, "VRE7_KPC_tree.pdf"), height=85, width=110)
plot(tree)
dev.off()

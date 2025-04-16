library(rlang)
library(ggplot2)
library(ggpubr)
source("../../public/public.R")
loadS = load("./NMDS_bac/smallbacplot.Rdata")
bacplot = eval(sym(loadS))
loadS = load("./NMDS_fgi/smallfgiplot.Rdata")
fgiplot = eval(sym(loadS))

ggpubr::ggarrange(bacplot$leaf,
                  bacplot$root,
                  bacplot$soil,
                  fgiplot$leaf,
                  fgiplot$root,
                  fgiplot$soil,
                  nrow = 6,ncol=1,labels = c('(b)','','','(c)'),
                  font.label=labelFont,align='hv')
ggsave(file="./both_arranged_plot.pdf",width=5,height=26.5,units="cm")


loadS = load("./NMDS_bac/00alldata.plot.Rdata")
bacAllPlot = eval(sym(loadS))
loadS = load("./NMDS_fgi//00alldata.plot.Rdata")
fgiAllPlot = eval(sym(loadS))
ggarrange(bacAllPlot,fgiAllPlot,ncol = 2,nrow=1,labels = c('(a)','(b)'),
          font.label = labelFont,common.legend = TRUE,align = 'hv')
ggsave(filename = "./allData.arranged.pdf",width = 18,height=10,unit="cm")

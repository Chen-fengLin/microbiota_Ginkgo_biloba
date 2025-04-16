source("../../public/public.R")
load("bac.19.plot.Rdata")
load("fgi.19.plot.Rdata")

library(ggpubr)
# ggpubr::ggarrange(bac.plot+theme(legend.position = "none"),
#                   fgi.plot+theme(legend.position = "none"),
#                   ncol=2,nrow=1,
#                   labels = c('(a)','(b)'),
#                   font.label = labelFont)
# 
# ggsave("bac+fgi.plot.pdf",width=20,height = 10,units = "cm")
# 


ggpubr::ggarrange(bac.plot,
                  fgi.plot,align = "hv",
                  ncol=2,nrow=1,
                  labels = c('(a)','(b)'),
                  font.label = labelFont)

ggsave("bac+fgi.plot.legend.pdf",width=20,height = 20,units = "cm")



#计算abundance图例
# 3 刻度 = 5.20155 pt
MIN_SCALE = 8.2
show_abund_tag=c(0.1,1,10,50)
data.frame(abund = show_abund_tag,
           `刻度` = log10(show_abund_tag*1e3/MIN_SCALE),
           pt = log10(show_abund_tag*1e3/MIN_SCALE)/3*5.20155)

# 60 刻度 = 8.1631 pt

show_contri_tag=c(20,40,60,80,100)

data.frame(contri = show_contri_tag,
           pt = show_contri_tag/60*8.1631)



#关闭
library(ggpubr)
source("../../public/public.R")
load("bac.29StatPlot.Rdata")
load("fgi.29StatPlot.Rdata")

ggpubr::ggarrange(bac.hos.plot+theme(axis.title.y = element_text(size = .fontsizePT)),
                  bac.sexDiff+theme(axis.title.y = element_text(size = .fontsizePT)),
                  fgi.hos.plot+theme(axis.title.y = element_text(size = .fontsizePT)),
                  fgi.sexDiff+theme(axis.title.y = element_text(size = .fontsizePT))+
                    scale_y_continuous(name="Sex-induced difference (%)",expand = c(0,0),
                                                 limits = function(o){max(abs(o))*c(-1,1)},
                                                 breaks = function(l){
                                                   if(l[2]<1){c(-0.5,0,0.5)}else{
                                                     if(l[2]<2){c(-1,0,1)}else{
                                                       c(-2,0,2)
                                                     }
                                                   }
                                                 }),
                  nrow=2,ncol=2,
                  align = "hv",
                  labels = c("(c)","(d)","(e)","(f)"),
                  font.label = labelFont,widths = c(2,2.9))

ggsave(file="bac+fgi.stat.pdf",width = 119/2,height = 23*2,units = "mm")

ggpubr::ggarrange(line.bac.plot,
                  line.fgi.plot,
                  nrow=2,ncol=1,align = "hv",
                  labels = c('(a)','(b)'),font.label = labelFont)

ggsave(filename = "../../02Ranalysis/63iCAMPCorr/iCAMP2BC.line.pdf",
       width = 9,height = 8,units = "cm")


ggpubr::ggarrange(bar.bac.plot,
                  bar.fgi.plot,
                  nrow = 2,ncol = 1,align = "hv",
                  labels = c('(a)','(b)'),
          common.legend = T,legend = "bottom",font.label = labelFont)
ggsave("../../02Ranalysis/34iCAMPamongNichceStage/bac.fgi.arrangeProcess.pdf",width = 9,height = 8,units = 'cm')


ggpubr::ggarrange(bar.bac.plot,
                  line.bac.plot,
                  bar.fgi.plot,
                  line.fgi.plot,
                  nrow = 2,ncol = 2,align = "hv",
                  labels = c('(a)','(c)','(b)','(d)'),
                  common.legend = T,legend = "bottom",font.label = labelFont)
ggsave("../../02Ranalysis/63iCAMPCorr/all.pdf",width = 19.5,height = 8.2,units = 'cm')


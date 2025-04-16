.linewidthPT = 1
.linewidthBoldPT = 1.5
.fontsizeTitlePT = 7.5
.fontsizePT = 6
.fontsizeSmallPT = 6
.symbolSizePT = 3

LINEPT = ggplot2::.pt*72.27/96

realLineWidth = function(linewidth){
  return(linewidth / LINEPT)
}

.linewidth = realLineWidth(.linewidthPT)
.linewidthBold = realLineWidth(.linewidthBoldPT)
.fontsizeTitle = .fontsizeTitlePT
.fontsize = .fontsizePT
.fontsizeSmall = .fontsizeSmallPT
.symbolSize = realLineWidth(.symbolSizePT)

leafColor = '#2ba62b'
rootColor = '#772ba6'
soilColor = 'brown'
colorList = c(leafColor,rootColor,soilColor)

S1Shape = 21
S2Shape = 24
S3shape = 22
shapeList = c(S1Shape,S2Shape,S3shape)
shapeList2 = c(S1Shape,S2Shape,S3shape,16,17,15)

labelFont = list(size=.fontsizeTitle,color="black",face="bold")

genelateBWtheme <- function(linewidthPT=.linewidthPT,linewidthBoldPT = .linewidthBoldPT,
                            fontsizeTitlePT = .fontsizeTitlePT,fontsizePT = .fontsizePT){
  linewidth = realLineWidth(linewidthPT)
  linewidthBold = realLineWidth(linewidthBoldPT)
  fontsizeTitle = fontsizeTitlePT
  fontsize = fontsizePT
  theme(axis.text = element_text(colour = 'black',size = fontsize),
        axis.title = element_text(colour = 'black',size=fontsizeTitle),
        legend.text = element_text(colour = 'black',size = fontsize),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_line(linewidth = linewidth),
        panel.border = element_rect(linewidth = linewidthBold,fill = "transparent",colour = "#000"),
        panel.background = element_blank(),plot.background = element_blank())
}

# genelateClassictheme <- function(linewidthPT=.linewidthPT,linewidthBoldPT = .linewidthBoldPT,
#                             fontsizeTitlePT = .fontsizeTitlePT,fontsizePT = .fontsizePT){
#   linewidth = realLineWidth(linewidthPT)
#   linewidthBold = realLineWidth(linewidthBoldPT)
#   fontsizeTitle = fontsizeTitlePT
#   fontsize = fontsizePT
#   theme_bw()+theme(axis.text = element_text(colour = 'black',size = fontsize),
#         axis.title = element_text(colour = 'black',size=fontsizeTitle),
#         legend.text = element_text(colour = 'black',size = fontsize),
#         legend.background = element_blank(),
#         legend.title = element_blank(),
#         axis.ticks = element_line(linewidth = linewidth))
# }


fmtMeta = function(sampleMeta){
  sampleMeta$sex = factor(sampleMeta$sex,levels = c('F','M'),ordered = TRUE)
  sampleMeta$period = factor(sampleMeta$period,levels = c('S1','S2','S3'),ordered = TRUE)
  sampleMeta$type = factor(sampleMeta$type,levels = c('L','R','S'),ordered = TRUE)
  sampleMeta$sexByName = factor(sampleMeta$sexByName,levels = c('Female','Male'),ordered = TRUE)
  sampleMeta$typeByName = factor(sampleMeta$typeByName,levels = c('Leaf','Root','Rhizosphere soil'),ordered = TRUE)
  sampleMeta$periodSex = factor(sampleMeta$periodSex,levels = c('S1F','S2F','S3F','S1M','S2M','S3M'),ordered = TRUE)
  sampleMeta$periodSexByName = factor(sampleMeta$periodSexByName,levels = c('S1 Female','S2 Female','S3 Female','S1 Male','S2 Male','S3 Male'),ordered = TRUE)
  rownames(sampleMeta) <- sampleMeta$sampleID
  sampleMeta
}

star = function(pvalue){
  ifelse(pvalue<0.001,'***',
         ifelse(pvalue<0.01,'**',
                ifelse(pvalue<0.05,'*',NA)))
}

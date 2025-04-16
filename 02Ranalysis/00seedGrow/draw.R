seed_data = read.table("seed.txt",header = T);

cal_summary <- function(dataset,vec,name){
  mean = aggregate(vec, list(seed_data$date,seed_data$tree),mean);
  sd = aggregate(vec, list(seed_data$date,seed_data$tree),sd);
  n = aggregate(vec, list(seed_data$date,seed_data$tree),length);
  summary.res = data.frame(
    date=mean$Group.1,
    tree=mean$Group.2,
    mean=mean$x,
    se=sd$x/sqrt(n$x),
    type=name
  );
  summary.res$date = as.Date(as.character(summary.res$date),"%Y%m%d")
  return(summary.res)
}

summary.length = cal_summary(seed_data,seed_data$length,'length (mm)')
summary.diameter = cal_summary(seed_data,seed_data$diameter,'diameter (mm)')
summary.fw = cal_summary(seed_data,seed_data$fw,'fresh weight (g)')
summary.dw = cal_summary(seed_data,seed_data$dw,'dry weight (g)')

summary.all = rbind(summary.length,summary.diameter,summary.fw,summary.dw)
summary.all$type = factor(summary.all$type,levels = c("dry weight (g)","fresh weight (g)","diameter (mm)","length (mm)"))

library(ggplot2)
lct <- Sys.getlocale("LC_TIME")  
#备份本地默认日期显示格式

Sys.setlocale("LC_TIME", "C")    
#指定标准日期显示格式

source("../../public/public.R")
plotColor = "#00ab00"
ggplot(data=summary.all,mapping = aes(x=date,y=mean))+
  geom_point(mapping=aes(),shape=21,color=plotColor,show.legend = F,size=.symbolSize/2)+
  geom_line(mapping=aes(group=tree),color=plotColor,show.legend = F,linewidth=.linewidth,alpha=0.6)+
  geom_errorbar(mapping=aes(x=date,ymin=mean-se,ymax=mean+se),
                color=plotColor,width=.symbolSize,show.legend=F,linewidth=.linewidth)+
  scale_y_continuous(position = "right",expand = expansion(c(0,0.05),c(0,0)))+
  scale_x_date(
    limits = c(as.Date("2023-04-01"),as.Date("2023-11-01")),
    breaks = c(as.Date("2023-04-01"),as.Date("2023-04-15"),as.Date("2023-05-01"),as.Date("2023-05-16"),as.Date("2023-06-01"),as.Date("2023-06-15"),as.Date("2023-07-01"),as.Date("2023-07-16"),as.Date("2023-08-01"),as.Date("2023-08-16"),as.Date("2023-09-01"),as.Date("2023-09-15"),as.Date("2023-10-01"),as.Date("2023-10-16"),as.Date("2023-11-01")),
    labels = c("","Apr","","May","","Jun","","Jul","","Aug","","Sep","","Oct","")
  )+
  facet_wrap(~type,ncol=1,scales="free_y",strip.position = "left")+
  genelateBWtheme()+
  theme(strip.background = element_blank(),
        strip.text = element_text(colour = 'black',size = .fontsize),
        axis.title = element_blank(),
        axis.ticks.x = element_line(colour = c("#000000",NA)))+
  geom_vline(xintercept = as.Date(c("2023-04-05","2023-06-14","2023-09-17")),lwd=.linewidthBold,color="#ff6000",alpha=0.6)

ggsave(filename = "./seedGrowth.pdf",width = 16,height = 16,units = 'cm')


Sys.setlocale("LC_TIME",lct)
#这一句是恢复默认系统日期显示格式

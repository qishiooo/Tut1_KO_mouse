library(ggplot2)

dat  <- read.delim('DSG.txt')
dat$PSI <- abs(dat$dPSI)
dat$Type <- factor(dat$Type, levels=c('AE', 'IR', 'A5SS', 'A3SS'))
for(group in unique(dat$Group)){
	dat2 <- dat[dat$Group==group,]
	p <- ggplot(dat2, aes(x=Type, y=PSI, color=Type))+
		geom_jitter(width=0.2)+
		theme_bw()+
		labs(x='', y='| delta PSI |')
	ggsave(paste0(group,'.dPSI.pdf'))

}

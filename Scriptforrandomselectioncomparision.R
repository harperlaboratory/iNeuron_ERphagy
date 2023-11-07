#Script for comparing random selection of proteins to organelle protein list of interest

library(ggplot2)
library(reshape2)
library(dplyr)

# Load data
# Input file must contain a column with the Gene names, a column with log2 fold-changes for each experiment and
# a column with the organelle/compartment annotation
file_Expression = "C:/Users/migue/Desktop/USA/people/Melissa/231020_DistrStats/additional sets/SuppData_CALCOCO1_allData.csv"
t = read.csv(file_Expression, header=T, stringsAsFactors=F, skipNul = TRUE)

# Adjust column names according to the input file
t3 = melt(t[, grepl("Gene.Symbol|log2FC_CALCOCO1KOd12.WTd12|log2FC_ATG12KOd12.WTd12|Compartment_ss",colnames(t))])
colnames(t3)[2] = "Organelle"

#Add all proteins for random selection
t3.all = t3
t3.all$Organelle = "All"
t3 = rbind(t3,t3.all)

# Select organelles/groups of interest
t3 = t3[!is.na(t3$Organelle),]
t3 = t3[grepl("Golgi|ER|All", t3$Organelle),]
t3 = t3[!is.na(t3$value),] # Remove NAs

# Generate violin plot
ggplot(data=t3, aes(x = variable, y = value, color=variable, fill=variable))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75),alpha=0.8, position = position_dodge(),color="black")+theme_classic()+geom_hline(yintercept = 0,linetype="dashed")+
  scale_fill_manual(values=c("maroon","#5EC962","#21918C","#3B528B","#440154"))+
  scale_color_manual(values=c("maroon","#5EC962","#21918C","#3B528B","#440154"))+
  facet_wrap(~Organelle)+
  labs(x="",y=bquote(""~Log[2]~"(KO / WT)"))+
  theme(legend.title = element_blank(),
        legend.text = element_text(size=16),
        legend.position = "top",
        axis.title = element_text(size=18),
        axis.text.x = element_text(size = 18,hjust=1,vjust = 0.2,angle=90),
        axis.text.y = element_text(size = 16),
        axis.ticks.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=16))+
  guides(fill = guide_legend(nrow = 1))

#Random protein selection per organelle and experiment, 100 times
j = 0
rndm.t = data.frame()
rndm.t.100 = data.frame()
while (j<100) {
  j = j+1
  rndm.t = data.frame()
  for (i in unique(t3$Organelle)) {
    tmp = t3[t3$Organelle == "All",] %>% group_by(variable) %>% slice_sample(n=nrow(t3[t3$Organelle == i,]) / length(unique(t3$variable)), replace = T)
    tmp$Rndm = i
    rndm.t = rbind(rndm.t,tmp)
  }
  rndm.t.100 = rbind(rndm.t.100, rndm.t)
}

#Subset random data keeping distribution to match organelle/group size
rndm.bin <- rndm.t.100 %>% #Break distribution
  group_by(Rndm, variable) %>%
  mutate(point_bin = cut_number(value, 100)) %>% #Adjust number of bins depending on organelle/group size
  mutate(binsize = length(unique(point_bin)))%>%
  mutate(groupsize = n())

rndm.bin <- rndm.bin %>% #Calculate probability each bin
  group_by(Rndm, variable, point_bin) %>%
  mutate(prob = (n() / groupsize)) %>%
  mutate(sample.size = round( prob*groupsize/100, 0))

rndm.bin.strat = rndm.bin %>% #Match organelle/group size keeping distribution
  group_by(Rndm, variable, point_bin) %>%
  sample_n(size=unique(sample.size), replace = T)

#Kolmogorov-Smirnov Test
rndm.bin.strat = as.data.frame(rndm.bin.strat)

t.stats = data.frame()
for (i in unique(t3$Organelle)) {
  for (j in unique(t3$variable)) {
    
    a = ks.test(t3[t3$Organelle == i & t3$variable == j, "value"],
                rndm.bin.strat[rndm.bin.strat$Rndm == i & rndm.bin.strat$variable == j, "value"])$p.val
    t.stats = rbind(t.stats,c(i,j,a))
  }
  
}
colnames(t.stats)= c("Organelle", "KO", "KS p-value")

# Generate density plots
for (i in unique(t3$Organelle)) {
  for (j in unique(t3$variable)) {
    
    plot(density(t3[t3$Organelle == i & t3$variable == j, "value"]),col = "orange",
         axes = T, xlab = bquote(""~Log[2]~"(KO / WT)"), ylab = "Density", xlim = c(-1, 1.5), lwd = 2, xaxs = "i",
         ylim = c(0, max(c(density(t3[t3$Organelle == i & t3$variable == j, "value"])$y,density(rndm.bin.strat[rndm.bin.strat$Rndm == i & rndm.bin.strat$variable == j, "value"])$y))),
         main = paste(i, " - ", j,"\n p-value", round(as.numeric(t.stats[t.stats$Organelle == i & t.stats$KO == j, "KS p-value"]),3)))
    lines(density(rndm.bin.strat[rndm.bin.strat$Rndm == i & rndm.bin.strat$variable == j, "value"]),
          col = "black",lwd = 2)
  }
  
}

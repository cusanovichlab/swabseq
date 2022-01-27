source("swab_plot_helper.R")
counts = read.table("RNA_vs_virus.data.txt")
counts = counts[,-2]
metadata = read.table("RNA_vs_virus.metadata.txt",header=T,sep="\t")
colnames(metadata)[1] = "Well_BC"
head(match(counts$V1,metadata$Well_BC))
test = cbind(counts,metadata[match(counts$V1,metadata$Well_BC),])

colnames(test)[1:3] = c("Centroid","target","Count")
test$Sample_ID = gsub("-1","",test$Sample_ID)
test$target[test$target == 0] = "S2"
test$target[test$target == 1] = "spike"
test$target[test$target == 4] = "RPP30"

pdf("Sequencingheatmap.pdf",height=5,width=10)
test %>%
  distinct(Sample_ID, Centroid, Count)  %>%
  count(Sample_ID, wt=Count, name = 'Well_Total') %>%
  separate(Sample_ID, into = c('Sample_Plate', 'Well'), sep = '-', remove=F) %>%
  mutate(
    Row = factor(str_sub(Well, 1, 1), levels = rev(LETTERS[1:16])),
    Col = str_sub(Well, 2)
  ) %>%
  ggplot(aes(x=Col, y=Row, fill=log10(Well_Total))) +
  geom_raster() +
  coord_equal() +
  scale_fill_viridis_c(option = 'plasma')

dev.off()

test.wide = test[test$target == "S2",]
spikesub = test[test$target == "spike",]
test.wide$Spike_Count = 0
test.wide$spike_Count[match(spikesub$Centroid,test.wide$Centroid)] = spikesub$Count
rppsub = test[test$target == "RPP30",]
test.wide$RPP30_Count = 0
test.wide$RPP30_Count[match(rppsub$Centroid,test.wide$Centroid)] = rppsub$Count
test.wide$S2_Spike_Ratio = (test.wide$Count+1)/(test.wide$Spike_Count + 1)
test.wide = test.wide[,c(1,5:8,15,9,11,3,16:18)]
colnames(test.wide)[1] = "Index"
colnames(test.wide)[2] = "Index1"
colnames(test.wide)[3] = "Index2"
colnames(test.wide)[9] = "S2_Count"

pdf("RNAvsVirus.pdf",height=8,width=8)
test.wide %>%
  mutate(RNA_virus_copies = if_else(RNA_virus_copies == 0, 3.162278, as.numeric(RNA_virus_copies))) %>%
  ggplot(aes(x=RNA_virus_copies, y=(S2_Count+1)/(Spike_Count+1), color=Plate_config)) +
  geom_quasirandom(alpha=1, aes(color=Plate_config),size=3) +
  geom_smooth(data=subset(test.wide,as.numeric(RNA_virus_copies) > 0),method=lm,se=FALSE) +
  geom_vline(xintercept = 5.623413) +
  scale_color_manual(values=c('#440154FF','#55C667FF')) +
  scale_x_log10(breaks = c(10^(-1:4)), labels = c(0,10^(0:4))) +
  scale_y_log10(limits = c(0.001,160)) +
  annotation_logticks() +
  theme(legend.position='none')
test.wide %>%
  mutate(RNA_virus_copies = if_else(RNA_virus_copies == 0, 3.162278, as.numeric(RNA_virus_copies))) %>%
  ggplot(aes(x=RNA_virus_copies, y=(S2_Count+1)/(Spike_Count+1), color=Plate_config)) +
  geom_quasirandom(alpha=1, aes(color=Plate_config),size=3) +
  geom_smooth(data=subset(test.wide,as.numeric(RNA_virus_copies) > 0),method=lm,se=FALSE) +
  geom_vline(xintercept = 5.623413) +
  scale_color_manual(values=c('#440154FF','#55C667FF')) +
  scale_x_log10(breaks = c(10^(-1:4)), labels = c(0,10^(0:4))) +
  scale_y_log10(limits = c(0.001,160)) +
  annotation_logticks()

test.wide %>%
  mutate(RNA_virus_copies = if_else(RNA_virus_copies == 0, 3.162278, as.numeric(RNA_virus_copies))) %>%
  ggplot(aes(x=RNA_virus_copies, y=(S2_Count+1)/(Spike_Count+1))) +
  geom_quasirandom(alpha=1, aes(color="indianred"),size=3) +
  geom_vline(xintercept = 5.623413) +
  scale_x_log10(breaks = c(10^(-1:4)), labels = c(0,10^(0:4))) +
  scale_y_log10(limits = c(0.001,160)) +
  annotation_logticks() +
  theme(legend.position='none')

dev.off()

write.table(test.wide,"RNA_vs_virus_wide_mat.txt",row.names=F,quote=F,sep="\t")

rpps = which(test.wide$RPP30_Count >= 1)
s2s = which(test.wide$S2_Count+test.wide$Spike_Count >= 500)
nots2s = which(test.wide$S2_Count+test.wide$Spike_Count < 500)
inconclusive = setdiff(nots2s,rpps)
pdf("ReadDepthFilters.pdf")
plot(log10(test.wide$RPP30_Count+1),log10(test.wide$S2_Count+test.wide$Spike_Count+1),pch=20,
     main = paste0(nrow(test.wide)," total, ",
                   length(s2s)," pass S2+Spike filter\n",
                   length(rpps)," pass RPP30 filter (",length(intersect(rpps,s2s))," pass both)\n",
                   length(inconclusive)," inconclusive"),xlab="Log10(RPP30+1)",
     ylab="Log10(S2+Spike+1)")
abline(v=log10(2),col="red",lty="dashed")
abline(h=log10(501),col="red",lty="dashed")
dev.off()

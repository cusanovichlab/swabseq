source("swab_plot_helper.R")
counts = read.table("Gargle_samples.data.txt")
counts = counts[,-2]
metadata = read.table("Gargle_samples.metadata.txt",header=T,sep="\t")
colnames(metadata)[1] = "Well_BC"
head(match(counts$V1,metadata$Well_BC))
test = cbind(counts,metadata[match(counts$V1,metadata$Well_BC),])

colnames(test)[1:3] = c("Centroid","target","Count")
twos = grep("-2",test$Sample_ID)
if(length(twos) > 0){
  test = test[-grep("-2",test$Sample_ID),]
}
test$Sample_ID = gsub("-1","",test$Sample_ID)
test$target[test$target == 0] = "S2"
test$target[test$target == 1] = "spike"
test$target[test$target == 4] = "RPP30"

pdf("Sequencingheatmap.pdf",height=10,width=20)
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

test.wide = test[test$target == "spike",]
s2sub = test[test$target == "S2",]
test.wide$S2_Count = 0
test.wide$S2_Count[match(s2sub$Centroid,test.wide$Centroid)] = s2sub$Count
spikesub = test[test$target == "spike",]
test.wide$Spike_Count = 0
test.wide$Spike_Count[match(spikesub$Centroid,test.wide$Centroid)] = spikesub$Count
rppsub = test[test$target == "RPP30",]
test.wide$RPP30_Count = 0
test.wide$RPP30_Count[match(rppsub$Centroid,test.wide$Centroid)] = rppsub$Count
test.wide$S2_Spike_Ratio = (test.wide$S2_Count+1)/(test.wide$Spike_Count + 1)
test.wide = test.wide[,c(1,5:8,15,10,11,16:21)]
colnames(test.wide)[1] = "Index"
colnames(test.wide)[2] = "Index1"
colnames(test.wide)[3] = "Index2"
test.wide$Ct_N1[grep("Undetermined",test.wide$Ct_N1)] = 41

test.wide.rep1 = test.wide[grep("rep1",test.wide$Replicates),]
test.wide.rep2 = test.wide[grep("rep2",test.wide$Replicates),]
test.wide.reps = test.wide.rep1
test.wide.reps$S2_Count_rep2 = 0
test.wide.reps$Spike_Count_rep2 = 0
test.wide.reps$RPP30_Count_rep2 = 0
test.wide.reps$S2_Spike_Ratio_rep2 = 0
test.wide.reps[order(test.wide.reps$Plate_samples),15:18] = test.wide.rep2[order(test.wide.rep2$Plate_samples),11:14]

pdf("Gargle_samples_swab_vs_qpcr.pdf",height=8,width=8)
test.wide.reps %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=S2_Spike_Ratio,color=Plate_ID)) +
  geom_point(size=5) +
  scale_color_manual(values='dodgerblue2') +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b") +
  geom_smooth(data=subset(test.wide.reps,as.numeric(Ct_N1) < 41 & S2_Spike_Ratio >= 0.002),method=lm,se=FALSE,color="black") +
  geom_vline(xintercept = 40) +
  geom_hline(yintercept = 0.002) +
  theme(legend.position='none')
test.wide.reps %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=S2_Spike_Ratio_rep2,color=Plate_ID)) +
  geom_point(size=5) +
  scale_color_manual(values='cadetblue2') +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b") +
  geom_smooth(data=subset(test.wide.reps,as.numeric(Ct_N1) < 41 & S2_Spike_Ratio_rep2 >= 0.002),method=lm,se=FALSE,color="black") +
  geom_vline(xintercept = 40) +
  geom_hline(yintercept = 0.002) +
  theme(legend.position='none')
test.wide.reps %>%
  ggplot(aes(x=S2_Spike_Ratio, y=S2_Spike_Ratio_rep2,color=Plate_ID)) +
  geom_point(size=5) +
  scale_color_manual(values='tomato') +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "bl") +
  geom_vline(xintercept = 0.002) +
  geom_hline(yintercept = 0.002) +
  geom_abline(color="black") +
  theme(legend.position='none')

dev.off()

tester1 = test.wide.reps[which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio >= 0.002),]
cor(as.numeric(tester1$Ct_N1),as.numeric(log10(tester1$S2_Spike_Ratio)))^2
cor.test(as.numeric(tester1$Ct_N1),as.numeric(log10(tester1$S2_Spike_Ratio)))$p.value

tester2 = test.wide.reps[which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002),]
cor(as.numeric(tester2$Ct_N1),as.numeric(log10(tester2$S2_Spike_Ratio_rep2)))^2
cor.test(as.numeric(tester2$Ct_N1),as.numeric(tester2$S2_Spike_Ratio_rep2))$p.value

tester3 = test.wide.reps[which(test.wide.reps$S2_Spike_Ratio >= 0.002 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002),]
cor(log10(tester3$S2_Spike_Ratio),as.numeric(log10(tester3$S2_Spike_Ratio_rep2)))^2
cor.test(as.numeric(tester3$S2_Spike_Ratio),as.numeric(tester3$S2_Spike_Ratio_rep2))$p.value

#squares
#N1/rep1 pos/pos
length(which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio >= 0.002))
#N1/rep1 neg/pos
length(which(as.numeric(test.wide.reps$Ct_N1) == 41 & test.wide.reps$S2_Spike_Ratio >= 0.002))
#N1/rep1 pos/neg
length(which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio < 0.002))
#N1/rep1 neg/neg
length(which(as.numeric(test.wide.reps$Ct_N1) == 41 & test.wide.reps$S2_Spike_Ratio < 0.002))

#N1/rep2 pos/pos
length(which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002))
#N1/rep2 neg/pos
length(which(as.numeric(test.wide.reps$Ct_N1) == 41 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002))
#N1/rep2 pos/neg
length(which(as.numeric(test.wide.reps$Ct_N1) < 41 & test.wide.reps$S2_Spike_Ratio_rep2 < 0.002))
#N1/rep2 neg/neg
length(which(as.numeric(test.wide.reps$Ct_N1) == 41 & test.wide.reps$S2_Spike_Ratio_rep2 < 0.002))

#rep1/rep2 pos/pos
length(which(test.wide.reps$S2_Spike_Ratio >= 0.002 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002))
#rep1/rep2 neg/pos
length(which(as.numeric(test.wide.reps$S2_Spike_Ratio) < 0.002 & test.wide.reps$S2_Spike_Ratio_rep2 >= 0.002))
#rep1/rep2 pos/neg
length(which(as.numeric(test.wide.reps$S2_Spike_Ratio) >= 0.002 & test.wide.reps$S2_Spike_Ratio_rep2 < 0.002))
#rep1/rep2 neg/neg
length(which(as.numeric(test.wide.reps$S2_Spike_Ratio) < 0.002 & test.wide.reps$S2_Spike_Ratio_rep2 < 0.002))

test.wide$Ct_N1[test.wide$Ct_N1 == 41] = "Undetermined"
test.wide$Ct_RPP30[test.wide$Ct_RPP30 == 41] = "Undetermined"
write.table(test.wide,"Gargle_samples_wide_mat.txt",row.names=F,quote=F,sep="\t")

rpps = which(test.wide$RPP30_Count >= 1)
s2s = which(test.wide$S2_Count+test.wide$Spike_Count >= 500)
nots2s = which(test.wide$S2_Count+test.wide$Spike_Count < 500)
inconclusive = setdiff(nots2s,rpps)
pdf("ReadDepthFilters2.pdf")
plot(log10(test.wide$RPP30_Count+1),log10(test.wide$S2_Count+test.wide$Spike_Count+1),pch=20,
     main = paste0(nrow(test.wide)," total, ",
                   length(s2s)," pass S2+Spike filter\n",
                   length(rpps)," pass RPP30 filter (",length(intersect(rpps,s2s))," pass both)\n",
                   length(inconclusive)," inconclusive"),xlab="Log10(RPP30+1)",
     ylab="Log10(S2+Spike+1)")
abline(v=log10(2),col="red",lty="dashed")
abline(h=log10(501),col="red",lty="dashed")
dev.off()

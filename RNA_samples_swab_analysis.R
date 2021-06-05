source("swab_plot_helper.R")
counts = read.table("RNA_samples.data.txt")
counts = counts[,-2]
metadata = read.table("RNA_samples.metadata.txt",header=T)
colnames(metadata)[1] = "Well_BC"
head(match(counts$V1,metadata$Well_BC))
metadata$Ct_N1[metadata$Ct_N1 == "Undetermined"] = 46
metadata$Ct_N2[metadata$Ct_N2 == "Undetermined"] = 46
metadata$Ct_RPP30[metadata$Ct_RPP30 == "Undetermined"] = 46
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
  facet_wrap(~Sample_Plate) +
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
test.wide = test.wide[,c(1,5,6,7,8,14,10,15:21)]
colnames(test.wide)[1] = "Index"
colnames(test.wide)[2] = "Index1"
colnames(test.wide)[3] = "Index2"

test.wide.subjects = test.wide[grep("ASYMPT|SYMPT",test.wide$Plate_config),]
test.wide.subjects = test.wide.subjects[(test.wide.subjects$S2_Count + test.wide.subjects$Spike_Count) >= 500,]
test.wide.controls = test.wide[-grep("ASYMPT|SYMPT",test.wide$Plate_config),]
test.wide.controls = test.wide.controls[(test.wide.controls$S2_Count + test.wide.controls$Spike_Count) >= 500,]
pdf("swab_vs_qpcr_samples_only2.pdf",height=8,width=8)
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=S2_Spike_Ratio,color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N2), y=S2_Spike_Ratio,color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_RPP30), y=S2_Spike_Ratio,color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_RPP30), y=(RPP30_Count+1)/(Spike_Count+1),color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=as.numeric(Ct_N2),color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(base = 2, sides = "bl")
test.wide.subjects %>%
  ggplot(aes(y=S2_Spike_Ratio, x=(RPP30_Count+1)/(Spike_Count+1),color=Plate_ID)) +
  geom_point(size=3) +
  scale_color_manual(values=c('#FDE725FF','#55C667FF', '#2D708EFF','#440154FF')) +
  scale_y_log10() +
  scale_x_log10() +
  annotation_logticks(sides = "bl")
test.wide.subjects %>%
  ggplot(aes(y=S2_Spike_Ratio, x=(RPP30_Count+1)/(Spike_Count+1),color=log10(Spike_Count+1))) +
  geom_point(size=3) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_viridis_c(option = 'plasma') +
  annotation_logticks(sides = "bl")

test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=S2_Spike_Ratio)) +
  geom_point(color="mediumorchid",size=3) +
  geom_smooth(data=subset(test.wide.subjects,as.numeric(Ct_N1) < 45),method=lm,se=FALSE,color="black") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_RPP30), y=S2_Spike_Ratio)) +
  geom_point(color="mediumseagreen",size=3) +
  geom_smooth(data=subset(test.wide.subjects,as.numeric(Ct_RPP30) < 45),method=lm,se=FALSE,color="black") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_RPP30), y=(RPP30_Count+1)/(Spike_Count+1))) +
  geom_point(color="mediumseagreen",size=3) +
  geom_smooth(data=subset(test.wide.subjects,as.numeric(Ct_RPP30) < 45),method=lm,se=FALSE,color="black") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=as.numeric(Ct_N2))) +
  geom_point(color="dodgerblue2",size=3) +
  geom_smooth(data=subset(test.wide.subjects,as.numeric(Ct_N1) < 45 & as.numeric(Ct_N2)<45),method=lm,se=FALSE,color="black") +
  scale_y_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  annotation_logticks(base = 2, sides = "bl")
dev.off()

#The Ct values provided for N1 and N2 were found to be in high correlation
#(R2 = 0.95, p-value = 4.2x10-38, Fig. 2B)
tester = test.wide.subjects[which(as.numeric(test.wide.subjects$Ct_N1) < 45 & test.wide.subjects$Ct_N2 < 45),]
cor(as.numeric(tester$Ct_N1),as.numeric(tester$Ct_N2))^2
cor.test(as.numeric(tester$Ct_N1),as.numeric(tester$Ct_N2))$p.value
#qRT-PCR Ct values for N1 and Swab-Seq Log10-transformed S2/Spike-in ratios
#were highly correlated (R2 = 0.92 for samples detectable by both methods,
#p-vlaue = 1.0x10-34)
tester2 = test.wide.subjects[which(as.numeric(test.wide.subjects$Ct_N1) < 45 & test.wide.subjects$S2_Spike_Ratio >= 0.01),]
cor(as.numeric(tester2$Ct_N1),as.numeric(log10(tester2$S2_Spike_Ratio)))^2
cor.test(as.numeric(tester2$Ct_N1),as.numeric(log10(tester2$S2_Spike_Ratio)))$p.value
#Ct values for RPP30 (human control gene) were not correlated with the
#log-transformed Swab-Seq ratio (R2 = 0.002, p-value = 0.50, Fig. S2E)
tester3 = test.wide.subjects[which(as.numeric(test.wide.subjects$Ct_RPP30) < 45),]
cor(as.numeric(tester3$Ct_RPP30),as.numeric(log10(tester3$S2_Spike_Ratio)))^2
cor.test(as.numeric(tester3$Ct_RPP30),as.numeric(log10(tester3$S2_Spike_Ratio)))$p.value
#the correlation of RPP30 counts from the Swab-Seq assay with the Ct value
#of RPP30 was initially low (R2 = 0.21, p-value = 1.0x10-13, Fig. S2F)
cor(as.numeric(tester3$Ct_RPP30),as.numeric(log10((tester3$RPP30_Count+1)/(tester3$Spike_Count+1))))^2
cor.test(as.numeric(tester3$Ct_RPP30),as.numeric(log10((tester3$RPP30_Count+1)/(tester3$Spike_Count+1))))$p.value
#which we utilized as a normalizing factor for both S2 and RPP30 in the
#Swab-Seq data (R2 = 0.53, p-value = 3.2x10-35, Fig. S2G)
tester4 = test.wide.subjects[which(as.numeric(test.wide.subjects$Ct_RPP30) < 45 & test.wide.subjects$Spike_Count > 600),]
cor((tester4$RPP30_Count),as.numeric(log10((tester4$RPP30_Count+1)/(tester4$Spike_Count+1))))^2
cor.test((tester4$RPP30_Count),as.numeric(log10((tester4$RPP30_Count+1)/(tester4$Spike_Count+1))))$p.value

#identifying positive samples at this Ct threshold (65 for N1, 63 for N2,
#59 detected in common)
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46))
length(which(as.numeric(test.wide.subjects$Ct_N2) < 46))
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & as.numeric(test.wide.subjects$Ct_N2) < 46))
#Swab-Seq was able to detect 94% (61/65) of N1-positive samples (Ct < 45)
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46))
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$S2_Spike_Ratio >= 0.01))
#of which the N2 qRT-PCR primer set only detected 91% (59/65)
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$Ct_N2 < 46))
#Among 124 negative control wells, Swab-Seq gave a positive result for four
falsepos = test.wide.controls[-grep("SPOS", test.wide.controls$Plate_config),]
dim(falsepos)
length(which(falsepos$S2_Spike_Ratio >= 0.01))
length(which(falsepos$S2_Spike_Ratio >= 0.01 & as.numeric(falsepos$Ct_N1) < 46))

#squares
#N1/N2 pos/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$Ct_N2 < 46))
#N1/N2 neg/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$Ct_N2 < 46))
#N1/N2 pos/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$Ct_N2 == 46))
#N1/N2 neg/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$Ct_N2 == 46))

#N1/swab pos/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$S2_Spike_Ratio >= 0.01))
#N1/swab neg/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$S2_Spike_Ratio >= 0.01))
#N1/swab pos/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$S2_Spike_Ratio < 0.01))
#N1/swab neg/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$S2_Spike_Ratio < 0.01))

test.wide$Ct_N1[test.wide$Ct_N1 == 46] = "Undetermined"
test.wide$Ct_N2[test.wide$Ct_N2 == 46] = "Undetermined"
test.wide$Ct_RPP30[test.wide$Ct_RPP30 == 46] = "Undetermined"
write.table(test.wide,"RNA_samples_wide_mat.txt",row.names=F,quote=F,sep="\t")

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

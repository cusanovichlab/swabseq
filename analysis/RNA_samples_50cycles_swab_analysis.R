source("swab_plot_helper.R")
counts = read.table("RNA_samples_50cycles.data.txt")
counts = counts[,-2]
metadata = read.table("RNA_samples_50cycles.metadata.txt",header=T,sep="\t")
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

test.wide.40 = read.table("RNA_samples_wide_mat.txt",header=T)
test.wide.40$Ct_N1[test.wide.40$Ct_N1 == "Undetermined"] = 46
test.wide.40$Ct_N2[test.wide.40$Ct_N2 == "Undetermined"] = 46
test.wide.40$Ct_RPP30[test.wide.40$Ct_RPP30 == "Undetermined"] = 46
subjects.common = intersect(test.wide.40$Plate_config,test.wide$Plate_config)
test.wide.40.common = test.wide.40[match(subjects.common,test.wide.40$Plate_config),]
test.wide.50.common = test.wide[match(subjects.common,test.wide$Plate_config),]
test.wide.both = cbind(test.wide.50.common,test.wide.40.common[,11:14])
colnames(test.wide.both)[15:18] = c("S2_Count_40","Spike_Count_40","RPP30_Count_40","S2_Spike_Ratio_40")
test.wide.both = test.wide.both[which(test.wide.both$S2_Count+test.wide.both$Spike_Count >= 500 & test.wide.both$S2_Count_40+test.wide.both$Spike_Count_40 >= 500),]
test.wide.both.subjects = test.wide.both[grep("ASYMPT|SYMPT",test.wide.both$Plate_config),]
test.wide.both.controls = test.wide.both[-grep("ASYMPT|SYMPT",test.wide.both$Plate_config),]

pdf("40vs50.pdf")
test.wide.both.subjects %>%
  ggplot(aes(y=S2_Spike_Ratio, x=S2_Spike_Ratio_40,color=as.numeric(Ct_N1))) +
  geom_point(size=3) +
  scale_y_log10() +
  scale_x_log10() +
  geom_abline() +
  geom_hline(yintercept = 0.002) +
  geom_vline(xintercept = 0.01) +
  scale_color_viridis_c(option = 'plasma') +
  annotation_logticks(sides = "bl")

dev.off()

test.wide.subjects = test.wide[grep("ASYMPT|SYMPT",test.wide$Plate_config),]
test.wide.subjects = test.wide.subjects[(test.wide.subjects$S2_Count + test.wide.subjects$Spike_Count) >= 500,]
test.wide.controls = test.wide[-grep("ASYMPT|SYMPT",test.wide$Plate_config),]
test.wide.controls = test.wide.controls[(test.wide.controls$S2_Count + test.wide.controls$Spike_Count) >= 500,]
pdf("swab50_vs_qpcr_samples_only.pdf",height=8,width=8)
test.wide.subjects %>%
  ggplot(aes(x=as.numeric(Ct_N1), y=S2_Spike_Ratio)) +
  geom_point(color="orange",size=3) +
  geom_smooth(data=subset(test.wide.subjects,as.numeric(Ct_N1) < 45),method=lm,se=FALSE,color="black") +
  scale_y_log10() +
  scale_x_continuous(breaks = c(10,15,20,25,30,35,40,45)) +
  geom_hline(yintercept=0.002) +
  annotation_logticks(sides = "l") +
  annotation_logticks(base = 2, sides = "b")
dev.off()

#N1-positives with an R2  = 0.91
tester2 = test.wide.subjects[which(as.numeric(test.wide.subjects$Ct_N1) < 45 & test.wide.subjects$S2_Spike_Ratio >= 0.002),]
cor(as.numeric(tester2$Ct_N1),as.numeric(log10(tester2$S2_Spike_Ratio)))^2
cor.test(as.numeric(tester2$Ct_N1),as.numeric(log10(tester2$S2_Spike_Ratio)))$p.value
#Swab-Seq40 and Swab-Seq50 (from four for Swab-Seq40 to 10 for Swab-Seq50)
test.wide.40.subjects = test.wide.40[grep("ASYMPT|SYMPT",test.wide.40$Plate_config),]
test.wide.40.subjects = test.wide.40.subjects[which(test.wide.40.subjects$S2_Count+test.wide.40.subjects$Spike_Count >= 500),]
length(which(test.wide.subjects$Ct_N1 == 46 & test.wide.subjects$S2_Spike_Ratio >= 0.002))
length(which(test.wide.40.subjects$Ct_N1 == 46 & test.wide.40.subjects$S2_Spike_Ratio >= 0.01))
#Comparing the two Swab-Seq tests directly on subjects, the R2 = 0.96
#(p-value = 5.2x10-40) for samples positive by both
tester_both = test.wide.both.subjects[which(test.wide.both.subjects$S2_Spike_Ratio >= 0.002 & test.wide.both.subjects$S2_Spike_Ratio_40 >= 0.01),]
tester_both = tester_both[which(tester_both$S2_Count + tester_both$Spike_Count >= 500 & tester_both$S2_Count_40 + tester_both$Spike_Count_40 >= 500),]
cor(log10(tester_both$S2_Spike_Ratio),log10(tester_both$S2_Spike_Ratio_40))^2
cor.test(log10(tester_both$S2_Spike_Ratio),log10(tester_both$S2_Spike_Ratio_40))$p.value
#Swab-Seq50 detected 88% (57/65) of the samples detected as positive with Swab-Seq40
length(which(test.wide.both.subjects$S2_Spike_Ratio >= 0.002 & test.wide.both.subjects$S2_Spike_Ratio_40 >= 0.01))
length(which(test.wide.both.subjects$S2_Spike_Ratio_40 >= 0.01))
#Swab-Seq40 identified eight positives that were not corroborated by Swab-Seq50
#(four of which were also positive by qRT-PCR)
length(which(test.wide.both.subjects$S2_Spike_Ratio_40 >= 0.01 & test.wide.both.subjects$S2_Spike_Ratio < 0.002))
length(which(test.wide.both.subjects$S2_Spike_Ratio_40 >= 0.01 & test.wide.both.subjects$S2_Spike_Ratio < 0.002 & test.wide.both.subjects$Ct_N1 < 46))
#while Swab-Seq50 identified 12 that were not detected with Swab-Seq40 (two of
#which were also qRT-PCR positive)
length(which(test.wide.both.subjects$S2_Spike_Ratio_40 < 0.01 & test.wide.both.subjects$S2_Spike_Ratio >= 0.002))
length(which(test.wide.both.subjects$S2_Spike_Ratio_40 < 0.01 & test.wide.both.subjects$S2_Spike_Ratio >= 0.002 & test.wide.both.subjects$Ct_N1 < 46))
#11 negative control samples were detected as positive
test.wide.controls.neg = test.wide.controls[-grep("SPOS",test.wide.controls$Plate_config),]
length(which(test.wide.controls.neg$S2_Spike_Ratio >= 0.002))
#a potential false positive rate as high as 7% (8/118)
length(which(test.wide.controls.neg$S2_Spike_Ratio >= 0.002 & test.wide.controls.neg$Ct_N1 == 46))
nrow(test.wide.controls.neg) - length(which(test.wide.controls.neg$S2_Spike_Ratio >= 0.002 & test.wide.controls.neg$Ct_N1 < 46))

#square
#N1/swab50 pos/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$S2_Spike_Ratio >= 0.002))
#N1/swab50 neg/pos
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$S2_Spike_Ratio >= 0.002))
#N1/swab50 pos/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) < 46 & test.wide.subjects$S2_Spike_Ratio < 0.002))
#N1/swab50 neg/neg
length(which(as.numeric(test.wide.subjects$Ct_N1) == 46 & test.wide.subjects$S2_Spike_Ratio < 0.002))

test.wide$Ct_N1[test.wide$Ct_N1 == 46] = "Undetermined"
test.wide$Ct_N2[test.wide$Ct_N2 == 46] = "Undetermined"
test.wide$Ct_RPP30[test.wide$Ct_RPP30 == 46] = "Undetermined"
write.table(test.wide,"RNA_samples_50cycles_wide_mat.txt",row.names=F,quote=F,sep="\t")

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

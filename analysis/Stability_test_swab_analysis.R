source("swab_plot_helper.R")
counts = read.table("Stability_test.data.txt")
counts = counts[,-2]
metadata = read.table("Stability_test.metadata.txt",header=T,sep="\t")
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
test.wide = test.wide[,c(1,5:8,18,9,11:14,19:22)]
colnames(test.wide)[1] = "Index"
colnames(test.wide)[2] = "Index1"
colnames(test.wide)[3] = "Index2"
test.wide$Storage_duration[test.wide$Storage_duration == "0 day"] = 0
test.wide$Storage_duration[test.wide$Storage_duration == "1 day"] = 1
test.wide$Storage_duration[test.wide$Storage_duration == "7 days"] = 7
test.wide$Storage_duration[test.wide$ATCC_Virus_copies == 0] = -1.75

test.wide.panels = test.wide[-grep("None",test.wide$Storage_condition),]
test.wide.zeroes = test.wide[grep("None",test.wide$Storage_condition),]
test.wide.freezer.zeroes = test.wide.zeroes
test.wide.freezer.zeroes$Storage_condition = "freezer"
test.wide.rt.zeroes = test.wide.zeroes
test.wide.rt.zeroes$Storage_condition = "RT"
test.wide.panels = rbind(test.wide.panels,test.wide.freezer.zeroes)
test.wide.panels = rbind(test.wide.panels,test.wide.rt.zeroes)
test.wide.panels$Saliva_vs_Gargle = factor(test.wide.panels$Saliva_vs_Gargle,levels=c("Saliva","Gargle"))

pdf("TimeCourse.pdf",height=8,width=8)
test.wide.panels %>%
  ggplot(aes(x=as.numeric(Storage_duration), y=S2_Spike_Ratio)) +
  geom_quasirandom(alpha=1, aes(color=Subjects),size=3) +
  geom_smooth(data=subset(test.wide.panels,as.numeric(ATCC_Virus_copies) > 1),method=lm,se=FALSE,color="gray") +
  scale_shape_manual(values=c(8, 19, 1))+
  scale_color_manual(values=c("dodgerblue2","orange","indianred")) +
  geom_vline(xintercept = -1) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0,7,1)) +
  facet_wrap(~Saliva_vs_Gargle + Storage_condition) +
  annotation_logticks(sides = "l") +
  theme(legend.position='none')
test.wide.panels %>%
  ggplot(aes(x=as.numeric(Storage_duration), y=S2_Spike_Ratio)) +
  geom_quasirandom(alpha=1, aes(color=Subjects),size=3) +
  geom_smooth(data=subset(test.wide.panels,as.numeric(ATCC_Virus_copies) > 1),method=lm,se=FALSE,color="gray") +
  scale_shape_manual(values=c(8, 19, 1))+
  scale_color_manual(values=c("dodgerblue2","orange","indianred")) +
  geom_vline(xintercept = -1) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(0,7,1)) +
  facet_wrap(~Saliva_vs_Gargle + Storage_condition) +
  annotation_logticks(sides = "l")

dev.off()

write.table(test.wide,"Stability_test_wide_mat.txt",row.names=F,quote=F,sep="\t")

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

test.wide.nonzero = test.wide.panels[test.wide.panels$ATCC_Virus_copies != 0,]
test.wide.gargle = test.wide.nonzero[test.wide.nonzero$Saliva_vs_Gargle == "Gargle",]
test.wide.saliva = test.wide.nonzero[test.wide.nonzero$Saliva_vs_Gargle == "Saliva",]
test.wide.gargle.rt = test.wide.gargle[test.wide.gargle$Storage_condition == "RT",]
test.wide.saliva.rt = test.wide.saliva[test.wide.saliva$Storage_condition == "RT",]
test.wide.gargle.freezer = test.wide.gargle[test.wide.gargle$Storage_condition == "freezer",]
test.wide.saliva.freezer = test.wide.saliva[test.wide.saliva$Storage_condition == "freezer",]

lrter = function(object,dfs=2){
  A = logLik(lm(object$S2_Spike_Ratio ~ object$Subjects))
  B = logLik(lm(object$S2_Spike_Ratio ~ as.numeric(object$Storage_duration) + object$Subjects))
  teststat <- -2 * (as.numeric(A)-as.numeric(B))
  p.val <- pchisq(teststat, df = dfs, lower.tail = FALSE)
  return(p.val)
}

lrter(test.wide.saliva.freezer)
#[1] 0.5195341

lrter(test.wide.saliva.rt)
#[1] 0.9881141

lrter(test.wide.gargle.freezer)
#[1] 0.09418925

lrter(test.wide.gargle.rt)
#[1] 0.6468406
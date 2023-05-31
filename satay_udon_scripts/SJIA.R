library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

groupsk10 = read.table("/Users/kairaveethakkar/Downloads/Revised-26-sample/ICGS-NMF-k10-ward-rho0.3-c10/FinalGroups-CellTypesFull_forFishers.txt",sep = "\t",stringsAsFactors = F,header = F)
groupsk15 = read.table("/Users/kairaveethakkar/Downloads/Revised-26-sample/ICGS-NMF-k15-ward-rho0.3-c12/FinalGroups-CellTypesFull_forFishers.txt",sep = "\t",stringsAsFactors = F,header = F)
groupsk20 = read.table("/Users/kairaveethakkar/Downloads/Revised-26-sample/ICGS-NMF-k20-ward-rho0.3-c14/FinalGroups-CellTypesFull_forFishers.txt",sep = "\t",stringsAsFactors = F,header = F)
groupsk30 = read.table("/Users/kairaveethakkar/Downloads/Revised-26-sample/ICGS-NMF-k30-ward-rho0.3-c17/FinalGroups-CellTypesFull_forFishers.txt",sep = "\t",stringsAsFactors = F,header = F)

final_list_k10_active = fisherstest(groups_df = groupsk10,diseasesubtype = "Active",p_val = 0.1)
final_list_k10_inactive = fisherstest(groups_df = groupsk10,diseasesubtype = "Inactive",p_val = 0.1)
final_list_k10_mas = fisherstest(groups_df = groupsk10,diseasesubtype = "MAS",p_val = 0.1)
final_list_k10_lungdisease = fisherstest(groups_df = groupsk10,diseasesubtype = "LungDisease",p_val = 0.1)

final_list_k20_active = fisherstest(groups_df = groupsk20,diseasesubtype = "Active",p_val = 0.1)
final_list_k20_inactive = fisherstest(groups_df = groupsk20,diseasesubtype = "Inactive",p_val = 0.1)
final_list_k20_mas = fisherstest(groups_df = groupsk20,diseasesubtype = "MAS",p_val = 0.1)
final_list_k20_lungdisease = fisherstest(groups_df = groupsk20,diseasesubtype = "LungDisease",p_val = 0.1)

final_list_k15_active = fisherstest(groups_df = groupsk15,diseasesubtype = "Active",p_val = 0.1)
final_list_k15_inactive = fisherstest(groups_df = groupsk15,diseasesubtype = "Inactive",p_val = 0.1)
final_list_k15_mas = fisherstest(groups_df = groupsk15,diseasesubtype = "MAS",p_val = 0.1)
final_list_k15_lungdisease = fisherstest(groups_df = groupsk15,diseasesubtype = "LungDisease",p_val = 0.1)

final_list_k30_active = fisherstest(groups_df = groupsk30,diseasesubtype = "Active",p_val = 0.1)
final_list_k30_inactive = fisherstest(groups_df = groupsk30,diseasesubtype = "Inactive",p_val = 0.1)
final_list_k30_mas = fisherstest(groups_df = groupsk30,diseasesubtype = "MAS",p_val = 0.1)
final_list_k30_lungdisease = fisherstest(groups_df = groupsk30,diseasesubtype = "LungDisease",p_val = 0.1)

pm = final_list_k10_mas$pval 
p_lgm = pm < 0.1
p_lgm = p_lgm*1

pa = final_list_k10_active$pval 
p_lga = pa < 0.1
p_lga = p_lga*1

pi = final_list_k10_inactive$pval 
p_lgi = pi < 0.1
p_lgi = p_lgi*1

pld = final_list_k10_lungdisease$pval 
p_lgld = pld < 0.1
p_lgld = p_lgld*1


orm = final_list_k10_mas$OR
orm[orm=="Inf"] = NA

ora = final_list_k10_active$OR 
ora[ora=="Inf"] = NA

ori = final_list_k10_inactive$OR
ori[ori=="Inf"] = NA

orld = final_list_k10_lungdisease$OR
orld[orld=="Inf"] = NA


plot1=pheatmap::pheatmap(p_lga,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Active_PVal")


plot2=pheatmap::pheatmap(p_lgi,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Inactive_PVal")


plot3=pheatmap::pheatmap(p_lgld,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "LungDisease_PVal")


plot4=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "MAS_PVal")


plot5=pheatmap::pheatmap(ora,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Active_ORVal")


plot6=pheatmap::pheatmap(ori,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Inactive_ORVal")


plot7=pheatmap::pheatmap(orld,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "LungDisease_ORVal")


plot8=pheatmap::pheatmap(orm,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "MAS_ORVal")


plot_list=list()
plot_list[['plot1']]=plot1[[4]]
plot_list[['plot2']]=plot2[[4]]
plot_list[['plot3']]=plot3[[4]]
plot_list[['plot4']]=plot4[[4]]
plot_list[['plot5']]=plot5[[4]]
plot_list[['plot6']]=plot6[[4]]
plot_list[['plot7']]=plot7[[4]]
plot_list[['plot8']]=plot8[[4]]


pdf("K10__Version2_26Samples_FishersExactTest.pdf",width = 25,height = 11)
x = grid.arrange(grobs=plot_list, ncol=4,nrow=2)
show(x)

dev.off()




pm = final_list_k15_mas$pval 
p_lgm = pm < 0.1
p_lgm = p_lgm*1

pa = final_list_k15_active$pval 
p_lga = pa < 0.1
p_lga = p_lga*1

pi = final_list_k15_inactive$pval 
p_lgi = pi < 0.1
p_lgi = p_lgi*1

pld = final_list_k15_lungdisease$pval 
p_lgld = pld < 0.1
p_lgld = p_lgld*1


orm = final_list_k15_mas$OR
orm[orm=="Inf"] = NA

ora = final_list_k15_active$OR 
ora[ora=="Inf"] = NA

ori = final_list_k15_inactive$OR
ori[ori=="Inf"] = NA

orld = final_list_k15_lungdisease$OR
orld[orld=="Inf"] = NA


plot1=pheatmap::pheatmap(p_lga,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Active_PVal")


plot2=pheatmap::pheatmap(p_lgi,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Inactive_PVal")


plot3=pheatmap::pheatmap(p_lgld,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "LungDisease_PVal")


plot4=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "MAS_PVal")


plot5=pheatmap::pheatmap(ora,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Active_ORVal")


plot6=pheatmap::pheatmap(ori,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Inactive_ORVal")


plot7=pheatmap::pheatmap(orld,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "LungDisease_ORVal")


plot8=pheatmap::pheatmap(orm,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "MAS_ORVal")


plot_list=list()
plot_list[['plot1']]=plot1[[4]]
plot_list[['plot2']]=plot2[[4]]
plot_list[['plot3']]=plot3[[4]]
plot_list[['plot4']]=plot4[[4]]
plot_list[['plot5']]=plot5[[4]]
plot_list[['plot6']]=plot6[[4]]
plot_list[['plot7']]=plot7[[4]]
plot_list[['plot8']]=plot8[[4]]


pdf("k15__Version2_26Samples_FishersExactTest.pdf",width = 25,height = 11)
x = grid.arrange(grobs=plot_list, ncol=4,nrow=2)
show(x)

dev.off()




pm = final_list_k20_mas$pval 
p_lgm = pm < 0.1
p_lgm = p_lgm*1

pa = final_list_k20_active$pval 
p_lga = pa < 0.1
p_lga = p_lga*1

pi = final_list_k20_inactive$pval 
p_lgi = pi < 0.1
p_lgi = p_lgi*1

pld = final_list_k20_lungdisease$pval 
p_lgld = pld < 0.1
p_lgld = p_lgld*1


orm = final_list_k20_mas$OR
orm[orm=="Inf"] = NA

ora = final_list_k20_active$OR 
ora[ora=="Inf"] = NA

ori = final_list_k20_inactive$OR
ori[ori=="Inf"] = NA

orld = final_list_k20_lungdisease$OR
orld[orld=="Inf"] = NA


plot1=pheatmap::pheatmap(p_lga,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Active_PVal")


plot2=pheatmap::pheatmap(p_lgi,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Inactive_PVal")


plot3=pheatmap::pheatmap(p_lgld,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "LungDisease_PVal")


plot4=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "MAS_PVal")


plot5=pheatmap::pheatmap(ora,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Active_ORVal")


plot6=pheatmap::pheatmap(ori,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Inactive_ORVal")


plot7=pheatmap::pheatmap(orld,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "LungDisease_ORVal")


plot8=pheatmap::pheatmap(orm,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "MAS_ORVal")


plot_list=list()
plot_list[['plot1']]=plot1[[4]]
plot_list[['plot2']]=plot2[[4]]
plot_list[['plot3']]=plot3[[4]]
plot_list[['plot4']]=plot4[[4]]
plot_list[['plot5']]=plot5[[4]]
plot_list[['plot6']]=plot6[[4]]
plot_list[['plot7']]=plot7[[4]]
plot_list[['plot8']]=plot8[[4]]


pdf("k20__Version2_26Samples_FishersExactTest.pdf",width = 25,height = 11)
x = grid.arrange(grobs=plot_list, ncol=4,nrow=2)
show(x)

dev.off()


pm = final_list_k30_mas$pval 
p_lgm = pm < 0.1
p_lgm = p_lgm*1

pa = final_list_k30_active$pval 
p_lga = pa < 0.1
p_lga = p_lga*1

pi = final_list_k30_inactive$pval 
p_lgi = pi < 0.1
p_lgi = p_lgi*1

pld = final_list_k30_lungdisease$pval 
p_lgld = pld < 0.1
p_lgld = p_lgld*1


orm = final_list_k30_mas$OR
orm[orm=="Inf"] = NA

ora = final_list_k30_active$OR 
ora[ora=="Inf"] = NA

ori = final_list_k30_inactive$OR
ori[ori=="Inf"] = NA

orld = final_list_k30_lungdisease$OR
orld[orld=="Inf"] = NA


plot1=pheatmap::pheatmap(p_lga,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Active_PVal")


plot2=pheatmap::pheatmap(p_lgi,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "Inactive_PVal")


plot3=pheatmap::pheatmap(p_lgld,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "LungDisease_PVal")


plot4=pheatmap::pheatmap(p_lgm,border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,color = c("black","red"),main = "MAS_PVal")


plot5=pheatmap::pheatmap(ora,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Active_ORVal")


plot6=pheatmap::pheatmap(ori,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "Inactive_ORVal")


plot7=pheatmap::pheatmap(orld,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "LungDisease_ORVal")


plot8=pheatmap::pheatmap(orm,color = colorRampPalette(c("blue", "black", "yellow"))(50),border_color = "white",na_col = "white",cluster_rows = F,cluster_cols = F,cellwidth = 10,cellheight = 10,main = "MAS_ORVal")


plot_list=list()
plot_list[['plot1']]=plot1[[4]]
plot_list[['plot2']]=plot2[[4]]
plot_list[['plot3']]=plot3[[4]]
plot_list[['plot4']]=plot4[[4]]
plot_list[['plot5']]=plot5[[4]]
plot_list[['plot6']]=plot6[[4]]
plot_list[['plot7']]=plot7[[4]]
plot_list[['plot8']]=plot8[[4]]


pdf("k30__Version2_26Samples_FishersExactTest.pdf",width = 25,height = 11)
x = grid.arrange(grobs=plot_list, ncol=4,nrow=2)
show(x)

dev.off()






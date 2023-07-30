#!/usr/bin/Rscript

library(devtools)
#install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
library(dplyr)
library(circlize)
library(RColorBrewer)
library(readxl)


########## VD #############
############
son <- read.csv('bbmap_VD_abundance.csv')
vd <- read_excel("~/paht/to/read_stats_taxonomy_VD.xlsx", 
sheet = 2)

rownames(son) <- son[,1]
son <- son[,-1]

#Grouping?
######################
rt <- as.matrix(replicate(11, 'Rotarix'))
tb <- as.matrix(replicate(11, 'Tb-Chen'))
au <- as.matrix(replicate(11, 'AU-1'))
group <- c(au, rt, tb)
son$Group=group


#another way
au <- son[son$Group == 'AU-1', ] 
rt <- son[son$Group == 'Rotarix', ] 
tb <- son[son$Group == 'Tb-Chen', ] 

son <- rbind(rt, tb, au)
################


#read count cut-off of 500
son[son<500 ] <- 0
son[nrow(son)+1, ] <- c(t(vd[,2]))
rownames(son)[34] <- 'QC_Reads'

son <- son[,1:52]

sonm <- as.matrix(son)
sonm <- apply(son, 2, function(x){x/max(x)*10000})
sonm <- sonm[c(-34),]


rownames(sonm) <- c('NSP1.2', 'NSP2.2', 'NSP3.2', 'NSP5.2', 'VP1.2', 'VP2.2', 'VP3.2','VP4.2', 'VP6.2', 'VP7.2',
                    'NSP4.2', 'VP4', 'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5', 'VP1', 'VP2', 'VP3', 'VP6', 
                    'VP7', 'NSp1', 'NSp2', 'NSp3', 'NSp4', 'NSp5', 'VP1.1', 'VP2.1', 'VP3.1','VP4.1', 'VP6.1', 'VP7.1')


#Normalization

slog <- log10(sonm)
#min value is determined by 
min(slog[slog != min(slog)])
slog[which(!is.finite(slog))] <- 0
slog[slog < 0] <- 0

#Add the trimmed read counts per sample again
slog <- rbind(slog, t(vd[,2]))

#remove columns with only zero values - MAYBE NOT DO - ONLY REMOVES ONLY AU-1 VP6
#s <- sonm[apply(sonm[,-1], 1, function(x) !all(x==0)),]

col_fun = colorRamp2(c(0, 200000, 400000, 600000, 800000, 1200000), c("white","bisque2", "burlywood2", "coral1", "brown3", "brown4"))
col_fun2 = colorRamp2(seq(min(slog), max(slog), length = 6), c("white","bisque2", "burlywood2", "coral1", "brown3", "brown4"), 
                      space = "RGB")

#check the hexa codes of the colors you want
brewer.pal(n = 8, name = "YlOrRd")
col_fun3 = colorRamp2(seq(min(slog1), max(slog1), length = 8), c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C" ,"#FC4E2A", "#E31A1C", "#B10026"))

#rowOrder <- rownames(slog1)
rowOrder <- c('VP7', 'VP4', 'VP6', 'VP1', 'VP2', 'VP3', 'NSP1', 'NSP2', 'NSP3', 'NSP4', 'NSP5',
                    'VP7.1', 'VP4.1', 'VP6.1', 'VP1.1', 'VP2.1', 'VP3.1', 'NSp1', 'NSp2', 'NSp3', 'NSp4', 'NSp5',
                    'VP7.2', 'VP4.2', 'VP6.2', 'VP1.2', 'VP2.2', 'VP3.2', 'NSP1.2', 'NSP2.2', 'NSP3.2', 'NSP4.2', 'NSP5.2')

#Add qc reads as barplot for the column annotation
col_ha <- HeatmapAnnotation(Total = anno_barplot(slog[34,]))

slog1 <- slog[-34,]

#VIRIDIS - magma
col_fun4 = colorRamp2(seq(min(slog1), max(slog1), length = 7), c("#FFFFCC", '#FCDF96', '#FA8657', '#E03B50', '#981D69', '#57096E', '#1E0848'))

## Vector to split the heatmap by rows
group1 <- c(rt, tb, au)
rowSplit <- group1

logAbundance=Heatmap(slog1, name='logAbundance', col = col_fun4, row_order = rowOrder, show_column_dend = F,
        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
        top_annotation = col_ha, row_names_side = "left", border=T, row_split = rev(rowSplit))

logAbundance

# Add seperating lines between genotype constellations
decorate_heatmap_body('logAbundance', {grid.lines(x=c(-0.2,0.99), y=c(0.33,0.33), gp = gpar(lty = 2, lwd = 2))})
decorate_heatmap_body('logAbundance', {grid.lines(x=c(-0.2,0.99), y=c(0.665,0.665), gp = gpar(lty = 2, lwd = 2))})

dev.off()

#deneme for splitting
slog1
slog11 = slog1[12:22,]
slog12 = slog1[23:33,]
slog13 = slog1[1:11,]

################

########## RVAB #########
#son1 <- read.csv('~/path/to/new_bbmap_RVAB_abundance.txt', sep = "\t")
son1 <- read.csv('~/paht/to/new_bbmap_RVAB_abundance1.csv')
son1 <- read.csv('RVAB_abundance.bwa_final.csv')

rownames(son1) <- son1[,1]
son1 <- son1[, c(-1)]

d <- read_excel("~/paht/to/RVAB_NGS1-2_readStats.xlsx", sheet = 3)
son <- son1[,c(-1, -104)]  #removed the length column for now


#Add the QC reads per sample at the end of each column and divide the read counts with them to get the 
#relative abundances

son[nrow(son)+1, ] <- t(d[,2])
rownames(son)[55] <- 'QC_Reads'

sonm <- as.matrix(son)

remove <- c('F02133', 'F05339', 'F09229', 'F09263', 'F09649')
sonx <- t(sonm)
sonx <- sonx[!rownames(sonx) %in% remove, ]
sonm <- t(sonx)

#read count cut-off - selected 100
sonm[sonm<1000 ] <- 0
sonm[sonm<500 ] <- 0
sonm[sonm<100 ] <- 0


sonm <- apply(sonm, 2, function(x){x/max(x)*100000})
sonm <- sonm[-55,]

#these are also removed automatically when column ordering by G type
#remove <- c('F03127', 'F05040', 'F06083', 'F09179')
remove <- c('F03127', 'F05040', 'F06083', 'F09179', 'F02133', 'F05339', 'F09229', 'F09263', 'F09649')
sonx <- t(sonm)
sonx <- sonx[!rownames(sonx) %in% remove, ]
sonm <- t(sonx)

low_ones <- sonx[c('F05339', 'F02133','F06184','F07027', 'F09155', 'F09229', 'F09263', 'F09649'), ]
low <- t(low_ones)
low1 <- t(low_ones)



#Normalization, tried scaling because that is preferred for not changing the data but changing the magnitudes
#but gives negative values

slog <- log10(sonm)
#min value is determined by 
min(slog[slog!=min(slog)])
slog[which(!is.finite(slog))] <- 0
slog[slog < 0] <- 0


#Seperate the columns (samples) according to G type - bwa
#####
gs <- as.data.frame(t(slog))
gsG1 <- subset(gs, gs[,16] > 1.8)
gsG2 <- subset(gs, gs[,18] > 0.8 )
gsG3 <- subset(gs, gs[,19] > 1.7 | gs[,20] > 2)
gsG4 <- subset(gs, gs[,21] > 0.7 )
gsG6_G8 <- subset(gs, gs[,22] > 3.5 )
gsG9 <- subset(gs, gs[,24] > 0.5 )
gsG12 <- subset(gs, gs[,17] > 1.3 )


ss <- rbind(gsG1, gsG2, gsG3, gsG4, gsG6_G8, gsG9, gsG12)
ss1 <- as.matrix(t(ss))
colOrder <- colnames(ss1)
########


#Seperate the columns (samples) according to G type - bbmap 0.95
gs <- as.data.frame(t(slog))
gsG1 <- subset(gs, gs[,16] > 1.8)
gsG2 <- subset(gs, gs[,18] > 0.7 )
gsG3 <- subset(gs, gs[,19] > 1.7 | gs[,20] > 2)
gsG4 <- subset(gs, gs[,21] > 0.69 )
gsG6_G8 <- subset(gs, gs[,22] > 3.5 )
gsG9 <- subset(gs, gs[,24] > 0.5 )
gsG12 <- subset(gs, gs[,17] > 1.3 )

#add <- gs[c('F06107', 'F09317'), ]
ss <- rbind(gsG1, gsG2, gsG3, gsG4, gsG6_G8, gsG9, gsG12)
ss1 <- as.matrix(t(ss))
colOrder <- colnames(ss1)

rownames(ss1)[16] <- 'G1'

#Add the trimmed read counts per sample again
#no
slog <- rbind(ss1, t(d[,2]))

#Scale
#sonm_scaled = t(apply(sonm, 1, scale))
#son_Scaled <- apply(son, 2, scale)


#Annotate the heatmap
#Column annotations are total # of trimmed reads/sample
c <- slog[45,]
slog1 <- slog[-45,]

#Row annotations are total read counts per segment/gene length - NOT DONE YET

#row_ha = rowAnnotation(Total = anno_barplot(d), show_annotation_name = c(Total = FALSE), gap = unit(1.5, "mm"))


col_ha <- HeatmapAnnotation(Total = anno_barplot(c), gap = unit(1.5, "mm"))
#red-yellow
col_fun3 = colorRamp2(seq(min(slog1), max(slog1), length = 8), c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C" ,"#FC4E2A", "#E31A1C", "#B10026"))
#viridis - magma? - purple&red
col_fun4 = colorRamp2(seq(min(ss1), max(ss1), length = 7), c("#FFFFCC", '#FCDF96', '#FA8657', '#E03B50', '#981D69', '#57096E', '#1E0848'))
#viridis - viridis - green&blue
col_fun4 <- colorRamp2(seq(min(ss1), max(ss1), length = 7), rev(c("#440154FF", "#443A83FF" ,"#31688EFF" ,"#21908CFF", "#35B779FF", "#8FD744FF", '#FFFFCC')))


#Heatmap with seriation but the genotypes are not all together still
#cluster_rows = as.dendrogram(o1[[1]]),
# top_annotation = col_ha
Heatmap(sonm, name='Abundances', col = col_fun3, row_order = get_order(o1), column_order = get_order(o2), 
        row_km = 3 , row_title = c('other', 'DS-1', 'Wa'))


#Genotype constellation row order
###Final!!!!

rowOrder <- c('G3', 'G12', 'G4', 'G9', 'G1','P8', 'I1', 'R1','R1-M0094', 'C1', 'M1', 'A1', 'N1', 'N1-SS454', 'N1-JES11', 'N1-Nov09',
              'N1-MRC', 'T1', 'T1.1', 'E1', 'E1-ESP455','E1-LL', 'H1', 'G3.1', 'P6', 'G2', 'P4', 'I2', 'R2', 'R2-557' ,'C2', 'M2', 
              'M2-3203WC', 'A2', 'N2', 'T2', 'E2','E2-MWI', 'E2-MRC', 'E2-RV09', 'H2', 'G8', 'G6', 'P14', 'I2-b', 'R2-b', 'C2-b', 
              'M2-b', 'A11', 'A3', 'N2-b','T6', 'E2-b', 'H3')




#save pdf 
logAbundance <- Heatmap(ss1, name='logAbundance', col = col_fun4, show_parent_dend_line = FALSE,
                     #column_km = 2, column_title = c('Wa-like', 'DS-like'),
                     row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                     border=T, row_order=rowOrder, row_names_side = "left",
                     column_order = colOrder)

#top_annotation = col_ha
logAbundance







#Extract the clustered column order from the heatmap so that we can  change
#the place of F03011 - or others if necessary 
#You have to first have a heatmap to extract the row order, and then add the colOrder vector back
#to the heatmap and run it again

ht <- draw(logAbundance)
ord <- as.vector(column_order(ht))

samples <- t(t(colnames(slog1)))
rownames(samples) <- seq(1, 102, by=1) 

ab <- as.matrix(c(ord$`1`,ord$`2`))
colOrder <- as.matrix(samples[order(match(rownames(samples), ab[,1]))])

#Add F03011 to the end
colOrd <- c(colOrder[c(1:11, 13:102),], 'F03011')

logAbundance <- Heatmap(slog1, name='logAbundance', col = col_fun4, show_parent_dend_line = FALSE,
                        #column_km = 2, column_title = c('Wa-like', 'DS-like'),
                        row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), column_names_rot = 60,
                        border=T, row_order=rowOrder, row_names_side = "left", top_annotation = col_ha,
                        column_order = colOrd)
                        
logAbundance


# Add seperating lines between Wa and DS-1 like genotype constellations
decorate_heatmap_body('logAbundance', {grid.lines(x=c(0,1), y=c(0.27,0.27), gp = gpar(lty = 5, lwd = 3, col='darkcyan'))})
decorate_heatmap_body('logAbundance', {grid.lines(x=c(0,1), y=c(0.615,0.615), gp = gpar(lty = 5, lwd = 3, col='darkcyan'))})

dev.off()



#Subset heatmaps - olmadı sonra 
split_1_2_3 <- c("first","first","first", "first","first","first","first", "first", "first", "first", "first","first","first","first","first","first","first", 
                              "second","second","second","second","second","second","second","second","second","second","second","second","second","second",
                              "third", "third", "third", "third", "third", "third", "third", "third", "third", "third", "third", "third", "third")

#Add to heatmap fnctn: row_split = split_1_2_3, row_title = NULL, row_gap = unit(1.5, "mm")


ab[1:17,]
ab[18:32,]
ab[33:44,]




genes <- t(t(rownames(slog1)))
rownames(genes) <- seq(1, 44, by=1) 
ht <- draw(logAbundance)
ord <- as.vector(row_order(ht))
ab <- as.matrix(c(ord$first,ord$second, ord$third))

roooowOrder <- as.matrix(genes[order(match(genes[,1], ab[,1]))])







logAbundance[1:17,]
logAbundance[18:32,]
logAbundance[33:44,]



#!/usr/bin/Rscript

#Plotting the map of Belgium and showing the locations of VD sample collection

library(raster)
library(mapproj)
library(scales)
library(clipr)
library(sf)

#post codes are copied/read from excel
#postcds <- read_excel('postcodes_VD_samples.xlsx')
postcds <- read_clip_tbl()
vd <- rep('VD', 75)
post <- cbind(postcds, vd)
colnames(post) <- c('Post_code', 'Group')

# read all commune post codes of Belgium
setwd("~/Documents/Bioinformatics/R_Stuff")
all_commune <- read.csv('bel_postcodes.csv')
all_commune <- all_commune[, c(1,3,4)]
all_commune <- all_commune[complete.cases(all_commune), ]

rotacord <- merge(post, all_commune, by="Post_code",all.x=F, all.y=F)

#get all unique rows (communes) to get the coordinates
rotacord_dedup <- rotacord[!duplicated(rotacord$Post_code), ]
rotacord_dedup <- rotacord_dedup[, -c(2)]

# plot the map of Belgium
bel <- getData('GADM', country='BEL', level=3)
plot(bel, col="lightgrey", border="white", axes=F)
points(rotacord_dedup$lon, rotacord_dedup$lat, pch = 16, col = alpha("navyblue", 0.7))
#box()
#Image deconvolution analysis
setwd("Your Directory for Images of X Stain")  
#Load Libraries
library(EBImage)
library(CRImage)
library(stringr)
library(reticulate)
#Import python libraries
ski<- import('skimage')
ski_color <- import('skimage.color')

#Set Paths for data output (Please Adjust These Based on Your Directory) **Please note this changes based on the stain***
segmentpath <- "~/Your Path/stain/segmentations"
featurespath <- "~/Your Path/stain/features"
deconpath <- "~/Your Path/stain/decon/"
overlaypath <- "~/Your Path/stain/overlay/"

#Read in images
my.list <- list.files(pattern = ".png") #Please email authors is you require image tiles
#Note: some images were saves as tif thus you may need to run another deconvolution on that set

###Before Deconvoluting Make Sure you run the mask_to_segmented functions
#Analysis
for (c in my.list){
  mask_to_segmentated(c)
}

#Function for png images
mask_to_segmentated <- function(i){
  print(i)
  #Read in image
  img <- readImage(i)
  oligo <- img[,,1:3]
  
  #Split the H and DAB stains
  ihc_hed <- ski_color$rgb2hed(oligo)
  null <- matrix(0, dim(ihc_hed)[1], dim(ihc_hed)[2])
  ihc_h <- ski_color$hed2rgb(abind(ihc_hed[,,1], null, null, along = 3))
  ihc_e <- ski_color$hed2rgb(abind(ihc_hed[,,2], null, null, along = 3))
  
  #Save
  im <- ihc_h
  filename = paste0(deconpath,tools::file_path_sans_ext(i), "_ihc_h.png")
  writeImage(im[,,1], files = paste0(filename))
  
  im2 <- ihc_e
  filename = paste0(deconpath,tools::file_path_sans_ext(i), "_ihc_e.png")
  writeImage(im2[,,1], files = paste0(filename))
  
  #Segment
  test_median <- medianFilter((1-im[,,1]), 2)
  test_median <- test_median>otsu(test_median)
  
  test_median2 <- medianFilter((1-im2[,,1]), 2)
  test_median2 <- test_median2>otsu(test_median2)
  #display(test_median2) #I need to add a model noise filtration step
  
  #Save
  filename = paste0(tools::file_path_sans_ext(i), "_segment_h.png")
  writeImage(test_median, files = paste0(segmentpath,"/", filename))
  
  filename = paste0(tools::file_path_sans_ext(i), "_segment_e.png")
  writeImage(test_median2, files = paste0(segmentpath,"/", filename))
  
  #Compute features
  R <- oligo[,,1]
  G <- oligo[,,2]
  B <- oligo[,,3]
  
  #Correct into masks
  test_median <- bwlabel(test_median)
  test_median2 <- bwlabel(test_median2)
  
  #extraction of features from the objects (H Layer)
  crenabled_features_Rh <- computeFeatures(x = test_median, ref = R, xname = "R")
  crenabled_features_Gh <- computeFeatures(x = test_median, ref = G, xname = "G")
  crenabled_features_Bh <- computeFeatures(x = test_median, ref = B, xname = "B")
  
  crenabled_features_Rh <- as.data.frame(crenabled_features_Rh)
  crenabled_features_Gh <- as.data.frame(crenabled_features_Gh)
  crenabled_features_Bh <- as.data.frame(crenabled_features_Bh)
  
  final_df <- cbind(crenabled_features_Rh, crenabled_features_Gh[,3:89], crenabled_features_Bh[,3:89])
  segmented1 <- paintObjects(test_median, oligo, col = 'red', thick = TRUE)
  filename = paste0(tools::file_path_sans_ext(i), "_overlay_h.png")
  writeImage(segmented1, files = paste0(overlaypath,"/", filename))
  
  filename = paste0(tools::file_path_sans_ext(i), "_features_h.csv")
  write.csv(final_df, file = paste0(featurespath, "/", filename))
  
  #extraction of features from the objects (H Layer)
  crenabled_features_Rd <- computeFeatures(x = test_median2, ref = R, xname = "R")
  crenabled_features_Gd <- computeFeatures(x = test_median2, ref = G, xname = "G")
  crenabled_features_Bd <- computeFeatures(x = test_median2, ref = B, xname = "B")
  
  crenabled_features_Rd <- as.data.frame(crenabled_features_Rd)
  crenabled_features_Gd <- as.data.frame(crenabled_features_Gd)
  crenabled_features_Bd <- as.data.frame(crenabled_features_Bd)
  
  final_df <- cbind(crenabled_features_Rd, crenabled_features_Gd[,3:89], crenabled_features_Bd[,3:89])
  segmented2 <- paintObjects(test_median2, oligo, col = 'red', thick = TRUE)
  filename = paste0(tools::file_path_sans_ext(i), "_overlay_e.png")
  writeImage(segmented2, files = paste0(overlaypath,"/", filename))
  
  filename = paste0(tools::file_path_sans_ext(i), "_features_e.csv")
  write.csv(final_df, file = paste0(featurespath, "/", filename))
}

#Function for tif images
mask_to_segmentated_tif <- function(i){
  print(i)
  #Read in image
  img <- readTIFF(i)
  oligo <- img[,,1:3]
  
  #Split the H and DAB stains
  ihc_hed <- ski_color$rgb2hed(oligo)
  null <- matrix(0, dim(ihc_hed)[1], dim(ihc_hed)[2])
  ihc_h <- ski_color$hed2rgb(abind(ihc_hed[,,1], null, null, along = 3))
  ihc_d <- ski_color$hed2rgb(abind(ihc_hed[,,3], null, null, along = 3))
  
  #Save
  im <- ihc_h
  filename = paste0(deconpath,tools::file_path_sans_ext(i), "_ihc_h.png")
  writeImage(im[,,1], files = paste0(filename))
  
  im2 <- ihc_d
  filename = paste0(deconpath,tools::file_path_sans_ext(i), "_ihc_d.png")
  writeImage(im2[,,1], files = paste0(filename))
  
  #Segment
  test_median <- medianFilter((1-im[,,1]), 2)
  test_median <- test_median>otsu(test_median)
  
  test_median2 <- medianFilter((1-im2[,,1]), 2)
  test_median2 <- test_median2>otsu(test_median2)
  #display(test_median2) #I need to add a model noise filtration step
  
  #Save
  filename = paste0(tools::file_path_sans_ext(i), "_segment_h.png")
  writeImage(test_median, files = paste0(segmentpath,"/", filename))
  
  filename = paste0(tools::file_path_sans_ext(i), "_segment_d.png")
  writeImage(test_median2, files = paste0(segmentpath,"/", filename))
  
  #Compute features
  R <- oligo[,,1]
  G <- oligo[,,2]
  B <- oligo[,,3]
  
  #Correct into masks
  test_median <- bwlabel(test_median)
  test_median2 <- bwlabel(test_median2)
  
  #extraction of features from the objects (H Layer)
  crenabled_features_Rh <- computeFeatures(x = test_median, ref = R, xname = "R")
  crenabled_features_Gh <- computeFeatures(x = test_median, ref = G, xname = "G")
  crenabled_features_Bh <- computeFeatures(x = test_median, ref = B, xname = "B")
  
  crenabled_features_Rh <- as.data.frame(crenabled_features_Rh)
  crenabled_features_Gh <- as.data.frame(crenabled_features_Gh)
  crenabled_features_Bh <- as.data.frame(crenabled_features_Bh)
  
  final_df <- cbind(crenabled_features_Rh, crenabled_features_Gh[,3:89], crenabled_features_Bh[,3:89])
  filename = paste0(tools::file_path_sans_ext(i), "_features_h.csv")
  write.csv(final_df, file = paste0(featurespath, "/", filename))
  
  #extraction of features from the objects (H Layer)
  crenabled_features_Rd <- computeFeatures(x = test_median2, ref = R, xname = "R")
  crenabled_features_Gd <- computeFeatures(x = test_median2, ref = G, xname = "G")
  crenabled_features_Bd <- computeFeatures(x = test_median2, ref = B, xname = "B")
  
  crenabled_features_Rd <- as.data.frame(crenabled_features_Rd)
  crenabled_features_Gd <- as.data.frame(crenabled_features_Gd)
  crenabled_features_Bd <- as.data.frame(crenabled_features_Bd)
  
  final_df <- cbind(crenabled_features_Rd, crenabled_features_Gd[,3:89], crenabled_features_Bd[,3:89])
  filename = paste0(tools::file_path_sans_ext(i), "_features_d.csv")
  write.csv(final_df, file = paste0(featurespath, "/", filename))
}




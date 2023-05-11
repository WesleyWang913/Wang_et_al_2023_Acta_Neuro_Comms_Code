library(EBImage)

#Set path to unet segmented masks
setwd("Your Directory")

#Pull file list
p_files <- list.files(path = "UNET_Masks_Processes")
length(p_files) #5809
b_files <- list.files(path = "UNET_Masks_Bodies")
length(b_files) #5809

#We are making a function that reads in both mask txt files, converts to image, and then merges them
txt_to_mask <- function(n){
  mask1 <- read.delim(paste0("UNET_Masks_Processes/", p_files[n]), sep = " ", header = FALSE)
  mask2 <- read.delim(paste0("UNET_Masks_Bodies/", b_files[n]), sep = " ", header = FALSE)
  mask3 <- mask1 + mask2
  mask3 <- as.Image(as.matrix(mask3))
  mask3 = flip(rotate(mask3, -90))
  writeImage(mask3, paste0("Merged_masks/", gsub("\\..*","",p_files[n]) ,".png"))
}

#Convert masks
for (i in 1:length(p_files)){
  print(i)
  txt_to_mask(i)
}


library(ANTsR)
tfn = "MyImage.nii.gz"
img <- antsImageRead(tfn, 3, "float")
imgn3 <- antsImageClone(img)

##N3BiasFieldCorrection(list(img@dimension, img, imgn3, "4"))
##N3BiasFieldCorrection(list(img@dimension, imgn3, imgn3, "2"))

mask <- antsImageClone(img, "unsigned int")
threshold = 1.e-6
mask[img > threshold ] <- 1
mask[img <= threshold ] <- 0
antsImageWrite(mask,"mask.nii.gz")
segs <- Atropos(d = 3, a = imgn3, m = "[0.2,1x1x1]", c = "[5,0]", i = "kmeans[3]",
    x = mask)

antsImageWrite(segs$segmentation,"segmentation.nii.gz")




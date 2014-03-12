
# STAPLE of prostate peripheral zone as in Warfield, Zou, Wells IEEE TMI 2004.
# Data from Figure 7.
crlSTAPLE ../peripheral/seg001.nrrd ../peripheral/seg002.nrrd \
 ../peripheral/seg003.nrrd ../peripheral/seg004.nrrd \
 ../peripheral/seg005.nrrd \
  --compressOutput --outputImage weights.nrrd

crlIndexOfMaxComponent weights.nrrd maxlikelihoodclassification.nrrd

crlExtractSmallerImageFromImage -i weights.nrrd \
 -o peripheralzoneprobability.nrrd -l 1 -m 2 -a 2

# To visualize the result:
# crlViz -s maxlikelihoodclassification.nrrd -p peripheralzoneprobability.nrrd ../mri/prostatemri.mhd


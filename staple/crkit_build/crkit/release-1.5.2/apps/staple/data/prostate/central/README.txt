
# 2D segmentations of the prostate central zone
for i in ../combined/*.mhd; do
   crlRelabelImages $i $i 4 1 `basename $i .mhd`.nrrd 0
done


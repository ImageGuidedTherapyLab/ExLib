
# 2D segmentations of the prostate peripheral zone
for i in ../combined/*.mhd; do
   crlRelabelImages $i $i 5 1 `basename $i .mhd`.nrrd 0
done


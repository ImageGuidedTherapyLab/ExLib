
# Data for : Figure 3 Warfield, Zou, Wells IEEE TMI 2004

crlSTAPLE -c --outputImage weights.nrrd \
  sqc10e-001.nrrd sqc10e-002.nrrd sqc10e-003.nrrd

crlIndexOfMaxComponent weights.nrrd maxlikelihoodclassification.nrrd
crlExtractSmallerImageFromImage -i weights.nrrd -o foregroundprobability.nrrd \
  -l 6 -m 7 -a 3


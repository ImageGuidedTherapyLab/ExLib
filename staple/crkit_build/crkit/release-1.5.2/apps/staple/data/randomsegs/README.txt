
# Create a set of random segmentations from a known ground truth but
# specified error rates.
#  Run crlSTAPLE and estimate the reference standard and performance
#  parameters.

# Test 10 random images with p=0.80, q=0.75
# Test 10 random images with p=0.95, q=0.90
for i in 1 2 3 4 5 6 7 8 9 10; do
crlRandomSegmentationImageGenerator 256 256 1 0.95 0.90 \
   exp2${i}.nrrd
done

crlSTAPLE --outputImage weights.nrrd \
  --stationaryPrior "0.5 0.5" \
  exp21.nrrd exp22.nrrd exp23.nrrd exp24.nrrd exp25.nrrd \
  exp26.nrrd exp27.nrrd exp28.nrrd exp29.nrrd exp210.nrrd 
crlIndexOfMaxComponent weights.nrrd maxlikelihoodclassification.nrrd
crlExtractSmallerImageFromImage -i weights.nrrd -o foregroundprobability.nrrd \
  -l 1 -m 2 -a 3

crlMeanFieldMRF weights.nrrd automatic 0.00001 1 5 mfweights.nrrd
crlIndexOfMaxComponent mfweights.nrrd mfmlc.nrrd
crlViz -p foregroundprobability.nrrd -s maxlikelihoodclassification.nrrd \
 -s mfmlc.nrrd



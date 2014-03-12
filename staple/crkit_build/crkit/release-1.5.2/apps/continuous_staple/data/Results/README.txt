for ((i=1;${i}<25;i++))
do
	crlCastSymMatDoubleToFloat ../images/Tensors_${i}.nrrd ../images/TensorsFloat_${i}.nrrd
	crlTensorLog -i ../images/TensorsFloat_${i}.nrrd -o ../logimages/logTensors_${i}.nrrd
done

itkContinuousSTAPLE -p 4 -i 100 -t 0 -v 1 -b 0 -m ComputationMask.nrrd -o yourLogGroundTruth_STAPLE.nrrd -l listLogs.txt -O yourResultingParams_STAPLE.m

crlCastSymMatDoubleToFloat yourLogGroundTruth_STAPLE.nrrd yourLogGroundTruthFloat_STAPLE.nrrd
crlTensorExp yourLogGroundTruthFloat_STAPLE.nrrd yourGroundTruth_STAPLE.nrrd

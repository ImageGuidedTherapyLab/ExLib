
# Need to add some test data so as to be able to verify correct operation.

castsymmattodouble test-tensors.nrrd test-tensors-double.nrrd

tensorclean test-tensors-double.nrrd test-tensors-double-clean.nrrd
tensorlog test-tensors-double-clean.nrrd test-tensors-double-clean-log.nrrd
tensorexp test-tensors-double-clean.nrrd test-tensors-double-clean-exp.nrrd
tensorexp test-tensors-double-clean-log.nrrd test-tensors-double-clean-log-exp.nrrd

castsymmattofloat test-tensors-double-clean-log-exp.nrrd /tmp/t1-log-exp.nrrd
castsymmattofloat test-tensors-double-clean-log.nrrd /tmp/t2-log.nrrd



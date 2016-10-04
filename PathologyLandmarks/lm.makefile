ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 

#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# apply transformation
pathology.lmreg.nii.gz: landmarktransform.tfm pathology.nii.gz dce.nii.gz
	 $(ANTSAPPLYTRANSFORMSCMD) -d 3 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun itksnap -g $(word 3, $^)  -o $@

# compute lm transformation
landmarktransform.tfm: dceLM.nii.gz pathologyLM.nii.gz 
	 $(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert path to nifti and set spacing
pathology.nii.gz: pathology.tif 
	c3d -verbose -mcs $< -origin 0x0x0mm -spacing 0.01x0.01x0.01mm -omc $@

view:
	vglrun itksnap -g pathology.nii.gz -s pathologyLM.nii.gz
	vglrun itksnap -g dce.nii.gz -s dceLM.nii.gz
	echo Layer Inspector -> General -> Display Mode -> RGB

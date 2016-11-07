ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 
ITKSNAP  = vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
MYSQL = mysql 
CONVERTSVS = python ./convertsvs.py
DIMENSION = 3
ATROPOSCMD=$(ANTSPATH)/Atropos -d $(DIMENSION)  -c [3,0.0] -m [0.1,1x1x1] 
NTISSUE = 4

################
# Dependencies #
################
-include datalocation/dependencies
dependencies: 
	bash -i -c '$(MYSQL) --local-infile < ./hccpathdb.sql'
	bash -i -c '$(MYSQL) -sNre "call Metadata.HCCPathDBList();"  > datalocation/dependencies'

HENIFTI:=  $(subst PathHE.svs,PathHE.nii.gz,$(PathologyHE))
PIMONIFTI:=  $(subst PathPIMO.svs,PathPIMO.nii.gz,$(PathologyPimo))
HEGMM:=  $(subst PathHE.svs,PathHE.gmm.nii.gz,$(PathologyHE))
PIMOGMM:=  $(subst PathPIMO.svs,PathPIMO.gmm.nii.gz,$(PathologyPimo))
convert:    $(HENIFTI) $(PIMONIFTI)
gmm:    $(HEGMM) $(PIMOGMM)
lm:         $(T2LM)
transform:  $(UpdateTransform)

# debug
jobs:
	echo $(T2LM)

#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# initialize LM
%/RefImgLM.nii.gz: %/RefImg.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/PathHELM.nii.gz: %/PathHE.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/PathPIMOLM.nii.gz: %/PathPIMO.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 

# apply transformation
pathology.lmreg.nii.gz: landmarktransform.tfm pathology.nii.gz dce.nii.gz
	 $(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@

# update transform lm
%/updatetransform: %/Pathology/PathHELM.nii.gz %/Pathology/PathHE.nii.gz %/T2wReference/RefImgLM.nii.gz %/T2wReference/RefImg.hdr %/Pathology/PathPIMOLM.nii.gz %/Pathology/PathPIMO.nii.gz 
	-$(ITKSNAP) -l LMLabels.txt -s $(word 1, $^)  -g $(word 2, $^)  &  PIDPATH=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 3, $^)  -g $(word 4, $^)  &  PIDMRI=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 5, $^)  -g $(word 6, $^)  &  PIDPIMO=$$!; \
        zenity --info --title="OutputFile" --text="Tools -> Layer Inspector -> General -> Display Mode -> RGB $(word 4, $^)  "; \
        kill -9 $$PIDPATH;kill -9 $$PIDMRI;kill -9 $$PIDPIMO;

# compute lm transformation
%/t2HElmtransform.tfm:   %/T2wReference/RefImgLM.nii.gz %/Pathology/PathHELM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
%/t2PIMOlmtransform.tfm: %/T2wReference/RefImgLM.nii.gz %/Pathology/PathPIMOLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
# FIXME - need script to parse header
%/PathHE.nii.gz: %/PathHE.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@
%/PathPIMO.nii.gz: %/PathPIMO.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@

%/PathHE.mask.nii.gz: %/PathHE.nii.gz
	$(C3DEXE) $<  -threshold 0% 40% 1 0 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

%/PathHE.gmm.nii.gz: %/PathHE.nii.gz %/PathHE.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $*/PathHEred.nii.gz $*/PathHEgreen.nii.gz $*/PathHEblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUE)] -x $(word 2,$^) -a $*/PathHEred.nii.gz -a $*/PathHEgreen.nii.gz -a $*/PathHEblue.nii.gz   -o [$@,$(WORKDIR)/$*.gmmPOSTERIORS%d.nii.gz] 

view:
	vglrun itksnap -g pathology.nii.gz -s pathologyLM.nii.gz
	vglrun itksnap -g dce.nii.gz -s dceLM.nii.gz
	echo Tools -> Layer Inspector -> General -> Display Mode -> RGB

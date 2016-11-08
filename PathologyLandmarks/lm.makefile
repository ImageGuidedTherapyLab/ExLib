ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 
ITKSNAP  = vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
MYSQL = mysql 
CONVERTSVS = python ./convertsvs.py
DIMENSION = 3
NTISSUE = 4
OTBOFFSET  = 30
OTBRADIUS  = 50
ATROPOSCMD=$(ANTSPATH)/Atropos -d $(DIMENSION)  -c [3,0.0] -m [0.1,1x1x1] 
OTBTEXTURE=/rsrch2/ip/dtfuentes/github/ExLib/otbScalarImageTextures/otbScalarImageToTexturesFilter
CLUSTERDIR=/rsrch2/ip/dtfuentes/github/ExLib/PathologyLandmarks/Processed
WORKDIR=Processed

################
# Dependencies #
################
-include datalocation/dependencies
dependencies: 
	bash -i -c '$(MYSQL) --local-infile < ./hccpathdb.sql'
	bash -i -c '$(MYSQL) -sNre "call Metadata.HCCPathDBList();"  > datalocation/dependencies'

HENIFTI:=  $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.nii.gz,$(PathologyHE)))
PIMONIFTI:=$(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.nii.gz,$(PathologyPimo)))
HEGMM:=    $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.gmm.nii.gz,$(PathologyHE)))
PIMOGMM:=  $(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.gmm.nii.gz,$(PathologyPimo)))
HEOTB:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE120.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyHE)))
convert:   $(HENIFTI) $(PIMONIFTI)
gmm:       $(HEGMM) $(PIMOGMM)
lm:        $(T2LM)
transform: $(UpdateTransform)

# debug
jobs:
	echo $(T2LM)

#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# initialize LM
%/RefImgLM.nii.gz: %/RefImg.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/T1postLM.nii.gz: %/T1post.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/DCEavgLM.nii.gz: %/DCEavg.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/PathHELM.nii.gz: %/PathHE.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 
%/PathPIMOLM.nii.gz: %/PathPIMO.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 

# apply transformation
%/Pathology/PathHE.lmreg.nii.gz: %/t2HElmtransform.tfm %/Pathology/PathHE.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/Pathology/PathHE.gmm.lmreg.nii.gz: %/t2HElmtransform.tfm %/Pathology/PathHE.gmm.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/Pathology/PathPIMO.lmreg.nii.gz: %/t2PIMOlmtransform.tfm %/Pathology/PathPIMO.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/Pathology/PathPIMO.gmm.lmreg.nii.gz: %/t2PIMOlmtransform.tfm %/Pathology/PathPIMO.gmm.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/DCE/DCEavg.lmreg.nii.gz: %/t2dcelmtransform.tfm %/DCE/DCEavg.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/DCE/DCE.lmreg.nii.gz: %/t2dcelmtransform.tfm %/DCE/DCE.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/T1Post/T1post.lmreg.nii.gz: %/t2t1lmtransform.tfm %/T1Post/T1post.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@
%/T1Pre/T1pre.lmreg.nii.gz: %/t2t1lmtransform.tfm %/T1Pre/T1pre.nii.gz %/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun $(ITKSNAP) -g $(word 3, $^)  -o $@

# update transform lm
%/updatetransform: %/Pathology/PathHELM.nii.gz %/Pathology/PathHE.nii.gz %/T2wReference/RefImgLM.nii.gz %/T2wReference/RefImg.hdr %/Pathology/PathPIMOLM.nii.gz %/Pathology/PathPIMO.nii.gz %/T1post/T1postLM.nii.gz %/T1post/T1post.hdr %/T1pre/T1pre.hdr %/DCE/DCEavgLM.nii.gz %/DCE/DCEavg.hdr  %/DCE/DCE.hdr 
	-$(ITKSNAP) -l LMLabels.txt -s $(word 1, $^)  -g $(word 2, $^)  &  PIDPATH=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 3, $^)  -g $(word 4, $^)  &  PIDMRI=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 5, $^)  -g $(word 6, $^)  &  PIDPIMO=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 7, $^)  -g $(word 8, $^)  -o $(word 9, $^)  &  PIDT1=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 10, $^) -g $(word 11, $^) -o $(word 12, $^) &  PIDDCE=$$!; \
        zenity --info --title="OutputFile" --text="Tools -> Layer Inspector -> General -> Display Mode -> RGB $(word 4, $^)  "; \
        kill -9 $$PIDPATH;kill -9 $$PIDMRI;kill -9 $$PIDPIMO; kill -9 $$PIDT1; kill -9 $$PIDDCE

# compute lm transformation
%/t2HElmtransform.tfm:   %/T2wReference/RefImgLM.nii.gz %/Pathology/PathHELM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
%/t2PIMOlmtransform.tfm: %/T2wReference/RefImgLM.nii.gz %/Pathology/PathPIMOLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
%/t2t1lmtransform.tfm: %/T2wReference/RefImgLM.nii.gz %/T1post/T1postLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
%/t2dcelmtransform.tfm: %/T2wReference/RefImgLM.nii.gz %/DCE/DCEavgLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
%/PathHE.nii.gz: %/PathHE.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@
%/PathPIMO.nii.gz: %/PathPIMO.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@

%/PathHE.mask.nii.gz: %/PathHE.nii.gz
	$(C3DEXE) $<  -threshold 0% 40% 1 0 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

%/PathPIMO.mask.nii.gz: %/PathPIMO.nii.gz
	$(C3DEXE) $<  -threshold 0% 40% 1 0 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

%/PathHE.gmm.nii.gz: %/PathHE.nii.gz %/PathHE.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $*/PathHEred.nii.gz $*/PathHEgreen.nii.gz $*/PathHEblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUE)] -x $(word 2,$^) -a $*/PathHEred.nii.gz -a $*/PathHEgreen.nii.gz -a $*/PathHEblue.nii.gz   -o [$@,$*/PathHEgmmPOSTERIORS%d.nii.gz] 
	echo $(ITKSNAP) -s  $@ -g $<

# FIXME - push to cluster
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE000.   3 50     0            0            0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE001.   3 50     0            0        $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE010.   3 50     0        $(OTBOFFSET)     0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE011.   3 50     0        $(OTBOFFSET) $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE100.   3 50 $(OTBOFFSET)     0            0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE101.   3 50 $(OTBOFFSET)     0        $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE110.   3 50 $(OTBOFFSET) $(OTBOFFSET)     0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE111.   3 50 $(OTBOFFSET) $(OTBOFFSET) $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
$(WORKDIR)/%/PathHE120.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(DATADIR)/%/PathHE.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	$(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5
	$(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	$(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	$(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	$(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathHE.nii.gz -o $(WORKDIR)/$*/PathHE120.Entropy_$(OTBRADIUS).nii.gz

%/PathPIMO.gmm.nii.gz: %/PathPIMO.nii.gz %/PathPIMO.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $*/PathPIMOred.nii.gz $*/PathPIMOgreen.nii.gz $*/PathPIMOblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUE)] -x $(word 2,$^) -a $*/PathPIMOred.nii.gz -a $*/PathPIMOgreen.nii.gz -a $*/PathPIMOblue.nii.gz   -o [$@,$*/PathPIMOgmmPOSTERIORS%d.nii.gz] 
	echo $(ITKSNAP) -s  $@ -g $<

view:
	vglrun itksnap -g pathology.nii.gz -s pathologyLM.nii.gz
	vglrun itksnap -g dce.nii.gz -s dceLM.nii.gz
	echo Tools -> Layer Inspector -> General -> Display Mode -> RGB

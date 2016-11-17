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
HELMREG:=  $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.lmreg.nii.gz,$(PathologyHE)))
PIMOLMREG:=$(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.lmreg.nii.gz,$(PathologyPimo)))
HEGMM:=    $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.gmm.nii.gz,$(PathologyHE)))
PIMOGMM:=  $(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.gmm.nii.gz,$(PathologyPimo)))
HEOTB:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE120.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyHE)))
PIMOOTB:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO120.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyPimo)))
convert:   $(HENIFTI) $(PIMONIFTI)
gmm:       $(HEGMM) $(PIMOGMM)
lm:        $(HELMREG) $(PIMOLMREG)
transform: $(addprefix $(DATADIR)/,$(UpdateTransform))
otb: $(HEOTB) $(PIMOOTB)

# debug
jobs:
	@echo $(HEOTB)

#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# initialize LM
$(DATADIR)/%/RefImgLM.nii.gz: $(DATADIR)/%/RefImg.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
$(DATADIR)/%/T1postLM.nii.gz: $(DATADIR)/%/T1post.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
$(DATADIR)/%/DCEavgLM.nii.gz: $(DATADIR)/%/DCEavg.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 
$(DATADIR)/%/PathHELM.nii.gz: $(DATADIR)/%/PathHE.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 
$(DATADIR)/%/PathPIMOLM.nii.gz: $(DATADIR)/%/PathPIMO.nii.gz
	-$(C3DEXE) $< -scale 0. -type char $@ 

# apply transformation
$(DATADIR)/%/Pathology/PathHE.lmreg.nii.gz: $(DATADIR)/%/t2HElmtransform.tfm $(DATADIR)/%/Pathology/PathHE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	echo $(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEred.nii.gz    -o $(@D)/PathHEred.lmreg.nii.gz    -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEgreen.nii.gz  -o $(@D)/PathHEgreen.lmreg.nii.gz  -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEblue.nii.gz   -o $(@D)/PathHEblue.lmreg.nii.gz   -r $(word 3, $^)  -n Linear   -t $< 
	$(C3DEXE) $(@D)/PathHEred.lmreg.nii.gz    $(@D)/PathHEgreen.lmreg.nii.gz  $(@D)/PathHEblue.lmreg.nii.gz   -omc $@
$(DATADIR)/%/Pathology/PathHE.gmm.lmreg.nii.gz: $(DATADIR)/%/t2HElmtransform.tfm $(DATADIR)/%/Pathology/PathHE.gmm.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
$(DATADIR)/%/Pathology/PathPIMO.lmreg.nii.gz: $(DATADIR)/%/t2PIMOlmtransform.tfm $(DATADIR)/%/Pathology/PathPIMO.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	echo $(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOEred.nii.gz    -o $(@D)/PathPIMOEred.lmreg.nii.gz    -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOEgreen.nii.gz  -o $(@D)/PathPIMOEgreen.lmreg.nii.gz  -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOEblue.nii.gz   -o $(@D)/PathPIMOEblue.lmreg.nii.gz   -r $(word 3, $^)  -n Linear   -t $< 
	$(C3DEXE) $(@D)/PathPIMOEred.lmreg.nii.gz    $(@D)/PathPIMOEgreen.lmreg.nii.gz  $(@D)/PathPIMOEblue.lmreg.nii.gz   -omc $@
$(DATADIR)/%/Pathology/PathPIMO.gmm.lmreg.nii.gz: $(DATADIR)/%/t2PIMOlmtransform.tfm $(DATADIR)/%/Pathology/PathPIMO.gmm.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
$(DATADIR)/%/DCE/DCEavg.lmreg.nii.gz: $(DATADIR)/%/t2dcelmtransform.tfm $(DATADIR)/%/DCE/DCEavg.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(DATADIR)/%/DCE/DCE.lmreg.nii.gz: $(DATADIR)/%/t2dcelmtransform.tfm $(DATADIR)/%/DCE/DCE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(DATADIR)/%/T1Post/T1post.lmreg.nii.gz: $(DATADIR)/%/t2t1lmtransform.tfm $(DATADIR)/%/T1Post/T1post.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(DATADIR)/%/T1Pre/T1pre.lmreg.nii.gz: $(DATADIR)/%/t2t1lmtransform.tfm $(DATADIR)/%/T1Pre/T1pre.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 

# update transform lm
# FIXME - add texture to dependencies
$(DATADIR)/%/updatetransform: $(DATADIR)/%/Pathology/PathHELM.nii.gz $(DATADIR)/%/Pathology/PathHE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/T2wReference/RefImg.hdr $(DATADIR)/%/Pathology/PathPIMOLM.nii.gz $(DATADIR)/%/Pathology/PathPIMO.nii.gz $(DATADIR)/%/T1post/T1postLM.nii.gz $(DATADIR)/%/T1post/T1post.hdr $(DATADIR)/%/T1pre/T1pre.hdr $(DATADIR)/%/DCE/DCEavgLM.nii.gz $(DATADIR)/%/DCE/DCEavg.hdr  $(DATADIR)/%/DCE/DCE.hdr $(DATADIR)/%/Pathology/PathHE.lmreg.nii.gz  $(DATADIR)/%/Pathology/PathHE.gmm.nii.gz
	-$(ITKSNAP) -l LMLabels.txt -s $(word 1, $^)  -g $(word 2, $^) -o $(WORKDIR)/$*/Pathology/PathHE000.Entropy_$(OTBRADIUS).nii.gz &  PIDPATH=$$!; \
	$(ITKSNAP) -l LMLabels.txt                   -g $(word 2, $^)  -s $(word 14, $^) &  PIDGMM=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 3, $^)  -g $(word 4, $^)  -o $(word 13, $^)  &  PIDMRI=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 5, $^)  -g $(word 6, $^)  -o $(WORKDIR)/$*/Pathology/PathPIMO000.Entropy_$(OTBRADIUS).nii.gz &  PIDPIMO=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 7, $^)  -g $(word 8, $^)  -o $(word 9, $^)  &  PIDT1=$$!; \
        $(ITKSNAP) -l LMLabels.txt -s $(word 10, $^) -g $(word 11, $^) -o $(word 12, $^) &  PIDDCE=$$!; \
        zenity --info --title="OutputFile" --text="Tools -> Layer Inspector -> General -> Display Mode -> RGB $(word 4, $^)  "; \
        kill -9 $$PIDPATH;kill -9 $$PIDMRI;kill -9 $$PIDPIMO; kill -9 $$PIDT1; kill -9 $$PIDDCE

# compute lm transformation
$(DATADIR)/%/t2HElmtransform.tfm:   $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/Pathology/PathHELM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(DATADIR)/%/t2PIMOlmtransform.tfm: $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/Pathology/PathPIMOLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(DATADIR)/%/t2t1lmtransform.tfm: $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/T1post/T1postLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(DATADIR)/%/t2dcelmtransform.tfm: $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/DCE/DCEavgLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
$(DATADIR)/%/PathHE.nii.gz: $(DATADIR)/%/PathHE.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@
$(DATADIR)/%/PathPIMO.nii.gz: $(DATADIR)/%/PathPIMO.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@

$(DATADIR)/%/PathHE.mask.nii.gz: $(DATADIR)/%/PathHE.nii.gz
	$(C3DEXE) $<  -threshold 0% 40% 1 0 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathPIMO.mask.nii.gz: $(DATADIR)/%/PathPIMO.nii.gz
	$(C3DEXE) $<  -threshold 0% 40% 1 0 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathHE.gmm.nii.gz: $(DATADIR)/%/PathHE.nii.gz $(DATADIR)/%/PathHE.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $*/PathHEred.nii.gz $*/PathHEgreen.nii.gz $*/PathHEblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUE)] -x $(word 2,$^) -a $*/PathHEred.nii.gz -a $*/PathHEgreen.nii.gz -a $*/PathHEblue.nii.gz   -o [$@,$*/PathHEgmmPOSTERIORS%d.nii.gz] 
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathPIMO.gmm.nii.gz: $(DATADIR)/%/PathPIMO.nii.gz $(DATADIR)/%/PathPIMO.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $*/PathPIMOred.nii.gz $*/PathPIMOgreen.nii.gz $*/PathPIMOblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUE)] -x $(word 2,$^) -a $*/PathPIMOred.nii.gz -a $*/PathPIMOgreen.nii.gz -a $*/PathPIMOblue.nii.gz   -o [$@,$*/PathPIMOgmmPOSTERIORS%d.nii.gz] 
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
	$(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathHE.nii.gz -o $(WORKDIR)/$*/PathHE120.Entropy_$(OTBRADIUS).nii.gz

$(WORKDIR)/%/PathPIMO120.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(DATADIR)/%/PathPIMO.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	$(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathPIMO.nii.gz -o $(WORKDIR)/$*/PathPIMO120.Entropy_$(OTBRADIUS).nii.gz

view:
	vglrun itksnap -g pathology.nii.gz -s pathologyLM.nii.gz
	vglrun itksnap -g dce.nii.gz -s dceLM.nii.gz
	echo Tools -> Layer Inspector -> General -> Display Mode -> RGB

ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 
ITKSNAP  = vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
CONVERTSVS = python ./convertsvs.py
DIMENSION = 3
NTISSUEHE   = 3
NTISSUEPIMO = 2
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
	$(MYSQL) --local-infile < ./hccpathdb.sql
	$(MYSQL) -sNre "call HCCPath.HCCPathDBList();"  > datalocation/dependencies

HENIFTI:=  $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.nii.gz,$(PathologyHE)))
PIMONIFTI:=$(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.nii.gz,$(PathologyPimo)))
HELMREG:=  $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.lmreg.nii.gz,$(PathologyHE)))
PIMOLMREG:=$(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.lmreg.nii.gz,$(PathologyPimo)))
HEGMM:=    $(addprefix $(DATADIR)/,$(subst PathHE.svs,PathHE.gmm.nii.gz,$(PathologyHE)))
PIMOGMM:=  $(addprefix $(DATADIR)/,$(subst PathPIMO.svs,PathPIMO.gmm.nii.gz,$(PathologyPimo)))
HEOTB:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE000.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyHE)))
PIMOOTB:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO000.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyPimo)))
HESTAT:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE000.HaralickCorrelation_$(OTBRADIUS).sql,$(PathologyHE)))
PIMOSTAT:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO000.HaralickCorrelation_$(OTBRADIUS).sql,$(PathologyPimo)))
convert:   $(HENIFTI) $(PIMONIFTI)
gmm:       $(HEGMM) $(PIMOGMM)
lm:        $(HELMREG) $(PIMOLMREG)
transform: $(addprefix $(DATADIR)/,$(UpdateTransform))
otb: $(HEOTB) $(PIMOOTB)
stat: $(HESTAT) $(PIMOSTAT)

# debug
jobs:
	@echo $(HEOTB)

#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# initialize LM
$(DATADIR)/%LM.nii.gz: $(DATADIR)/%.hdr
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
$(DATADIR)/%/updatetransform: $(DATADIR)/%/Pathology/PathHELM.nii.gz $(DATADIR)/%/Pathology/PathHE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz $(DATADIR)/%/T2wReference/RefImg.hdr $(DATADIR)/%/Pathology/PathPIMOLM.nii.gz $(DATADIR)/%/Pathology/PathPIMO.nii.gz $(DATADIR)/%/T1post/T1postLM.nii.gz $(DATADIR)/%/T1post/T1post.hdr $(DATADIR)/%/T1pre/T1pre.hdr $(DATADIR)/%/DCE/DCEavgLM.nii.gz $(DATADIR)/%/DCE/DCEavg.hdr  $(DATADIR)/%/DCE/DCE.hdr $(DATADIR)/%/Pathology/PathHE.lmreg.nii.gz  $(DATADIR)/%/Pathology/PathHE.gmm.nii.gz  $(DATADIR)/%/Pathology/PathHE.mask.nii.gz $(DATADIR)/%/Pathology/PathPIMO.gmm.nii.gz  $(DATADIR)/%/Pathology/PathPIMO.mask.nii.gz
	echo $@
	-$(ITKSNAP) -l LMLabels.txt -g $(word  2, $^) -s $(word  1, $^) -o $(WORKDIR)/$*/Pathology/PathHE000.Entropy_$(OTBRADIUS).nii.gz &  PIDPATH=$$!;   \
	 $(ITKSNAP) -l LMLabels.txt -g $(word  2, $^) -s $(word 14, $^) -o $(word 15, $^) &  PIDHEGMM=$$!; \
	 $(ITKSNAP) -l LMLabels.txt -g $(word  6, $^) -s $(word 16, $^) -o $(word 17, $^) &  PIDPIMOGMM=$$!; \
         $(ITKSNAP) -l LMLabels.txt -g $(word  4, $^) -s $(word  3, $^) -o $(word 13, $^) &  PIDLMREG=$$!; \
         $(ITKSNAP) -l LMLabels.txt -g $(word  6, $^) -s $(word  5, $^) -o $(WORKDIR)/$*/Pathology/PathPIMO000.Entropy_$(OTBRADIUS).nii.gz &  PIDPIMO=$$!; \
         $(ITKSNAP) -l LMLabels.txt -g $(word  8, $^) -s $(word  7, $^) -o $(word 9, $^)  &  PIDT1=$$!;  \
         $(ITKSNAP) -l LMLabels.txt -g $(word 11, $^) -s $(word 10, $^) -o $(word 12, $^) &  PIDDCE=$$!; \
        zenity --info --title="OutputFile" --text="Tools -> Layer Inspector -> General -> Display Mode -> RGB $(word 4, $^)  "; \
        kill -9 $$PIDPATH;kill -9 $$PIDLMREG;kill -9 $$PIDPIMO; kill -9 $$PIDT1; kill -9 $$PIDDCE; kill -9 $$PIDHEGMM; kill -9 $$PIDPIMOGMM

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
	$(C3DEXE) -verbose -mcs $< -popas b -popas g -popas r -push b -push g -scale -1 -add -dup -multiply -sqrt -popas bg  -push b -push r  -scale -1 -add -dup -multiply -sqrt -push bg -add -popas rgb -push g -push r -scale -1 -add -dup -multiply -sqrt -push rgb -add -threshold 0 10 0 1 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathPIMO.mask.nii.gz: $(DATADIR)/%/PathPIMO.nii.gz
	$(C3DEXE) -verbose -mcs $< -popas b -popas g -popas r -push b -push g -scale -1 -add -dup -multiply -sqrt -popas bg  -push b -push r  -scale -1 -add -dup -multiply -sqrt -push bg -add -popas rgb -push g -push r -scale -1 -add -dup -multiply -sqrt -push rgb -add -threshold 0 10 0 1 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathHE.gmm.nii.gz: $(DATADIR)/%/PathHE.nii.gz $(DATADIR)/%/PathHE.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $(DATADIR)/$*/PathHEred.nii.gz $(DATADIR)/$*/PathHEgreen.nii.gz $(DATADIR)/$*/PathHEblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUEHE)] -x $(word 2,$^) -a $(DATADIR)/$*/PathHEred.nii.gz -a $(DATADIR)/$*/PathHEgreen.nii.gz -a $(DATADIR)/$*/PathHEblue.nii.gz   -o [$@,$(DATADIR)/$*/PathHEgmmPOSTERIORS%d.nii.gz] 
	echo $(ITKSNAP) -s  $@ -g $<

$(DATADIR)/%/PathPIMO.gmm.nii.gz: $(DATADIR)/%/PathPIMO.nii.gz $(DATADIR)/%/PathPIMO.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $(DATADIR)/$*/PathPIMOred.nii.gz $(DATADIR)/$*/PathPIMOgreen.nii.gz $(DATADIR)/$*/PathPIMOblue.nii.gz
	$(ATROPOSCMD) -i kmeans[$(NTISSUEPIMO)] -x $(word 2,$^) -a $(DATADIR)/$*/PathPIMOred.nii.gz -a $(DATADIR)/$*/PathPIMOgreen.nii.gz -a $(DATADIR)/$*/PathPIMOblue.nii.gz   -o [$@,$(DATADIR)/$*/PathPIMOgmmPOSTERIORS%d.nii.gz] 
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
# HE
$(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(DATADIR)/%/PathHE.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	$(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathHE.nii.gz -o $(WORKDIR)/$*/PathHE000.Entropy_$(OTBRADIUS).nii.gz

# push to database
$(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS).sql: $(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS)/lstat.csv
	$(MYSQLIMPORT) --replace --fields-terminated-by=',' --lines-terminated-by='\n' --ignore-lines 1 HCCPath $<

$(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS)/lstat.csv: $(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS).nii.gz $(DATADIR)/%/PathHELM.nii.gz
	mkdir -p $(@D)
	$(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt
	sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@

# PIMO
$(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(DATADIR)/%/PathPIMO.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	$(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathPIMO.nii.gz -o $(WORKDIR)/$*/PathPIMO000.Entropy_$(OTBRADIUS).nii.gz

# push to database
$(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS).sql: $(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS)/lstat.csv
	$(MYSQLIMPORT) --replace --fields-terminated-by=',' --lines-terminated-by='\n' --ignore-lines 1 HCCPath $<

$(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS)/lstat.csv: $(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS).nii.gz $(DATADIR)/%/PathPIMOLM.nii.gz
	mkdir -p $(@D)
	$(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt
	sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@


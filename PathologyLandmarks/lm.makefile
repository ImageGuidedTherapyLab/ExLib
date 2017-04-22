SHELL := /bin/bash
ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 
ANTSWarpImageMultiTransformCMD=$(ANTSPATH)/WarpImageMultiTransform 
ITKSNAP  = vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
CONVERTSVS = python ./convertsvs.py
DIMENSION = 3
OTBOFFSET  = 30
OTBRADIUS  = 20
ATROPOSCMD=$(ANTSPATH)/Atropos -d $(DIMENSION)  -c [3,0.0] -m [0.1,1x1x1] 
#OTBTEXTURE=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//otbScalarImageToTexturesFilter
OTBTEXTURE=/rsrch2/ip/dtfuentes/github/ExLib/otbScalarImageTextures/otbScalarImageToTexturesFilter
CLUSTERDIR=/rsrch2/ip/dtfuentes/github/ExLib/PathologyLandmarks/Processed
WORKDIR=Processed

################
# Dependencies #
################
-include $(PWD)/datalocation/dependencies
datalocation/dependencies: ./hccpathdb.sql ./datalocation/CorrelativePathFiles.csv
	$(MYSQL) --local-infile < $< 
	$(MYSQL) -sNre "call HCCPath.HCCPathDBList('$(word 2, $^)');"  > $@

HENIFTI:=  $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE.nii.gz,$(PathologyHE)))
PIMONIFTI:=$(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO.nii.gz,$(PathologyPimo)))
HELMREG:=  $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE.t2reg.nii.gz,$(PathologyHE)))
PIMOLMREG:=$(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO.t2reg.nii.gz,$(PathologyPimo)))
HEPIMOREG:=$(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO.hereg.nii.gz,$(PathologyPimo))) $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO.gmm.hereg.nii.gz,$(PathologyPimo)))
HEGMM:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE.gmm.nii.gz,$(PathologyHE)))
PIMOGMM:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO.gmm.nii.gz,$(PathologyPimo))) 
HEOTB:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE000.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyHE)))
PIMOOTB:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO000.HaralickCorrelation_$(OTBRADIUS).nii.gz,$(PathologyPimo)))

HESTAT:=   $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE000.HaralickCorrelation_$(OTBRADIUS).sql,$(PathologyHE))) $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHE000.Entropy_$(OTBRADIUS).sql,$(PathologyHE))) $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHELMdist.sql,$(PathologyHE)))
PIMOSTAT:= $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO000.HaralickCorrelation_$(OTBRADIUS).sql,$(PathologyPimo))) $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMO000.Entropy_$(OTBRADIUS).sql,$(PathologyPimo))) $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMOLMdist.sql,$(PathologyPimo)))
HEDIST:=    $(addprefix $(WORKDIR)/,$(subst PathHE.svs,PathHELMdist.nii.gz,$(PathologyHE)))
PIMODIST:=  $(addprefix $(WORKDIR)/,$(subst PathPIMO.svs,PathPIMOLMdist.nii.gz,$(PathologyPimo)))
convert:   $(HENIFTI) $(PIMONIFTI)
gmm:       $(HEGMM) $(PIMOGMM)
#lm:        $(HELMREG) $(PIMOLMREG) $(HEPIMOREG)
lm:        $(HEPIMOREG)
overlap:   $(addprefix $(WORKDIR)/,$(subst PathHE.svs,overlap.sql,$(PathologyHE)))
transform: $(addprefix $(WORKDIR)/,$(UpdateTransform))
otb: $(HEOTB) $(PIMOOTB)
stat: $(HESTAT) $(PIMOSTAT)
dist: $(HEDIST) $(PIMODIST)
imagestat:  $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T2WeightedReference)))       \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T2starMapOxygen)))           \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T2statMapMedicalAir)))       \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T2starMapAbsoluteChange)))   \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T2statMapPercentageChange))) \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(DCE)))                       \
            $(addprefix $(WORKDIR)/,$(subst DCE.hdr,DCEavg.sql,$(DCE)))              \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T1PreContrast)))             \
            $(addprefix $(WORKDIR)/,$(subst .hdr,.sql,$(T1PostContrast)))            

# debug
jobs:
	@echo $(HEOTB)

check:
	@$(foreach idfile,$(T2WeightedReference)      , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T2starMapOxygen)          , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T2statMapMedicalAir)      , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T2starMapAbsoluteChange)  , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T2statMapPercentageChange), if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(DCE)                      , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T1PreContrast)            , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(T1PostContrast)           , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(PathologyHE)              , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)
	@$(foreach idfile,$(PathologyPimo)            , if [ ! -f $(DATADIR)/$(idfile)  ]  ; then echo 'missing ' $(DATADIR)/$(idfile) ;fi;)


#https://www.gnu.org/software/make/manual/html_node/Special-Targets.html
# do not delete secondary files
.SECONDARY: 

# initialize LM
$(DATADIR)/%LM.nii.gz: $(DATADIR)/%.hdr
	-$(C3DEXE) $< -scale 0. -type char $@ 

$(WORKDIR)/%LMdist.nii.gz: $(DATADIR)/%LM.nii.gz
	$(C3DEXE) $< -thresh 4 4 1 0 -dilate 1 5x5x0vox -erode 1 5x5x0vox -sdt -o $@

# apply transformation
$(WORKDIR)/%/Pathology/PathHE.t2reg.nii.gz: $(WORKDIR)/%/t2HElmtransform.tfm $(WORKDIR)/%/Pathology/PathHE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEred.nii.gz    -o $(@D)/PathHEred.t2reg.nii.gz    -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEgreen.nii.gz  -o $(@D)/PathHEgreen.t2reg.nii.gz  -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathHEblue.nii.gz   -o $(@D)/PathHEblue.t2reg.nii.gz   -r $(word 3, $^)  -n Linear   -t $< 
	$(C3DEXE) $(@D)/PathHEred.t2reg.nii.gz    $(@D)/PathHEgreen.t2reg.nii.gz  $(@D)/PathHEblue.t2reg.nii.gz   -omc $@
$(WORKDIR)/%/Pathology/PathHE.gmm.t2reg.nii.gz: $(WORKDIR)/%/t2HElmtransform.tfm $(WORKDIR)/%/Pathology/PathHE.gmm.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
$(WORKDIR)/%/Pathology/PathPIMO.t2reg.nii.gz: $(WORKDIR)/%/t2PIMOlmtransform.tfm $(WORKDIR)/%/Pathology/PathPIMO.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOred.nii.gz    -o $(@D)/PathPIMOred.t2reg.nii.gz    -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOgreen.nii.gz  -o $(@D)/PathPIMOgreen.t2reg.nii.gz  -r $(word 3, $^)  -n Linear   -t $< 
	$(ANTSAPPLYTRANSFORMSCMD) -d 3      -i $(word 2, $(^D))/PathPIMOblue.nii.gz   -o $(@D)/PathPIMOblue.t2reg.nii.gz   -r $(word 3, $^)  -n Linear   -t $< 
	$(C3DEXE) $(@D)/PathPIMOred.t2reg.nii.gz    $(@D)/PathPIMOgreen.t2reg.nii.gz  $(@D)/PathPIMOblue.t2reg.nii.gz   -omc $@
$(WORKDIR)/%/Pathology/PathPIMO.gmm.t2reg.nii.gz: $(WORKDIR)/%/t2PIMOlmtransform.tfm $(WORKDIR)/%/Pathology/PathPIMO.gmm.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n NearestNeighbor   -t $<  --float 0 
$(WORKDIR)/%/DCE/DCEavg.t2reg.nii.gz: $(WORKDIR)/%/t2dcelmtransform.tfm $(WORKDIR)/%/DCE/DCEavg.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(WORKDIR)/%/DCE/DCE.t2reg.nii.gz: $(WORKDIR)/%/t2dcelmtransform.tfm $(WORKDIR)/%/DCE/DCE.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(WORKDIR)/%/T1Post/T1post.t2reg.nii.gz: $(WORKDIR)/%/t2t1lmtransform.tfm $(WORKDIR)/%/T1Post/T1post.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
$(WORKDIR)/%/T1Pre/T1pre.t2reg.nii.gz: $(WORKDIR)/%/t2t1lmtransform.tfm $(WORKDIR)/%/T1Pre/T1pre.nii.gz $(DATADIR)/%/T2wReference/RefImgLM.nii.gz
	$(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 

# update transform lm
# FIXME - add texture to dependencies
$(WORKDIR)/%/updatetransform: $(DATADIR)/%/Pathology/PathHELM.nii.gz                $(DATADIR)/%/T2wReference/RefImgLM.nii.gz             $(DATADIR)/%/Pathology/PathPIMOLM.nii.gz              $(DATADIR)/%/T1post/T1postLM.nii.gz                   $(DATADIR)/%/DCE/DCEavgLM.nii.gz                      $(DATADIR)/%/BOLD/T2mapAbsChange/T2AbsChangeLM.nii.gz 
	echo $@
	$(C3DEXE) $(WORKDIR)/$*/Pathology/PathHE.gmm.nii.gz -as AA  $(WORKDIR)/$*/Pathology/PathPIMO.gmm.hereg.nii.gz  -as BB -push AA -thresh 2 3 1 0 -push BB -thresh 1 1 1 0 -overlap 1 
	-$(ITKSNAP) -l LMLabels.txt -g $(WORKDIR)/$*/Pathology/PathHE.nii.gz             -s $(DATADIR)/$*/Pathology/PathHELM.nii.gz                -o $(WORKDIR)/$*/Pathology/PathHE000.Entropy_$(OTBRADIUS).nii.gz $(WORKDIR)/$*/Pathology/PathHELMdist.nii.gz &  PIDPATH=$$!;  \
         $(ITKSNAP) -l LMLabels.txt -g $(WORKDIR)/$*/Pathology/PathPIMO.nii.gz           -s $(DATADIR)/$*/Pathology/PathPIMOLM.nii.gz              -o $(WORKDIR)/$*/Pathology/PathPIMO000.Entropy_$(OTBRADIUS).nii.gz  $(WORKDIR)/$*/Pathology/PathPIMOLMdist.nii.gz &  PIDPIMO=$$!; \
	 $(ITKSNAP) -l LMLabels.txt -g $(WORKDIR)/$*/Pathology/PathHE.nii.gz             -s $(WORKDIR)/$*/Pathology/PathHE.gmm.nii.gz              -o $(WORKDIR)/$*/Pathology/PathPIMO.hereg.nii.gz         &  PIDHEGMM=$$!;  \
	 $(ITKSNAP) -l LMLabels.txt -g $(WORKDIR)/$*/Pathology/PathPIMO.hereg.nii.gz     -s $(WORKDIR)/$*/Pathology/PathPIMO.gmm.hereg.nii.gz      -o $(WORKDIR)/$*/Pathology/PathHE.nii.gz      &  PIDPIMOGMM=$$!;\
         $(ITKSNAP) -l LMLabels.txt -g $(DATADIR)/$*/T2wReference/RefImg.hdr             -s $(DATADIR)/$*/T2wReference/RefImgLM.nii.gz             -o $(WORKDIR)/$*/Pathology/PathHE.t2reg.nii.gz     $(DATADIR)/$*/BOLD/T2mapAbsChange/T2AbsChange.hdr  $(DATADIR)/$*/BOLD/T2mapMedAir/MedAirT2.hdr $(DATADIR)/$*/BOLD/T2mapOxy/OxyT2.hdr $(DATADIR)/$*/BOLD/T2mapPctChange/T2PctChange.hdr             &  PIDLMREG=$$!; \
         $(ITKSNAP) -l LMLabels.txt -g $(DATADIR)/$*/T1post/T1post.hdr                   -s $(DATADIR)/$*/T1post/T1postLM.nii.gz                   -o $(DATADIR)/$*/T1pre/T1pre.hdr               &  PIDT1=$$!;  \
         $(ITKSNAP) -l LMLabels.txt -g $(DATADIR)/$*/DCE/DCEavg.hdr                      -s $(DATADIR)/$*/DCE/DCEavgLM.nii.gz                      -o $(DATADIR)/$*/DCE/DCE.hdr                &  PIDDCE=$$!; \
         $(ITKSNAP) -l LMLabels.txt -g $(DATADIR)/$*/BOLD/T2mapAbsChange/T2AbsChange.hdr -s $(DATADIR)/$*/BOLD/T2mapAbsChange/T2AbsChangeLM.nii.gz -o $(DATADIR)/$*/BOLD/T2mapMedAir/MedAirT2.hdr $(DATADIR)/$*/BOLD/T2mapOxy/OxyT2.hdr $(DATADIR)/$*/BOLD/T2mapPctChange/T2PctChange.hdr  &  PIDBOLD=$$!; \
        zenity --info --title="OutputFile" --text="Tools -> Layer Inspector -> General -> Display Mode -> RGB   $*   "; \
        pkill -9 ITK-SNAP

# apply tranformation
$(WORKDIR)/%/Pathology/PathPIMO.hereg.nii.gz: $(WORKDIR)/%/HEPIMOlmtransform.txt  $(WORKDIR)/%/Pathology/PathPIMO.nii.gz  $(WORKDIR)/%/Pathology/PathHE.nii.gz 
	$(ANTSWarpImageMultiTransformCMD) 3 $(word 2, $(^D))/PathPIMOred.nii.gz       $(@D)/PathPIMOred.hereg.nii.gz    $< --use-NN -R  $(word 3, $^) 
	$(ANTSWarpImageMultiTransformCMD) 3 $(word 2, $(^D))/PathPIMOgreen.nii.gz     $(@D)/PathPIMOgreen.hereg.nii.gz  $< --use-NN -R  $(word 3, $^) 
	$(ANTSWarpImageMultiTransformCMD) 3 $(word 2, $(^D))/PathPIMOblue.nii.gz      $(@D)/PathPIMOblue.hereg.nii.gz   $< --use-NN -R  $(word 3, $^) 
	$(C3DEXE) $(@D)/PathPIMOred.hereg.nii.gz    $(@D)/PathPIMOgreen.hereg.nii.gz  $(@D)/PathPIMOblue.hereg.nii.gz   -omc $@
$(WORKDIR)/%/Pathology/PathPIMO.gmm.hereg.nii.gz: $(WORKDIR)/%/Pathology/PathHE.nii.gz $(WORKDIR)/%/Pathology/PathPIMO.gmm.nii.gz  $(WORKDIR)/%/HEPIMOlmtransform.txt
	$(ANTSWarpImageMultiTransformCMD) 3 $(word 2, $^)  $@ $(word 3,$^) --use-NN -R  $<

# affine registration of pathology mask
$(WORKDIR)/%/Pathology/HEPIMOantsintroAffine.txt: $(WORKDIR)/%/Pathology/PathHE.mask.nii.gz $(WORKDIR)/%/Pathology/PathPIMO.mask.nii.gz 
	cd $(@D); ../../../../antsIntroduction.sh -d 2 -i $(word 2,$(^F))  -r $(<F)   -o HEPIMOantsintro  -n 0 -s MI -t RA -m 30x90x20 > HEPIMOantsintro.log 2>&1
	cat $(@D)/HEPIMOantsintro.log
# preprocess landmarks
$(WORKDIR)/%/Pathology/PathHELM.nii.gz: $(DATADIR)/%/Pathology/PathHELM.nii.gz 
	$(C3DEXE) $< -replace 4 0 -o  $@
$(WORKDIR)/%/Pathology/PathPIMOLM.nii.gz: $(DATADIR)/%/Pathology/PathPIMOLM.nii.gz 
	$(C3DEXE) $< -replace 4 0 -o  $@
$(WORKDIR)/%/T2wReference/RefImgLM.nii.gz: $(WORKDIR)/%/T2wReference/RefImgLM.nii.gz
	$(C3DEXE) $< -replace 4 0 -o  $@
# compute lm transformation
$(WORKDIR)/%/HEPIMOlmtransform.txt: $(WORKDIR)/%/Pathology/PathHELM.nii.gz $(WORKDIR)/%/Pathology/PathPIMOLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(WORKDIR)/%/t2HElmtransform.txt:   $(WORKDIR)/%/T2wReference/RefImgLM.nii.gz $(WORKDIR)/%/Pathology/PathHELM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(WORKDIR)/%/t2PIMOlmtransform.txt: $(WORKDIR)/%/T2wReference/RefImgLM.nii.gz $(WORKDIR)/%/Pathology/PathPIMOLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(WORKDIR)/%/t2t1lmtransform.txt: $(WORKDIR)/%/T2wReference/RefImgLM.nii.gz $(WORKDIR)/%/T1post/T1postLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1
$(WORKDIR)/%/t2dcelmtransform.txt: $(WORKDIR)/%/T2wReference/RefImgLM.nii.gz $(WORKDIR)/%/DCE/DCEavgLM.nii.gz 
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
$(WORKDIR)/%/PathHE.nii.gz: $(DATADIR)/%/PathHE.svs
	if [ ! -f $@  ] ; then $(CONVERTSVS) --svsfile=$< --outimage=$@ ; fi
$(WORKDIR)/%/PathPIMO.nii.gz: $(DATADIR)/%/PathPIMO.svs
	if [ ! -f $@  ] ; then $(CONVERTSVS) --svsfile=$< --outimage=$@ ; fi

$(WORKDIR)/%/PathHE.mask.nii.gz: $(WORKDIR)/%/PathHE.nii.gz $(DATADIR)/%/PathHELM.nii.gz
	$(C3DEXE) -verbose -mcs $< -popas b -popas g -popas r -push b -push g -scale -1 -add -dup -multiply -sqrt -popas bg  -push b -push r  -scale -1 -add -dup -multiply -sqrt -push bg -add -popas rgb -push g -push r -scale -1 -add -dup -multiply -sqrt -push rgb -add -threshold 0 10 0 1 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0 $(word 2,$^) -thresh 4 4 1 0 -multiply -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(WORKDIR)/%/PathPIMO.mask.nii.gz: $(WORKDIR)/%/PathPIMO.nii.gz $(DATADIR)/%/PathPIMOLM.nii.gz
	$(C3DEXE) -verbose -mcs $< -popas b -popas g -popas r -push b -push g -scale -1 -add -dup -multiply -sqrt -popas bg  -push b -push r  -scale -1 -add -dup -multiply -sqrt -push bg -add -popas rgb -push g -push r -scale -1 -add -dup -multiply -sqrt -push rgb -add -threshold 0 10 0 1 -erode 1 8x8x0vox  -comp -threshold 1 3 1 0  $(word 2,$^) -thresh 4 4 1 0 -multiply -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(WORKDIR)/%/PathHE.gmm.nii.gz: $(WORKDIR)/%/PathHE.nii.gz $(WORKDIR)/%/PathHE.mask.nii.gz $(DATADIR)/%/PathHELM.nii.gz
	$(C3DEXE) -mcs $< -oo $(WORKDIR)/$*/PathHEred.nii.gz $(WORKDIR)/$*/PathHEgreen.nii.gz $(WORKDIR)/$*/PathHEblue.nii.gz
	$(ATROPOSCMD) -i kmeans[2] -x $(word 2,$^) -a $(WORKDIR)/$*/PathHEred.nii.gz -a $(WORKDIR)/$*/PathHEgreen.nii.gz -a $(WORKDIR)/$*/PathHEblue.nii.gz   -o [$(WORKDIR)/$*/PathHE.gmm02.nii.gz] 
	$(ATROPOSCMD) -i kmeans[3] -x $(word 2,$^) -a $(WORKDIR)/$*/PathHEred.nii.gz -a $(WORKDIR)/$*/PathHEgreen.nii.gz -a $(WORKDIR)/$*/PathHEblue.nii.gz   -o [$(WORKDIR)/$*/PathHE.gmm03.nii.gz] 
	$(ATROPOSCMD) -i kmeans[4] -x $(word 2,$^) -a $(WORKDIR)/$*/PathHEred.nii.gz -a $(WORKDIR)/$*/PathHEgreen.nii.gz -a $(WORKDIR)/$*/PathHEblue.nii.gz   -o [$(WORKDIR)/$*/PathHE.gmm04.nii.gz] 
	$(C3DEXE) $(WORKDIR)/$*/PathHE.gmm02.nii.gz $(word 3,$^) -thresh 4 4 1 0 -multiply -o $@ 
	echo $(ITKSNAP) -s  $@ -g $<

$(WORKDIR)/%/PathPIMO.gmm.nii.gz: $(WORKDIR)/%/PathPIMO.nii.gz $(WORKDIR)/%/PathPIMO.mask.nii.gz
	$(C3DEXE) -mcs $< -oo $(WORKDIR)/$*/PathPIMOred.nii.gz $(WORKDIR)/$*/PathPIMOgreen.nii.gz $(WORKDIR)/$*/PathPIMOblue.nii.gz
	$(ATROPOSCMD) -i kmeans[2] -x $(word 2,$^) -a $(WORKDIR)/$*/PathPIMOred.nii.gz -a $(WORKDIR)/$*/PathPIMOgreen.nii.gz -a $(WORKDIR)/$*/PathPIMOblue.nii.gz   -o [$(WORKDIR)/$*/PathPIMO.gmm02.nii.gz] 
	$(ATROPOSCMD) -i kmeans[3] -x $(word 2,$^) -a $(WORKDIR)/$*/PathPIMOred.nii.gz -a $(WORKDIR)/$*/PathPIMOgreen.nii.gz -a $(WORKDIR)/$*/PathPIMOblue.nii.gz   -o [$(WORKDIR)/$*/PathPIMO.gmm03.nii.gz] 
	$(ATROPOSCMD) -i kmeans[4] -x $(word 2,$^) -a $(WORKDIR)/$*/PathPIMOred.nii.gz -a $(WORKDIR)/$*/PathPIMOgreen.nii.gz -a $(WORKDIR)/$*/PathPIMOblue.nii.gz   -o [$(WORKDIR)/$*/PathPIMO.gmm04.nii.gz] 
	$(C3DEXE) $(WORKDIR)/$*/PathPIMO.gmm04.nii.gz -replace 4 2 3 2 -o $@
	echo $(ITKSNAP) -s  $@ -g $<

$(WORKDIR)/%/overlap.csv: $(WORKDIR)/%/PathHE.gmm.nii.gz  $(WORKDIR)/%/PathPIMO.gmm.hereg.nii.gz
	mkdir -p $(@D) 
	$(C3DEXE) $< -as AA  $(word 2,$^) -as BB -push AA -thresh 2 3 1 0 -push BB -thresh 1 1 1 0 -overlap 1 > $(@D)/overlap.txt
	sed "s/OVL: \([0-9]\),/\1,overlap.0\1.nii.gz,/g;s/OVL: 1\([0-9]\),/1\1,overlap.1\1.nii.gz,/g;s/^/$(firstword $(subst /, ,$*)),LABELSGMM.nii.gz,Truth.nii.gz,/g;"  $(@D)/overlap.txt | sed "1 i InstanceUID,FirstImage,SecondImage,LabelID,SegmentationID,MatchingFirst,MatchingSecond,SizeOverlap,DiceSimilarity,IntersectionRatio" > $@
## push to db
$(WORKDIR)/%/overlap.sql: $(WORKDIR)/%/overlap.csv
	$(MYSQLIMPORT) --replace --fields-terminated-by=',' --lines-terminated-by='\n' --ignore-lines 1 HCCPath $<

# FIXME - push to cluster
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE001.   3 50     0            0        $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE010.   3 50     0        $(OTBOFFSET)     0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE011.   3 50     0        $(OTBOFFSET) $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE100.   3 50 $(OTBOFFSET)     0            0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE101.   3 50 $(OTBOFFSET)     0        $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE110.   3 50 $(OTBOFFSET) $(OTBOFFSET)     0        0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] $(OTBTEXTURE) $(@D)/$(<F)  $(@D)/PathHE111.   3 50 $(OTBOFFSET) $(OTBOFFSET) $(OTBOFFSET) 0.5 3.5 > $(@D)/PathHE.otb.log 2>&1 '; rsync -avz $(CLUSTERDIR)/$(@D)/ $(@D)/ ; 
# HE
$(WORKDIR)/%/PathHE000.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(WORKDIR)/%/PathHE.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192]  $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5  > $(@D)/PathHE.otb.log 2>&1 '; 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathHE120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathHE.nii.gz -o $(WORKDIR)/$*/PathHE000.Entropy_$(OTBRADIUS).nii.gz

# PIMO
$(WORKDIR)/%/PathPIMO000.HaralickCorrelation_$(OTBRADIUS).nii.gz: $(WORKDIR)/%/PathPIMO.gmm.nii.gz
	mkdir -p $(CLUSTERDIR)/$*; rsync --exclude '*.svs' -avz $(<D)/ $(CLUSTERDIR)/$*/;  
	ssh dtfuentes@eagle 'bsub -J glcm -Ip -cwd $(CLUSTERDIR)/ -n 6 -q short -W 0:30 -M 8192 -R rusage[mem=8192] -R span[ptile=6] $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO000.   3 $(OTBRADIUS)     0             0            0     0.5 3.5  > $(@D)/PathPIMO.otb.log 2>&1 '; 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO010.   3 $(OTBRADIUS)     0         $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO110.   3 $(OTBRADIUS) $(OTBOFFSET)  $(OTBOFFSET)     0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO100.   3 $(OTBRADIUS) $(OTBOFFSET)      0            0     0.5 3.5 
	echo $(OTBTEXTURE) $(CLUSTERDIR)/$*/$(<F)  $(CLUSTERDIR)/$*/PathPIMO120.   3 $(OTBRADIUS) $(OTBOFFSET) -$(OTBOFFSET)     0     0.5 3.5 
	mkdir -p $(WORKDIR)/$*; rsync  -avz  $(CLUSTERDIR)/$*/ $(WORKDIR)/$*/;
	echo $(ITKSNAP) -s  $< -g $(WORKDIR)/$*/PathPIMO.nii.gz -o $(WORKDIR)/$*/PathPIMO000.Entropy_$(OTBRADIUS).nii.gz

LSTATCMD = mkdir -p $(@D); $(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt && sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@
#image parameters
$(WORKDIR)/%/T2wReference/RefImg/lstat.csv:                       $(DATADIR)/%/T2wReference/RefImg.hdr               $(DATADIR)/%/T2wReference/RefImgLM.nii.gz           
	$(LSTATCMD)
$(WORKDIR)/%/BOLD/T2mapAbsChange/T2AbsChange/lstat.csv:           $(DATADIR)/%/BOLD/T2mapAbsChange/T2AbsChange.hdr   $(DATADIR)/%/T2wReference/RefImgLM.nii.gz           
	$(LSTATCMD)
$(WORKDIR)/%/BOLD/T2mapMedAir/MedAirT2/lstat.csv:                 $(DATADIR)/%/BOLD/T2mapMedAir/MedAirT2.hdr         $(DATADIR)/%/T2wReference/RefImgLM.nii.gz           
	$(LSTATCMD)
$(WORKDIR)/%/BOLD/T2mapOxy/OxyT2/lstat.csv:                       $(DATADIR)/%/BOLD/T2mapOxy/OxyT2.hdr               $(DATADIR)/%/T2wReference/RefImgLM.nii.gz           
	$(LSTATCMD)
$(WORKDIR)/%/BOLD/T2mapPctChange/T2PctChange/lstat.csv:           $(DATADIR)/%/BOLD/T2mapPctChange/T2PctChange.hdr   $(DATADIR)/%/T2wReference/RefImgLM.nii.gz           
	$(LSTATCMD)
$(WORKDIR)/%/T1post/T1post/lstat.csv:                             $(DATADIR)/%/T1post/T1post.hdr                     $(DATADIR)/%/T1post/T1postLM.nii.gz
	$(LSTATCMD)
$(WORKDIR)/%/T1pre/T1pre/lstat.csv:                               $(DATADIR)/%/T1pre/T1pre.hdr                       $(DATADIR)/%/T1post/T1postLM.nii.gz
	$(LSTATCMD)
$(WORKDIR)/%/DCE/DCEavg/lstat.csv:                                $(DATADIR)/%/DCE/DCEavg.hdr                        $(DATADIR)/%/DCE/DCEavgLM.nii.gz 
	$(LSTATCMD)
$(WORKDIR)/%/DCE/DCE/lstat.csv:                                   $(DATADIR)/%/DCE/DCE.hdr                           $(DATADIR)/%/DCE/DCEavgLM.nii.gz  
	$(LSTATCMD)

# entropy
$(WORKDIR)/%000.Entropy_$(OTBRADIUS)/lstat.csv: $(WORKDIR)/%000.Entropy_$(OTBRADIUS).nii.gz $(DATADIR)/%LM.nii.gz
	mkdir -p $(@D)
	$(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt
	sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@

# HaralickCorrelation
$(WORKDIR)/%000.HaralickCorrelation_$(OTBRADIUS)/lstat.csv: $(WORKDIR)/%000.HaralickCorrelation_$(OTBRADIUS).nii.gz $(DATADIR)/%LM.nii.gz
	mkdir -p $(@D)
	$(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt
	sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@

# dist
$(WORKDIR)/%LMdist/lstat.csv: $(WORKDIR)/%LMdist.nii.gz  $(WORKDIR)/%.gmm.nii.gz
	mkdir -p $(@D)
	$(C3DEXE) $< $(word 2,$^) -lstat > $(@D).txt
	sed "s/^\s\+/$(word 1,$(subst /, ,$*)),$(word 2,$(^F)),$(<F),/g;s/\s\+/,/g;s/LabelID/SeriesInstanceUID,SegmentationID,FeatureID,LabelID/g;s/Vol(mm^3)/Vol.mm.3/g;s/Extent(Vox)/ExtentX,ExtentY,ExtentZ/g" $(@D).txt > $@

# push to database
$(WORKDIR)/%.sql: $(WORKDIR)/%/lstat.csv
	$(MYSQLIMPORT) --replace --fields-terminated-by=',' --lines-terminated-by='\n' --ignore-lines 1 HCCPath $<


ANTSAPPLYTRANSFORMSCMD=/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//antsApplyTransforms                       
ANTSLANDMARKCMD       =/opt/apps/ANTsR/dev//ANTsR_src/ANTsR/src/ANTS/ANTS-build//bin//ANTSUseLandmarkImagesToGetAffineTransform 
ITKSNAP  = vglrun /opt/apps/itksnap/itksnap-3.2.0-20141023-Linux-x86_64/bin/itksnap
C3DEXE=/rsrch2/ip/dtfuentes/bin/c3d
MYSQL = mysql 
CONVERTSVS = python ./convertsvs.py

################
# Dependencies #
################
-include datalocation/dependencies
dependencies: 
	$(MYSQL) -sNre "call Metadata.HCCPathDBList();"  > datalocation/dependencies

convert:    $(HENIFTI)
lm:         $(T2LM)
transform:  $(TransformT2HELM)

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

# apply transformation
pathology.lmreg.nii.gz: landmarktransform.tfm pathology.nii.gz dce.nii.gz
	 $(ANTSAPPLYTRANSFORMSCMD) -d 3 -e 1 -i $(word 2, $^)  -o $@ -r $(word 3, $^)  -n Linear   -t $<  --float 0 
	@echo vglrun itksnap -g $(word 3, $^)  -o $@

# compute lm transformation
%/t2HElmtransform.tfm: %/Pathology/PathHELM.nii.gz %/Pathology/PathHE.nii.gz %/T2wReference/RefImgLM.nii.gz %/T2wReference/RefImg.hdr
	vglrun itksnap -l LMLabels.txt -s $(word 1, $^)  -g $(word 2, $^)  &  PIDPATH=$$!; \
        vglrun itksnap -l LMLabels.txt -s $(word 3, $^)  -g $(word 4, $^)  &  PIDMRI=$$!; \
        zenity --info --title="OutputFile" --text="$(word 2, $^) $(word 4, $^)"; \
        kill -9 $$PIDPATH;kill -9 $$PIDMRI;
	$(ANTSLANDMARKCMD) $^  rigid $@  >> $(basename $@).log 2>&1

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
# FIXME - need script to parse header
%/PathHE.nii.gz: %/PathHE.svs
	$(CONVERTSVS) --svsfile=$< --outimage=$@
#aperio.MPP: '0.5027'
#44704x36680 [0,100 43823x36580] (240x240) J2K/KDU Q=70|AppMag = 20|StripeWidth = 2032|ScanScope ID = SS7063|Filename = 72348|Date = 09/07/16|Time = 15:17:46|Time Zone = GMT-07:00|User = 1a5d2360-04a8-4058-a773-a454a6657563|MPP = 0.5027|Left = 19.602404|Top = 21.335609|LineCameraSkew = 0.000917|LineAreaXOffset = 0.020941|LineAreaYOffset = 0.005602|Focus Offset = -0.000500|DSR ID = DCPWPAPERIO|ImageID = 72348|Exposure Time = 32|Exposure Scale = 0.000001|DisplayColor = 0|OriginalWidth = 44704|OriginalHeight = 36680|ICC Profile = ScanScope v1'
#openslide.layer[2].downsample: '16.003614116350867'
#openslide.layer[2].height: '2286'
#openslide.layer[2].width: '2738'
#44704x36680 [0,100 43823x36580] (240x240) J2K/KDU Q=70|AppMag = 20|StripeWidth = 2032|ScanScope ID = SS7063|Filename = 72348|Date = 09/07/16|Time = 15:17:46|Time Zone = GMT-07:00|User = 1a5d2360-04a8-4058-a773-a454a6657563|MPP = 0.5027|Left = 19.602404|Top = 21.335609|LineCameraSkew = 0.000917|LineAreaXOffset = 0.020941|LineAreaYOffset = 0.005602|Focus Offset = -0.000500|DSR ID = DCPWPAPERIO|ImageID = 72348|Exposure Time = 32|Exposure Scale = 0.000001|DisplayColor = 0|OriginalWidth = 44704|OriginalHeight = 36680|ICC Profile = ScanScope v1'

view:
	vglrun itksnap -g pathology.nii.gz -s pathologyLM.nii.gz
	vglrun itksnap -g dce.nii.gz -s dceLM.nii.gz
	echo Layer Inspector -> General -> Display Mode -> RGB

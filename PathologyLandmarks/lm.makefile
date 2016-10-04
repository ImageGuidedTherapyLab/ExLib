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

#Convert svs to nifti using open slide to read header
# http://openslide.org/    MPP = micron per pixel
# FIXME - need script to parse header
pathology.nii.gz: pathology.svs 
	openslide-show-properties $<   | grep -i 'layer\[2\|mpp'
	openslide-write-png $< 0 0 2 2738 2286 pathology.png
	# remove alpha channel from rgba image and write rgb nii.gz
	c3d -verbose -mcs pathology.png -pop -origin 0x0x0mm -spacing 0.008x0.008x0.008mm -omc $@
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

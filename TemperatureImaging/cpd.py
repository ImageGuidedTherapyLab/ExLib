import vtk
# echo vtk version info
print "using vtk version", vtk.vtkVersion.GetVTKVersion()
import vtk.util.numpy_support as vtkNumPy 
import numpy
import os
import sys
import time
import dicom
import scipy.io as scipyio
import csv

class RealTimeDicomFileRead:
  """ Base Class for realtime image header parsing...  """
  def __init__(self,rootDirectory,ExpectedFileSize,DefaultNstep):
    print " base class constructor called \n\n" 
    # dictionary key template = time echo slice type
    self.keyTemplate   = "%04d_%03d_%03d_%02d" 
    self.dataDirectory = rootDirectory
    self.FileSize = ExpectedFileSize
    self.DicomDataDictionary = {}
    self.FilesReadIn = set()
    self.NumTimeStep = DefaultNstep

  def GetHeaderInfo(self):
    """ get initial header info"""
    # get initial file data
    dcmHeaderFileName =  self.QueryDictionary(0,1,0,2)

    headerData = dicom.read_file( "%s/%s" % (self.dataDirectory,dcmHeaderFileName) )

    # get number of slices
    self.nslice = headerData[0x0021,0x104f].value
    # get number of echos
    self.NumberEcho  = headerData[0x0019,0x107e].value
    self.dimensions = [headerData.Rows,headerData.Rows,self.nslice,self.NumberEcho]
    self.FullSize   =  headerData.Rows*headerData.Rows*self.nslice*self.NumberEcho 
    spacing_mm = headerData.PixelSpacing
    spacing_mm.append(headerData.SpacingBetweenSlices) 
    origin_mm  = headerData.ImagePositionPatient
    #convert to meter
    self.spacing = [ 0.001 * float(dXi) for dXi in spacing_mm  ]
    self.origin  = [ 0.001 * float(Xi) for Xi in origin_mm   ]
    print self.dimensions, self.spacing, self.origin
    
    # FIXME should be negative but phase messed up somewhere
    alpha = +0.0097   
    # FIXME need to read in multiple echo times to process MFGRE should 
    # FIXME still be fine for 1st echo
    echoTime = float(headerData.EchoTime)
    imagFreq = float(headerData.ImagingFrequency)
    # temperature map factor
    self.tmap_factor = 1.0 / (2.0 * numpy.pi * imagFreq * alpha * echoTime * 1.e-3)

    # should be equal and imaginary data
    expectedNtime = int(headerData.ImagesinAcquisition)/self.nslice/self.NumberEcho/2
    if( self.NumTimeStep != expectedNtime ):
       print headerData.ImagesinAcquisition,self.nslice,self.NumberEcho
       raise RuntimeError("expecting %d total time points" % expectedNtime )
    # end GetHeaderInfo(self):
    return

  # return image data from raw file names
  def QueryDictionary(self,timeInstance,echo_id,slice_id,imagetype):
    "infinite loop until file is available and of proper size "
    localFileKey = self.keyTemplate % (timeInstance,echo_id,slice_id,imagetype) 
    while ( True ):
      try: 
        # get current list of ONLY files
        directoryList = os.listdir(self.dataDirectory)
        # check if this has already been read in
        filename = self.DicomDataDictionary[localFileKey][0]
        return filename 
      except OSError: 
        print "waiting for directory %s ... " % ( self.dataDirectory) 
        time.sleep(1)
      except KeyError: 
        print "waiting ... time %04d Echo %03d slice %03d type %02d" % (timeInstance,echo_id,slice_id,imagetype) 
        time.sleep(1)
        files = set( filter(lambda x:os.path.isfile("%s/%s" % (self.dataDirectory,x) ) ,directoryList) )
        ### filestmp = filter(lambda x:os.path.isfile("%s/%s" % (self.dataDirectory,x) ) ,directoryList) 
        ### files = set (filter(lambda x:int(x.split(".").pop()) < 50 ,filestmp ) )

        # we will only read files that have not been read
        FilesNotYetRead = files - self.FilesReadIn
        for filename in FilesNotYetRead :
          FullPathToFile = "%s/%s" %(self.dataDirectory,filename) 
          if ( os.path.getsize( FullPathToFile ) > self.FileSize ):
            print "found", FullPathToFile
            dcmimage = dicom.read_file( FullPathToFile )
            deltat = dcmimage[0x0019,0x105a].value/dcmimage[0x0019,0x10f2].value*1.e-6
            sliceIntID = int(abs(round((dcmimage.SliceLocation - dcmimage[0x0019,0x1019].value)/ dcmimage.SpacingBetweenSlices)))
            if(sliceIntID < 0 ) :
              print "SliceLocation", dcmimage.SliceLocation , "spacing between slices",dcmimage.SpacingBetweenSlices, "first scan location " , dcmimage[0x0019,0x1019].value
              raise RuntimeError("slice integer %d < 0 " % sliceIntID )
            rawdataNumber = dcmimage[0x0019,0x10a2].value
            numEchoes = dcmimage[0x0019,0x107e].value
            #check if default ntime not set
            if (self.NumTimeStep == None):
              try:
                self.NumTimeStep = dcmimage.NumberofTemporalPositions
              except AttributeError:
                raise RuntimeError("NumberofTemporalPositions Not found try setting directly\n\t\t--nstep=...")
            #compute timeIntID
            if ( numEchoes == 1 ) : 
               # for 1 echo assume CPD
               numberSlice = dcmimage[0x0021,0x104f].value
               timeIntID = int(dcmimage.InstanceNumber - 1)/int(numberSlice*2)
            elif( numEchoes > 1 ) :
               # for multiple echo assume MFGRE
               tmptimeID = dcmimage.InstanceNumber - 1 - self.NumTimeStep * numEchoes * sliceIntID * 2
               timeIntID = tmptimeID /numEchoes / 2
            else :
               raise RuntimeError("unknown sequence ")
            #error check
            if( timeIntID < 0  or timeIntID >= self.NumTimeStep ):
               print 'timeIntID', timeIntID ,"numEchoes ", numEchoes 
               print "InstanceNumber ", dcmimage.InstanceNumber, "NumberofTemporalPositions", self.NumTimeStep, "sliceIntID" ,sliceIntID 
               print "TriggerTime", dcmimage.TriggerTime, "deltat", deltat, "number slice ", dcmimage[0x0021,0x104f].value
               raise RuntimeError("time error: %d not \\notin [0,%d) " % (timeIntID,self.NumTimeStep) )
            datatype = int(dcmimage[0x0043,0x102f].value)
            #error check
            if ( datatype == 2 or datatype == 3 ) : 
              keyID = self.keyTemplate % ( timeIntID, int(dcmimage.EchoNumbers), 
                                 sliceIntID, datatype ) 
            else :
               raise RuntimeError("\n\n\t unknown datatype %d : expecting real and imaginary data" % datatype)
            #error check key
            if ( keyID in self.DicomDataDictionary) : 
               raise RuntimeError("\n\n\t duplicate keyID %s not allowed...error parsing header" % keyID )
            self.DicomDataDictionary[keyID]=(filename,dcmimage.EchoTime)
            print "deltat", deltat, "raw data", rawdataNumber, "key", keyID
            # not all headers have this ? 
            try:
              print "trigger ",dcmimage.TriggerTime,"temporal ID ",dcmimage.TemporalPositionIdentifier
            except AttributeError:
              pass
            # ensure we do not read this in again
            self.FilesReadIn.add( filename )
          else: 
            print "filesize too small", os.path.getsize( FullPathToFile )
        print "read in ", len(self.FilesReadIn), "files"

  # return image data from raw file names
  def GetRawDICOMData(self,idtime,outDirectoryID):
    
    # loop until files are ready to be read in
    realImageFilenames = []
    imagImageFilenames = []
    # FIXME: index fun, slice start from 0, echo start from 1
    for idslice in range(self.nslice):
      for idecho in range(1,self.NumberEcho+1):
        realImageFilenames.append( self.QueryDictionary( idtime,idecho,idslice,2 ) )  
        imagImageFilenames.append( self.QueryDictionary( idtime,idecho,idslice,3 ) )  
  
    #create local vars
    rootdir = self.dataDirectory
    dim = self.dimensions
    real_array=numpy.zeros(self.FullSize,dtype=numpy.float32)
    imag_array=numpy.zeros(self.FullSize,dtype=numpy.float32) 

    vtkAppendReal = vtk.vtkImageAppendComponents()
    vtkAppendImag = vtk.vtkImageAppendComponents()

    for idEchoLoc,(fileNameReal,fileNameImag) in enumerate(zip(realImageFilenames,imagImageFilenames)):
      # FIXME: index nightmare
      # FIXME: will be wrong for different ordering
      # arrange such that echo varies fast then x, then y, then z
      # example of how slicing behaves
      #>>> x = range(100)
      #>>> x
      #[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99]
      #>>> x[0:100:10]
      #[0, 10, 20, 30, 40, 50, 60, 70, 80, 90]
      #>>> x[1:100:10]
      #[1, 11, 21, 31, 41, 51, 61, 71, 81, 91]
      #>>> x[2:100:10]
      #[2, 12, 22, 32, 42, 52, 62, 72, 82, 92]
      #>>> x[3:100:10]
      #[3, 13, 23, 33, 43, 53, 63, 73, 83, 93]
      idEcho  = idEchoLoc % dim[3]
      idSlice = idEchoLoc / dim[3]
      beginIndex = dim[0]*dim[1]*dim[3]* idSlice   +idEcho
      finalIndex = dim[0]*dim[1]*dim[3]*(idSlice+1)
      stepIndex  = dim[3]
  
      ## realds = dicom.read_file( "%s/%s"%(rootdir,fileNameReal) )
      ## imagds = dicom.read_file( "%s/%s"%(rootdir,fileNameImag) )
      ## realsliceID = int( round((float(realds.SliceLocation) - float(realds[0x0019,0x1019].value))/ realds.SliceThickness))
      ## imagsliceID = int( round((float(imagds.SliceLocation) - float(imagds[0x0019,0x1019].value))/ imagds.SliceThickness))
  
      ## print "%03d echo %03d slice %03d slice [%d:%d:%d] %03d %s %03d %s "% (idEchoLoc,idEcho,idSlice,beginIndex,finalIndex,stepIndex,realsliceID,fileNameReal,imagsliceID,fileNameImag )
  
      vtkRealDcmReader = vtk.vtkDICOMImageReader()
      vtkRealDcmReader.SetFileName("%s/%s"%(rootdir,fileNameReal) )
      vtkRealDcmReader.Update()
      vtkRealData = vtk.vtkImageCast()
      vtkRealData.SetOutputScalarTypeToFloat()
      vtkRealData.SetInput( vtkRealDcmReader.GetOutput() )
      vtkRealData.Update( )
      real_image = vtkRealData.GetOutput().GetPointData() 
      real_array[ beginIndex: finalIndex : stepIndex ] = vtkNumPy.vtk_to_numpy(real_image.GetArray(0)) 
  
      vtkImagDcmReader = vtk.vtkDICOMImageReader()
      vtkImagDcmReader.SetFileName("%s/%s"%(rootdir,fileNameImag) )
      vtkImagDcmReader.Update()
      vtkImagData = vtk.vtkImageCast()
      vtkImagData.SetOutputScalarTypeToFloat()
      vtkImagData.SetInput( vtkImagDcmReader.GetOutput() )
      vtkImagData.Update( )
      imag_image = vtkImagData.GetOutput().GetPointData() 
      imag_array[ beginIndex: finalIndex : stepIndex ] = vtkNumPy.vtk_to_numpy(imag_image.GetArray(0)) 
  
      vtkAppendReal.SetInput( idEchoLoc ,vtkRealDcmReader.GetOutput() )
      vtkAppendImag.SetInput( idEchoLoc ,vtkImagDcmReader.GetOutput() )
      vtkAppendReal.Update( )
      vtkAppendImag.Update( )
  
    vtkRealDcmWriter = vtk.vtkDataSetWriter()
    vtkRealDcmWriter.SetFileName("Processed/%s/realrawdata.%04d.vtk" % (outDirectoryID,idtime) )
    vtkRealDcmWriter.SetInput(vtkAppendReal.GetOutput())
    vtkRealDcmWriter.Update()
  
    vtkImagDcmWriter = vtk.vtkDataSetWriter()
    vtkImagDcmWriter.SetFileName("Processed/%s/imagrawdata.%04d.vtk" % (outDirectoryID,idtime) )
    vtkImagDcmWriter.SetInput(vtkAppendImag.GetOutput())
    vtkImagDcmWriter.Update()
  
    # write numpy to disk in matlab
    echoTimes = []
    for idecho in range(1,self.NumberEcho+1):
       localKey = self.keyTemplate % ( idtime,idecho,0,2 )
       echoTimes.append(self.DicomDataDictionary[localKey][1])
    scipyio.savemat("Processed/%s/rawdata.%04d.mat"%(outDirectoryID,idtime), {'dimensions':dim,'echoTimes':echoTimes,'real':real_array,'imag':imag_array})
  
    # end GetRawDICOMData
    return (real_array,imag_array)

  # write a numpy data to disk in vtk format
  def ConvertNumpyVTKImage(self,NumpyImageData):
    # Create initial image
    dim = self.dimensions
    # imports raw data and stores it.
    dataImporter = vtk.vtkImageImport()
    # array is converted to a string of chars and imported.
    data_string = NumpyImageData.tostring()
    dataImporter.CopyImportVoidPointer(data_string, len(data_string))
    # The type of the newly imported data is set to unsigned char (uint8)
    dataImporter.SetDataScalarTypeToFloat()
    # Because the data that is imported only contains an intensity value (it isnt RGB-coded or someting similar), the importer
    # must be told this is the case.
    dataImporter.SetNumberOfScalarComponents(dim[3])
    # The following two functions describe how the data is stored and the dimensions of the array it is stored in. For this
    # simple case, all axes are of length 75 and begins with the first element. For other data, this is probably not the case.
    # I have to admit however, that I honestly dont know the difference between SetDataExtent() and SetWholeExtent() although
    # VTK complains if not both are used.
    dataImporter.SetDataExtent( 0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)
    dataImporter.SetWholeExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)
    dataImporter.SetDataSpacing( self.spacing )
    dataImporter.SetDataOrigin(  self.origin )
    dataImporter.Update()
    return dataImporter.GetOutput()
  
# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--dbfile", 
                  action="store", dest="dbfile", default=None,
                  help="[REQUIRED] full path to data directory", metavar="FILE")
parser.add_option("--echoID", 
                  action="store", dest="echoID", type="int", default=0,
                  help="[OPTIONAL] echo # to display for CPD", metavar="INT")
parser.add_option("--sliceID", 
                  action="store", dest="sliceID", type="int", default=None,
                  help="[OPTIONAL] slice to display", metavar="INT")
parser.add_option("--nstep", 
                  action="store", dest="nstep", type="int", default=None,
                  help="[OPTIONAL] # of expected time steps ", metavar="INT")
parser.add_option("--baseline", 
                  action="store", dest="baseline", type="float", default=0.0,
                  help="[OPTIONAL] initial temperature ", metavar="FLOAT")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")
(options, args) = parser.parse_args()


if (options.dbfile != None):
  
  with open('./mrtidb.txt', 'r') as csvfile:
      reader = csv.reader(csvfile, delimiter=',')
      rawdata = [row for row in reader]
  header = rawdata.pop(0)
  
  for patientid,datasetid,numfiles,numslice,datadirectory,description in rawdata:
    print "%s/%s" % (patientid,datasetid) , numfiles, datadirectory

    # instantiate helper class
    fileHelper = RealTimeDicomFileRead( datadirectory, 256*256,options.nstep  )
    outputDirID = "%s/%s" % (patientid,datasetid)
    os.system( "mkdir -p Processed/%s" % outputDirID )

    # Get Header data
    fileHelper.GetHeaderInfo( )
    print fileHelper.tmap_factor

    # display the center slice by default
    # FIXME should we display more than one slice ? 
    # FIXME should we display the max of all slice ? 
    DisplaySlice = fileHelper.nslice/2
    if ( options.sliceID!= None ) : 
       DisplaySlice = options.sliceID
    
    print fileHelper.NumTimeStep , int(numfiles)/2,int(numslice),int(fileHelper.NumberEcho)
    if ( fileHelper.NumTimeStep > int(numfiles)*int(numslice)*int(fileHelper.NumberEcho)/2  ):
       fileHelper.NumTimeStep = int(int(numfiles)*int(numslice)*int(fileHelper.NumberEcho)/2)
       print "updating number of time steps"

    print "# time steps %d, display slice %d " % ( fileHelper.NumTimeStep,DisplaySlice )
    

    deltat = 6.0
    pvd=open("Processed/%s/temperature.pvd" % outputDirID ,"w")
    pvd.write('<?xml version="1.0"?>\n')
    pvd.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">\n')
    pvd.write('  <Collection>\n')
    for idtime in range(fileHelper.NumTimeStep):
         pvd.write('   <DataSet timestep="%f" part="0" file="%s.%04d.vti"/>\n' % (idtime*deltat,"temperature",idtime) )
    pvd.write('  </Collection>\n')
    pvd.write('</VTKFile>\n')
    
    # create initial image as 1d array
    absTemp = numpy.zeros(fileHelper.FullSize,
                           dtype=numpy.float32) + options.baseline
    vtkTempImage = fileHelper.ConvertNumpyVTKImage(absTemp)
    vtkTempWriter = vtk.vtkXMLImageDataWriter()
    vtkTempWriter.SetFileName( "Processed/%s/temperature.%04d.vti" % (outputDirID,0))
    vtkTempWriter.SetInput( vtkTempImage )
    vtkTempWriter.Update()
    
    # create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
     
    # create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    
    # try image viewer for multiple windows
    imageViewer = vtk.vtkImageViewer2()

    # loop and compute tmap
    vtkPreviousImage = fileHelper.GetRawDICOMData( 0, outputDirID )
    # do not finish until all files processed
    print "writing timeID " 
    for idfile in range(fileHelper.NumTimeStep): 
      try:
        # get current data set
        vtkCurrent_Image = fileHelper.GetRawDICOMData( idfile, outputDirID )
    
        #  - \delta \theta = atan( conj(S^i) * S^{i+1} ) 
        #                  = atan2(Im,Re) 
        #                  = atan2( S^{i+1}_y S^i_x - S^{i+1}_x S^i_y ,
        #                           S^{i+1}_x S^i_x + S^{i+1}_y S^i_y ) 
        deltaTemp = fileHelper.tmap_factor * numpy.arctan2(
                              vtkPreviousImage[0] * vtkCurrent_Image[1] 
                            - vtkPreviousImage[1] * vtkCurrent_Image[0] ,
                              vtkPreviousImage[1] * vtkCurrent_Image[1] 
                            + vtkPreviousImage[0] * vtkCurrent_Image[0]  )
        absTemp  = absTemp + deltaTemp 
    
        # write numpy to disk in vtk format
        sys.stdout.write(" %d " % (idfile))
        sys.stdout.flush()
        
        vtkTempImage = fileHelper.ConvertNumpyVTKImage(absTemp)
        vtkTempWriter = vtk.vtkXMLImageDataWriter()
        vtkTempWriter.SetFileName( "Processed/%s/temperature.%04d.vti" % (outputDirID,idfile))
        vtkTempWriter.SetInput( vtkTempImage )
        vtkTempWriter.Update()
    
        # write numpy to disk in matlab
        scipyio.savemat("Processed/%s/temperature.%04d.mat"%(outputDirID,idfile), {'temp':absTemp})
    
        # update for next time step
        vtkPreviousImage = vtkCurrent_Image 
    
        # color table
        # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
        # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/ImageProcessing/Python/ImageSlicing.py
        hueLut = vtk.vtkLookupTable()
        hueLut.SetNumberOfColors (256)
        #FIXME: adjust here to change color  range
        hueLut.SetRange (-10.0, 30.0)  
        #hueLut.SetSaturationRange (0.0, 1.0)
        #hueLut.SetValueRange (0.0, 1.0)
        hueLut.SetHueRange (0.667, 0.0)
        hueLut.SetRampToLinear ()
        hueLut.Build()
    
        # colorbar
        # http://www.vtk.org/doc/release/5.8/html/c2_vtk_e_3.html#c2_vtk_e_vtkLookupTable
        scalarBar = vtk.vtkScalarBarActor()
        scalarBar.SetTitle("Temperature")
        scalarBar.SetNumberOfLabels(4)
        scalarBar.SetLookupTable(hueLut)
    
        # image viewer
        imageViewer.SetInput(vtkTempImage)
        imageViewer.SetSize(512,512)
        imageViewer.SetSlice(DisplaySlice)
        imageViewer.SetPosition(512,0)
        imageViewer.Render()

        # extract VOI to display
        extractVOI = vtk.vtkExtractVOI()
        extractVOI.SetInput(vtkTempImage) 
        extractVOI.SetVOI([0,fileHelper.dimensions[0],0,fileHelper.dimensions[1],DisplaySlice,DisplaySlice]) 
        extractVOI.Update()
        
        # mapper
        #mapper = vtk.vtkDataSetMapper()
        mapper = vtk.vtkImageMapToColors()
        mapper.SetInput(extractVOI.GetOutput())
        # set echo to display
        mapper.SetActiveComponent( options.echoID )
        mapper.SetLookupTable(hueLut)
    
        # actor
        actor = vtk.vtkImageActor()
        actor.SetInput(mapper.GetOutput())
         
        # assign actor to the renderer
        ren.AddActor(actor)
        ren.AddActor2D(scalarBar)
         
        # uncomment to enable user interface interactor
        #iren.Initialize()
        renWin.SetSize(512,512)
        renWin.Render()
        #iren.Start()

        # save as a movie for animation
        windowToImage = vtk.vtkWindowToImageFilter() 
        windowToImage.SetInput(renWin)
        windowToImage.Update()
        jpgWriter     = vtk.vtkJPEGWriter() 
        jpgWriter.SetFileName( "Processed/%s/temperature.%04d.jpg" % (outputDirID,idfile))
        #jpgWriter.SetInput(extractVOI.GetOutput())
        jpgWriter.SetInput(windowToImage.GetOutput())
        jpgWriter.Write()

      except KeyboardInterrupt:
        print "moving onto next dataset... "
        break
        #reset reference phase
        ## print "reseting base phase image at time ", idfile
        ## time.sleep(1)
        ## vtkPreviousImage = fileHelper.GetRawDICOMData( idfile, outputDirID )
        ## absTemp = numpy.zeros(fileHelper.FullSize,
        ##                       dtype=numpy.float32) + options.baseline
    # save gif animations
    print "saving animations... "
    os.system( "convert -delay 30 -resize 50%% -loop 0 Processed/%s/temperature.*.jpg Processed/%s/temperature.gif  " % (outputDirID,outputDirID) ) 
else:
  parser.print_help()
  print options

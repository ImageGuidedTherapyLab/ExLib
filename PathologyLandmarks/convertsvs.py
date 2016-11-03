import subprocess
import os

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option( "--svsfile",
                  action="store", dest="svsfile", default=None,
                  help="FILE containing image info", metavar="FILE")
parser.add_option( "--layer",
                  action="store", dest="layer", default=2,
                  help="image layer write out ", metavar="INT")
parser.add_option( "--outimage",
                  action="store", dest="outimage", default=None,
                  help="FILE to write out ", metavar="basename FILE")
(options, args) = parser.parse_args()

c3dexe = '/rsrch2/ip/dtfuentes/bin/c3d'

if (options.svsfile != None and options.svsfile.split('.').pop() == 'svs' and options.outimage != None):
     getHeaderCmd = 'openslide-show-properties %s ' % (options.svsfile)
     print getHeaderCmd
     headerProcess = subprocess.Popen(getHeaderCmd ,shell=True,stdout=subprocess.PIPE )
     while ( headerProcess.poll() == None ):
        pass
     layerinfo = dict([tuple(line.strip('\n').split(":",1)) for line in headerProcess.stdout if (line[:15] =='openslide.layer' or line[:10] =='aperio.MPP')])
     print layerinfo

     #centroid = eval( rawcentroidinfo.replace('CENTROID_VOX',''))
     #print centroid , int(centroid [2])
     downsample = float(layerinfo['openslide.layer[%d].downsample' % options.layer].replace("'",'')  )
     micronperpixel = float(layerinfo['aperio.MPP' ].replace("'",'')  )
     width      = int(  layerinfo['openslide.layer[%d].width'  % options.layer].replace("'",'') )
     height     = int(  layerinfo['openslide.layer[%d].height' % options.layer].replace("'",'') )
     spacing = downsample * micronperpixel  * .001
     pngfile = options.svsfile.replace('.svs','.png')
     pngcmd = 'openslide-write-png %s 0 0 %d %d %d %s ' %(options.svsfile,options.layer,width,height,pngfile )
     print pngcmd
     os.system(pngcmd)
     nifticmd = '%s -verbose -mcs %s -pop -origin 0x0x0mm -spacing %fx%fx1mm -omc %s' %(c3dexe,pngfile,spacing ,spacing  ,options.outimage)
     print nifticmd 
     os.system(nifticmd)

     #try:
     #except Exception as excp:
     #  print excp
else:
  parser.print_help()

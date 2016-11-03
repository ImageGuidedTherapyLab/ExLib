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
     layerinfo = dict([tuple(line.strip('\n').split(":",1)) for line in headerProcess.stdout if line[:15] =='openslide.layer'])
     print layerinfo 
     print layerinfo['openslide.layer[2].height'] 
     print layerinfo['openslide.layer[2].width'] 
     print layerinfo['openslide.layer[2].downsample'] 

     #centroid = eval( rawcentroidinfo.replace('CENTROID_VOX',''))
     #print centroid , int(centroid [2])
     width = 10
     height = 100
     pngfile = options.svsfile.replace('.svs','.png')
     pngcmd = 'openslide-write-png %s 0 0 %d %d %d %s ' %(options.svsfile,options.layer,width,height,pngfile )
     print pngcmd
     #os.system(pngcmd)
     nifticmd = '%s -verbose -mcs %s -pop -origin 0x0x0mm -spacing 0.008x0.008x0.008mm -omc %s' %(c3dexe,pngfile ,options.outimage)
     print nifticmd 

     #try:
     #except Exception as excp:
     #  print excp
else:
  parser.print_help()

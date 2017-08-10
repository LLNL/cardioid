#! /usr/bin/python

# This script creates .ddcMD visualization files (can be opened in Visit) from full resolution anatomy#* files (or non-anatomy files that are ddcMD-ish style).  It may create multiple .ddcMD files, in which case, each one will have to be openend in Visit one by one to visualize entire heart anatomy.  Companion python script "loadddcMDsInVisit.py" aids in automating importing .ddcMD files into Visit.

import sys
import glob
import os
import math
import subprocess

if len(sys.argv) < 3:
  print '\n\nERROR, usage:  anatomy2ddcMDvis [anatomy files base name (also works for non-anatomy files so long as files are ddcMD-ish format)] [maximum number of files you want assigned to each .ddcMD visualization file]\n\n*(smaller number gives more .ddcMD files, but each one takes less time to visualize in Visit, which avoids intrinsic 30 minute execution timeout in Visit)\n\n'
  print 'e.g., anatomy2ddcMDvis anatomy 1\n\n'
  print 'e.g., anatomy2ddcMDvis domains 1\n\n'
  print 'Note:  using 1 for max. number of files per .ddcMD visualization file gives best performance in Visit'
  exit(1)

fb = sys.argv[1]
nbatch = int(sys.argv[2])
 
files = glob.glob(fb+'#*')
print '\nFiles that I will process are as follows: \n'
print files

# Make first directory that will hold nbatch number of files
dirname = 'dir0'
print '\nMaking directory '+dirname+' and corresponding ddcMD header file '+dirname+'.ddcMD\n'
os.system('mkdir '+dirname)
fout = open(dirname+'.ddcMD','w')
fout.write('ddcMD FILEHEADER {files = '+fb+';}')
fout.close()

# Gather header data from existing fb#000000
print 'Gathering header data from existing '+fb+'#000000, and copying its data only (not header) to '+dirname+'/'+fb+'#000000'
myf = fb+'#000000'
fin = open(myf,'r')
fout = open(dirname+'/'+fb+'#000000','a')
ct = -1000000
for line in fin:
 line0 = line[:len(line)-1]
 linesp = line0.split(';')
 for f in linesp:
  line = f
  if ct<2:
    line=line+';\n'
  else:
    line=line+'\n'
  if line.find('datatype') >= 0:
     if line.find('VAR') >= 0:
      print '\n\nERROR:  datatype must be FIXRECORDASCII, which means all records must have the same line length, defined by "lrec" in header\n\n'
      exit()
     else:
      datatype = line
  if line.find('lrec') >= 0:
    lsplrec = line.split('=')
    lsplrec[1] = lsplrec[1].strip()
    lsp1 = lsplrec[1].strip(';')
    lrec = int(lsp1)
  if line.find('nfields') >= 0:
    nfields = line
  if line.find('endian') >= 0:
    endian = line
  if line.find('field_names') >= 0:
    fieldnames = line
  if line.find('field_types') >= 0:
    fieldtypes = line
  if line.find('field_units') >= 0:
    fieldunits = line
  if line.find('nx') >= 0:
    nx = line
  if line.find('ny') >= 0:
    ny = line
  if line.find('nz') >= 0:
    nz = line
  
  
  if line.find('}') >= 0:
    ct = 0

  if ct >= 2:
    if (len(line) != lrec):
      print 'ERROR:  length of line\n\n'+line+'\nis '+str(len(line))+', which is not what header indicates (lrec = '+str(lrec)+')'
      exit()
    fout.write(line)

  ct += 1
fin.close()
fout.close()

ct = 1

if nbatch <= 1:
    print '\nMaking modifications to '+dirname+'/'+fb+'#000000'+', so copied it to '+dirname+'/'+fb+'#000000'+'ORIG before making changes\n'
    os.system('mv '+dirname+'/'+fb+'#000000'+' '+dirname+'/'+fb+'#000000'+'ORIG')
    print '\nReading '+dirname+'/'+fb+'#000000ORIG and writing it to the file'+dirname+'/'+fb+'000000\n'
    print '\nSwitching nx and nz, as Visit expects them to be switched to what we have from these '+fb+' files\n'
    fin = open(dirname+'/'+fb+'#000000'+'ORIG','r')

    fout = open(dirname+'/'+fb+'#000000','a')
    nrecord = os.popen('cat '+dirname+'/'+fb+'* | wc -l').read()
    nxsp = nx.split('=')
    nzsp = nz.split('=')
    nxtmp = nxsp[1]
    nxsp[1] = nzsp[1]
    nzsp[1] = nxtmp
    nxnew = nxsp[0]+'='+nxsp[1]
    nznew = nzsp[0]+'='+nzsp[1]
    if len(files) <= 1:
      nfiles = 1
    else:
      nfiles = int(ct - math.floor((ct-1)/nbatch)*nbatch)
    fout.write('stateVariable FILEHEADER {\n   create_time = Fri Oct 21 08:56:35 2016;\n   user = jpc;   host = vulcanio50-ib0;\n   exe_version = cardioid-bgq-spi r2371M;   srcpath = /g/g19/jpc/CARDIOID/EP_r2371/trunk/src;\n   compile_time = Sep 29 2016 17:07:47;\n'+datatype+'nfiles = '+str(nfiles)+';\n   nrecord = '+str(nrecord).rstrip()+';\n'+nfields+'   lrec = '+str(lrec)+';\n'+endian+'   time = 2080.000000;\n   loop = 208000;\n'+fieldnames+fieldtypes+fieldunits+nxnew+ny+nznew+'}\n\n')

    for line in fin:
      fout.write(line)

    fin.close()
    fout.close()
    print 'Removing '+dirname+'/'+fb+'#000000ORIG'
    os.system('rm '+dirname+'/'+fb+'#000000ORIG') 
    
    if (ct % nbatch == 0) and (len(files)>1):
      dirname = 'dir'+str(int(math.floor(ct/nbatch)))
      print '\nMaking directory '+dirname+' and corresponding ddcMD header file '+dirname+'.ddcMD\n'
      os.system('mkdir '+dirname)
      fout = open(dirname+'.ddcMD','w')
      fout.write('ddcMD FILEHEADER {files = '+fb+';}')
      fout.close()
 
ct = 2

for myf in files:
 if len(files) <= 1 or myf != fb+'#000000':
#  print 'ct is '+str(ct)
  if len(files) > 1 and ct >= 2:
    newf = fb+'#%06d' % int(ct - 1 - math.floor((ct-1)/nbatch)*nbatch)
    print '\nCopying file '+myf+' to '+dirname+'/'+newf+'\n'
    os.system('cp -p '+myf+' '+dirname+'/'+newf)
  if (len(files)>1) and (ct%nbatch == 0 or (ct%nbatch != 0 and ct >= len(files))):
    print '\nMaking modifications to '+dirname+'/'+fb+'#000000'+', so copied it to '+dirname+'/'+fb+'#000000'+'ORIG before making changes\n'
    os.system('mv '+dirname+'/'+fb+'#000000'+' '+dirname+'/'+fb+'#000000'+'ORIG')
    print '\nReading '+dirname+'/'+fb+'#000000ORIG and writing it to the file'+dirname+'/'+fb+'000000\n'
    print '\nSwitching nx and nz, as Visit expects them to be switched to what we have from these '+fb+' files\n'
    fin = open(dirname+'/'+fb+'#000000'+'ORIG','r')
 
    fout = open(dirname+'/'+fb+'#000000','a')
    nrecord = os.popen('cat '+dirname+'/'+fb+'* | wc -l').read()
    nxsp = nx.split('=')
    nzsp = nz.split('=')
    nxtmp = nxsp[1]
    nxsp[1] = nzsp[1]
    nzsp[1] = nxtmp
    nxnew = nxsp[0]+'='+nxsp[1]
    nznew = nzsp[0]+'='+nzsp[1]
    if len(files) <= 1:
      nfiles = 1
    else:
      nfiles = int(ct - math.floor((ct-1)/nbatch)*nbatch)
    fout.write('stateVariable FILEHEADER {\n   create_time = Fri Oct 21 08:56:35 2016;\n   user = jpc;   host = vulcanio50-ib0;\n   exe_version = cardioid-bgq-spi r2371M;   srcpath = /g/g19/jpc/CARDIOID/EP_r2371/trunk/src;\n   compile_time = Sep 29 2016 17:07:47;\n'+datatype+'nfiles = '+str(nfiles)+';\n   nrecord = '+str(nrecord).rstrip()+';\n'+nfields+'   lrec = '+str(lrec)+';\n'+endian+'   time = 2080.000000;\n   loop = 208000;\n'+fieldnames+fieldtypes+fieldunits+nxnew+ny+nznew+'}\n\n')
    
    for line in fin:
      fout.write(line)
    
    fin.close()
    fout.close()
    print 'Removing '+dirname+'/'+fb+'#000000ORIG'
    os.system('rm '+dirname+'/'+fb+'#000000ORIG')

    if (ct % nbatch == 0 and ct < len(files)):
      dirname = 'dir'+str(int(math.floor(ct/nbatch)))
      print '\nMaking directory '+dirname+' and corresponding ddcMD header file '+dirname+'.ddcMD\n'
      os.system('mkdir '+dirname)
      fout = open(dirname+'.ddcMD','w')
      fout.write('ddcMD FILEHEADER {files = '+fb+';}')
      fout.close()
          
  ct += 1  

print '\n\n\n\n###### Instructions for visualization in Visit:\n'

print '1)  Start a session of Visit, reserve 30 nodes on surface\n'
print '2)  Open each *.ddcMD file that this script creates, one by one, in Visit to creat full heart visualization.  For each *.ddcMD file, You will see a big cube appear in visualization window, refer to next step...\n'
print '3)  Unfortunately, visit assigns color value 0 to all non-tissue voxels, so we need to threshold out the tissue.  First set up an "expression" in Visit to reassign value of 0 to something nonsensical, e.g., if(eq(cellType,0),-1000,cellType)\n'
print '4)  Threshold out the -1000 values, set minimum threshold to minimum value of data\n'
print '5)  Set Pseudocolor minimum value to minimum value of data\n'
print 'NOTE:   If number of files per .ddcMD visualization file is set to 1 in running this script, each .ddcMD files takes about 1 minute to visualize in Visit.  Setting number of files per .ddcMD visualization file to 1 is best and yields fastest drawing of viz in Visit\n'
print 'BONUS:  the script loadddcMDsInVisit.py accomplishes steps 2 and  4-5 if you run it in the Visit GUI, but you will need to write in your expression for step 3, and of course, edit the parameters in the script to fit your own needs'

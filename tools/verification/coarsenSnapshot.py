#generate coarsen snapshot by keeping (1/n)**3 of original data
#usage:
#python coarsenSnapshot.py snapshot.xxxxxxxxxxxx n
#create snapshot.xxxxxxxxxxxx_coarsen
import sys, string,os

file_types=['cellViz','anatomy','coarsened_anatomy','state','activationTime']

nx=-1
ny=-1
nz=-1
dx=-1
dy=-1
dz=-1
lrec=-1
filetype=''

#######################################################
def isCoarseNode(ix,iy,iz,factor):
  factor2=factor/2
  if( ix%factor!=factor2 ):
    return False
  if( iy%factor!=factor2 ):
    return False
  if( iz%factor!=factor2 ):
    return False
  return True

#def isCoarseNode(ix,iy,iz,factor):
#  if( ix%factor==1 and iy%factor==1 and iz%factor==1 ):
#     return True
#  else:
#    return False

#######################################################
def read_header(L):
  print 'read_header...'
  global nx,ny,nz,lrec,dx,dy,dz,filetype,factor
  oldline=''
  flag = 1
  nrecords = 0
  
  for line in L: ## loop over lines of file
    if 'FILEHEADER' in line:
      flag=0
    if flag==0:
      if '}' in line:
        flag=1
      if ';' in line:
        line=oldline+line
        oldline=''
        words=string.split(line,';')
      else:
        oldline=oldline+line
        continue
      newline=string.join(words)
      words=string.split(newline,'=')
      newline=string.join(words)
      words=string.split(newline)
      if words[0] in file_types:
        filetype=words[0]
      if 'nx' in words:
        for i in range(len(words)):
          if words[i]=='nx':
            nx=eval(words[i+1])
      if 'ny' in words:
        for i in range(len(words)):
          if words[i]=='ny':
            ny=eval(words[i+1])
      if 'nz' in words:
        for i in range(len(words)):
          if words[i]=='nz':
            nz=eval(words[i+1])
      if 'dx' in words:
        for i in range(len(words)):
          if words[i]=='dx':
            dx=eval(words[i+1])
      if 'dy' in words:
        for i in range(len(words)):
          if words[i]=='dy':
            dy=eval(words[i+1])
      if 'dz' in words:
        for i in range(len(words)):
          if words[i]=='dz':
            dz=eval(words[i+1])
      if 'lrec' in words:
        for i in range(len(words)):
          if words[i]=='lrec':
            lrec=eval(words[i+1])
      if 'h' in words:
        nx=eval(words[1])
        ny=eval(words[5])
        nz=eval(words[9])
    else:
      words=string.split(line)
      if len(words)>0:
        ix=0
        iy=0
        iz=0
        if filetype=='anatomy' or filetype=='coarsened_anatomy' or filetype=='state':
          gid=eval(words[0])
          iz=gid/(ny*nx)
          iy=(gid-iz*ny*nx)/nx
          ix=gid%nx
        elif filetype=='cellViz' or filetype=='activationTime':
          ix=eval(words[0])
          iy=eval(words[1])
          iz=eval(words[2])
        else:
          print 'file type=',filetype
          print 'error in determining ix,iy,iz...'
          return -1
        if( isCoarseNode(ix,iy,iz,factor) ):
          nrecords=nrecords+1
  return nrecords

#######################################################
def read_and_print_newfile(L,ofilename):
  print 'New file:',ofilename
  ofile=open(ofilename,'w')
  global nx,ny,nz,lrec,dx,dy,dz,filetype,nrecords,factor
  oldline=''
  flag = 1
  count=0

  for line in L: ## loop over lines of file
    if 'FILEHEADER' in line:
      flag=0
    if flag==0: #print header
      if ';' in line or '{' in line or '}' in line:
        line=oldline+line
        oldline=''
      else:
        oldline=oldline+line
        continue
      words=string.split(line,';')
      newline=string.join(words)
      splitwords=string.split(newline,'=')
      newline=string.join(splitwords)
      splitwords=string.split(newline)
      if words[0] in file_types:
        filetype=words[0]
      if 'nx' in splitwords:
        s='  nx = '+`ncx`+'; ny = '+`ncy`+'; nz ='+`ncz`+';'+'\n'
        ofile.write(s)
      elif 'h' in splitwords:
        nxstring=`ncx`
        nystring=`ncy`
        nzstring=`ncz`
        s='  h = '+nxstring.rjust(4)+'  0  '+'  0'
        ofile.write(s)
        ofile.write('\n')
        s='        0    '+nystring.rjust(4)+'  0'
        ofile.write(s)
        ofile.write('\n')
        s='        0    0    '+nzstring.rjust(4)+';'
        ofile.write(s)
        ofile.write('\n')
      elif 'dx' in splitwords:
        s='  dx = '+`dx`+'; dy = '+`dy`+'; dz ='+`dz`+';'+'\n'
        ofile.write(s)
      elif 'nrecord' in splitwords or 'nrecords' in splitwords:
        for i in range(len(splitwords)):
          if splitwords[i]=='nrecord' or splitwords[i]=='nrecords':
            splitwords[i+1]=`nrecords`
            newword=string.join(splitwords[i:i+2],' = ')
        for i in range(len(words)):
          if 'nrecord' in words[i]:
            if i==0:
              words[i]='  '+newword
            else:
              words[i]=newword
        newline=string.join(words,';')
        ofile.write(newline)
      elif '}' in splitwords:
        flag=1
        ofile.write(line)
      else:
        ofile.write(line)
    else: #print data
      #sys.exit()
      words=string.split(line)
      if len(words)>0:
        if filetype=='anatomy' or filetype=='coarsened_anatomy' or filetype=='state':
          gid=eval(words[0])
          iz=gid/(ny*nx)
          iy=(gid-iz*ny*nx)/nx
          ix=gid%nx
          if( isCoarseNode(ix,iy,iz,factor) ):
            gid=ncx*ncy*(iz/factor)+ncx*(iy/factor)+ix/factor
            newline=`gid`+' '+string.join(words[1:])
            if lrec>0:
              newline=newline.ljust(lrec-1)
            ofile.write(newline)
            ofile.write('\n')
            count=count+1
        elif filetype=='cellViz' or filetype=='activationTime':
          ix=eval(words[0])
          iy=eval(words[1])
          iz=eval(words[2])
          if( isCoarseNode(ix,iy,iz,factor) ):
            ixstring=`ix/factor`
            iystring=`iy/factor`
            izstring=`iz/factor`
            newline=ixstring.rjust(5)+' '+iystring.rjust(5)+' '+izstring.rjust(5)
            if filetype=='cellViz':
              newline=newline+' '+words[3].rjust(4)+' '+words[4].rjust(8)+' '+words[5].rjust(18)
            elif filetype=='activationTime':
              newline=newline+' '+words[3].rjust(18)
            if lrec>0:
              newline=newline.ljust(lrec-1)
            ofile.write(newline)
            ofile.write('\n')
            count=count+1
      else:
        ofile.write(line)
  return count

#######################################################

filesdir=sys.argv[1]
factor  =eval(sys.argv[2])
inputs=os.listdir(filesdir)
filenametypes=[]
for filename in inputs:
  if '#000000' in filename:
    filetype=filename.split('#')
    filenametypes.append(filetype[0])

print 'File types: ',filenametypes

inputs = [filename for filename in inputs if 'profile' not in filename]
inputs.sort()

#print inputs

#read files and determine nrecords, nx, ny, nz, ...
for ftype in filenametypes:
  filenames = [filename for filename in inputs if ftype in filename]
  nrecords = 0
  for filename in filenames:
    if '#000000' in filename:
      newfilename=filesdir+'/'+filename
      print 'Open file:',newfilename
      ifile=open(newfilename,'r')
      L=ifile.readlines()
      nnewrecords=read_header(L)
      if nnewrecords<0:
        break
      nrecords=nrecords+nnewrecords
      print newfilename,', nrecords=',nrecords
  for filename in filenames:
    if '#' in filename and '#000000' not in filename:
      newfilename=filesdir+'/'+filename
      print newfilename
      ifile=open(newfilename,'r')
      L=ifile.readlines()
      nnewrecords=read_header(L)
      if nnewrecords<0:
        break
      nrecords=nrecords+nnewrecords
      print newfilename,', nrecords=',nrecords

  ncx=nx/factor
  ncy=ny/factor
  ncz=nz/factor
  dx=dx*factor
  dy=dy*factor
  dz=dz*factor
  print 'Coarsened mesh: nx = ',ncx,'; ny = ',ncy,'; nz =',ncz
  print 'Read ',nrecords,' records...'
  if nx<1:
    print 'Could not determine nx. Stop'
    sys.exit()
  if ny<1:
    print 'Could not determine ny. Stop'
    sys.exit()
  if nz<1:
    print 'Could not determine nz. Stop'
    sys.exit()

  #write coarsened files...
  coarsedir=filesdir+'_coarse'+`factor`
  if not os.path.isdir(coarsedir):
    os.mkdir(coarsedir)
  nwrecords = 0
  for filename in filenames:
    if '#' in filename:
      ifile=open(filesdir+'/'+filename,'r')
      L=ifile.readlines()
      ofilename=coarsedir+'/'+filename
      nwrecords=nwrecords+read_and_print_newfile(L,ofilename)

  print 'Wrote ',nwrecords,' records...'

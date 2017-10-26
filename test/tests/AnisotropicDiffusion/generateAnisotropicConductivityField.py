#usage: python generateAnisotropicConductivityField.py nx ny nz

import sys, string
from math import sin, cos, pi

nx=eval(sys.argv[1])
ny=eval(sys.argv[2])
nz=eval(sys.argv[3])

eta1=1.
eta2=1.
eta3=1.

alpha=1.22
beta =1.33
gamma=1.46

sigma1=[]
sigma1.append( 0.13) #11
sigma1.append(-0.24) #12
sigma1.append( 0.11) #13
sigma1.append( 1.00) #22
sigma1.append( 0.20) #23
sigma1.append( 1.00) #33

sigma2=[]
sigma2.append( 1.00) #11
sigma2.append(-0.13) #12
sigma2.append( 0.20) #13
sigma2.append( 0.24) #22
sigma2.append(-0.11) #23
sigma2.append( 1.00) #33

sigma3=[]
sigma3.append( 1.00) #11
sigma3.append( 0.20) #12
sigma3.append( 0.13) #13
sigma3.append( 1.00) #22
sigma3.append(-0.24) #23
sigma3.append( 0.11) #33

dx=pi/float(nx)
dy=pi/float(ny)
dz=pi/float(nz)

print 'anatomy FILEHEADER{'
print '  datatype = VARRECORDASCII;'
print '  nfiles = 1;'
print '  nrecords = ',nx*ny*nz,';'
print '  nfields = 8;'
print '  field_names = gid cellType sigma11 sigma12 sigma13 sigma22 sigma23 sigma33;'
print '  field_types = u u f f f f f f;'
print '  nx = ',nx,'; ny = ',ny,'; nz = ',nz,';'
print '  field_units = 1 1 mS/mm mS/mm mS/mm mS/mm mS/mm mS/mm;'
print '}'
print ''

#loop over mesh
for k in range(nz):
  z=(float(k)+0.5)*dz
  for j in range(ny):
    y=(float(j)+0.5)*dy
    for i in range(nx):
      x=(float(i)+0.5)*dx
      
      gid=k*nx*ny+j*nx+i
      print gid,
      print '101',

      #print 6 components of tensor
      for l in range(6):
        print eta1*sin(z)*sin(y)*sigma1[l] \
             +eta2*sin(z)*sin(x)*sigma2[l] \
             +eta3*sin(y)*sin(x)*sigma3[l],
      print ''
      

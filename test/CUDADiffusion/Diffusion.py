import numpy as np
import re

N = 35

mV = np.ndarray(shape=(N,N,N), dtype=float, order='F')
for i in xrange(N):
  for j in xrange(N):
    for k in xrange(N):
      mV[i,j,k] = np.sin(i+j+k)

XX = np.ndarray(shape=(N-1,N,N), dtype=float, order='F')
XY = np.ndarray(shape=(N-1,N,N), dtype=float, order='F')
XZ = np.ndarray(shape=(N-1,N,N), dtype=float, order='F')
Xm = np.ndarray(shape=(N-1,N,N), dtype=float, order='F')

YX = np.ndarray(shape=(N,N-1,N), dtype=float, order='F')
YY = np.ndarray(shape=(N,N-1,N), dtype=float, order='F')
YZ = np.ndarray(shape=(N,N-1,N), dtype=float, order='F')
Ym = np.ndarray(shape=(N,N-1,N), dtype=float, order='F')

ZX = np.ndarray(shape=(N,N,N-1), dtype=float, order='F')
ZY = np.ndarray(shape=(N,N,N-1), dtype=float, order='F')
ZZ = np.ndarray(shape=(N,N,N-1), dtype=float, order='F')
Zm = np.ndarray(shape=(N,N,N-1), dtype=float, order='F')

for i in xrange(N-1):
  XX[i,:,:] =  mV[i+1,:,:] - mV[i,:,:]
  Xm[i,:,:] = (mV[i+1,:,:] + mV[i,:,:])/2.0

for i in xrange(1,N-1):
  XY[:,i,:] = (Xm[:,i+1,:] - Xm[:,i-1,:])/2.0
  XZ[:,:,i] = (Xm[:,:,i+1] - Xm[:,:,i-1])/2.0


for i in xrange(N-1):
  YY[:,i,:] =  mV[:,i+1,:] - mV[:,i,:]
  Ym[:,i,:] = (mV[:,i+1,:] + mV[:,i,:])/2.0

for i in xrange(1,N-1):
  YX[i,:,:] = (Ym[i+1,:,:] - Ym[i-1,:,:])/2.0
  YZ[:,:,i] = (Ym[:,:,i+1] - Ym[:,:,i-1])/2.0


for i in xrange(N-1):
  ZZ[:,:,i] =  mV[:,:,i+1] - mV[:,:,i]
  Zm[:,:,i] = (mV[:,:,i+1] + mV[:,:,i])/2.0

for i in xrange(1,N-1):
  ZX[i,:,:] = (Zm[i+1,:,:] - Zm[i-1,:,:])/2.0
  ZY[:,i,:] = (Zm[:,i+1,:] - Zm[:,i-1,:])/2.0

inF = open('10457.out','r')
progx = re.compile('^x=(\d+) y=(\d+) z=(\d+) xflux=\(([\d\.]+),([\d\.]+),([\d\.]+)\)')
progy = re.compile('^x=(\d+) y=(\d+) z=(\d+) yflux=\(([\d\.]+),([\d\.]+),([\d\.]+)\)')
progz = re.compile('^x=(\d+) y=(\d+) z=(\d+) zflux=\(([\d\.]+),([\d\.]+),([\d\.]+)\)')

total=0;

for line in inF:
  m = progx.search(line)
  if m:
     total=total+1
     x = m.group(1)
     y = m.group(2)
     z = m.group(3)
     dx = m.group(4)
     dy = m.group(5)
     dz = m.group(6)
     if (abs(float(dx)-XX[x,y,z])>0.00001) or (abs(float(dy)-XY[x,y,z])>0.00001) or (abs(float(dz)-XZ[x,y,z])>0.00001):
       print 'fluxx',x,y,z,dx,dy,dz,XX[x,y,z],XY[x,y,z],XZ[x,y,z]
  else :
     m=progy.search(line)
     if m:
       total=total+1
       x = m.group(1)
       y = m.group(2)
       z = m.group(3)
       dx = m.group(4)
       dy = m.group(5)
       dz = m.group(6)
       if (abs(float(dx)-YX[x,y,z])>0.00001) or (abs(float(dy)-YY[x,y,z])>0.00001) or (abs(float(dz)-YZ[x,y,z])>0.00001):
         print 'fluxy',x,y,z,dx,dy,dz,YX[x,y,z],YY[x,y,z],YZ[x,y,z]
       else :
         m=progz.search(line)
         if m:
           total=total+1
           x = m.group(1)
           y = m.group(2)
           z = m.group(3)
           dx = m.group(4)
           dy = m.group(5)
           dz = m.group(6)
           if (abs(float(dx)-ZX[x,y,z])>0.00001) or (abs(float(dy)-ZY[x,y,z])>0.00001) or (abs(float(dz)-ZZ[x,y,z])>0.00001):
             print 'fluxz',x,y,z,dx,dy,dz,ZX[x,y,z],ZY[x,y,z],ZZ[x,y,z]

inF.close()

print 'comparisons:',total


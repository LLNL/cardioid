import numpy as np
import re

N = 40

mV = np.zeros(shape=(N+2,N+2,N+2), dtype=float)
#cell size larger to make indexing easier
#cell[0:-1,0:-1,0:-1] recovers (N+2)*(N+2)*(N+2) array
cell = np.zeros(shape=(N+3,N+3,N+3), dtype=int)

sigX  = np.zeros(shape=(3,N+2,N+2,N+2), dtype=float)
sigY  = np.zeros(shape=(3,N+2,N+2,N+2), dtype=float)
sigZ  = np.zeros(shape=(3,N+2,N+2,N+2), dtype=float)
sigXe = np.zeros(shape=(3,N+2,N,N), dtype=int)
sigYe = np.zeros(shape=(3,N,N+2,N), dtype=int)
sigZe = np.zeros(shape=(3,N,N,N+2), dtype=int)

#inF = open('snapshot.000000000001/state#000000','r')
inF = open('state#000000','r')
Vp = re.compile('(\d+) (\d+)')
for line in inF:
  m = Vp.search(line)
  if m:
     x = m.group(1)
     v = m.group(2)
     kk = int(x)/N/N + 1
     jj = (int(x)/N)%N + 1
     ii = (int(x))%N + 1
#     print 'gid=%s %d %d %d V=%s' % (x,ii,jj,kk,v)
     cell[ii,jj,kk]=1
     mV[ii,jj,kk]=float(v)
inF.close()


Sscale = 0.001*0.2
#inF = open('snapshot.initial/anatomy#000000','r')
inF = open('anatomy#000000','r')
Vp = re.compile('(\d+) (\d+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+)')
for line in inF:
  m = Vp.search(line)
  if m:
     x = m.group(1)
     kk = int(x)/N/N + 1
     jj = (int(x)/N)%N + 1
     ii = (int(x))%N + 1
     for y in xrange(3):
       sigX[y,ii,jj,kk]=float(m.group(y+3))*Sscale
     sigY[0,ii,jj,kk]=sigX[1,ii,jj,kk]
     for y in xrange(2):
       sigY[y+1,ii,jj,kk]=float(m.group(y+6))*Sscale
     sigZ[0,ii,jj,kk]=sigX[2,ii,jj,kk]
     sigZ[1,ii,jj,kk]=sigY[2,ii,jj,kk]
     sigZ[2,ii,jj,kk]=float(m.group(8))*Sscale

inF.close()

#define conductivity
shift=[[[0,0,0],[1,0,0]]]
shift.append([[0,0,0],[1,0,0], [0,1,0],[1,1,0], [0,-1,0],[1,-1,0]])
shift.append([[0,0,0],[1,0,0], [0,0,1],[1,0,1], [0,0,-1],[1,0,-1]])
for i,sh in enumerate(shift):
  for x in sh:
    dx = x[0] - 1; dy = x[1]; dz = x[2]
    print "%d %d %d %d" % (5+dx,7+dy,10+dz,cell[5+dx,7+dy,10+dz])

for i,sh in enumerate(shift):
  for x in sh:
    dx = x[0] - 1; dy = x[1]; dz = x[2]
    sigXe[i,1:,:,:] += cell[(1+dx):(-1+dx),(1+dy):(-2+dy),(1+dz):(-2+dz)]
    sigYe[(i+1)%3,:,1:,:] += cell[(1+dz):(-2+dz),(1+dx):(-1+dx),(1+dy):(-2+dy)]
    sigZe[(i+2)%3,:,:,1:] += cell[(1+dy):(-2+dy),(1+dz):(-2+dz),(1+dx):(-1+dx)]

cnum=[2,6,6]
for i,c in enumerate(cnum):
  check = lambda x : x==c
  sigXe[i,:,:,:]       = check(sigXe[i,:,:,:])
  sigYe[(i+1)%3,:,:,:] = check(sigYe[(i+1)%3,:,:,:])
  sigZe[(i+2)%3,:,:,:] = check(sigZe[(i+2)%3,:,:,:])

sigX = sigX[:,:,1:-1,1:-1] * sigXe
sigY = sigY[:,1:-1,:,1:-1] * sigYe
sigZ = sigZ[:,1:-1,1:-1,:] * sigZe

#inF = open('snapshot.initial/anatomy#000000','r')
inF = open('anatomy#000000','r')
Vp = re.compile('(\d+) (\d+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+) +([-\d\.]+)')
for line in inF:
  m = Vp.search(line)
  if m:
     x = m.group(1)
     kk = int(x)/N/N
     jj = (int(x)/N)%N
     ii = (int(x))%N
#     print '%d %d %d %s' % (ii+1,jj+1,kk+1,x)
#     print sigX[0:3,ii+1,jj,kk]
#     print sigY[0:3,ii,jj+1,kk]
#     print sigZ[0:3,ii,jj,kk+1]
inF.close()

#for i in xrange(N):
#  for j in xrange(N):
#    for k in xrange(N):
#      mV[i,j,k] = np.sin(i+j+k)

XX =   mV[1:,:,:] - mV[:-1,:,:]
Xm = ( mV[1:,:,:] + mV[:-1,:,:] )/2.0
XY = (Xm[:,2:,:] - Xm[:,:-2,:])/2.0
XZ = (Xm[:,:,2:] - Xm[:,:,:-2])/2.0

YY =   mV[:,1:,:] - mV[:,:-1,:]
Ym =  (mV[:,1:,:] + mV[:,:-1,:])/2.0
YX = (Ym[2:,:,:] - Ym[:-2,:,:])/2.0
YZ = (Ym[:,:,2:] - Ym[:,:,:-2])/2.0

ZZ =   mV[:,:,1:] - mV[:,:,:-1]
Zm =  (mV[:,:,1:] + mV[:,:,:-1])/2.0
ZX = (Zm[2:,:,:] - Zm[:-2,:,:])/2.0
ZY = (Zm[:,2:,:] - Zm[:,:-2,:])/2.0

#for x in [XX,Xm,XY,XZ,YY,Ym,YX,YZ,ZZ,Zm,ZX,ZY]:
#  print(x.shape) 

fluxX =(sigX[0,1:,:,:]*XX[:,1:-1,1:-1] + 
        sigX[1,1:,:,:]*XY[:,:,1:-1] +
        sigX[2,1:,:,:]*XZ[:,1:-1,:] )

fluxY =(sigY[0,:,1:,:]*YX[:,:,1:-1] + 
        sigY[1,:,1:,:]*YY[1:-1,:,1:-1] +
        sigY[2,:,1:,:]*YZ[1:-1,:,:] )

fluxZ =(sigZ[0,:,:,1:]*ZX[:,1:-1,:] + 
        sigZ[1,:,:,1:]*ZY[1:-1,:,:] +
        sigZ[2,:,:,1:]*ZZ[1:-1,1:-1,:] )

dV =(fluxX[1:,:,:] - fluxX[:-1,:,:] +
     fluxY[:,1:,:] - fluxY[:,:-1,:] +
     fluxZ[:,:,1:] - fluxZ[:,:,:-1] )

print '2 7 10 =',sigX[:,3,6,9]
print '2 7 10 =',sigY[:,1,8,9]
print '2 7 10 =',sigZ[:,1,6,11]
print 'X =',XX[2,7,10],XY[2,6,10],XZ[2,7,9]
print 'Y =',YX[1,7,10],YY[2,7,10],YZ[2,7,9]
print 'Z =',ZX[1,7,10],ZY[2,6,10],ZZ[2,7,10]
print 'FX=',fluxX[2,6,9],fluxX[1,6,9]
print 'FY=',fluxY[1,7,9],fluxY[1,6,9]
print 'FZ=',fluxZ[1,6,10],fluxZ[1,6,9]

print 'YX =',mV[1,7,10],mV[1,8,10],mV[3,7,10],mV[3,8,10]
print 'dV=',dV[1,6,9]
print dV.shape

#inF = open('snapshot.000000000001/state#000000','r')
inF = open('state#000000','r')
Vp = re.compile('(\d+) (\d+)')
for line in inF:
  m = Vp.search(line)
  if m:
     x = m.group(1)
     v = m.group(2)
     kk = int(x)/N/N
     jj = (int(x)/N)%N
     ii = (int(x))%N
     print '(%d,%d,%d) %s %e' %(ii+1,jj+1,kk+1,x,float(v)+dV[ii,jj,kk])
inF.close()
#for ii in xrange(N):
#  for jj in xrange(N):
#   for kk in xrange(N):
#     cid = ii + N*(jj + N*(kk))
##     if dV[ii,jj,kk]>0.:
##       print cid,' ',mV[ii,jj,kk]+dV[ii,jj,kk]
#     print "%d %d %d dV=%f %f" %(ii+1,jj+1,kk+1,dV[ii,jj,kk],dV[ii,jj,kk]+mV[ii+1,jj+1,kk+1])


exit()

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


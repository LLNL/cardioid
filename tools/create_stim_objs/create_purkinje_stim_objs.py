# This script takes as input a Paraview state file (.pvsm) with box stimuli in it (place box stimuli using Paraview GUI, and save it as a state file within the GUI), and produces an object.stim.data with STIMULUS objects in it that are compatible with object.data STIMULUS objects.  This saves the trouble of having to manually type in box stimulus data into object.data after positioning them in Paraview, which can take a really long time if there are many box stimuli.It assumes there is a left and a right branch of stimuli.  Box stimuli in Paraview must use the naming convention that an "L" or a "R" in the name means box is either in the Left or Right branch respectively, and boxes within each branch must be numbered in with subsequently larger numbers, e.g., BoxL0 BoxL1 BoxL2 etc.

import re	# import regular expressions
import bisect	# import array bisection algorithm, which uses binary search to insert item into array
import math	# import math
import sys	# import system-specific params. and fcns

# Stimulus parameters, feel free to alter
COND_VEL = 3.0	# conduction velocity along purkinje fibers in mm/ms, (distance btw/ stimuli) / (COND_VEL) = tStart (tStart is stimulus delay start time)
vStim = -72	# object.data parameter stimulus magnitude associated with stimulus in mV/ms
dur = 1		# object.data parameter stimulus duration in ms
per = 1000	# object.data parameter stimulus period in ms
t0 = 0		# object.data parameter minimum start time of stimulus and simulation
tf = 100000	# object.data parameter maximum end time of stimulus and simulation
timeDelay = 1	# do delay stimuli firing based on their location, just fire all at beginning of each heart beat
num_branches = 3	# number of branches, for example, left and right would be 2 branches
durrerDelay = 1	# set stimulus time delays according to Durrer 1970


# Debug parameters, feel free to alter
print_debug = 1	# 1 if want to print debugging info. to std out (messier looking), 0 if not print debuggin info. (cleaner looking)

# User defined functions
def get_ind_first_digit(mystr):
  ind = 0
  for i in mystr:
    if i.isdigit():
      return ind
    ind+=1;

# Check no inconsitencies in parameter choices
if (durrerDelay and not timeDelay):
  timeDelay = 1

# Open Paraview state file (.pvsm file)
try:
	print '\nOpening Paraview state file %s for reading\n' % sys.argv[1]
except:
	print '\nFATAL ERROR:  User must provide a Paraview save state file (.pvsm file) on command line that has box stimuli specified in it according to proper Paraview box stimuli syntax.\nFurthermore, this code is meant to be used assuming there are 2 branches of Purkinje fiber branches going down the septum and up the heart walls, thus,  box stimuli in Paraview must be specified using the naming scheme BoxL# and BoxR# for left and right branches, where # are subsequently larger numbers as traverse respective purkinje fiber branches (e.g., BoxR0 BoxR1 etc))'
	exit()
f = open(sys.argv[1],'r')


# Open and read the the Paraview state file to gather the names of the box stimuli and Paraview-assigned ids of the box stimuli
# Variables in program, probably should not alter unless expert user
listen = 0		# trigger for reading the .pvsm file, if set to 1 then program will search each subsequent line for important box stimuli information until set to 0 again
boxnumlist = []		# list of # of each box stimuli, just the number part, e.g. 0 1 2 etc.
boxfullnamelist = []	# list of full names of box stimuli, e.g. BoxR0 BoxR1 etc.
boxidlist = []		# list of box ids, which are assigned by Paraview to identify boxes instead of using the name of boxes.  
for ii in range(1,num_branches):	# initialize all branches
  boxnumlist.append([])
  boxfullnamelist.append([])
  boxidlist.append([])

for line in f:
	if line.find('<ProxyCollection name="sources">') != -1:
		listen = 1					# start listening for box stim. names!
	
	if listen:
		boxnameind1 = line.find('name="Box')
		if boxnameind1 != -1:
			boxnameind2 = line.find('"',boxnameind1+6)
			name = line[boxnameind1+6:boxnameind2]	# the full name of current box stimulus
			boxidind1 = line.find('id')
			boxidind2 = line.find('"',boxidind1+4)
			boxid = line[boxidind1+4:boxidind2]	# the id of current box stimulus
			indFirstDigit = get_ind_first_digit(name)
			print 'Index of first digit in string %s is %d' % (name,indFirstDigit)
			print 'boxnumlist is ',boxnumlist
			if name.find('R') != -1:		# is box on left or right Purkinje branch?, store on different rows of boxnumlist array depending on which branch stim. is in
				boxlistrow = 1
			elif name.find('S') != -1:
				boxlistrow = 2
			else:
				boxlistrow = 0
			try:					# do binary search to find index to insert new box name in numerical order
				nameind = bisect.bisect_right(boxnumlist[boxlistrow],float(name[indFirstDigit:]))
				boxnumlist[boxlistrow].insert(nameind,float(name[indFirstDigit:]))			# insert box sequential order # (which should be part of full box stim. name) in boxnumlist, keeping list in sequential order by box #
				print 'boxnumlist is', boxnumlist[boxlistrow]
				boxfullnamelist[boxlistrow].insert(nameind,name)			# insert box full name in boxfullnamelist, keeping list in sequential order by box #
				boxidlist[boxlistrow].insert(nameind,boxid)				# insert box id in boxidlist, keeping list in sequential order by box #
			except IndexError:			# if binary search fails, that means the current row of the list does not exist, must make it and then just append item to it
				boxnumlist.append([])							# append sublist to boxnumlist
				boxfullnamelist.append([])						# append sublist to boxfullnamelist
				boxidlist.append([])							# append sublist to boxidlist
				nameind = bisect.bisect_right(boxnumlist[boxlistrow],int(name[indFirstDigit:]))	# now can do binary search on newly created sublist
				boxnumlist[boxlistrow].insert(nameind,int(name[indFirstDigit:]))			# insert box sequential order # (which should be part of full box stim. name) in boxnumlist, keeping list in sequential order by box #
				boxfullnamelist[boxlistrow].insert(nameind,name)			# insert box full name in boxfullnamelist, keeping list in sequential order by box #
				boxidlist[boxlistrow].insert(nameind,boxid)				# insert box id in boxidlist, keeping list in sequential order by box #
			
	if line.find('</ProxyCollection>') != -1 and listen:
		listen = 0					# stop listening for box stim. names!
f.close()

if print_debug:
	print 'Box stim. full names are %s' % boxfullnamelist	
	print 'Box stim. numbers are %s' % boxnumlist
	print 'Box stim. IDs are %s\n' % map(str,boxidlist)

# Check for duplicate box stimuli names
rep_boxes = []
for branch in range(0,len(boxfullnamelist)):
  for box in range(0,len(boxfullnamelist[branch])):
    mybox = boxfullnamelist[branch][box]
    for box2 in boxfullnamelist[branch][box+1:]:
      if (box2 == mybox):
	rep_boxes.append(mybox)
    for branch2 in range(branch+1,len(boxfullnamelist)):
      for box2 in range(0,len(boxfullnamelist[branch2])):
	if (box2 == mybox):
	  rep_boxes.append(mybox)    
if rep_boxes:
  print "\n\nFATAL ERROR:  The following box stimuli names are repeated, every name must be unique, check your .pvsm file contains no duplicate box stimuli names\n\n",rep_boxes
  exit(1)


# Now that have names, numbers, and ids of all box stimuli, open the Paraview state file again and find the the values of center and x,y,z lengths for each of the box stimuli
# Variables in program, probably should not alter unless expert user
f = open(sys.argv[1],'r')
listen = 0
clist = [[[]]]
llist = [[[]]]
listen_c = 0
listen_l = 0

for line in f:
	if line.find('CubeSource') != -1:
		listen = 1						# found a box stimuli entry, so start listening
		boxidind1 = line.find('id')
		boxidind2 = line.find('"',boxidind1+4)
		myboxid = line[boxidind1+4:boxidind2]
		for row in range(0,len(boxidlist)):			# find row and column where current box id exists in boxidlist
			try:
				col = boxidlist[row].index(myboxid)
				break
			except ValueError:
				continue
		if print_debug:
			print 'Found boxid %s at row %d and col %d in boxidlist...' % (myboxid,row,col)
			
	
	if listen:
		if line.find('Property name="Center"') != -1:
			listen_c = 1;						# listen to find center of box stimuli
		if (line.find('XLength') != -1) or (line.find('YLength') != -1) or (line.find('ZLength') != -1):
			listen_l = 1;						# listen to find x,y,z length of box stimuli

	if listen_c:								# get center of current box stimulus
		valind1 = line.find('value')
		if valind1 != -1:
			mynum = re.findall('\d+',line[valind1:])		# get the x,y, or z component of the coordinate
			mynumf = map(float, mynum)
			while True:
				try:						# try insert the x,y, or z component of coordinate into row,col of clist
					clist[row][col].append(mynumf[0])
					break
				except IndexError:
					try:					# if insertion at row,col fails, then may mean we need to append new sublist for coordinate of current box stimuli
						clist[row][0]
						clist[row].append([])
					except IndexError:			# if appending the sublist above fails, then means we need to append a new sublist representing a new Purkinje branch (right or left)
						clist.append([[]])
						continue

	
	if listen_l:								# get x,y,z lengths of current box stimulus
		valind1 = line.find('value')
		if valind1 != -1:
			mynum = re.findall('\d+',line[valind1:])		# get the x,y, or z length of box stimulus
			mynumf = map(float, mynum)
			while True:
				try:						# try insert the x,y, or z component of length into row,col of llist
					llist[row][col].append(mynumf[0])
					break
				except IndexError:
					try:					# if insertion at row,col fails, then may mean we need to append new sublist for length of current box stimuli
						llist[row][0]
						llist[row].append([])
					except IndexError:			# if appending the sublist above fails, then means we need to append a new sublist representing a new Purkinje branch (right or left)
						llist.append([[]])
						continue
	
	if line.find('</Property>') != -1:					# stop listening for coordinates of box center and x,y,z lengths of box!
		listen_c = 0;
		listen_l = 0;
	if line.find('</Proxy>') != -1:						# stop listening to anything!
		listen = 0	
f.close()


if print_debug:
	print '\nclist is %s\n' % ' '.join(map(str,clist))
	print 'llist is %s\n' % ' '.join(map(str,llist))

# Check if any stimuli overlap, that is, stimulate the same region of tissue.  If so, warn the user, because this could cause instabilities in the model if stimulate any one region too much (instabilities in cell level model, at the very least)
wmsg = "";
for branch in range(0,len(clist)):
  for box in range(0,len(clist[branch])):
    wmsgT = "";
    xc = clist[branch][box][0]; yc = clist[branch][box][1]; zc = clist[branch][box][2]
    xl = llist[branch][box][0]; yl = llist[branch][box][1]; zl = llist[branch][box][2]
    for box2 in range(box+1,len(clist[branch])):
       xc2 = clist[branch][box2][0]; yc2 = clist[branch][box2][1]; zc2 = clist[branch][box2][2]
       xl2 = llist[branch][box2][0]; yl2 = llist[branch][box2][1]; zl2 = llist[branch][box2][2]
       if ( (abs(xc-xc2) - abs(0.5*xl + 0.5*xl2)) <= 1.00001  and (abs(yc-yc2) - abs(0.5*yl + 0.5*yl2)) <= 1.00001 and (abs(zc-zc2) - abs(0.5*zl + 0.5*zl2)) <= 1.00001 ):
	   wmsgT = wmsgT + "%s, " % (boxfullnamelist[branch][box2])
    for branch2 in range(branch+1,len(clist)):
      for box2 in range(0,len(clist[branch2])):
         xc2 = clist[branch2][box2][0]; yc2 = clist[branch2][box2][1]; zc2 = clist[branch2][box2][2]
         xl2 = llist[branch2][box2][0]; yl2 = llist[branch2][box2][1]; zl2 = llist[branch2][box2][2]
         if ( (abs(xc-xc2) - abs(0.5*xl + 0.5*xl2)) <= 1.00001  and (abs(yc-yc2) - abs(0.5*yl + 0.5*yl2)) <= 1.00001 and (abs(zc-zc2) - abs(0.5*zl + 0.5*zl2)) <= 1.00001 ):
	     wmsgT = wmsgT + "%s, " % (boxfullnamelist[branch2][box2])
    if wmsgT:
      wmsgT = (boxfullnamelist[branch][box]) + " overlaps with boxes: " + wmsgT
      wmsg += wmsgT + "\n\n"
if wmsg:
  wmsg = "!!!!!!!Warning, the following box stimuli overlap and will stimulate same piece of tissue, potentially leading to numerical instabilities if cumulative stimulus too high:\n\n" + wmsg


# Open object.data to get dx, dy, and dz, so can later calculate distances in mm between stimuli, and then find tStart for each stimulus using this distance and conduction velocity in Purkinje fibers
listen = 0
f = open('object.data','r')
for line in f:
	if line.find('ANATOMY') != -1:
		listen = 1;
	
	if listen:
		dxind1 = line.find('dx')
		if dxind1 != -1:
			mynum = re.findall('\d+.\d+',line[dxind1:])
			dx = map(float, mynum)
		dyind1 = line.find('dy')
		if dyind1 != -1:
			mynum = re.findall('\d+.\d+',line[dyind1:])
			dy = map(float, mynum)
		dzind1 = line.find('dz')
		if dzind1 != -1:
			mynum = re.findall('\d+.\d+',line[dzind1:])
			dz = map(float, mynum)
f.close()			

if print_debug:
	print 'dx = %f, dy = %f, dz = %f' % (dx[0],dy[0],dz[0])

# Calculate distance between successive box stimuli
D=[[]]
for i in range(0,len(clist)):				# loop over Purkinje branches of clist
	for j in range(1,len(clist[i])):		# loop over box stimuli of each branch of clist
		while True:
			try:				# try adding current value of distance to list D
				D[i].append( math.sqrt( ((clist[i][j][0]-clist[i][j-1][0])*dx[0])**2 + ((clist[i][j][1]-clist[i][j-1][1])*dy[0])**2 + ((clist[i][j][2]-clist[i][j-1][2])*dz[0])**2 ) )
				break
			except IndexError:		# if there is an IndexError, then probably need to append sublist to current row (branch) in list D
				try:
					D[i][0]
					D[i].append([])
				except IndexError:	# if cannot append sublist to current row (branch) in list D, then need to append a sublist representing a new Purkinje branch (left or right)
					D.append([])
					continue
if print_debug:
	print '\nD is %s' % ' '.join(map(str,D))

# From list D and conduction velocity (COND_VEL) within Purkinje fiber, find the time delay for each stimulus (tStart in object.data STIMULUS object)
tstart = [[]]
for i in range(0,len(D)):
	for j in range(0,len(D[i])):
		while True:
			try:				# try appending current value of tStart to list tstart
				tstart[i].append( D[i][j] / COND_VEL )
				if j>0:
					tstart[i][j] = tstart[i][j] + tstart[i][j-1]
				break
			except IndexError:		# if appending fails, then may need to append a new sublist for current value of tStart
				try:
					tstart[i][0]
					tstart[i].append([])
				except IndexError:	# if previous append fails, then need to append new sublist representing a new Purkinje branch (left or right)
					tstart.append([])
					continue
if print_debug:
	print '\ntstart is %s' % ' '.join(map(str,tstart))
		
# Open object.stim.data and write to it all box stimuli with proper parameters, in correct format for a STIMULUS object in object.data 
f = open('object.stim.data','w')
for i in range(0,len(boxfullnamelist)):
	for j in range(0,len(boxfullnamelist[i])):
		f.write('%s STIMULUS\n{\n' % boxfullnamelist[i][j])
		f.write('   method = box;\n')
		f.write('   xMin = %f; xMax = %f;\n' % (clist[i][j][0]-llist[i][j][0]/2,clist[i][j][0]+llist[i][j][0]/2))
		f.write('   yMin = %f; yMax = %f;\n' % (clist[i][j][1]-llist[i][j][1]/2,clist[i][j][1]+llist[i][j][1]/2))
		f.write('   zMin = %f; zMax = %f;\n' % (clist[i][j][2]-llist[i][j][2]/2,clist[i][j][2]+llist[i][j][2]/2))
		f.write('   vStim = %f mV/ms;\n' % vStim) 
		if (j == 0 and not durrerDelay) or (durrerDelay and boxfullnamelist[i][j].find("Zero") != -1) or not timeDelay:					# if first box stimulus in Purkinje branch, or Durrer level 0, or no time delays, then set tStart = 0
			f.write('   tStart = 0;\n')
		else:																		# else, tStart comes from the tstart list which was made above, unless using Durrer timing
		        if (durrerDelay and boxfullnamelist[i][j].find("One") != -1):
			  f.write('   tStart = 5;\n')						
			if (not durrerDelay):
			  f.write('   tStart = %f;\n' % tstart[i][j-1])
		f.write('   duration = %f;\n' % dur)
		f.write('   period = %f;\n' % per)
		f.write('   t0 = %f;\n' % t0)
		f.write('   tf = %f;\n}\n\n' % tf)		
f.close()	


# Print directions so user knows how to incorporate data generated in this script in Cardioid
print '\n*******Two-Step directions to use these box stimulations in Cardioid*******\n1)  Copy and paste the below box stimuli names into object.data under SIMULATE object'
myboxarr=[]
for mybox in boxfullnamelist:
	myboxstr = ' '.join(map(str,mybox))
	print '%s' % myboxstr,
print '\n\n2) Append object.data with the STIMULUS objects in the newly created object.stim.data (by copying from object.stim.data and pasting into object.data\n'
print wmsg

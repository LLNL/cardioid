'''
Created on 09/01/2013

@author: butler
'''

class Creator():
    '''
    This Class recreates an arbitary object.data file
    '''


    def __init__(self):
        self.loop = 0;
        self.maxLoop = 2320000
        self.dt_ms = "0.01 ms"
        self.time = 0
        self.printRate = 100;
        self.checkRanges = 0;
        self.snapshotRate = 500;
        self.parallelDiffusionReaction = 1
        self.nDiffusionCores = 2
        self.checkpointRate = 2000000
        self.snapshotCellList = "sensor.txt"
        self.sensor = {}
        self.stimulus = {}
        self.stimulus_index = 1;
        # There are parameters we are unlikely to change and only ahve one object
        

    def Write_file(self, file_name=object.data):
        self.fd = open(file_name,'w')
        self.fd.write("stimulus STIMULUS\n")
        self.fd.write("{\n")
        self.fd.write("   anatomy = anat;\n")
        self.fd.write("   decomposition = decomp;\n")
        self.fd.write("   diffusion = diff;\n")
        self.fd.write("   reaction = react;\n")
        self.fd.write(("   loop = " + repr(self.loop) + ";\n"))
        self.fd.write(("   maxLoop = " + repr(self.maxLoop) + ";\n"))
        self.fd.write(("   dt = " + self.dt_ms + ";\n"))
        self.fd.write(("   time = " + repr(self.time) + ";\n"))
        self.fd.write(("   printRate = " + repr(self.printRate) + ";\n"))
        self.fd.write(("   checkRanges = " + repr(self.checkRanges) + ";\n"))
        self.fd.write(("   snapshotRate = " + repr(self.snapshotRate) + ";\n"))
        self.fd.write(("   parallelDiffusionReaction = " + repr(self.parallelDiffusionReaction) + ";\n"))
        self.fd.write(("   nDiffusionCores = " + repr(self.nDiffusionCores) + ";\n"))
        self.fd.write(("   checkpointRate = " + repr(self.checkpointRate) + ";\n"))
        self.fd.write(("   snapshotCellList = " + self.snapshotCellList + ";\n"))
        self.fd.write(("   stimulus = "))
        for stimuli in self.stimulus.keys():
            self.fd.write(stimuli)
            self.fd.write(" ")
        self.fd.write(";\n")
        self.fd.write(("   sensor = "))
        for sensor in self.sensor.keys():
            self.fd.write(sensor)
            
        self.fd.write("}\n")
        self.fd.write("}\n")
        
        self.fd.write(self.anatomy)
        self.fd.write(self.reaction)
        self.fd.write(self.conductivity)
        self.fd.write(self.diffusion)
        self.fd.write(self.decomp)
        
        
    def create_pioAnatomy(self,dx,dy,dz):
        self.anatomy
    
    def createConducivity(self):
        self.conductivity = "cond CONDUCTIVITY\n"
        self.conductivity.append("{\n")
        self.conductivity.append("method = pio;\n")
        self.conductivity.append("}\n")
        
    def createReaction(self):
        self.reaction = "\
react REACTION\n\
{\n\
   method = TT06Opt;\n\
   tolerance = 0.0001;\n\
   mod = 1;\n\
  fitFile = fit.data;\n\
   fastGate = 1;\n\
   fastNonGate = 1;\n\
   cellTypes = endo mid epi;\n\
}\n\
endo CELLTYPE { clone=endoRRG;}\n\
mid CELLTYPE { clone=midRRG;}\n\
epi CELLTYPE { clone=epiRRG;}\n\
\n"


    
    def createDiffusion(self):
        self.diffusion = "\
diff DIFFUSION\n\
{\n\
   method = FGR;\n\
   diffusionScale = 1000 mm^3/mF;\n\
}\n"

    def create_grid_decomposition(self,X,Y,Z):
        self.decomp = "decomp DECOMPOSITION\n"
        self.decomp.append("{\n")
        self.decomp.append("   method = grid;\n")
        self.decomp.append(("   nx = " + repr(X) + ";\n"))
        self.decomp.append(("   ny = " + repr(Y) + ";\n"))
        self.decomp.append(("   nz = " + repr(Z) + ";\n"))
        self.decomp.append(("   minimize = volume;\n"))
        self.decomp.append(("   threshold = 4;\n"))
        self.decomp.append(("   printStats = 0;\n"))
        self.decomp.append(("   visgrid = 0;\n"))
        self.decomp.append(("}\n"))


    def add_Ca_sensor(self,rate=self.snapshotRate):
        sentext = "ca SENSOR\n"
        sentext.append("{ \n")
        sentext.append("   method = averageCa\n;")
        sentext.append(("   cellList = " + self.snapshotCellList + ";\n"))
        sentext.append(("   printRate = " + repr(self.snapshotRate) + ";\n"))
        sentext.append(("}\n"))
        self.sensor['ca'] = sentext
    
    def add_dv_sensor(self,rate=self.snapshotRate,filename="maxDV.dat"):
        sentext = "dv SENSOR\n"
        sentext.append("{ \n")
        sentext.append("   method = maxDV\n;")
        sentext.append(("   filename = " + self.filename + ";\n"))
        sentext.append(("   printRate = " + repr(self.snapshotRate) + ";\n"))
        sentext.append(("}\n"))
        self.sensor['dv'] = sentext
        
    def add_gradient_sensor(self,radius,rate=self.snapshotRate):
        sentext = "grad SENSOR\n"
        sentext.append("{ \n")
        sentext.append("   method = gradientVoronoiCoarsening\n;")
        sentext.append(("   maxDistance = " + repr(radius) + " mm;\n"))
        sentext.append(("   cellList = " + self.snapshotCellList + ";\n"))
        sentext.append(("   evalRate = " + repr(self.snapshotRate) + ";\n"))
        sentext.append(("   printRate = " + repr(self.snapshotRate) + ";\n"))
        sentext.append(("}\n"))
        self.sensor['grad'] = sentext

    def add_box_periodic_stimulus(self,x_min, x_max, y_min, y_max, z_min, z_max, period, t0=0, tf= 0):
        stimname = "stim_" + repr(self.stimulus_index)
        self.stimulus_index = self.stimulus_index + 1;
        stimulus = stimname + " STIMULUS\n"
        stimulus.append("{\n")
        stimulus.append("    method = box;\n")
        stimulus.append(("    xMin = " + repr(coords[0]) + ";\n"))
        stimulus.append(("    xMax = " + repr(coords[1]) + ";\n"))
        stimulus.append(("    yMin = " + repr(coords[2]) + ";\n"))
        stimulus.append(("    yMax = " + repr(coords[3]) + ";\n"))
        stimulus.append(("    zMin = " + repr(coords[4]) + ";\n"))
        stimulus.append(("    zMax = " + repr(coords[5]) + ";\n")))
        stimulus.append(("    vStim = -36 mV/ms;\n"))
        stimulus.append(("    tStart = 0;\n"))
        stimulus.append(("    duration = 1;\n"))
        stimulus.append(("    period = " + repr(period) + ";\n"))
        stimulus.append(("    t0 = " + repr(t0) + ";\n"))
        if tf == 0:
            tf = self.maxLoop
        stimulus.append(("    tf = " + repr(tf) + ";\n"))
        stimulus.append("}\n")
        self.stimuli[stimname] = stimulus
        
    def add_box_reent_stimulus(self,coords, tStart):
        stimname = "stim_" + repr(self.stimulus_index)
        self.stimulus_index = self.stimulus_index + 1;
        stimulus = stimname + " STIMULUS\n"
        stimulus.append("{\n")
        stimulus.append("    method = box;\n")
        stimulus.append(("    xMin = " + repr(coords[0]) + ";\n"))
        stimulus.append(("    xMax = " + repr(coords[1]) + ";\n"))
        stimulus.append(("    yMin = " + repr(coords[2]) + ";\n"))
        stimulus.append(("    yMax = " + repr(coords[3]) + ";\n"))
        stimulus.append(("    zMin = " + repr(coords[4]) + ";\n"))
        stimulus.append(("    zMax = " + repr(coords[5]) + ";\n"))
        stimulus.append(("    vStim = -36 mV/ms;\n"))
        stimulus.append(("    tStart = " + repr(tStart) + ";\n"))
        stimulus.append(("    duration = 1;\n"))
        stimulus.append(("    period = +  20000000000;\n"))
        stimulus.append("}\n")
        self.stimuli[stimname] = stimulus
        
        
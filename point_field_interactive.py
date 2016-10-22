########################################################################
#
# point_field.py - a 3D correlated random field generator
#
# generate a correlated random field by sequential addition of points
#
########################################################################

from numpy import *
from scipy.spatial import *
from scipy.stats import norm
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from functools import partial
from Tkinter import *

#########################################
#
# classes
#
#########################################

class Params:
    def __init__(self):
        # read in basic constraints on model from file to populate basic model parameter set
        line_input = []        
        input_file = open('params.txt','r')
        for line in input_file: line_input.append(line.split())
        input_file.close()
        self.num_extra_seeds = int(line_input[0][1])
        self.max_pts = int(line_input[1][1])
        self.r_search = float(line_input[2][1])        
        self.ref_dist = float(line_input[3][1])
        self.min_value = float(line_input[4][1])
        self.max_value = float(line_input[5][1])
        self.exp_gen = float(line_input[6][1])
        self.epsilon = float(line_input[7][1])
        self.Interface()
    def Button_click(self,entry):
        self.num_extra_seeds = int(entry[0].get())
        self.max_pts = int(entry[1].get())
        self.r_search = float(entry[2].get())        
        self.ref_dist = float(entry[3].get())
        self.min_value = float(entry[4].get())
        self.max_value = float(entry[5].get())
        self.exp_gen = float(entry[6].get())
        self.epsilon = float(entry[7].get())        
        self.WriteParams()
    def Interface(self):
        # create a data input window for modifying parameter values
        root = Tk()
        container = Frame(root)
        container.grid()
        Label(container,text='Model Parameters', font=('Courier',14)).grid()
        labels = ['No. of extra seeds','Max. no. of spawn points','Search radius','Reference distance for spawned points','Minimum parameter value in domain','Maximum parameter value in domain','Inverse distance exponent','Smoothing distance']
        param_list = [self.num_extra_seeds,self.max_pts,self.r_search,self.ref_dist,self.min_value,self.max_value,self.exp_gen,self.epsilon]
        entry = []
        for i in xrange(len(labels)):
            Label(container,text=labels[i]).grid(row=i+1,sticky=W)
            entry.append(Entry(container))
            entry[i].grid(row=i+1,column=1)
            entry[i].insert(0,str(param_list[i]))
        btn = Button(container,text='UPDATE',command=lambda:self.Button_click(entry)).grid(column=1)
        root.mainloop()
    def WriteParams(self):
        # write present value set to text file
        output_file = open('params.txt','w')
        output_file.writelines(['extra_seeds','\t',str(self.num_extra_seeds),'\n'])
        output_file.writelines(['max_pts','\t',str(self.max_pts),'\n'])
        output_file.writelines(['r_search','\t',str(self.r_search),'\n'])
        output_file.writelines(['ref_dist','\t',str(self.ref_dist),'\n'])
        output_file.writelines(['min_value','\t',str(self.min_value),'\n'])
        output_file.writelines(['max_value','\t',str(self.max_value),'\n'])
        output_file.writelines(['exp_gen','\t',str(self.exp_gen),'\n'])
        output_file.writelines(['smooth_fct','\t',str(self.epsilon),'\n'])
        output_file.close()

class Grid:
    def __init__(self):
        # read in grid constraints
        line_input = []   
        input_file = open('domain.txt','r')
        for line in input_file: line_input.append(line.split())
        input_file.close()
        self.start = array([float(line_input[1][1]),float(line_input[1][2]),float(line_input[1][3])])
        self.end = array([float(line_input[2][1]),float(line_input[2][2]),float(line_input[2][3])])
        self.N = array([int(line_input[3][1]),int(line_input[3][2]),int(line_input[3][3])])
        self.aniso = array([float(line_input[4][1]),float(line_input[4][2]),float(line_input[4][3])])
        self.slope = array([float(line_input[5][1]),float(line_input[5][2])])
        # introduce interface; opportunity to modify grid parameters
        root = Tk()
        container = Frame(root)
        container.grid()
        Label(container,text='Grid Attributes', font=('Courier',14)).grid()
        column_labels = ['X','Y','Z']
        row_labels = ['Start','End','N','Anisotropy (0-1)','Slope']
        for i in xrange(3): Label(container,text=column_labels[i]).grid(row=1,column=i+1,sticky=W)
        for i in xrange(5): Label(container,text=row_labels[i]).grid(row=i+2,sticky=W)
        entry_start = []
        entry_end = []
        entry_N = []
        entry_aniso = []
        entry_slope = []
        for i in xrange(3):
            entry_start.append(Entry(container))
            entry_start[i].grid(row=2,column=i+1)
            entry_start[i].insert(0,str(self.start[i]))
        for i in xrange(3):
            entry_end.append(Entry(container))
            entry_end[i].grid(row=3,column=i+1)
            entry_end[i].insert(0,str(self.end[i]))
        for i in xrange(3):
            entry_N.append(Entry(container))
            entry_N[i].grid(row=4,column=i+1)
            entry_N[i].insert(0,str(self.N[i]))
        for i in xrange(3):
            entry_aniso.append(Entry(container))
            entry_aniso[i].grid(row=5,column=i+1)
            entry_aniso[i].insert(0,str(self.aniso[i]))
        for i in xrange(2):
            entry_slope.append(Entry(container))
            entry_slope[i].grid(row=6,column=i+1)
            entry_slope[i].insert(0,str(self.slope[i]))
        btn = Button(container,text='UPDATE',command=lambda:self.Button_click(entry_start,entry_end,entry_N,entry_aniso,entry_slope)).grid(column=2)
        root.mainloop()
        # assign additional attributes
        self.tensor = 1.0/self.aniso        
        self.dl = (self.end - self.start)/self.N
        self.xg = arange(self.start[0],self.end[0],self.dl[0]) + 0.5*self.dl[0]
        self.yg = arange(self.start[1],self.end[1],self.dl[1]) + 0.5*self.dl[1]
        self.zg = arange(self.start[2],self.end[2],self.dl[2]) + 0.5*self.dl[2]        
        self.X,self.Y,self.Z = meshgrid(self.xg,self.yg,self.zg,indexing='ij',)
        self.grid = array([self.X.flatten(),self.Y.flatten(),self.Z.flatten()]).T
        self.values = zeros(self.N.prod(),float)                        # placeholder
    def Button_click(self,entry_start,entry_end,entry_N,entry_aniso,entry_slope):
        self.start = array([float(entry_start[0].get()),float(entry_start[1].get()),float(entry_start[2].get())])
        self.end = array([float(entry_end[0].get()),float(entry_end[1].get()),float(entry_end[2].get())])
        self.N = array([int(entry_N[0].get()),int(entry_N[1].get()),int(entry_N[2].get())])
        self.aniso = array([float(entry_aniso[0].get()),float(entry_aniso[1].get()),float(entry_aniso[2].get())])
        self.slope = array([float(entry_slope[0].get()),float(entry_slope[1].get())])
        self.WriteParams()
    def WriteParams(self):
        # write present value set to text file
        output_file = open('domain.txt','w')
        output_file.writelines(['\t','X','\t','Y','\t','Z','\n'])
        output_file.writelines(['start','\t',str(self.start[0]),'\t',str(self.start[1]),'\t',str(self.start[2]),'\n'])
        output_file.writelines(['end','\t',str(self.end[0]),'\t',str(self.end[1]),'\t',str(self.end[2]),'\n'])
        output_file.writelines(['N','\t',str(self.N[0]),'\t',str(self.N[1]),'\t',str(self.N[2]),'\n'])
        output_file.writelines(['anisotrophy_(0-1)','\t',str(self.aniso[0]),'\t',str(self.aniso[1]),'\t',str(self.aniso[2]),'\n'])
        output_file.writelines(['slope','\t',str(self.slope[0]),'\t',str(self.slope[1]),'\n'])
        output_file.close()
    def ApplySlope(self):
        # alter grid by applying slope vector to z-values (done post-interpolation to avoid problems with anisotropy, etc.)
        self.Z += (self.X - self.start[0])*self.slope[0] + (self.Y - self.start[1])*self.slope[1]
        self.grid = array([self.X.flatten(),self.Y.flatten(),self.Z.flatten()]).T        
    def InterpGrid(self,pts,v):
        # interpolate field values at grid points
        self.values = griddata(pts*self.tensor, v, self.grid*self.tensor, method='nearest', fill_value=nan)
    def Nodes(self,log_flag):
        # create MODFLOW import-ready grid files ...
        ig = arange(1,self.N[0]+1,1)        # column index
        jg = arange(self.N[1],0,-1)         # row index
        kg = arange(self.N[2],0,-1)         # layer index
        I,J,K = meshgrid(ig,jg,kg,indexing='ij',)
        node_grid = array([K.flatten(),J.flatten(),I.flatten()]).T
        z_bottom = self.Z.flatten() - 0.5*self.dl[2]           # note that 'z' is at the cell center-point
        z_top = self.Z.flatten() + 0.5*self.dl[2]
        if log_flag:
            WriteOutput(node_grid,10.**self.values,'val_mf.txt',0)
        else:
            WriteOutput(node_grid,self.values,'val_mf.txt',0)
        WriteOutput(node_grid,z_bottom,'zbottom_mf.txt',0)
        WriteOutput(node_grid,z_top,'ztop_mf.txt',0)
    def Stretch(self,sigma_new,min_value,max_value):
        # stretch or compress normal distribution of field values
        nu = self.values.mean()
        sigma_0 = self.values.std()                                     # grid (initial) standard deviation (typically smaller/tighter)
        cum_dist = norm.cdf(self.values,nu,sigma_0)                     # match grid values to corresponding cumulative distributions points
        grid_raw = norm.ppf(cum_dist,nu,sigma_new)                      # map cumulative distribution points to new grid values
        grid_bottom = where(grid_raw >= min_value,grid_raw,min_value)
        self.values = where(grid_bottom <= max_value,grid_bottom,max_value)
    def PlotSlice(self):
        # for 2-D realizations, generate color-map plot ...
        if self.N[0] == 1:
            # y-z 2D problem
            b = self.values.reshape(self.N[1],self.N[2])
            x_label = 'y'
            y_label = 'z'
            Contour(array(self.yg),array(self.zg),b.T,x_label,y_label)
        if self.N[1] == 1:
            # x-z 2D problem
            b = self.values.reshape(self.N[0],self.N[2])
            x_label = 'x'
            y_label = 'z'
            Contour(array(self.xg),array(self.zg),b.T,x_label,y_label)
        if self.N[2] == 1:
            # x-y 2D problem
            b = self.values.reshape(self.N[0],self.N[1])
            x_label = 'x'
            y_label = 'y'
            Contour(array(self.xg),array(self.yg),b.T,x_label,y_label)


#########################################
#
# support functions
#
#########################################

        
def NormVector(x,y,z):
    tot = abs(x)+abs(y)+abs(z)
    return array([x/tot,y/tot,z/tot]).T

def WriteOutput(points,values,file_name,header_flag=1):
    output_file = open(file_name,'w')
    if header_flag:
        line_out = ['x','\t','y','\t','z','\t','value','\n']
        output_file.writelines(line_out)
    for i in xrange(len(values)):
        line_out = []
        line_out.append(str(points[i,0]))
        line_out.append('\t')
        line_out.append(str(points[i,1]))
        line_out.append('\t')
        line_out.append(str(points[i,2]))
        line_out.append('\t')
        line_out.append(str(values[i]))
        line_out.append('\n')
        output_file.writelines(line_out)
    output_file.close()

def InvDistSquared(p,pts,v_pts,tensor,k,epsilon):
    # interpolate field value at p(x,y,z) via inverse-distance-squared weighting interpolation
    # p and pts are M x 3 numpy arrays; v_pts is a 1-D numpy array of length M; tensor is a 1 x 3 numpy array
    nv = NormVector(p[0] - pts[:,0],p[1] - pts[:,1],p[2] - pts[:,2])    # normalize the point direction with reference to axes
    f = sqrt(1.0/sum(nv**2/tensor**2,axis=1))                                # find the distance distoration factor from the anisotropy tensor
    mapfunc=partial(distance.euclidean,v=p)
    d = map(mapfunc,pts)
    h = sqrt((f*array(d))**2. + (f*epsilon)**2.)
    a = sum(v_pts/h**k)
    b = sum(1.0/h**k)
    return a/b

def ReadSeeds(params,grid):
    # read in initial point set (used as seeds)
    input_file = open('seeds.txt','r')
    pts = []
    values = []
    i = 0
    for line in input_file:
        if i:                       # don't parse header
            line_input = line.split()
            x = float(line_input[0])
            y = float(line_input[1])
            z = float(line_input[2])
            v = float(line_input[3])
            pts.append([x,y,z])
            values.append(v)  
        i += 1
    # add additional seed points
    for i in xrange(params.num_extra_seeds):
        x = random.uniform(grid.start[0],grid.end[0])
        y = random.uniform(grid.start[1],grid.end[1])
        z = random.uniform(grid.start[2],grid.end[2])
        v = random.uniform(params.min_value,params.max_value)        
        pts.append([x,y,z])
        values.append(v) 
    input_file.close()
    return array(pts),array(values)

def Contour(U,V,M,x_label,y_label):
    # contour distribution, depending on geometry
    plt.pcolor(U,V,M,cmap=cm.RdBu)      
    plt.colorbar()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()

def SetStats(v,grid,params):
    # assess histograms; modify if requested
    root = Tk()
    container = Frame(root)
    container.grid()
    Label(container,text='Value Distributions', font=('Courier',14)).grid()
    column_labels = ['Points','Grid']
    row_labels = ['Mean','Std. dev.','Minimum','Maximum']
    sets = [v,grid.values]
    for i in xrange(4): Label(container,text=row_labels[i]).grid(row=i+2,sticky=W)
    for i in xrange(2):
        Label(container,text=column_labels[i]).grid(row=1,column=i+1,sticky=W)
        Label(container,text=str(sets[i].mean())).grid(row=2,column=i+1,sticky=W)
        Label(container,text=str(sets[i].std())).grid(row=3,column=i+1,sticky=W)
        Label(container,text=str(sets[i].min())).grid(row=4,column=i+1,sticky=W)
        Label(container,text=str(sets[i].max())).grid(row=5,column=i+1,sticky=W)
    Label(container,text='New standard deviation >').grid(row=6,column=0,sticky=W)
    sigma_entry = Entry(container)
    sigma_entry.grid(row=6,column=1)
    sigma_entry.insert(0,str(grid.values.std()))
    btn = Button(container,text='STRETCH HISTOGRAM',command=lambda:grid.Stretch(float(sigma_entry.get()),params.min_value,params.max_value)).grid(row=6,column=2)
    root.mainloop()
    return grid

def SpecialOutput(grid):
    # write to output files formatted for MODFLOW (ModelMuse)
    log_out = 0             # flag to un-log transform K-values (default = not)
    root = Tk()
    container = Frame(root)
    container.grid()
    Label(container,text='Write MODFLOW-Formatted I-J-K Files', font=('Courier',14)).grid()
    log_out = IntVar() 
    Checkbutton(container,text='Un-log transform K-values',variable=log_out).grid(sticky=W)
    btn = Button(container,text='WRITE SPECIAL OUTPUT',command=lambda:grid.Nodes(log_out)).grid()
    root.mainloop()

def Bin(A,n):
    # divide array A into bins by equal-interval method; return bin index numbers
    bins = linspace(min(A), max(A), n, endpoint=False)
    top = bins[1:]
    bottom = bins[:-1]
    return digitize(A, bins)

def ManageGroups(grid, n):
    # instructions for button in CreateGroups window
    group_indices = Bin(grid.values, n)
    WriteOutput(grid.grid, group_indices, 'group_distrib.txt')

def CreateGroups(grid):
    # create bins for assigning properties and write to output file
    root = Tk()
    container = Frame(root)
    container.grid()
    Label(container,text='Create Property Group Indices', font=('Courier',14)).grid()
    Label(container,text='Number of groups').grid(row=1, column=0)
    ngroups_entry = Entry(container)
    ngroups_entry.grid(row=1, column=1)
    ngroups_entry.insert(0,'5')
    btn = Button(container,text='WRITE GROUPS OUTPUT',command=lambda: ManageGroups(grid, int(ngroups_entry.get()))).grid()
    root.mainloop()


#########################################
#
# main script
#
#########################################


def Points():

    # read/modify model parameters
    params = Params()
    #num_extra_seeds,max_pts,r_search,ref_dist,min_value,max_value,exp_gen,epsilon = ReadModelParams()
    print 'Read model parameters.'

    # read grid settings
    grid = Grid()
    print 'Read grid settings.'

    # read seed points, if used
    pts,v = ReadSeeds(params,grid)    
    tree = KDTree(pts)
    print 'Read and supplemented seed points.'

    # populate point set by bootstrapping
    print 'Spawning points ...'

    while len(pts) < params.max_pts:

        # generate new point location
        xp = random.uniform(grid.start[0],grid.end[0])
        yp = random.uniform(grid.start[1],grid.end[1])
        zp = random.uniform(grid.start[2],grid.end[2])        
        p = array([xp,yp,zp])

        # search for points within r_search
        near_point_set = list(tree.query_ball_point(p,params.r_search))

        # extract those points falling within the search radius into a collapse matrix
        near_pts = pts[near_point_set]
        near_vals = v[near_point_set]

        # assign a value associated with the new (x,y,z) location
        if len(near_pts)==0:
            # assign random value from a uniform distribution
            p_value = random.uniform(params.min_value,params.max_value)
        else:
            # create a virtual random point on boundary of cylindrical search zone of radius ref_dist, assign random value to it, and then process along with rest of data set
            theta = random.uniform(0.0,2*pi)
            dxp = params.ref_dist * cos(theta)
            dyp = params.ref_dist * sin(theta)            
            zp = random.uniform(grid.start[2],grid.end[2]) 
            vp = array([xp+dxp,yp+dyp,zp])
            r_value = random.uniform(params.min_value,params.max_value)       # random value for point (some influence from other points, so glaring outlier less likely)
            near_pts_v = concatenate((near_pts,[vp]),axis=0)
            near_vals_v = concatenate((near_vals,[r_value]))
            p_value = InvDistSquared(p,near_pts_v,near_vals_v,grid.tensor,params.exp_gen,params.epsilon)           

        # add new point to location and value arrays
        pts = concatenate((pts,[p]),axis=0)
        v = concatenate((v,[p_value]))

    # write point set to output file
    WriteOutput(pts,v,'point_set.txt')
    print '\t... wrote point set to point_set file.'

    # interpolate point set across grid
    print 'Gridding ...'
    grid.InterpGrid(pts,v)

    # summarize set statistics and stretch grid histogram, if requested
    grid = SetStats(v,grid,params)

    # process output
    grid.ApplySlope()                                   # alter z-values, post-interpolation, to account for slope
    WriteOutput(grid.grid,grid.values,'grid_out.txt')   # values output
    CreateGroups(grid)                                  # create property group distributions output, if indicated
    SpecialOutput(grid)                                 # check if special output (e.g., MODFLOW input) files are to be written
    print 'Wrote grid point values to output file(s), as appropriate.'
    grid.PlotSlice()
    print 'Done.'


########################################################
#    
# run script
#
########################################################

Points()

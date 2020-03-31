import numpy as np
from ipywidgets import widgets
import ase
import ase.io as aseio
import ase.neighborlist as asenl
from itertools import permutations
from tqdm.notebook import tqdm, trange
import nglview

import warnings
warnings.simplefilter('ignore')

# Class to define and operate on environments
class Environments:
    def __init__(self):
        self.myindex = -1
        self.indeces = np.array([],dtype=np.int)
        self.delta = np.array([])
        self.distance =  np.array([])
        
# Class to find environments
class EnvironmentFinder:
    
    def __init__(self):
        self.fractionalFlag=False
        self.uniqueFlag=True
        self.fastFlag=True
        self.allEnvs=np.ndarray((0,),dtype=np.object)
        self.uniqueEnvs = np.ndarray((0,),dtype=np.object)
        self.tolerance=0
        #env=Environments()
        #env.delta.shape = (0,3)
        #uniqueEnvs = np.append(uniqueEnvs,env)
        self.box_layout = widgets.Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='solid',
                    width='90%')

    def plotConf(self):
        v = nglview.show_ase(self.conf)
        v.add_representation("unitcell")
        boxy=widgets.Box([v],layout=self.box_layout)
        display(boxy)

    def chooseConfiguration(self,Configuration):
        filename=Configuration
        self.conf = aseio.read(filename)
        self.atom_types = np.unique(self.conf.get_chemical_symbols())

    def chooseAndPlotConfiguration(self,Configuration):
        self.chooseConfiguration(Configuration)
        self.plotConf()
        
    def calculateEnvironments(self,lista,listb,cutoff):
        self.allEnvs = np.ndarray((0,),dtype=np.object)
        # Find all environments
        # Loop over input particles of type atom_type_1:
        myindeces=lista 
        nl = asenl.NeighborList((cutoff/2.0)*np.ones(self.conf.get_number_of_atoms()),skin=0.1,bothways=True)
        nl.update(self.conf)
        mat = self.conf.get_cell()
        for index in myindeces:
            neighbors, offsets = nl.get_neighbors(index)
            # Iterate over the neighbors of the current particle:
            Environment = Environments()
            Environment.myindex = index
            Environment.delta.shape = (0,3)
            for neigh, offset in zip(neighbors, offsets):
                if (neigh in listb): # Only if neighbor is in listb
                    delta = self.conf.positions[neigh] + np.dot(offset, mat) - self.conf.positions[index]
                    distance=np.linalg.norm(delta)
                    if (distance>1.e-10):
                        if (not self.fractionalFlag):
                            Environment.delta = np.append(Environment.delta,np.array([[delta[0], delta[1], delta[2]]])*0.1,axis=0)
                        else:
                            Environment.delta = np.append(Environment.delta,np.array([[delta[0]/mat[0][0], delta[1]/mat[1][1], delta[2]/mat[2][2]]]),axis=0)
                        Environment.distance = np.append(Environment.distance,distance)
                        Environment.indeces = np.append(Environment.indeces,neigh)
            self.allEnvs = np.append(self.allEnvs,Environment)

    def CalculateUniqueEnvironments(self):
        num_of_templates=self.allEnvs.shape[0]
        flag_unique=np.ones(num_of_templates)
        # Disregard empty environments
        for i in range(num_of_templates):
            if (self.allEnvs[i].indeces.shape[0]==0):
                flag_unique[i]=0
        for i in trange(num_of_templates,desc="Environment 1", leave=False):
            if ( self.allEnvs[i].indeces.shape[0]>0 and flag_unique[i]==1):
                for j in trange(i+1,num_of_templates,desc="Environment 2", leave=False):
                    # Not empty,  unique, and same number of neighbors
                    if ( self.allEnvs[j].indeces.shape[0]>0 and flag_unique[j]==1 and self.allEnvs[i].indeces.shape[0]==self.allEnvs[j].indeces.shape[0]):
                        # Compare against all permutations
                        #for perm in tqdm(permutations(np.arange(self.allEnvs[j].indeces.shape[0])),total=np.math.factorial(self.allEnvs[j].indeces.shape[0]),desc="Permutations of environment", leave=None):
                        for perm in permutations(np.arange(self.allEnvs[j].indeces.shape[0])):
                            p=np.array(perm)
                            # If both have the same vectors, then the second environment is not unique
                            if (np.count_nonzero(np.isclose(self.allEnvs[i].delta,self.allEnvs[j].delta[np.array(p),:],atol=self.tolerance))==self.allEnvs[i].indeces.shape[0]*3 ):
                                flag_unique[j]=0
                                break
        self.uniqueEnvs = np.ndarray((0,),dtype=np.object)
        for i in range(num_of_templates):
            if (flag_unique[i]==1):
                self.uniqueEnvs = np.append(self.uniqueEnvs,self.allEnvs[i])

    def CalculateUniqueEnvironmentsFast(self):
        # This is a fast but sometimes wrong algorithm to compare environments
        num_of_templates=self.allEnvs.shape[0]
        # Sort according to x+y+z,z,y,x
        for i in range(num_of_templates):
            sortindeces=np.argsort(self.allEnvs[i].delta[:,0]+self.allEnvs[i].delta[:,1]+self.allEnvs[i].delta[:,2],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,2],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,1],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,0],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
        flag_unique=np.ones(num_of_templates)
        # Disregard empty environments
        for i in range(num_of_templates):
            if (self.allEnvs[i].indeces.shape[0]==0):
                flag_unique[i]=0
        for i in range(num_of_templates):
            # Disregard this environment if it is empty or already not unique
            if ( self.allEnvs[i].indeces.shape[0]>0 and flag_unique[i]==1):
                for j in range(i+1,num_of_templates):
                    # Not empty,  unique, and same number of neighbors
                    if ( self.allEnvs[j].indeces.shape[0]>0 and flag_unique[j]==1 and self.allEnvs[i].indeces.shape[0]==self.allEnvs[j].indeces.shape[0]):
                        # If both have the same vectors, then the second environment is not unique
                        if (np.count_nonzero(np.isclose(self.allEnvs[i].delta,self.allEnvs[j].delta,atol=self.tolerance))==self.allEnvs[i].indeces.shape[0]*3 ):
                            flag_unique[j]=0
        self.uniqueEnvs = np.ndarray((0,),dtype=np.object)
        for i in range(num_of_templates):
            if (flag_unique[i]==1):
                self.uniqueEnvs = np.append(self.uniqueEnvs,self.allEnvs[i])

    def printEnvironmentSummaryInfo(self,Environments):
        num_of_templates=Environments.shape[0]
        print("Found " + str(num_of_templates) + " environments")

    def printEnvironments(self,Environments):
        # Print C++ syntax or plumed input syntax
        num_of_templates=Environments.shape[0]
        for i in range(num_of_templates):
            env_atom_types = np.asarray(self.conf.get_chemical_symbols())[Environments[i].indeces.astype(int)]
            env_positions = Environments[i].delta*10 
            env = ase.Atoms(env_atom_types,env_positions)
            aseio.write('-', env, format='proteindatabank')

                
    def plotEnv(self,myEnvironment): 
        env_atom_types = np.asarray(self.conf.get_chemical_symbols())[myEnvironment.indeces.astype(int)] #,np.array('C')] #Template[env_number].myindex)]
        env_atom_types = np.append(env_atom_types,np.array('Si'))
        env_positions = np.vstack((myEnvironment.delta*10,np.array([0,0,0]) ) ) 
        env = ase.Atoms(env_atom_types,env_positions)
        v = nglview.show_ase(env)
        boxy=widgets.Box([v],layout=self.box_layout)
        display(boxy)

    def chooseEnvPlotAll(self,number):
        if (self.allEnvs.shape[0]>0 and ((number-1)<self.allEnvs.shape[0])):
            self.plotEnv(self.allEnvs[number-1])
        else:
            print("Error")
    
    def chooseEnvPlotUnique(self,number):
        if (self.uniqueEnvs.shape[0]>0):
            self.plotEnv(self.uniqueEnvs[number-1])
        
    def chooseEnvPlot(self,number):
        if (self.uniqueFlag and self.uniqueEnvs.shape[0]>0):
            self.chooseEnvPlotUnique(number)
        elif (not(self.uniqueFlag) and self.allEnvs.shape[0]>0):
            self.chooseEnvPlotAll(number)
        else:
            print("Error")

    def calculateEnvironmentsType(self,atom_type_1,atom_type_2,cutoff,tolerance):
        self.tolerance=tolerance
        chemical_symbols = np.asarray(self.conf.get_chemical_symbols())
        lista=np.argwhere(chemical_symbols == atom_type_1).flatten()
        listb=np.argwhere(chemical_symbols == atom_type_2).flatten()
        self.calculateEnvironments(lista,listb,cutoff)
        if (self.uniqueFlag and not(self.fastFlag)):
            self.CalculateUniqueEnvironments()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        elif (self.uniqueFlag and self.fastFlag):
            self.CalculateUniqueEnvironmentsFast()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        else:
            self.printEnvironmentSummaryInfo(self.globalTemplate)

    def calculateEnvironmentsString(self,listastring,listbstring,cutoff,tolerance):
        self.tolerance=tolerance
        lista=np.fromstring(listastring, dtype=int, sep=',')-1
        listb=np.fromstring(listbstring, dtype=int, sep=',')-1
        self.calculateEnvironments(lista,listb,cutoff)
        if (self.uniqueFlag and not(self.fastFlag)):
            self.CalculateUniqueEnvironments()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        elif (self.uniqueFlag and self.fastFlag):
            self.CalculateUniqueEnvironmentsFast()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        else:
            self.printEnvironmentSummaryInfo(self.allEnvs)
    
    def calculateEnvironmentsMinMaxStride(self,mina,maxa,stridea,minb,maxb,strideb,cutoff,tolerance):
        self.tolerance=tolerance
        lista=np.arange(int(mina),int(maxa)+1,int(stridea))-1
        listb=np.arange(int(minb),int(maxb)+1,int(strideb))-1
        #print(lista,listb)
        self.calculateEnvironments(lista,listb,cutoff)
        if (self.uniqueFlag and not(self.fastFlag)):
            self.CalculateUniqueEnvironments()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        elif (self.uniqueFlag and self.fastFlag):
            self.CalculateUniqueEnvironmentsFast()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        else:
            self.printEnvironmentSummaryInfo(self.allEnvs)
    

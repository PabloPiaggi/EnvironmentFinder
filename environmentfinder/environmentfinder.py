import numpy as np
from ipywidgets import widgets
import ase
import ase.io
import ase.neighborlist
from itertools import permutations
from tqdm.notebook import tqdm, trange
import nglview
from IPython.display import display, FileLink
from zipfile import ZipFile
import os
import glob

import warnings
warnings.simplefilter('ignore')

# Class to define and operate on environments
class Environments:
    def __init__(self):
        self.myindex = -1
        self.indeces = np.array([],dtype=int)
        self.delta = np.array([])
        self.distance =  np.array([])
        # degeneracy defaults to 1 and is different from 1 for unique environments
        self.degeneracy = 1


# Class to find environments
class EnvironmentFinder:

    def __init__(self):
        # self.fractionalFlag is True if the coordinates are fractional
        self.fractionalFlag=False
        # uniqueFlags is True if unique environments are to be calculated
        self.uniqueFlag=True
        # self.fastFlag is True if the fast algorithm is to be employed
        self.fastFlag=True
        # self.allEnvs stores all calculated environments
        self.allEnvs=np.ndarray((0,),dtype=object)
        # self.uniqueEnvs stores the unique environments
        self.uniqueEnvs = np.ndarray((0,),dtype=object)
        # self.tolerance is the tolerance to determine if two elements of the
        # distance vectors are the same
        self.tolerance=0
        # self.atom_types contains the unique atom types
        self.atom_types = np.empty(1)
        # Max number of neighbors in environment
        self.max_neighbors=1
        self.max_neighbors_flag=False
        #env=Environments()
        #env.delta.shape = (0,3)
        #uniqueEnvs = np.append(uniqueEnvs,env)
        # self.box_layout contains options to display a box
        self.box_layout = widgets.Layout(display='flex',
                    flex_flow='column',
                    align_items='stretch',
                    border='solid',
                    width='90%')

    def plotConf(self):
        """ Plot configuration using nglview """

        v = nglview.show_ase(self.conf)
        v.add_representation("unitcell")
        boxy=widgets.Box([v],layout=self.box_layout)
        display(boxy)

    def chooseConfiguration(self,filename):
        """ Choose configuration by filename

        Args:
            filename (string): the file name

        Returns:
            no value
        """

        self.conf = ase.io.read(filename)
        self.atom_types = np.unique(self.conf.get_chemical_symbols())
        self.atom_types = np.append(self.atom_types,"Any")

    def chooseAndPlotConfiguration(self,filename):
        """ Choose and plot configuration by filename

        Args:
            filename (string): the file name

        Returns:
            no value
        """

        self.chooseConfiguration(filename)
        self.plotConf()

    def setMaximumNumberOfNeighbors(self,max_neighbors):
        self.max_neighbors_flag = True
        self.max_neighbors = max_neighbors

    def calculateEnvironments(self,lista,listb,cutoff):
        """ Find environments

        This function calculates all environments centered at the atom in
        lista using atoms in listb as neighbors if they are within the cutoff.

        Args:
            lista (numpy vector): list of central atoms
            listb (numpy vector): list of neighboring atoms
            cutoff (float): cut off for the calculation of neighbors

        Returns:
            no value
        """

        self.allEnvs = np.ndarray((0,),dtype=object)
        # Find all environments
        # Loop over input particles of type atom_type_1:
        myindeces=lista
        nl = ase.neighborlist.NeighborList((cutoff/2.0)*np.ones(self.conf.get_number_of_atoms()),skin=0.1,bothways=True)
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
            # Use maximum number of neighbors
            if (self.max_neighbors_flag==True):
                 # Sort
                 sortindeces=np.argsort(Environment.distance,kind='stable')
                 Environment.delta =    Environment.delta[sortindeces,:]
                 Environment.distance = Environment.distance[sortindeces]
                 Environment.indeces =  Environment.indeces[sortindeces]
                 # Choose first max_neighbors
                 Environment.delta =    Environment.delta[:self.max_neighbors,:]
                 Environment.distance = Environment.distance[:self.max_neighbors]
                 Environment.indeces =  Environment.indeces[:self.max_neighbors]
            self.allEnvs = np.append(self.allEnvs,Environment)


    def CalculateUniqueEnvironments(self):
        """ Calculate unique Environments

        This algorithm calculates unique environments by comparing all possible
        permutations of the environments. However, it discards entire sets of
        permutations if it is found that there is no possibility of match.

        The number of operations scales roughly as as N^2*(\sum_{i=1}^M M+1-i)
        where N is the total number of environments and M is the number of i
        atoms in the environments.

        Args:
            none

        Returns:
            no value
        """

        num_of_templates=self.allEnvs.shape[0]
        flag_unique=np.ones(num_of_templates)
        same_as=np.linspace(0,num_of_templates-1,num_of_templates)
        # Disregard empty environments
        for i in range(num_of_templates):
            if (self.allEnvs[i].indeces.shape[0]==0):
                flag_unique[i]=0
        for i in trange(num_of_templates,desc="Progress", leave=False):
            if ( self.allEnvs[i].indeces.shape[0]>0 and flag_unique[i]==1):
                for j in range(i+1,num_of_templates):
                    # Not empty,  unique, and same number of neighbors
                    if ( self.allEnvs[j].indeces.shape[0]>0 and flag_unique[j]==1 and self.allEnvs[i].indeces.shape[0]==self.allEnvs[j].indeces.shape[0]):
                        atom_match = np.zeros(self.allEnvs[j].indeces.shape[0])
                        for k1 in range(self.allEnvs[i].indeces.shape[0]):
                            for k2 in range(self.allEnvs[j].indeces.shape[0]):
                                if (np.count_nonzero(np.isclose(self.allEnvs[i].delta[k1,:],self.allEnvs[j].delta[k2,:],atol=self.tolerance))==3 ):
                                    #print(i,j,self.allEnvs[i].delta[k1,:],self.allEnvs[j].delta[k2,:],np.isclose(self.allEnvs[i].delta[k1,:],self.allEnvs[j].delta[k2,:],atol=self.tolerance))
                                    atom_match[k1] = 1
                                    break
                            # If no match was found for a given atom the whole self.allEnvs[j] must be discarded
                            if (atom_match[k1]==0):
                                break
                        # If a match was found for every atom in self.allEnvs[i]
                        #print(i,j,atom_match)
                        if np.all(atom_match==np.ones(atom_match.shape[0])):
                            flag_unique[j]=0
                            same_as[j] = i
                            #break
        self.uniqueEnvs = np.ndarray((0,),dtype=object)
        for i in range(num_of_templates):
            if (flag_unique[i]==1):
                self.uniqueEnvs = np.append(self.uniqueEnvs,self.allEnvs[i])
                self.uniqueEnvs[-1].degeneracy = np.sum(np.ones(num_of_templates)[same_as==i])

    def CalculateUniqueEnvironmentsOld(self):
        """ Calculate unique Environments

        This algorithm to calculate unique environments compares all possible
        permutations of the environments. It is thus always correct but often
        very slow.

        The number of operations scales roughly as as N^2 * M! where
        N is the total number of environments and M is the number of atoms in
        the environments.

        CalculateUniqueEnvironmentsFast is a faster alternative to this
        function.

        Args:
            none

        Returns:
            no value
        """

        num_of_templates=self.allEnvs.shape[0]
        flag_unique=np.ones(num_of_templates)
        same_as=np.linspace(0,num_of_templates-1,num_of_templates)
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
                                same_as[j] = i
                                break
        self.uniqueEnvs = np.ndarray((0,),dtype=object)
        for i in range(num_of_templates):
            if (flag_unique[i]==1):
                self.uniqueEnvs = np.append(self.uniqueEnvs,self.allEnvs[i])
                self.uniqueEnvs[-1].degeneracy = np.sum(np.ones(num_of_templates)[same_as==i])

    def CalculateUniqueEnvironmentsFast(self):
        """ Calculate unique Environments

        This algorithm to calculate unique environments first sorts the
        distance vectors (central atom to neighbor vectots) and then compares
        them. The algorithm is fast but it can be wrong in some pathological
        cases.

        The number of operations scales roughly as as N^2 * M * log(M)
        where N is the total number of environments and M is the number of
        atoms in the environments.

        CalculateUniqueEnvironments is a much slower alternative to this
        function but it is always correct.

        Args:
            none

        Returns:
            no value
        """

        # This is a fast but sometimes wrong algorithm to compare environments
        num_of_templates=self.allEnvs.shape[0]
        # Sort according to x+y+z,z,y,x
        for i in range(num_of_templates):
            sortindeces=np.argsort(self.allEnvs[i].delta[:,0]+self.allEnvs[i].delta[:,1]+self.allEnvs[i].delta[:,2],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            self.allEnvs[i].distance = self.allEnvs[i].distance[sortindeces]
            self.allEnvs[i].indeces = self.allEnvs[i].indeces[sortindeces]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,2],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            self.allEnvs[i].distance = self.allEnvs[i].distance[sortindeces]
            self.allEnvs[i].indeces = self.allEnvs[i].indeces[sortindeces]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,1],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            self.allEnvs[i].distance = self.allEnvs[i].distance[sortindeces]
            self.allEnvs[i].indeces = self.allEnvs[i].indeces[sortindeces]
            sortindeces=np.argsort(self.allEnvs[i].delta[:,0],kind='stable')
            self.allEnvs[i].delta = self.allEnvs[i].delta[sortindeces,:]
            self.allEnvs[i].distance = self.allEnvs[i].distance[sortindeces]
            self.allEnvs[i].indeces = self.allEnvs[i].indeces[sortindeces]
        flag_unique=np.ones(num_of_templates)
        same_as=np.linspace(0,num_of_templates-1,num_of_templates)
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
                            same_as[j] = i
        self.uniqueEnvs = np.ndarray((0,),dtype=object)
        for i in range(num_of_templates):
            if (flag_unique[i]==1):
                self.uniqueEnvs = np.append(self.uniqueEnvs,self.allEnvs[i])
                self.uniqueEnvs[-1].degeneracy = np.sum(np.ones(num_of_templates)[same_as==i])

    def printEnvironmentSummaryInfo(self,Environments):
        num_of_templates=Environments.shape[0]
        if (num_of_templates>0):
            avg_num_neighbors = 0.0
            for i in range(num_of_templates):
                avg_num_neighbors += Environments[i].indeces.shape[0]
            avg_num_neighbors /= num_of_templates
        else:
            avg_num_neighbors = 0
        print("Found " + str(num_of_templates) + " environments each with " + str(int(avg_num_neighbors))  + " neighbors on average")
        if (num_of_templates>0):
            for i in range(num_of_templates):
                print("Environment " + str(int(i+1)) + ": degeneracy = " + str(int(Environments[i].degeneracy)) + " - number of neighbors = ", str(int(Environments[i].indeces.shape[0])) )

    def printEnvironments(self,Environments):
        num_of_templates=Environments.shape[0]
        for i in range(num_of_templates):
            env_atom_types = np.asarray(self.conf.get_chemical_symbols())[Environments[i].indeces.astype(int)]
            env_positions = Environments[i].delta*10
            env = ase.Atoms(env_atom_types,env_positions)
            ase.io.write('-', env, format='proteindatabank')


    def printEnvironmentsToFile(self,Environments):
        num_of_templates=Environments.shape[0]
        for i in range(num_of_templates):
            env_atom_types = np.asarray(self.conf.get_chemical_symbols())[Environments[i].indeces.astype(int)]
            env_positions = Environments[i].delta*10
            env = ase.Atoms(env_atom_types,env_positions)
            fileName="env" + str(i+1) + ".pdb"
            ase.io.write(fileName, env, format='proteindatabank')

    def printEnvironmentsToZipFile(self,Environments):
        num_of_templates=Environments.shape[0]
        os.mkdir("Download")
        zipObj = ZipFile('Download/download.zip', 'w')
        for i in range(num_of_templates):
            env_atom_types = np.asarray(self.conf.get_chemical_symbols())[Environments[i].indeces.astype(int)]
            env_positions = Environments[i].delta*10
            env = ase.Atoms(env_atom_types,env_positions)
            fileName="Download/env" + str(i+1) + ".pdb"
            ase.io.write(fileName, env, format='proteindatabank')
            zipObj.write(fileName)
        zipObj.close()
        local_file = FileLink('Download/download.zip', result_html_prefix="Click here to download the environments in PDB format: ")
        display(local_file)
        fileList = glob.glob('Download/env*')
        for filePath in fileList:
            os.remove(filePath)
    

    def plotEnv(self,myEnvironment):
        env_atom_types = np.asarray(self.conf.get_chemical_symbols())[myEnvironment.indeces.astype(int)] #,np.array('C')] #Template[env_number].myindex)]
        env_atom_types = np.append(env_atom_types,np.asarray(self.conf.get_chemical_symbols())[myEnvironment.myindex])
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
            print("Degeneracy of the environment is " + str(int(self.uniqueEnvs[number-1].degeneracy)) )
            print("Number of neighbors in the environment is " + str(int(self.uniqueEnvs[number-1].indeces.shape[0])) )

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
        if (atom_type_1=="Any"):
            lista=np.linspace(0,chemical_symbols.shape[0]-1,chemical_symbols.shape[0]).astype(int)
        if (atom_type_2=="Any"):
            listb=np.linspace(0,chemical_symbols.shape[0]-1,chemical_symbols.shape[0]).astype(int)
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
        self.calculateEnvironments(lista,listb,cutoff)
        if (self.uniqueFlag and not(self.fastFlag)):
            self.CalculateUniqueEnvironments()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        elif (self.uniqueFlag and self.fastFlag):
            self.CalculateUniqueEnvironmentsFast()
            self.printEnvironmentSummaryInfo(self.uniqueEnvs)
        else:
            self.printEnvironmentSummaryInfo(self.allEnvs)

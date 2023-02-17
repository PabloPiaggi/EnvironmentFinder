import numpy as np
import os

from environmentfinder import EnvironmentFinder
# Define instance of class EnvironmentFinder
ef = EnvironmentFinder()
filename='../Examples/IceIh.pdb'
ef.chooseConfiguration(filename)
ef.setMaximumNumberOfNeighbors(16)
ef.calculateEnvironmentsType('O','O',5.5,0.25)
ef.printEnvironmentsToFile(ef.uniqueEnvs)


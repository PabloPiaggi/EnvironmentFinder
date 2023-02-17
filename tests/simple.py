import numpy as np
import os

from environmentfinder import EnvironmentFinder
# Define instance of class EnvironmentFinder
ef = EnvironmentFinder()
filename='../Examples/IceIh.pdb'
ef.chooseConfiguration(filename)
ef.calculateEnvironmentsType('O','O',3.0,0.02)
ef.printEnvironments(ef.uniqueEnvs)

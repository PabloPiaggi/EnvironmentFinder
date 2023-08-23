import numpy as np
import os

from environmentfinder import EnvironmentFinder
# Define instance of class EnvironmentFinder
ef = EnvironmentFinder()
ef.verboseFlag=False
filename='environmentfinder/examples/IceIh.pdb'
ef.chooseConfiguration(filename)
ef.calculateEnvironmentsType('O','O',3.0,0.02)

# Test number of unique environments
def test_env_num():
   assert ef.uniqueEnvs.shape[0] == 4, "Should be 4"

def test_env_num_neighbors():
   assert ef.uniqueEnvs[0].indeces.shape[0] == 4, "Should be 4"
   assert ef.uniqueEnvs[1].indeces.shape[0] == 4, "Should be 4"
   assert ef.uniqueEnvs[2].indeces.shape[0] == 4, "Should be 4"
   assert ef.uniqueEnvs[3].indeces.shape[0] == 4, "Should be 4"

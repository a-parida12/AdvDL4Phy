import tensorflow as tf
import numpy as np

sys.path.append("/home/a_parida/MantaFlow/manta/tensorflow/tools/")
import uniio
import glob
basePath = 'data/'

trainingEpochs = 10000
batchSize      = 10
inSize         = 64 * 64 * 1
######################################################
# Ex 2.1 – Saving and Loading Training Data
#####################################################

vel_files= glob.glob('./data/**/*.uni', recursive=True)
# load data
velocities = []

for uniPath in vel_files:
    header, content = uniio.readUni(uniPath)# returns [Z,Y,X,C] np array
    h = header['dimX']
    w  = header['dimY']
    arr = content[:, ::-1, :, :]  # reverse order of Y axis
    arr = np.reshape(arr, [w, h, 3])
    arr=arr[:,:,1:]# discard Z from [Z,Y,X]
    velocities.append( arr )
    
loadNum = len(densities)
if loadNum<200:
	print("Error - use at least two full sims, generate data by running 'manta ./manta_genSimSimple.py' a few times..."); exit(1)

velocities = np.reshape( velocities, (len(velocities), 64,64,2) )

print("Read uni files, total data " + format(velocities.shape) )
valiSize = max(100, int(loadNum * 0.1)) # at least 1 full sim...
valiData = velocities[loadNum-valiSize:loadNum,:]
velocities  = velocities[0:loadNum-valiSize,:]
print("Split into %d training and %d validation samples" % (velocities.shape[0], valiData.shape[0]) )
loadNum = velocities.shape[0]

############################################################

##############################################################
# Ex 2.2– First Network Architecture
##############################################################


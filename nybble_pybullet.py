import pybullet as p
import time
import pybullet_data
import numpy as np


physicsClient = p.connect(p.GUI) #, options="--width=1920 --height=1080 --mp4=\"test.mp4\" --mp4fps=30") #or p.DIRECT for non-graphical version
p.configureDebugVisualizer(p.COV_ENABLE_GUI, 0)
p.setAdditionalSearchPath(pybullet_data.getDataPath()) #optionally

p.setGravity(0,0,-9.81)
planeId = p.loadURDF("plane.urdf")
startPos = [0,0,0.07]
startOrientation = p.getQuaternionFromEuler([0,0,0])
boxId = p.loadURDF("nybble.urdf",startPos, startOrientation)
jointIds = []
paramIds = []


for j in range(p.getNumJoints(boxId)):
  info = p.getJointInfo(boxId, j)
  print(info)
  jointName = info[1]
  jointType = info[2]
  if (jointType == p.JOINT_PRISMATIC or jointType == p.JOINT_REVOLUTE):
    jointIds.append(j)
    paramIds.append(p.addUserDebugParameter(jointName.decode("utf-8")))


targetPos = np.ones(len(paramIds))
p.setJointMotorControlArray(boxId, jointIds, p.POSITION_CONTROL, targetPos)


gait = np.array([
    [26,  66, -66, -23,   3,  29, -29,  -5],
    [29,  70, -66, -19,   3,  21, -28,  -8],
    [33,  73, -66, -15,   1,  10, -27, -11],
    [38,  73, -66, -11,   0,  -1, -26, -13],
    [42,  68, -66,  -8,  -1, -11, -24, -16],
    [46,  57, -66,  -5,  -1, -14, -21, -19],
    [50,  44, -66,  -2,  -1, -15, -18, -22],
    [53,  29, -66,   0,   0,  -9, -16, -24],
    [56,  15, -65,   2,   1,   2, -13, -26],
    [59,   5, -64,   4,   3,  13, -10, -27],
    [61,  -2, -63,   5,   5,  23,  -7, -28],
    [63,  -5, -61,   5,   8,  29,  -5, -29],
    [64,  -4, -57,   1,  11,  28,  -3, -21],
    [65,  -3, -54,  -7,  13,  27,  -1, -10],
    [66,  -2, -52, -18,  16,  26,   0,   1],
    [66,   0, -49, -32,  19,  24,   1,  11],
    [66,   3, -45, -48,  22,  21,   1,  16],
    [66,   5, -41, -61,  24,  18,   1,  15],
    [66,   9, -37, -67,  26,  16,   0,   7],
    [66,  12, -32, -73,  27,  13,  -1,  -2],
    [66,  16, -28, -72,  28,  10,  -3, -13],
    [66,  20, -23, -69,  29,   7,  -5, -23],
    [70,  26, -19, -66,  21,   3,  -8, -29],
    [73,  29, -15, -66,  10,   3, -11, -28],
    [73,  33, -11, -66,  -1,   1, -13, -27],
    [68,  38,  -8, -66, -11,   0, -16, -26],
    [57,  42,  -5, -66, -14,  -1, -19, -24],
    [44,  46,  -2, -66, -15,  -1, -22, -21],
    [29,  50,   0, -66,  -9,  -1, -24, -18],
    [15,  53,   2, -66,   2,   0, -26, -16],
    [ 5,  56,   4, -65,  13,   1, -27, -13],
    [-2,  59,   5, -64,  23,   3, -28, -10],
    [-5,  61,   5, -63,  29,   5, -29,  -7],
    [-4,  63,   1, -61,  28,   8, -21,  -5],
    [-3,  64,  -7, -57,  27,  11, -10,  -3],
    [-2,  65, -18, -54,  26,  13,   1,  -1],
    [ 0,  66, -32, -52,  24,  16,  11,   0],
    [ 3,  66, -48, -49,  21,  19,  16,   1],
    [ 5,  66, -61, -45,  18,  22,  15,   1],
    [ 9,  66, -67, -41,  16,  24,   7,   1],
    [12,  66, -73, -37,  13,  26,  -2,   0],
    [16,  66, -72, -32,  10,  27, -13,  -1],
    [20,  66, -69, -28,   7,  28, -23,  -3]

])

gait = gait[:,[0,4,1,5,2,6,3,7]]

while(True):
    for i in range (np.size(gait, 0)):
        p.stepSimulation()
        targetPos = np.deg2rad(gait[i,:]) 
        p.setJointMotorControlArray(boxId, jointIds, p.POSITION_CONTROL, targetPos)
        time.sleep(5./240.)

cubePos, cubeOrn = p.getBasePositionAndOrientation(boxId)
print(cubePos,cubeOrn)
p.disconnect()


"""
Custom Static Optimization written on the basis of the official OpenSim example.
 Starting from experimental marker data (in .trc format) the optimal
 muscle activations are found that can reproduce the motion, solving:

   min   sum (w_i * a_i^2)+ sum (w_j * c_j^2)
   a,c    i                  j

   s.t.        0<=a_i<=1            for muscle activations
              -l<=c_j<=l            for coordinateActuators controls (if present)
        acc_{j,FD} = acc_{j,data}   constraint on accelerations


 The code is written specifically to consider a thoracoscapular shoulder
 model that has been already scaled to the biometrics of the subject of
 interest, and the data are loaded assuming to be considering a particular
 dataset (Waterloo's one). However, this script can be generalized to
 consider other models and data without changing its main structure.

 Please check the SETTINGS section, where you are asked if you want to
 include an external force (if so, specify the magnitude, the direction
 with respect to the subject, the body on which it is applied).

 author: Italo Belli (i.belli@tudelft.nl) May 2022
"""
import sys
from pathlib import Path
import opensim as osim
import numpy as np
import matplotlib.pyplot as plt
import time
from dataclasses import dataclass, field
from typing import Any, List
import scipy.optimize as sopt
import DataProcessing as dp
import Utilities_OSIM as ut_osim
import rref
from scipy import sparse
import osqp

# defining objective and constraint functions for the minimization problem --------------------------------------------
def obj_funct(x, coeff):
    epsilon = 0.01
    c = 10
    w = np.concatenate((np.ones([1, coeff]), epsilon*np.ones([1, 8]), c*np.ones([1, 9])), axis=1)
    cost = w.dot(np.square(x))
    return cost
    # return cost.reshape((1,))
    
def jntrxncon_linForce(x, params):
    directionVector = params.directionVector
    maxAngle = params.maxAngle
    A_f = params.A_f
    F_r0 = params.F_r0
    x = x.reshape((len(x), 1))
    # computing the reaction force vector at the given joint
    force_vec = np.matmul(A_f, np.atleast_2d(x)) + F_r0

    # evaluating the relative angle between the reaction force and Vec_H2GH
    angle = np.dot(directionVector, force_vec) / (np.linalg.norm(directionVector)*np.linalg.norm(force_vec))
    cosTheta = np.maximum(np.minimum(angle, 1), -1)
    rel_angle = np.real(np.rad2deg(np.arccos(cosTheta)))

    # value of the constraint violation
    c = np.square(rel_angle/maxAngle)-1
    ceq = 0
    return c, ceq

# Settings from the user (check every time!) ---------------------------------------------------------------------------
subject_considered = 'Ajay2019_'

dynamic_activation_bounds = True
apply_external_force = False
force_value = 0
force_direction = 'up'
body = 'hand'

motion_file_name = 'abd21.mot'
print_flag = False
withviz = False

weight_abd = 0.0001
weight_elev = 0.0001
weight_up_rot = 0.0001
weigth_wing = 0.0001

# Set the correct paths -------------------------------------------------------
path_to_repo = Path('C:/Users/italobelli/Desktop/GitHub/PTbot_officialCBL')
sys.path.append(path_to_repo)

# Load model and get initial state ---------------------------------------------
# Select the model you want to use, among the ones available in the
# OpenSim Models folder of this repository

modelFile = 'C:/Users/italobelli/Desktop/GitHub/PTbot_officialCBL/Personal_Results/TSM_Ajay2019_2kgWeight.osim'
model = osim.Model(modelFile)

# Select the trc file to be considered for Static Optimization
trcPath = 'C:/Users/italobelli/Desktop/Ajay old studies/FrontiersPubMaterials-latest/ThoracoscapularShoulderPaperMaterials/ThoracoscapularShoulderPaperMaterials/ExperimentalData/Markers/'
trc_file_used = 'ABD21.trc'

markersExp, timesExp, labelsExp, unitsExp = dp.readTRC(trcPath + trc_file_used)
start_time = timesExp[0]
end_time = timesExp[-1]

if unitsExp == 'mm':
    markersExp = markersExp/1000
    unitsExp = 'm'

frequency_trc_data = 1/(timesExp[1]-timesExp[0])

# Getting quantities about Glenohumeral joint --------------------------------------
model_temp = model.clone()    # create a temporary copy of the model, to be used in the tool
state = model_temp.initSystem()
maxAngle = ut_osim.get_glenoid_status(model_temp, state)[0]

# get the glenohumeral joint
alljoints = model.getJointSet()
glen = alljoints.get('GlenoHumeral')

#  getting the values of default scapula coordinate -------------------------------
# we get the values of the coordinates describing the scapula position from
# the general model in default pose
scapula_abd = model.getJointSet().get(2).get_coordinates(0)
scapula_ele = model.getJointSet().get(2).get_coordinates(1)
scapula_urt = model.getJointSet().get(2).get_coordinates(2)
scapula_wng = model.getJointSet().get(2).get_coordinates(3)

default_sa = scapula_abd.get_default_value()
default_se = scapula_ele.get_default_value()
default_su = scapula_urt.get_default_value()
default_sw = scapula_wng.get_default_value()

# Performing IK -------------------------------------------------------------------
# perform IK on the basis of marker data to retrieve the motion file for
# the coordinates of the model
ikSetupFile = str(path_to_repo) + '/Code/Scale Model/Utilities/setup_files/Standford2019_experiments/IKSetup_2019.xml'
ikTool = osim.InverseKinematicsTool(ikSetupFile)
ikTool.setMarkerDataFileName(trcPath + trc_file_used)
ikTool.setOutputMotionFileName(str(path_to_repo) + '/Personal_Results/' + motion_file_name)
ikTool.set_report_marker_locations(True)
ikTool.setStartTime(start_time)
ikTool.setEndTime(end_time)
ikTool.setModel(model_temp)

# set the reference values for the scapula coordinates (last 4 tasks)
num_IK_tasks = ikTool.getIKTaskSet().getSize()

# set the weight of each coordinate in the tracking tasks
ikTool.getIKTaskSet().get(num_IK_tasks-4).setWeight(weight_abd)
ikTool.getIKTaskSet().get(num_IK_tasks-3).setWeight(weight_elev)
ikTool.getIKTaskSet().get(num_IK_tasks-2).setWeight(weight_up_rot)
ikTool.getIKTaskSet().get(num_IK_tasks-1).setWeight(weigth_wing)

# set also the values here
osim.IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet().get(num_IK_tasks-4)).setValue(default_sa)
osim.IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet().get(num_IK_tasks-3)).setValue(default_se)
osim.IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet().get(num_IK_tasks-2)).setValue(default_su)
osim.IKCoordinateTask.safeDownCast(ikTool.getIKTaskSet().get(num_IK_tasks-1)).setValue(default_sw)
ikTool.printToXML('SO_autogenerated_IK_setup.xml')

ikTool.run()

del model_temp    # delete the temporary copy. It is just useful to avoid initializing/modifying the main model at this stage

# getting the kinematic data that we need ------------------------------------------
# Use the loadFilterCropArray() function provided by OpenSim Tutorial to load the
# coordinate kinematic and generalized force data into MATLAB arrays. This
# function also filters and crops the loaded array based on its two input
# arguments (more details in loadFilterCropArray.m).
lowpassFreq = 3.0   # Hz
timeRange = [start_time, end_time]

# get the coordinates from the output of the IK in rad
[coordinates, coordNames] = dp.loadFilterCropArray(str(path_to_repo) + '/Personal_Results/' + motion_file_name, lowpassFreq, timeRange)[0:2]
coordinates[:, 0:3] = np.deg2rad(coordinates[:, 0:3])
coordinates[:, 6:17] = np.deg2rad(coordinates[:, 6:17])

# get the velocities for each joint in rad/s
time_step_trc = timesExp[1]-timesExp[0]
speeds = np.zeros(coordinates.shape)
for i in range(0, coordinates.shape[1]-1):
    speeds[:, i] = np.gradient(coordinates[:, i], time_step_trc)
speedNames = coordNames

# get the accelerations for each coordinate in rad/s^2
accelerations = np.zeros(coordinates.shape)
for i in range(0, coordinates.shape[1]-1):
    accelerations[:,i] = np.gradient(speeds[:,i], time_step_trc)
accNames = speedNames

# check the values of joint states, speeds and accelerations
fig, axs = plt.subplots(4, 4)
count = 0
if print_flag:
    for i in range(0, 4):
        for j in range(0, 4):
            axs[i, j].plot(coordinates[:, count])
            axs[i, j].plot(speeds[:, count])
            axs[i, j].plot(accelerations[:, count])
            axs[i, j].set_title(coordNames[count])
            count = count + 1
            plt.legend(['coords', 'speeds', 'accs'])
    plt.show()

# Apply the external force to the model (if needed) -------------------------------
if apply_external_force:
    # decode the direction of the force
    if force_direction == 'right':
        force_direction_OS = '-Z'
    elif force_direction == 'left':
        force_direction_OS = '+Z'
    elif force_direction == 'up':
        force_direction_OS = '-Y'
    elif force_direction == 'down':
        force_direction_OS = '+Y'
    elif force_direction == 'push':
        force_direction_OS = '-Z'
    elif force_direction == 'pull':
        force_direction_OS = '+Z'

    # getting markers on the hand body
    # mcp5_traj = markersExp(:, 52: 54);
    # mcp2_traj = markersExp(:, 55: 57)

    # defining the application point as the center of the two markers on the hand
    # application_point = (mcp2_traj + mcp5_traj) / 2
    application_point = markersExp[:, 3:6]

    # generate the external force file
    external_force = ut_osim.generate_external_force(force_value, force_direction_OS, application_point, body, frequency_trc_data, 'TSM_external_force')

    # apply the external force to the model
    model.addForce(external_force)

    # reset the data source file for the force applied
    external_force_storage = osim.Storage('TSM_external_force.mot', False)
    n_forces = model.getForceSet().getSize()
    force_in_model = model.getForceSet().get(n_forces - 1)
    osim.ExternalForce.safeDownCast(force_in_model).setDataSource(external_force_storage)

    model.finalizeConnections()


# Store max isometric force values and disable muscle dynamics
muscles = model.getMuscles()
numMuscles = muscles.getSize()

muscleNames = []

for i in range(0, numMuscles):
    # Downcast base muscle to Millard2012EquilibriumMuscle
    muscle = osim.Millard2012EquilibriumMuscle.safeDownCast(muscles.get(i))
    muscleNames.append(muscle.getName())
    muscle.set_ignore_tendon_compliance(True)
    muscle.set_ignore_activation_dynamics(True)

if withviz:
    model.setUseVisualizer(True)

# Update the system to include any muscle modeling changes
state = model.initSystem()

# Get coordinate actuators
allActs = model.getActuators()
num_acts = allActs.getSize()
acts = []
act_names = []

# get all actuators and override actuation for the muscles only
for i in range(0, num_acts):
    act_names.append(allActs.get(i).getName())
    acts.append(osim.ScalarActuator.safeDownCast(allActs.get(i)))
    if i <= numMuscles:
        acts[i].overrideActuation(state, True)

# Perform static optimization ------------------------------------------------------
# We use FMINCON to solve the static optimization problem at selected time points.
# We set the 'timeInterval' variable to select the time points to be included in the
# optimization. For example, if set to 10, every 10th time point is selected. A
# time interval of 1 will select all available time points.
timeInterval = 10
time_step_SO = time_step_trc * timeInterval

# Update data arrays based on the time interval.
N = coordinates.shape[0]
coordinates = coordinates[0::timeInterval, :]
speeds = speeds[0::timeInterval, :]
accelerations = accelerations[0::timeInterval, :]
numTimePoints = coordinates.shape[0]

# SOLVER OPTION
opts = {'maxiter': 10000, 'ftol': 1E-3, 'disp': True, 'eps': 1E-6}
# opts = {'disp': True}
#            

# Construct initial guess and bounds arrays
numCoords = coordNames.shape[0]
numCoordActs = num_acts-numMuscles
k = 600
t_act = 0.01
t_deact = 0.04

w = np.concatenate((np.ones((1, numMuscles)), 0.01*np.ones((1, 8)), 10*np.ones((1, 9))), axis=1)
P = sparse.csc_matrix(2*np.diagflat(w))
lb = np.concatenate((np.zeros((1, numMuscles)), -k*np.ones((1, numCoordActs))), axis=1)
ub = np.concatenate((np.ones((1, numMuscles)), k*np.ones((1, numCoordActs))), axis=1)
x0 = np.concatenate((0.1*np.ones((1, numMuscles)), 1e-5*np.ones((1, numCoordActs))), axis=1).reshape(50, )

# Pre-allocate arrays to be filled in the optimization loop
fl = np.zeros((1, numMuscles))
fv = np.zeros((1, numMuscles))
fp = np.zeros((1, numMuscles))
cosPenn = np.zeros((1, numMuscles))
Fmax = np.zeros((1, numMuscles))
A_eq_acc = np.zeros((numCoords, num_acts))
A_eq_force = np.zeros((3, num_acts))
xsol = np.zeros((numTimePoints,))
simulatedAccelerations = np.zeros((numTimePoints, coordNames.shape[0]))
optimizationStatus = np.zeros((numTimePoints, 1))
norm_fv_in_ground = np.zeros((numTimePoints, 3))
norm_fv_rotated = np.zeros((numTimePoints, 3))
rel_angle = np.zeros((numTimePoints, 1))

# get model quantities we still need
coords = model.getCoordinateSet()

t_start = time.time()
solution = []
# enter in the optimization loop
for time_instant in range(0, numTimePoints):
    print('Optimizing...time step',  time_instant+1, '/', numTimePoints)

    # Loop through model coordinates to set coordinate values and speeds.
    # We set all coordinates to make sure we have the correct kinematic state
    # when compute muscle multipliers and moment arms.
    for j in range(0, coordNames.shape[0]):
        coord = coords.get(j)
        coord.setValue(state, coordinates[time_instant, j], False)
        coord.setSpeedValue(state, speeds[time_instant, j])

    # equilibrate the muscles to make them start in the correct state
    model.realizeVelocity(state)
    model.equilibrateMuscles(state)

    modelControls = model.getControls(state)

    # Populate the muscle multiplier arrays.To do this, we must realize the
    # system to the Velocity stage
    for k in range(0, numMuscles):
        muscle = osim.Millard2012EquilibriumMuscle.safeDownCast(muscles.get(k))
        fl[0, k] = muscle.getActiveForceLengthMultiplier(state)
        fv[0, k] = muscle.getForceVelocityMultiplier(state)
        fp[0, k] = muscle.getPassiveForceMultiplier(state)
        cosPenn[0, k] = muscle.getCosPennationAngle(state)
        Fmax[0, k] = muscle.getMaxIsometricForce()

    # get the vector Vac_H2GC between humeral head and the glenoid center it is expressed in the ground fire
    Vec_H2GC = ut_osim.get_glenoid_status(model, state)[1]

    # store the values of active and passive maximum force in the current configuration
    AMuscForce = np.transpose(fl*fv*Fmax*cosPenn)
    PMuscForce = np.transpose(Fmax*fp*cosPenn)

    # create a data class containing parameters to feed to the findInducedAccelerations function
    @dataclass
    class params_indAcc:
        acts: list = field(default_factory=list)
        model: Any = model
        state: Any = state
        AMuscForce: np.array = AMuscForce
        PMuscForce: np.array = PMuscForce
        coords: np.array = coords
        coordNames: np.array = coordNames
        muscles: Any = muscles
        useMuscles: bool = 1
        useControls: bool = 1
        modelControls: Any = modelControls
    params_indAcc.acts = acts

    # evaluate the matrices that allow the acceleration constraint to be written as a linear one
    q_ddot_0 = ut_osim.findInducedAccelerations(np.zeros((1, num_acts))[0], params_indAcc)
    delQ_delX = np.eye(num_acts)

    for k in range(0, num_acts):
        incrementalForceAccel = ut_osim.findInducedAccelerations(delQ_delX[k, :], params_indAcc)
        kthColumn = incrementalForceAccel - q_ddot_0
        A_eq_acc[:, k] = kthColumn[:, 0]

    Beq = np.transpose(np.atleast_2d(accelerations[time_instant, :])) - q_ddot_0

    # evaluate the matrices that allow the reaction force to be expressed as linear
    @dataclass
    class params_force:
        acts: list = field(default_factory=list)
        model: Any = model
        state: Any = state
        muscles: Any = muscles
        AMuscForce: np.array = AMuscForce
        PMuscForce: np.array = PMuscForce
        glen: Any = glen
        useControls: bool = 1
        modelControls: Any = modelControls
    params_force.acts = acts

    # Compute linearized expression for the joint reaction force, and equivalent force matrix to be used in the
    # optimization
    [F_r0, moment_r0] = ut_osim.findReactionForceAndMomentGH(np.zeros((1, num_acts))[0], params_force)

    for k in range(0, num_acts):
        [F_rk, moment_rk] = ut_osim.findReactionForceAndMomentGH(delQ_delX[k, :], params_force)
        kthColumn = F_rk-F_r0
        A_eq_force[:, k] = kthColumn

    @dataclass
    class params_nlc:
        directionVector: Any = Vec_H2GC
        F_r0: np.array = F_r0.reshape([3, 1])
        A_f: Any = A_eq_force
        maxAngle: float = maxAngle

    aux = (numMuscles)
    x0 = x0.reshape((len(np.transpose(x0)),))
    Beq = Beq.reshape((len(Beq),))


    # x_bounds = np.concatenate((np.transpose(lb), np.transpose(ub)), axis=1)
    # x_bounds_tuple = tuple(map(tuple, x_bounds))
    #
    # # reduce matrix A_eq_acc
    # A_eq_acc_reg, indexes, rows_rearranged = rref.rref(A_eq_acc)
    # rank = len(indexes)
    # A_eq_acc_reg = A_eq_acc_reg[0:rank, :]
    # Beq_reg = Beq[rows_rearranged]
    #
    # # set up the constraints
    # lc = sopt.LinearConstraint(A_eq_acc_reg, lb=Beq_reg * 0.95, ub=Beq_reg * 1.05,
    #                            keep_feasible=False)  # linear constraint on the accelerations
    # # nlc = sopt.NonlinearConstraint(jntrxncon_linForce(params_nlc), -np.inf, 0.0)
    #
    # nlc_ib = sopt.NonlinearConstraint(lambda x: jntrxncon_linForce(x, params_nlc), -1.0, 0.0)
    #
    # cons = [lc]
    # # cons = ({'type': 'eq', 'fun': lambda x_sym: np.matmul(A_eq_acc, x_sym) - Beq})
    # result = sopt.minimize(obj_funct, x0, args=aux, method='SLSQP', constraints=cons, bounds=x_bounds_tuple,
    #                        options=opts)
    #
    # x = result.x
    # solution.append(x)

    A_eq_osqp = sparse.csc_matrix(np.concatenate((A_eq_acc, np.eye(50)), axis=0))
    lb_osqp = np.concatenate((np.transpose(np.atleast_2d(Beq)), np.transpose(lb)), axis=0)
    ub_osqp = np.concatenate((np.transpose(np.atleast_2d(Beq)), np.transpose(ub)), axis=0)

    prob = osqp.OSQP()
    prob.setup(P=P, A=A_eq_osqp, l=lb_osqp, u=ub_osqp, verbose=True,
               adaptive_rho_interval=50)  # , adaptive_rho_fraction = 0)
    res = prob.solve()
    x = res.x

    # dynamically update the upper and lower bounds for the activations
    if dynamic_activation_bounds:
        for k in range(1, numMuscles):
            lb[:, k][0] = np.max((x[k]-x[k]*(0.5+1.5*x[k])*time_step_SO/t_deact, 0))
            ub[:, k][0] = np.min((x[k] + (1-x[k]) * time_step_SO / (t_act * (0.5 + 1.5*x[k])), 1))

    # retrieve the optimal accelerations
    simulatedAccelerations[time_instant, :] = ut_osim.findInducedAccelerations(x, params_indAcc).reshape(numCoords, )

    # retrieve the position of the joint reaction force on the approximated glenoid
    # (computing the reacton force vector at the given joint)
    force_vec = np.matmul(A_eq_force, x) + F_r0

    # evaluate the relative angle between the reaction force and Vec_H2GC
    cosTheta = np.max((np.min((np.dot(Vec_H2GC, force_vec)/(np.linalg.norm(Vec_H2GC)*np.linalg.norm(force_vec)), 1)), -1))
    rel_angle[time_instant] = np.arccos(cosTheta)*180/np.pi

    # evaluate the position on the glenoid where the reaction force is exerted
    norm_Vec_H2GC = Vec_H2GC/np.linalg.norm(Vec_H2GC)
    norm_fv_in_ground[time_instant, :] = force_vec/np.linalg.norm(force_vec)

    beta_angle = np.arctan(norm_Vec_H2GC[2]/norm_Vec_H2GC[0])
    alpha_angle = np.arctan(norm_Vec_H2GC[2]/(np.sin(beta_angle)*norm_Vec_H2GC[1]))

    Ry = np.array([np.cos(beta_angle), 0.0, np.sin(beta_angle), 0.0, 1.0, 0.0, -np.sin(beta_angle), 0.0, np.cos(beta_angle)]).reshape((3, 3))
    Rz = np.array([np.cos(alpha_angle), -np.sin(alpha_angle), 0.0, np.sin(alpha_angle), np.cos(alpha_angle), 0.0, 0.0, 0.0, 1.0]).reshape((3, 3))

    norm_fv_rotated[time_instant, :] = np.transpose(np.matmul(Rz, np.matmul(Ry, norm_fv_in_ground[time_instant, :])))

    if withviz:
        model.getVisualizer().show(state)

tOptim = time.time() - start_time

solution = np.concatenate(solution, axis=0).reshape(numTimePoints, num_acts)
# if needed plot results...

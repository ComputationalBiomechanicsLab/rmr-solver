import opensim as osim
import numpy as np
import math
import DataProcessing as dp

def get_glenoid_status(model, state):
    """
    This function returns parameters describing the status of the
    glenohumeral joint in the thoracoscapular model.
    
    INPUT
    :param model: opensim thoracoscapular model, that must be already provided
                  with markers on the glenoid center, humerus head and glenoid edge
                  (in this order, and they should be the last ones in the markerset)
    :param state: state of the model
    
    OUTPUT
    :return angle: the maximum angle representing the cone in which the reaction forces must
                   be contained is returned
    :return Vec_H2GC: 3D vector defined between the humeral head center (origin)
                      and the glenoid center. It is expressed in the ground frame
    """

    mkrs = model.getMarkerSet()
    nmkrs = mkrs.getSize()

    # manually hardcode the markers that we want(last 3 in the MarkerSet)
    G_Cent = mkrs.get(nmkrs - 3)
    HH_Cent = mkrs.get(nmkrs - 2)
    G_Edge = mkrs.get(nmkrs - 1)

    # get the location in ground of the 3 markers
    G_Cent_Loc = G_Cent.getLocationInGround(state).to_numpy()
    HH_Cent_Loc = HH_Cent.getLocationInGround(state).to_numpy()
    G_Edge_Loc = G_Edge.getLocationInGround(state).to_numpy()

    # define the vector from the glenoid center to the humerus head
    Vec_H2GC = G_Cent_Loc - HH_Cent_Loc

    # define the vector from the glenoid edge to the humerus head
    Vec_H2GE = G_Edge_Loc - HH_Cent_Loc

    # get the cosine of the angle between the two vectors
    CosTheta = max(min(np.dot(Vec_H2GC, Vec_H2GE) / (np.linalg.norm(Vec_H2GC) * np.linalg.norm(Vec_H2GE)), 1), -1)

    # find the maximum angle to be returned
    angle_rad = math.acos(CosTheta)
    angle = angle_rad * 180/math.pi

    # get additional informations about the glenohumeral joint
    Glenoid_Rad = np.linalg.norm(G_Cent_Loc-G_Edge_Loc)
    Head2Glen_dist = np.linalg.norm(HH_Cent_Loc-G_Cent_Loc)

    return angle, Vec_H2GC


def generate_external_force(force_magnitude, force_direction_in_ground, application_point_trajectory_in_ground, body_of_application, frequency, file_name):
    """
    This function allows to create and save an ExternalForce object given:
    :param force_magnitude: a (constant) value of force in N to be applied
    :param force_direction_in_ground: direction of the force expressed in the
           ground reference frame, with sign (string like '+X', '-Z', ...)
    :param application_point_trajectory_in_ground: array of Nx3, describing 3D
          trajectory (with N points) of the point at which the force is applied
    :param body_of_application: string with the name of the body to which the force
           will be applied
    :param frequency: frequency of sampling of the force value (in Hz)
    :param loads_name: name of the ExternalLoads object created
    :param file_name: name of the file to be created -> file_name.mot will store
           the data of the force, and file_name.xml will contain the
           external force object.

    The output is:
    :return external_loads: the ExternalForce object created

    Also, an XML file (saved in the working directory, as file_name.xml) is
    created, together with the storage file (file_name.mot).

    NOTE: as for now, you need to set again the Storage object that the
    forces references to (since the pointer gets lost otherwise):
    reset the data source file for the force applie - use the following lines 
    as an example of what needs to be done:
        external_force_storage = Storage('file_name.mot', false);
        n_forces = model.getForceSet().getSize();
        force_in_model = model.getForceSet().get(n_forces-1);
        ExternalForce.safeDownCast(force_in_model).setDataSource(external_force_storage);

    author: Italo Belli (i.belli@tudelft.nl) 2022
    """

    # Create the.mot file where forces and moments are stored
    # get the length of the trajectory
    len_trj = application_point_trajectory_in_ground.shape[0]

    # create the timestamps for the experiment
    time_stamps = list(range(0, len_trj)) * 1 / frequency
    time_stamps = time_stamps.reshape(time_stamps.shape[0], 1)

    # get the force direction in ground
    force_sign = force_direction_in_ground[0]
    force_axis = force_direction_in_ground[1]

    if force_sign == '-':
        force_magnitude = -force_magnitude

    force_vector = np.ones([len_trj, 1]) * force_magnitude
    force_moments_matrix = np.zeros([len_trj, 6])

    if force_axis == 'X':
        force_moments_matrix[:, 0] = force_vector[:, 0]
    elif force_axis == 'Y':
        force_moments_matrix[:, 1] = force_vector[:, 0]
    elif force_axis == 'Z':
        force_moments_matrix[:, 2] = force_vector[:, 0]

    data = np.hstack((np.array(time_stamps), force_moments_matrix, application_point_trajectory_in_ground))

    # create the motion file to be used in the ExternalForce object
    columnNames = ["F_x", "F_y", "F_z", "M_x", "M_y", "M_z", "p_x", "p_y", "p_z"]

    dp.writeMotStoData(data, 1, columnNames, file_name)

    mot_file_name = file_name + '.mot'

    # creating storage object
    force_storage = osim.Storage(mot_file_name, False)
    force_storage.setName(file_name)
    force_storage.printToFile(str(mot_file_name), 1/frequency, 'w')

    # create the external force object
    external_force = osim.ExternalForce()
    external_force.setName('Force')
    external_force.set_applied_to_body(body_of_application)
    external_force.set_force_expressed_in_body('ground')
    external_force.set_point_expressed_in_body('ground')
    external_force.set_force_identifier('F_')
    external_force.set_point_identifier('p_')
    external_force.set_torque_identifier('M_')
    external_force.set_data_source_name(mot_file_name)
    external_force.setDataSource(force_storage)

    xml_file_name = file_name + '.xml'
    external_force.printToXML(xml_file_name)

    # returning the external force object. This is not the best way to do it
    # since the pointer to the storage object gets lost and needs to be reset
    return external_force


def findInducedAccelerations(x, params):
    """
    This function returns the simulated accelerations for each coordinate of
    a Model, given the model state and the forces exerted by the Muscles (and
    CoordinateActuators). Most of the necessary inputs are passed as a single
    parameter struct, with the following fields:
    * model: the OpenSim Model considered
    * state: the State of the model
    * coords: the CoordinateSet of the model (retrieved as model.getCoordinateSet())
    * coordNames: the names of the occordinates coming from the model
    * acts: actuators in the model, where the i-th element is retrieved as
            acts(i+1) = ScalarActuator.safeDownCast(model.getActuators().get(i));
    * muscles: muscles present in the model, muscles = model.getMuscles()
    * useMuscles: flag to indicate whether the muscles are used (1) or they are ignored in the computation of the
                  accelerations (0)
    * AMuscForce: vector containing the muscle multipliers to find the muscle force dependent from the activation level
                  (for each muscle)
    * PMuscForce: vector containing the muscle multipliers to find the passive muscle force (for each muscle)
    * useControls: flag to indicate whether to use .setControls (1) for the CoordinateActuators, or to overwrite their
                   actuation with .setOverrideActuation (0). In the second case, remeber that the actuation should be
                   overrideable!
    * modelControls: model controls (model.getControls(state)) used if the
                 corresponding flag is set to 1.

    :param x: vector containing the muscle activations (for each muscle) and the controls for the CoordinateActuators (for each one of them)
    :param params: collecting all the parameters presented above
    :return simulatedAccelerations: accelerations for each of the coordinates are returned.

    author: Sagar Joshi (s.d.joshi@tudelft.nl) 2022
    """

    model = params.model
    state = params.state
    coords = params.coords
    coordNames = params.coordNames
    acts = params.acts
    muscles = params.muscles
    nMuscles = muscles.getSize()
    useMuscles = params.useMuscles
    useControls = params.useControls

    # initialize the muscles to produce the required forces
    if useMuscles:
        AMuscForce = params.AMuscForce
        PMuscForce = params.PMuscForce
        MuscleForces = AMuscForce*np.transpose(np.atleast_2d(x[0:nMuscles])) + PMuscForce
        for k in range(0, nMuscles):
            muscle = muscles.get(k)
            muscle.setOverrideActuation(state, MuscleForces[k, 0])

    # initialize the CoordinateActuators to produce the required effect
    if not(useControls):
        for k in range(nMuscles+1, np.size(acts)):
            acts[k].setOverrideActuation(state, x[k])
    else:
        modelControls = params.modelControls
        for k in range(nMuscles+1, np.size(acts)):
            acts[k].setControls(osim.Vector(1, x[k]), modelControls)
        model.realizeVelocity(state)
        model.setControls(state, modelControls)

    # realize the model to the acceleration stage
    model.realizeAcceleration(state)

    # retrieve the simulated accelerations for each coordinate
    simulatedAccelerations = np.zeros((np.size(coordNames), 1))

    for j in range(0, np.size(coordNames)):
        coord = coords.get(coordNames[j])
        simulatedAccelerations[j, 0] = coord.getAccelerationValue(state)

    return simulatedAccelerations


def findReactionForceAndMomentGH(x, params):
    """"
    This function returns the joint reaction force and moment at the gleno-humeral
    (GH) joint, given the model state and the forces exerted by the Muscles (and
    CoordindateActuators). Most of the necessary inputs are passed as a single
    parameter struct, with the followign fields:

    * model: the OpenSim Model considered;
    * state: the State of the model
    * acts: actuators in the model, where the i-th element is retrieved as acts(i+1) = ScalarActuator.safeDownCast(model.getActuators().get(i));
    * muscles: muscles present in the model, muscles = model.getMuscles()
    * AMuscForce: vector containing the muscle multipliers to find the muscle force dependent from the activation level (for each muscle)
    * PMuscForce: vector containing the muscle multipliers to find the passive muscle force (for each muscle)
    * glen: glenohumeral joint (model.getJointSet.get('GlenoHumeral');
    * useControls: flag to indicate whether to use .setControls (1) for the CoordinateActuators, or to overwrite their
                   actuation with .setOverrideActuation (0). In the second case, remember that the actuation should be overrideable!
    * modelControls: model controls (model.getControls(state)) used if the corresponding flag is set to 1.

    Overall, the input to the function are:
    :param x: vector containing the muscle activations (for each muscle) and the controls for the CoordinateActuators
              (for each one of them)
    :param params: collecting all the parameters presented above
    :return force, moment: reaction force and moment at the (GH) joint

    author: Italo Belli (i.belli@tudelft.nl) 2022
    """
    model = params.model
    state = params.state
    acts = params.acts
    muscles = params.muscles
    AMuscForce = params.AMuscForce
    PMuscForce = params.PMuscForce
    glen = params.glen
    useControls = params.useControls

    nMuscles = muscles.getSize()

    # initialize the muscles to produce the required forces
    MuscleForces = AMuscForce*np.transpose(np.atleast_2d(x[0:nMuscles])) + PMuscForce
    for k in range(0, nMuscles):
        muscle = muscles.get(k)
        muscle.setOverrideActuation(state, MuscleForces[k, 0])

    # initialize the CoordinateActuators to produce the required effect
    if not (useControls):
        for k in range(nMuscles + 1, np.size(acts)):
            acts[k].setOverrideActuation(state, x[k])
    else:
        modelControls = params.modelControls
        for k in range(nMuscles + 1, np.size(acts)):
            acts[k].setControls(osim.Vector(1, x[k]), modelControls)
        model.realizeVelocity(state)
        model.setControls(state, modelControls)

    # realize the model to the acceleration stage
    model.realizeAcceleration(state)

    # get moment and force at the GlenoHumeral joint
    moment = glen.calcReactionOnParentExpressedInGround(state).get(0).to_numpy()
    force = glen.calcReactionOnParentExpressedInGround(state).get(1).to_numpy()

    return force, moment

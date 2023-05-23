function [ModelOut] = UpdateAC(ModelIn)
% Function that updates the points the AC contstraint in scaled versions of
% the Thoracoscapular Shoulder Model as the native OpenSim Scale Tool does
% not adjust this in the scaling process
%
% INPUTS:
% model_in_name: Name of the .osim model that has been scaled using the
%                OpenSim scale tool that needs adjustments to the AC
%                constraint. Make sure this model is present in the current
%                directory.
%
% OUTPUTS:
% model_out_name: Name of the output model in which the points in the AC
%                 constraint have been adjusted. Model is saved
%                 automatically in the current directory.


% import opensim libraries
import org.opensim.modeling.*

% get constraints
AC = PointConstraint.safeDownCast(ModelIn.getConstraintSet().get('AC'));

% get PhysicalOffsetFrames in clavicle and scapula
ClavFrame = PhysicalOffsetFrame().safeDownCast(ModelIn.getBodySet().get('clavicle').getComponent('cp_in_clavicle'));
ScapFrame = PhysicalOffsetFrame().safeDownCast(ModelIn.getBodySet().get('scapula').getComponent('cp_in_scapula'));

% get Updated values from Offset Frames
ACclav = ClavFrame.get_translation();
ACscap = ScapFrame.get_translation();

% Set updated values
AC.set_location_body_1(ACclav);
AC.set_location_body_2(ACscap);

% Output updated model
ModelOut = ModelIn;

end
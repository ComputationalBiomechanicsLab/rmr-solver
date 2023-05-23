function writeMotStoData(Data,optionType, ColumnNames, FileName) 

% Adila Quesadilla Papaya
% Last Updated: November 5th 2010
% 
% This function takes the input: nRows, nColumns, range and Data
% and outputs the .mot file to be used as an example input motion file by
% OpenSim for the purpose of running the various toolboxes.
% 
%% INPUTS
% nRows = number of rows
% nColumns = number of columns
% range = [tstart tend]
% Data = Joint Coordinate in units of degrees 
% optionType: set to 1 for a .mot output file, set to 2 for a .sto output
% file.
% optionVel: set to 1 if we want that data included, else set to 0
%% OUTPUTS
% MotionFile.mot or MotionFile.sto motion file containing the Joint Coordinates/Angles to be
% used by OpenSim
% You can call the Motion File whatever you want -- just change the
% .../MotionFile_C_.mot or .../MotionFile_C_.sto part 

% Note: an easy way to get the path name, is to just copy and paste (or,
% drag and drop) a file from the same folder into the Matlab window. 
nRows = size(Data,1);
nColumns = size(Data,2);
if optionType == 1 %.mot file
    fid = fopen(strcat(FileName,'.mot'),'w');
    fprintf(fid, strcat(FileName, '\n'));
    fprintf(fid,'version=1 \n');
    fprintf(fid,'nRows=%d\n',nRows);
    fprintf(fid,'nColumns=%d\n',nColumns);
    fprintf(fid,'inDegrees=yes \n');
    fprintf(fid,'endheader\n');
    fprintf(fid,'time\t');
end

if optionType == 2 %.sto file
    fid = fopen(strcat(FileName,'.sto'),'w');
    fprintf(fid, strcat(FileName, '\n'));
    fprintf(fid,'version=1 \n');
    fprintf(fid,'nRows=%d\n',nRows);
    fprintf(fid,'nColumns=%d\n',nColumns);
    fprintf(fid,'inDegrees=yes \n');
    fprintf(fid,'endheader\n');
    fprintf(fid,'time\t');
end


%%%%%%%% Code to generate the first line detailing the DOF %%%%%%%%
for i = 1:(nColumns-1)
    fprintf(fid,ColumnNames{i});
    fprintf(fid,' \t ');
end

% fprintf(fid,'pelvic_tilt \t Abs_FE \t Abs_LB \t Abs_AR \t L5_S1_FE \t L5_S1_LB \t L5_S1_AR \t L4_L5_FE \t L4_L5_LB \t L4_L5_AR \t L3_L4_FE \t L3_L4_LB \t L3_L4_AR \t L2_L3_FE \t L2_L3_LB \t L2_L3_AR \t L1_L2_FE \t L1_L2_LB \t L1_L2_AR \t T12_L1_FE \t T12_L1_LB \t T12_L1_AR \t T11_T12_FE \t T11_T12_LB \t T11_T12_AR \t T10_T11_FE \t T10_T11_LB \t T10_T11_AR \t T9_T10_FE \t T9_T10_LB \t T9_T10_AR \t T8_T9_FE \t T8_T9_LB \t T8_T9_AR \t T7_T8_FE \t T7_T8_LB \t T7_T8_AR \t T6_T7_FE \t T6_T7_LB \t T6_T7_AR \t T5_T6_FE \t T5_T6_LB \t T5_T6_AR \t T4_T5_FE \t T4_T5_LB \t T4_T5_AR \t T3_T4_FE \t T3_T4_LB \t T3_T4_AR \t T2_T3_FE \t T2_T3_LB \t T2_T3_AR \t T1_T2_FE \t T1_T2_LB \t T1_T2_AR \t T1_head_neck_FE \t T1_head_neck_LB \t T1_head_neck_AR \t T12_r12R_X \t T12_r12R_Y \t T12_r12R_Z \t T11_r11R_X \t T11_r11R_Y \t T11_r11R_Z \t T10_r10R_X \t T10_r10R_Y \t T10_r10R_Z \t T9_r9R_X \t T9_r9R_Y \t T9_r9R_Z \t T8_r8R_X \t T8_r8R_Y \t T8_r8R_Z \t T7_r7R_X \t T7_r7R_Y \t T7_r7R_Z \t T6_r6R_X \t T6_r6R_Y \t T6_r6R_Z \t T5_r5R_X \t T5_r5R_Y \t T5_r5R_Z \t T4_r4R_X \t T4_r4R_Y \t T4_r4R_Z \t T3_r3R_X \t T3_r3R_Y \t T3_r3R_Z \t T2_r2R_X \t T2_r2R_Y \t T2_r2R_Z \t T1_r1R_X \t T1_r1R_Y \t T1_r1R_Z \t T12_r12L_X \t T12_r12L_Y \t T12_r12L_Z \t T11_r11L_X \t T11_r11L_Y \t T11_r11L_Z \t T10_r10L_X \t T10_r10L_Y \t T10_r10L_Z \t T9_r9L_X \t T9_r9L_Y \t T9_r9L_Z \t T8_r8L_X \t T8_r8L_Y \t T8_r8L_Z \t T7_r7L_X \t T7_r7L_Y \t T7_r7L_Z \t T6_r6L_X \t T6_r6L_Y \t T6_r6L_Z \t T5_r5L_X \t T5_r5L_Y \t T5_r5L_Z \t T4_r4L_X \t T4_r4L_Y \t T4_r4L_Z \t T3_r3L_X \t T3_r3L_Y \t T3_r3L_Z \t T2_r2L_X \t T2_r2L_Y \t T2_r2L_Z \t T1_r1L_X \t T1_r1L_Y \t T1_r1L_Z \t SternumX \t SternumY \t SternumZ \t SternumRotX \t SternumRotY \t SternumRotZ \t shoulder_elv_r \t shoulder_rot_r \t elv_angle_r \t elbow_flexion_r \t pro_sup_r \t wrist_dev_r \t wrist_flex_r \t shoulder_elv_l \t shoulder_rot_l \t elv_angle_l \t elbow_flexion_l \t pro_sup_l \t wrist_dev_l \t wrist_flex_l');
fprintf(fid, '\n');
for i = 1:nRows
    fprintf(fid,'%f',Data(i,1));
    fprintf(fid,'\t %f',Data(i,2:end));
    fprintf(fid, '\n');
end
fprintf(fid, '\n');


%If you need the velocities as well 
% if optionVel==1 
%         fprintf(fid,'flex_extension \t lat_bending \t axial_rotation\t');
%         fprintf(fid,'flex_extension_u \t lat_bending_u \t axial_rotation_u\t'); 	
%         fprintf(fid,'\n')
%         fprintf(fid,'\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n',Data);
% end

fclose(fid);

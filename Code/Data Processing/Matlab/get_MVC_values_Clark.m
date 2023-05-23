function array_MVC = get_MVC_values_Clark(path_to_CSV_subject_folder)
% This function can be used to compute the MVC values of a subject present
% in the Clark's dataset. It accepts as input:
% * path_to_CSV_subject_folder = the path to the folder containing all the
%     CVS files for a given subject in the dataset from Clark's. Such a
%     folder should be located in 'PTbot\DataShare\dataset_Clark\SMAP\SMAP 4
%     RAW'
% It returns as output:
% * array_MVC = it is a 14x1 array, containing the MVC value of each of the
%     muscles present in the experiments.
%
% Currently, this function is harcoded to consider a specific naming of the
% files in the dataset. It should be changed accordingly in case the naming
% of the files are changed. Could be improved with a for loop...

% The muscles that are present in the dataset are:
% 'Anterior Deltoid', 'Middle Deltoid', 'Posterior Deltoid', 'Biceps', 
% 'Triceps', 'Infra, 'Supra', 'PecMajorC', 'PecMajorS', 'Lats', 'Sert', 
% 'TrapLow', 'TrapMid', 'TrapUp'.
%
% The naming of the files produced after the MVC experiments are,
% correspondingly: 'ADEL*', 'MDEL*', 'PDEL*', 'BICE*', 'TRIC*', 'INFR*', 
% 'SUPR*', 'PECC*', 'PECS*', 'LATS*', 'SERA*', 'LTRA*', 'MTRA*', 'UTRA*'

array_MVC = zeros(14,1);

% 1. 'Anterior Deltoid' - 'ADEL*'
search_string = fullfile(path_to_CSV_subject_folder, 'ADEL*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Anterior Deltoid');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(1) = value_MVC;

% 2. 'Middle Deltoid' - 'MDEL*'
search_string = fullfile(path_to_CSV_subject_folder, 'MDEL*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Middle Deltoid');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(2) = value_MVC;

% 3. 'Posterior Deltoid' - 'PDEL*'
search_string = fullfile(path_to_CSV_subject_folder, 'PDEL*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Posterior Deltoid');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(3) = value_MVC;

% 4. 'Biceps' - 'BICE*'
search_string = fullfile(path_to_CSV_subject_folder, 'BICE*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Biceps');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(4) = value_MVC;

% 5. 'Triceps' - 'TRIC*'
search_string = fullfile(path_to_CSV_subject_folder, 'TRIC*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Triceps');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(5) = value_MVC;

% 6. 'Infra' - 'INFR*'
search_string = fullfile(path_to_CSV_subject_folder, 'INFR*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Infra');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(6) = value_MVC;

% 7. 'Supra' - 'SUPR*'
search_string = fullfile(path_to_CSV_subject_folder, 'SUPR*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Supra');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(7) = value_MVC;

% 8. 'PecMajorC' - 'PECC*'
search_string = fullfile(path_to_CSV_subject_folder, 'PECC*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'PecMajorC');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(8) = value_MVC;

% 9. 'PecMajorS' - 'PECS*'
search_string = fullfile(path_to_CSV_subject_folder, 'PECS*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'PecMajorS');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(9) = value_MVC;

% 10. 'Lats' - 'LATS*'
search_string = fullfile(path_to_CSV_subject_folder, 'LATS*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Lats');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(10) = value_MVC;

% 11. 'Sert' - 'SERA*'
search_string = fullfile(path_to_CSV_subject_folder, 'SERA*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'Sert');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(11) = value_MVC;

% 12. 'TrapLow' - 'LTRA*'
search_string = fullfile(path_to_CSV_subject_folder, 'LTRA*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'TrapLow');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(12) = value_MVC;

% 13. 'TrapMid' - 'MTRA*'
search_string = fullfile(path_to_CSV_subject_folder, 'MTRA*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'TrapMid');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(13) = value_MVC;


% 14. 'TrapUp' - 'UTRA*'
search_string = fullfile(path_to_CSV_subject_folder, 'UTRA*.CSV');
list_of_files = dir(search_string);

value_MVC = 0;

for k = 1: length(list_of_files)
    baseFileName = list_of_files(k).name;
    fullFileName = fullfile(list_of_files(k).folder, baseFileName); 
    new_MVC = get_MVC_value_per_muscle_Clark(fullFileName, 'TrapUp');
    value_MVC = max(value_MVC, new_MVC);
end

array_MVC(14) = value_MVC;

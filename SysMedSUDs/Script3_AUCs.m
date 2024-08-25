%% AUCs global measures

% Define the number of participants
numParticipants = numel(global_measures);

% Define unique groups from participantGroups
uniqueGroups = unique(participantGroups);

% Initialize variables to store AUC values for each measure and each group
Cglob_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
Tglob_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
Eglob_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
Pglob_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);

% Initialize variables to store all AUC values without separating into groups
Cglob_all = [];
Tglob_all = [];
Eglob_all = [];
Pglob_all = [];

% Compute AUC for each participant and each measure separated by group
for p = 1:numParticipants
    % Extract measures for participant p
    Clustering = global_measures{p}.Clustering;
    Transitivity = global_measures{p}.Transitivity;
    Efficiency = global_measures{p}.Efficiency;
    PathLength = global_measures{p}.PathLength;
    
    % Find the group name for participant p
    groupName = participantGroups{p};
    
    % Compute AUC values for each measure
    AUC_Clustering = trapz(Clustering);
    AUC_Transitivity = trapz(Transitivity);
    AUC_Efficiency = trapz(Efficiency);
    AUC_PathLength = trapz(PathLength);
    
    % Store AUC values in the appropriate group structure
    Cglob_AUC.(groupName) = [Cglob_AUC.(groupName); AUC_Clustering];
    Tglob_AUC.(groupName) = [Tglob_AUC.(groupName); AUC_Transitivity];
    Eglob_AUC.(groupName) = [Eglob_AUC.(groupName); AUC_Efficiency];
    Pglob_AUC.(groupName) = [Pglob_AUC.(groupName); AUC_PathLength];
    
    % Store AUC values in the overall structure
    Cglob_all = [Cglob_all; AUC_Clustering];
    Tglob_all = [Tglob_all; AUC_Transitivity];
    Eglob_all = [Eglob_all; AUC_Efficiency];
    Pglob_all = [Pglob_all; AUC_PathLength];
end

% Save AUC results for each measure separately
save('Cglob_AUC.mat', 'Cglob_AUC');
save('Tglob_AUC.mat', 'Tglob_AUC');
save('Eglob_AUC.mat', 'Eglob_AUC');
save('Pglob_AUC.mat', 'Pglob_AUC');

% Save AUC results without separating into groups
save('Cglob_all.mat', 'Cglob_all');
save('Tglob_all.mat', 'Tglob_all');
save('Eglob_all.mat', 'Eglob_all');
save('Pglob_all.mat', 'Pglob_all');

disp('AUC calculation for global measures separated by group completed and saved into separate files.');

%% AUC local measures separated by group and saved into separate files

% Define the number of participants
numParticipants = numel(local_measures);
numNodes = size(local_measures{1}.Clustering, 1);

% Define unique groups from participantGroups
uniqueGroups = unique(participantGroups);

% Initialize structures to store AUC values for each measure and each group
Cloc_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
Eloc_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
BCloc_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
ECloc_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);
Dloc_AUC = struct('Cannabis', [], 'Opiates', [], 'Nicotine', [], 'Controls', []);

% Initialize variables to store all AUC values without separating into groups
Cloc_all = [];
Eloc_all = [];
BCloc_all = [];
ECloc_all = [];
Dloc_all = [];

% Compute AUC for each participant and each local measure separated by group
for p = 1:numParticipants
    % Extract measures for participant p
    Clustering = local_measures{p}.Clustering;
    Efficiency = local_measures{p}.Efficiency;
    Betweenness = local_measures{p}.Betweenness;
    EigenvectorCentrality = local_measures{p}.EigenvectorCentrality;
    Degrees = local_measures{p}.Degrees;
    
    % Find the group name for participant p
    groupName = participantGroups{p};
    
    % Compute AUC values for each measure and each node
    AUC_Clustering = trapz(Clustering, 2);
    AUC_Efficiency = trapz(Efficiency, 2);
    AUC_Betweenness = trapz(Betweenness, 2);
    AUC_EigenvectorCentrality = trapz(EigenvectorCentrality, 2);
    AUC_Degrees = trapz(Degrees, 2);
    
    % Store AUC values in the appropriate group structure
    Cloc_AUC.(groupName) = [Cloc_AUC.(groupName); AUC_Clustering'];
    Eloc_AUC.(groupName) = [Eloc_AUC.(groupName); AUC_Efficiency'];
    BCloc_AUC.(groupName) = [BCloc_AUC.(groupName); AUC_Betweenness'];
    ECloc_AUC.(groupName) = [ECloc_AUC.(groupName); AUC_EigenvectorCentrality'];
    Dloc_AUC.(groupName) = [Dloc_AUC.(groupName); AUC_Degrees'];
    
    % Store AUC values in the overall structure
    Cloc_all = [Cloc_all; AUC_Clustering'];
    Eloc_all = [Eloc_all; AUC_Efficiency'];
    BCloc_all = [BCloc_all; AUC_Betweenness'];
    ECloc_all = [ECloc_all; AUC_EigenvectorCentrality'];
    Dloc_all = [Dloc_all; AUC_Degrees'];
end

% Save AUC results for each measure separately
save('Cloc_AUC.mat', 'Cloc_AUC');
save('Eloc_AUC.mat', 'Eloc_AUC');
save('BCloc_AUC.mat', 'BCloc_AUC');
save('ECloc_AUC.mat', 'ECloc_AUC');
save('Dloc_AUC.mat', 'Dloc_AUC');

% Save AUC results without separating into groups
save('Cloc_all.mat', 'Cloc_all');
save('Eloc_all.mat', 'Eloc_all');
save('BCloc_all.mat', 'BCloc_all');
save('ECloc_all.mat', 'ECloc_all');
save('Dloc_all.mat', 'Dloc_all');

disp('AUC calculation for local measures separated by group completed and saved into separate files.');

%% Clear
clear AUC_Betweenness AUC_Clustering AUC_Efficiency AUC_PathLength AUC_Betweenness AUC_EigenvectorCentrality AUC_Degrees AUC_Transitivity 

%% Load packages and data

% Define the path to the saved .mat file
dataFile = 'C:\Users\luisa\MATLAB\Projects\SysMedSUDs\new\allCorrelationMatrices_new.mat';

% Load the data from the .mat file
load(dataFile);

% Setup paths
path_BCT = 'C:\Users\luisa\MATLAB\packages\BCT';
path_GraphVar = 'C:\Users\luisa\MATLAB\packages\GraphVar_2.03a';
path_ComAlg = 'C:\Users\luisa\MATLAB\packages\communityalg-master';
addpath(genpath(path_BCT));
addpath(genpath(path_GraphVar));
addpath(genpath(path_ComAlg));


%% Thresholded Dataset
% Define the range of proportional thresholds
threshold_percentages = 0.1:0.01:0.5; 

% Initialize a structure to store thresholded matrices
thresholded_dataset = struct();

% Loop through each participant in the dataset
numParticipants = size(allCorrMatrices, 3);
for k = 1:numParticipants
    corrMat = allCorrMatrices(:, :, k); % Get the connectivity matrix for the current participant
    
    % Initialize a structure to store thresholded matrices for the current participant
    participant_name = ['Participant_' num2str(k)];
    thresholded_dataset.(participant_name) = struct();
    
    % Loop through each threshold percentage
    for i = 1:length(threshold_percentages)
        p = threshold_percentages(i); % Get the current threshold percentage
        
        % Apply proportional thresholding
        W_thresholded = threshold_proportional(corrMat, p);
        
        % Store the thresholded matrix in the participant's struct
        threshold_name = ['Threshold_' num2str(i)];
        thresholded_dataset.(participant_name).(threshold_name) = W_thresholded;
        
        % Display or use the thresholded matrix W_thresholded
        disp(['Thresholded matrix for ' participant_name ' at p = ' num2str(p)]);
    end
end

% Save the thresholded dataset
outputFile = 'C:\Users\luisa\MATLAB\Projects\SysMedSUDs\new\thresholded_allCorrMatrices.mat';
save(outputFile, 'thresholded_dataset');
disp('All thresholded data saved successfully.');

clear corrMat

%% Calculate graph measures for each participant

% Calculate the number of nodes from allCorrMatrices
numNodes = size(allCorrMatrices, 1); 
numParticipants = size(allCorrMatrices, 3);  
numThresholds = 41; 

% Initialize structures to store results
global_measures = cell(numParticipants, 1);
local_measures = cell(numParticipants, 1);

% Loop through each participant
for p = 1:numParticipants
    % Initialize arrays to store results for this participant
    Cglob_group = zeros(1, numThresholds);
    Tglob_group = zeros(1, numThresholds);
    Eglob_group = zeros(1, numThresholds);
    Pglob_group = zeros(1, numThresholds);
    
    Cloc_group = zeros(numNodes, numThresholds);
    Eloc_group = zeros(numNodes, numThresholds);
    BCloc_group = zeros(numNodes, numThresholds);
    ECloc_group = zeros(numNodes, numThresholds);
    Dloc_group = zeros(numNodes, numThresholds);
    
    % Loop through each threshold for this participant
    for t = 1:numThresholds
        % Retrieve correlation matrix for this threshold and participant
        corr_matrix_thresh_bin = thresholded_dataset.(['Participant_' num2str(p)]).(['Threshold_' num2str(t)]);
        
        % Global measures
        Cglob_group(t) = clusterMean_bu(corr_matrix_thresh_bin);
        Tglob_group(t) = transitivity_bu(corr_matrix_thresh_bin);
        Eglob_group(t) = efficiency_bin(corr_matrix_thresh_bin);
        Pglob_group(t) = charpath_B(corr_matrix_thresh_bin);
        
        % Local measures
        Cloc_group(:, t) = clustering_coef_bu(corr_matrix_thresh_bin);
        Eloc_group(:, t) = efficiency_local_bin(corr_matrix_thresh_bin);
        BCloc_group(:, t) = betweenness_bin(corr_matrix_thresh_bin);
        ECloc_group(:, t) = eigenvector_centrality_und(corr_matrix_thresh_bin);
         Dloc_group(:, t) = degrees_und(corr_matrix_thresh_bin);
    end
    
    % Store results for this participant in cell arrays
    global_measures{p} = struct(...
        'Clustering', Cglob_group, ...
        'Transitivity', Tglob_group, ...
        'Efficiency', Eglob_group, ...
        'PathLength', Pglob_group);
    
    local_measures{p} = struct(...
        'Clustering', Cloc_group, ...
        'Efficiency', Eloc_group, ...
        'Betweenness', BCloc_group, ...
        'EigenvectorCentrality', ECloc_group, ...
        'Degrees', Dloc_group);
end

% Save global and local measures to files
save('global_measures.mat', 'global_measures');
save('local_measures.mat', 'local_measures');

% Save numNodes for reference
save('numNodes.mat', 'numNodes');
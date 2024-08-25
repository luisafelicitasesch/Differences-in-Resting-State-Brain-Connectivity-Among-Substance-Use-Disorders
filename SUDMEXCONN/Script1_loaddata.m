%% Load and Standardize Correlation Matrices
% Define the directory containing the data files
dataDir = 'C:\Users\luisa\MATLAB\Projects\sudmex_conn\corrmatfinal';

% Get a list of all .mat files in the directory
dataFiles = dir(fullfile(dataDir, '*.mat'));

% Check if dataFiles is not empty
if isempty(dataFiles)
    error('No .mat files found in the specified directory: %s', dataDir);
end

% Loop through each file and load the data to get dimensions
sampleData = load(fullfile(dataDir, dataFiles(1).name));
if isfield(sampleData, 'CorrMat')
    [numNodes, ~] = size(sampleData.CorrMat);
else
    error('The field "CorrMat" is not found in the sample file: %s', dataFiles(1).name);
end

% Initialize storage for data and participant names
numParticipants = length(dataFiles);
allCorrMatrices = zeros(numNodes, numNodes, numParticipants);  % 3D array to store data matrices
participantNames = cell(numParticipants, 1);  % Cell array to store participant names or IDs

% Load and standardize the correlation matrices
for i = 1:numParticipants
    % Load the current .mat file
    filePath = fullfile(dataDir, dataFiles(i).name);
    data = load(filePath);
    
    % Check if the file contains the expected field
    if isfield(data, 'CorrMat')
        corrMat = data.CorrMat;
    else
        error('The field "CorrMat" is not found in file %s', dataFiles(i).name);
    end
    
    % Standardize the correlation matrix
    corrMat_mean = mean(corrMat(:));
    corrMat_std = std(corrMat(:));
    
    % Check if the standard deviation is zero
    if corrMat_std == 0
        warning('Standard deviation is zero for participant %s. Skipping standardization.', dataFiles(i).name);
        corrMat_standardized = corrMat;
    else
        corrMat_standardized = (corrMat - corrMat_mean) / corrMat_std;
    end
    
    % Set the diagonal elements to zero
    corrMat_standardized(1:numNodes+1:end) = 0;
    
    % Store the data in the allCorrMatrices 3D array
    allCorrMatrices(:, :, i) = corrMat_standardized;
    
    % Extract and store participant name or ID from the filename
    [~, filename, ~] = fileparts(dataFiles(i).name);
    participantNames{i} = filename;  % Store participant name or ID
    
    % Debug: Display the processed participant ID
    fprintf('Processed participant: %s\n', filename);
end

% Save the data in the specified output file
outputFile = 'C:\Users\luisa\MATLAB\Projects\sudmex_conn\allCorrelationMatrices_new.mat';
save(outputFile, 'allCorrMatrices', 'participantNames');
disp('All data saved successfully.');

%% Check if corrmatrices are complete

% Assuming allCorrMatrices already exists in the workspace
[numNodes, ~, numParticipants] = size(allCorrMatrices);

% Check if all matrices have values for all 246x246 fields
missingValues = false(numParticipants, 1);
for i = 1:numParticipants
    if any(isnan(allCorrMatrices(:, :, i)), 'all')
        missingValues(i) = true;
    end
end

% Display participants with missing values
if any(missingValues)
    fprintf('The following participants have missing values:\n');
    disp(find(missingValues));
else
    disp('All participants have complete 246x246 matrices.');
end


%% Assigning new groups
% Assuming 'remainingParticipants' table already exists in the workspace with columns:
% 'participant_id', 'group', and 'tobclastyear'

% Initialize the vector to store group codes
groupCodes = zeros(height(participants), 1);

% Loop through each participant and assign group codes
for i = 1:height(participants)
    group = string(participants.group(i)); % Convert group to string
    if group == "1"
        if participants.tobclastyear(i) == 1
            groupCodes(i) = 2; % Group 1, smoked in the last year
        else
            groupCodes(i) = 1; % Group 1, did not smoke in the last year
        end
    elseif group == "2"
        groupCodes(i) = 3; % Group 2
    end
end

% Display the group codes
disp(groupCodes);

% Optionally, you can save the group codes vector if needed
% save('groupCodes.mat', 'groupCodes');

% Initialize statistics structure
groupStats = struct();

% List of variables to calculate statistics for
variables = {'age', 'income', 'educ_yr', 'sex'};

% Loop through each group and calculate statistics for each variable
for group = 1:3
    % Filter participants by group
    groupIdx = (groupCodes == group);
    groupData = participants(groupIdx, :);
    
    % Initialize statistics for the current group
    groupStats(group).group = group;
    
    for varIdx = 1:length(variables)
        varName = variables{varIdx};
        if ismember(varName, groupData.Properties.VariableNames)
            data = groupData.(varName);
            
            % Calculate statistics
            groupStats(group).(varName).mean = mean(data, 'omitnan');
            groupStats(group).(varName).min = min(data, [], 'omitnan');
            groupStats(group).(varName).max = max(data, [], 'omitnan');
            groupStats(group).(varName).std = std(data, 'omitnan');
            groupStats(group).(varName).median = median(data, 'omitnan');
        else
            warning('Variable %s not found in the participants table.', varName);
        end
    end
end

% Display the statistics
for group = 1:3
    disp(['Statistics for Group ' num2str(group) ':']);
    disp(groupStats(group));
end

% Optionally, save the group statistics to a MAT-file
save('groupStats.mat', 'groupStats');
disp('Group statistics saved successfully.');


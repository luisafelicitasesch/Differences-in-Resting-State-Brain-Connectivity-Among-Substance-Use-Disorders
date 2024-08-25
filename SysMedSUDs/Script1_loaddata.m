%% Loading Data and Standardizing Matrices

% Define the directory containing the data files
dataDir = 'C:\Users\luisa\MATLAB\Projects\SysMedSUDs\new\alldata';

% Get a list of all .mat files in the directory
dataFiles = dir(fullfile(dataDir, '*.mat'));

% Load the first file to get the number of nodes
firstFile = load(fullfile(dataDir, dataFiles(1).name));
firstFieldName = fieldnames(firstFile);
sampleData = firstFile.(firstFieldName{1});
[numNodes, ~] = size(sampleData);

% Initialize storage for data and group labels
numParticipants = length(dataFiles);
allCorrMatrices = zeros(numNodes, numNodes, numParticipants);  
participantNames = cell(numParticipants, 1);  
participantGroups = cell(numParticipants, 1);  

% Define group labels and numbers
groupLabelMap = containers.Map({'Opiates', 'Cannabis', 'Nicotine', 'Controls', 'Unknown'}, [1, 2, 3, 4, 0]);

validFilesCount = 0;
for i = 1:numParticipants
    % Load the current .mat file
    filePath = fullfile(dataDir, dataFiles(i).name);
    fileData = load(filePath);
    fieldName = fieldnames(fileData);
    corrMat = fileData.(fieldName{1});
    
    % Check if the matrix is square
    [rows, cols] = size(corrMat);
    if rows ~= cols
        error('The file %s does not contain a square matrix.', dataFiles(i).name);
    end
    
    % Z-standardize the correlation matrix
    corrMat_mean = mean(corrMat(:));
    corrMat_std = std(corrMat(:));
    corrMat = (corrMat - corrMat_mean) / corrMat_std;
    
    % Set the diagonal elements to zero
    corrMat(1:numNodes+1:end) = 0;
    
    % Store the data in the allCorrMatrices 3D array
    validFilesCount = validFilesCount + 1;
    allCorrMatrices(:, :, validFilesCount) = corrMat;
    
    % Extract participant name or ID from the filename
    [~, filename, ~] = fileparts(dataFiles(i).name);
    participantNames{validFilesCount} = filename;  
    
    % Determine participant group based on filename or other criteria
    if contains(filename, 'BO') || contains(filename, 'MO')
        participantGroups{validFilesCount} = 'Opiates';
    elseif contains(filename, 'BC') || contains(filename, 'MC')
        participantGroups{validFilesCount} = 'Cannabis';
    elseif contains(filename, 'BR') || contains(filename, 'MR')
        participantGroups{validFilesCount} = 'Nicotine';
    elseif contains(filename, 'BG') || contains(filename, 'MG')
        participantGroups{validFilesCount} = 'Controls';
    else
        participantGroups{validFilesCount} = 'Unknown';  
    end
end

% Truncate the arrays to the valid files count
allCorrMatrices = allCorrMatrices(:, :, 1:validFilesCount);
participantNames = participantNames(1:validFilesCount);
participantGroups = participantGroups(1:validFilesCount);


%% Validation of Correlation Matrices: Size and NaN Value Chcks

% Initialize variables for storing results
numParticipants = size(allCorrMatrices, 3);
matrixSizeCheck = true(numParticipants, 1);  
nanCheck = true(numParticipants, 1);    

% Loop through each participant's correlation matrix
for i = 1:numParticipants
    % Extract the current correlation matrix
    corrMat = allCorrMatrices(:, :, i);
    
    % Check if the matrix is square
    [rows, cols] = size(corrMat);
    if rows ~= cols
        matrixSizeCheck(i) = false;
    end
    
    % Check for any NaN values
    if any(isnan(corrMat(:)))
        nanCheck(i) = false;
    end
end

% Report results
if all(matrixSizeCheck)
    disp('All matrices are square.');
else
    disp('Some matrices are not square.');
    disp(find(~matrixSizeCheck));  
end

if all(nanCheck)
    disp('No missing (NaN) values found in any matrices.');
else
    disp('Some matrices contain missing (NaN) values.');
    disp(find(~nanCheck));  
end


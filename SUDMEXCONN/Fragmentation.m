graphvar_path = "C:\Users\luisa\MATLAB\packages\GraphVar_2.03a"
addpath(genpath(graphvar_path))

%% Initialize a Structure and Loop Through Participants for Fragmentation Analysis

% Initialize a structure to store fragmentation info
fragmentation_info = struct();

% Number of participants
numParticipants = length(fieldnames(thresholded_dataset));

% Loop through each participant
for participant = 1:numParticipants
    participantField = sprintf('Participant_%d', participant);
    % Get the thresholded matrices for the current participant
    participantMatrices = thresholded_dataset.(participantField);
    
    % Get the number of thresholds
    numThresholds = length(fieldnames(participantMatrices));
    
    % Loop through each thresholded matrix
    for threshold = 1:numThresholds
        thresholdField = sprintf('Threshold_%d', threshold);
        adjMatrix = participantMatrices.(thresholdField);
        
        % Apply CheckFrag to determine if the network is fragmented
        isFragmented = CheckFrag(adjMatrix);
        
        % Convert adjacency matrix to a graph object
        G = graph(adjMatrix);

        % Find connected components
        bins = conncomp(G);

        % Calculate the size of each component
        componentSizes = histcounts(bins, 'BinMethod', 'integers');
        
        % Identify the largest connected component
        [maxSize, maxIdx] = max(componentSizes);
        
        % Store the fragmentation information
        fragmentation_info.(participantField).(thresholdField).isFragmented = isFragmented;
        fragmentation_info.(participantField).(thresholdField).largestComponentSize = maxSize;
        fragmentation_info.(participantField).(thresholdField).componentSizes = componentSizes;
        fragmentation_info.(participantField).(thresholdField).largestComponent = find(bins == maxIdx);
        
        % Check if the largest component size is less than 90% of the total nodes
        if maxSize < 0.9 * 247
            fprintf('Participant %d, Threshold %d: Largest component size is %d, indicating significant fragmentation.\n', ...
                    participant, threshold, maxSize);
        end
    end
end

% Save the fragmentation information for further analysis
save('fragmentation_info.mat', 'fragmentation_info');
disp('Fragmentation information saved successfully.');

%% Load Fragmentation Data and Analyze Significant Fragmentation Across Thresholds

% Load the fragmentation information
load('fragmentation_info.mat');

% Number of participants and thresholds
numParticipants = length(fieldnames(thresholded_dataset));
numThresholds = length(fieldnames(thresholded_dataset.Participant_1));

% Initialize an array to count the number of participants with significant fragmentation at each threshold
significantFragmentationCount = zeros(1, numThresholds);

% Loop through each participant and threshold to count significant fragmentation
for participant = 1:numParticipants
    participantField = sprintf('Participant_%d', participant);
    
    for threshold = 1:numThresholds
        thresholdField = sprintf('Threshold_%d', threshold);
        
        % Get the size of the largest component for the current participant and threshold
        largestComponentSize = fragmentation_info.(participantField).(thresholdField).largestComponentSize;
        
        % Check if the largest component size is less than 222
        if largestComponentSize < 222
            significantFragmentationCount(threshold) = significantFragmentationCount(threshold) + 1;
        end
    end
end

% Plot the count of participants with significant fragmentation at each threshold
thresholdLevels = 1:numThresholds;
figure;
bar(thresholdLevels, significantFragmentationCount);
xlabel('Threshold Level');
ylabel('Number of Participants with Significant Fragmentation');
title('Number of Participants with Largest Component Size < 222 at Each Threshold');
grid on;

% Save the plot as an image
saveas(gcf, 'fragmentation_plot.png');
disp('Fragmentation plot saved successfully.');


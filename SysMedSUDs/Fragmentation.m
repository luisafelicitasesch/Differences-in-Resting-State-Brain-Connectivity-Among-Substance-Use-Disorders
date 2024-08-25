%% Calculate largest components under 90% fragementation
% Initialize a structure to store fragmentation info
fragmentation_info_new = struct();

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
        fragmentation_info_new.(participantField).(thresholdField).isFragmented = isFragmented;
        fragmentation_info_new.(participantField).(thresholdField).largestComponentSize = maxSize;
        fragmentation_info_new.(participantField).(thresholdField).componentSizes = componentSizes;
        fragmentation_info_new.(participantField).(thresholdField).largestComponent = find(bins == maxIdx);
        
        % Check if the largest component size is less than 90% of the total nodes
        if maxSize < 0.9 * 247
            fprintf('Participant %d, Threshold %d: Largest component size is %d, indicating significant fragmentation.\n', ...
                    participant, threshold, maxSize);
        end
    end
end

% Save the fragmentation information for further analysis
save('fragmentation_info_new.mat', 'fragmentation_info_new');
disp('Fragmentation information saved successfully.');

%% Plot largest components under 90% fragementation
% Number of participants and thresholds
numParticipants = length(fieldnames(thresholded_dataset));
numThresholds = length(fieldnames(thresholded_dataset.Participant_1));

% Initialize an array to count the number of participants with significant fragmentation at each threshold
significantFragmentationCount_new = zeros(1, numThresholds);

% Loop through each participant and threshold to count significant fragmentation
for participant = 1:numParticipants
    participantField = sprintf('Participant_%d', participant);
    
    for threshold = 1:numThresholds
        thresholdField = sprintf('Threshold_%d', threshold);
        
        % Get the size of the largest component for the current participant and threshold
        largestComponentSize = fragmentation_info_new.(participantField).(thresholdField).largestComponentSize;
        
        % Check if the largest component size is less than 222
        if largestComponentSize < 222
            significantFragmentationCount_new(threshold) = significantFragmentationCount_new(threshold) + 1;
        end
    end
end

% Plot the count of participants with significant fragmentation at each threshold
thresholdLevels = 1:numThresholds;

figure;
bar(thresholdLevels, significantFragmentationCount_new, 'FaceColor', [0.6 0.8 1.0]);

% Customize the plot appearance
xlabel('Threshold Level', 'FontSize', 10);
ylabel('Number of Participants with Significant Fragmentation', 'FontSize', 10);
title('Number of Participants with Largest Component Size < 222 at Each Threshold', 'FontSize', 12);

% Remove the grid and the top and right lines
ax = gca;
ax.GridLineStyle = 'none';
ax.Box = 'off'; % Removes the top and right lines

% Adjust font size
ax.XAxis.FontSize = 8;
ax.YAxis.FontSize = 8;

% Save the plot as an image
saveas(gcf, 'fragmentation_plot_new.png');
disp('Fragmentation plot saved successfully.');

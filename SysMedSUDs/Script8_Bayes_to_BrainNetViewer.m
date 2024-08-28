% Load the Brainnetome Atlas subregion data
% Adjust the path to where you have saved the Excel file
filename = % add path to BNA subregions file from Brainnetome website;
% Read the Excel file
T = readtable(filename);

% Initialize cell arrays to store the combined data
combined_coords = {};
combined_labels = {};
combined_ids = {};

% Loop through each row of the table
for i = 1:height(T)
    % Extract the left hemisphere data
    lh_coords_str = T{i, 'lh_MNI_X_Y_Z_'};
    lh_coords = str2num(lh_coords_str{1}); % Convert string to numeric array
    lh_label = strrep(T{i, 'LeftAndRightHemisphere'}{1}, '(R)', ''); % Remove (R)
    lh_label = [lh_label, '_L']; % Append _L
    lh_label = regexprep(lh_label, '^"|"$', ''); % Remove leading and trailing quotes
    lh_id = T{i, 'LabelID_L'};
    
    % Extract the right hemisphere data
    rh_coords_str = T{i, 'rh_MNI_X_Y_Z_'};
    rh_coords = str2num(rh_coords_str{1}); % Convert string to numeric array
    rh_label = strrep(T{i, 'LeftAndRightHemisphere'}{1}, '(R)', ''); % Remove (R)
    rh_label = [rh_label, '_R']; % Append _R
    rh_label = regexprep(rh_label, '^"|"$', ''); % Remove leading and trailing quotes
    rh_id = T{i, 'LabelID_R'};
    
    % Combine the data
    combined_coords{end+1, 1} = lh_coords; %#ok<*SAGROW>
    combined_labels{end+1, 1} = lh_label;
    combined_ids{end+1, 1} = lh_id;
    combined_coords{end+1, 1} = rh_coords;
    combined_labels{end+1, 1} = rh_label;
    combined_ids{end+1, 1} = rh_id;
end

% Convert cell arrays to table
coords_table = cell2mat(combined_coords);
labels_table = string(combined_labels);
ids_table = cell2mat(combined_ids);

% Add two columns of zeros
zero_column = zeros(size(coords_table, 1), 2);

% Combine all data into the final table
final_table = table(coords_table(:,1), coords_table(:,2), coords_table(:,3), ...
    zero_column(:,1), zero_column(:,2), labels_table, ids_table, ...
    'VariableNames', {'x', 'y', 'z', 'zero1', 'zero2', 'label', 'id'});

% Save the final table as an ASCII text file with tab delimiter
writetable(final_table, 'brainnetome_combined_with_ids.node', 'FileType', 'text', 'Delimiter', '\t');

disp('Node file created successfully.');

%%
% Perform Bayesian comparisons for each measure
measures = fieldnames(residuals_by_group);
measure_ids = struct('Dloc', 1, 'BCloc', 2, 'ECloc', 3, 'Cloc', 4, 'Eloc', 5, 'Ploc', 6);

% Define comparisons
comparisons = {
    'Controls vs Nicotine'; 
    'Controls vs Opiates'; 
    'Controls vs Cannabis'; 
    'Nicotine vs Opiates'; 
    'Nicotine vs Cannabis'; 
    'Opiates vs Cannabis'
};

% Define the log10 threshold for significance
log_bf_threshold = 1;  % This corresponds to a Bayes factor threshold of 10

all_sig_regions = [];

for i = 1:numel(measures)
    measure = measures{i};
    measure_id = measure_ids.(measure);
    posthoc_results = [];
    data = residuals_by_group.(measure);
    
    % Check the number of columns in each group and find the maximum
    num_cols = max([size(data.Controls, 2), size(data.Nicotine, 2), size(data.Opiates, 2), size(data.Cannabis, 2)]);
    
    % Pad the data with NaNs to ensure consistent dimensions
    data.Controls(:, end+1:num_cols) = NaN;
    data.Nicotine(:, end+1:num_cols) = NaN;
    data.Opiates(:, end+1:num_cols) = NaN;
    data.Cannabis(:, end+1:num_cols) = NaN;
    
    % Loop through each region
    for region = 1:size(data.Controls, 1)
        % Extract data for the current region, ignoring NaNs
        data_controls = data.Controls(region, :);
        data_nicotine = data.Nicotine(region, :);
        data_opiates = data.Opiates(region, :);
        data_cannabis = data.Cannabis(region, :);
        
        % Remove NaNs
        data_controls = data_controls(~isnan(data_controls));
        data_nicotine = data_nicotine(~isnan(data_nicotine));
        data_opiates = data_opiates(~isnan(data_opiates));
        data_cannabis = data_cannabis(~isnan(data_cannabis));
        
        % Perform pairwise Bayesian t-tests
        bf_CvsNic = bf.ttest2(data_controls, data_nicotine);
        bf_CvsOpi = bf.ttest2(data_controls, data_opiates);
        bf_CvsCan = bf.ttest2(data_controls, data_cannabis);
        bf_NicvsOpi = bf.ttest2(data_nicotine, data_opiates);
        bf_NicvsCan = bf.ttest2(data_nicotine, data_cannabis);
        bf_OpivsCan = bf.ttest2(data_opiates, data_cannabis);
        
        % Log10 transform the Bayes factors
        log_bf_values = log10([bf_CvsNic; bf_CvsOpi; bf_CvsCan; bf_NicvsOpi; bf_NicvsCan; bf_OpivsCan]);
        
        % Store results in the table
        region_results = table(repmat(region, 6, 1), comparisons, log_bf_values, ...
                               'VariableNames', {'Region', 'Comparison', 'log_bayes_factor'});
        posthoc_results = [posthoc_results; region_results];
    end
    
    % Find significant regions
    significant_regions.(measure) = find_significant_regions_bayes(posthoc_results, log_bf_threshold);
    
    % Add significant regions to the all_sig_regions array
    sig_regions = significant_regions.(measure);
    for j = 1:size(sig_regions, 1)
        region = sig_regions{j, 1};
        comparison = sig_regions{j, 2};
        log_bf = sig_regions{j, 3};
        id = final_table.id(region);
        coords = final_table{region, {'x', 'y', 'z'}};
        label = final_table.label(region);
        all_sig_regions = [all_sig_regions; {coords(1), coords(2), coords(3), measure_id, log_bf, label, comparison}];
    end
end

% Create a table for the significant regions
sig_regions_table = cell2table(all_sig_regions, 'VariableNames', {'x', 'y', 'z', 'measure', 'log_bayes_factor', 'label', 'comparison'});

% Save the table as an ASCII text file with tab delimiter and .node suffix
writetable(sig_regions_table, 'significant_regions_combined.node', 'FileType', 'text', 'Delimiter', '\t');

disp('Significant regions with IDs saved to node file successfully.');

%%
function sig_regions = find_significant_regions_bayes(posthoc_results, log_bf_threshold)
    unique_regions = unique(posthoc_results.Region);
    sig_regions = {};
    for region = unique_regions'
        log_bf_values = posthoc_results.log_bayes_factor(posthoc_results.Region == region);
        comparisons = posthoc_results.Comparison(posthoc_results.Region == region);
        for i = 1:length(log_bf_values)
            if log_bf_values(i) > log_bf_threshold
                sig_regions = [sig_regions; {region, comparisons{i}, log_bf_values(i)}];
            end
        end
    end
end

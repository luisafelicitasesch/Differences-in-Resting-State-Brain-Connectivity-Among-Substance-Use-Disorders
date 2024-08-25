%%
path_bayes = 'C:\Users\luisa\MATLAB\packages\bayesFactor-master';
addpath(genpath(path_bayes));

% Potentially exclude GraphVar from path

%% Bayesian Analysis of Group Comparisons Across Brain Regions: Identifying Significant Differences
% Define comparisons
comparisons = {
    'Controls vs Nicotine'; 
    'Controls vs Opiates'; 
    'Controls vs Cannabis'; 
    'Nicotine vs Opiates'; 
    'Nicotine vs Cannabis'; 
    'Opiates vs Cannabis'
};

% Initialize structures to store significant regions for each measure
significant_regions = struct();
significant_regions.Dloc = {};
significant_regions.BCloc = {};
significant_regions.ECloc = {};
significant_regions.Cloc = {};
significant_regions.Eloc = {};

% Structure to hold residuals by group
residuals_by_group = struct();
residuals_by_group.Dloc = Dloc_residuals_by_group;
residuals_by_group.BCloc = BCloc_residuals_by_group;
residuals_by_group.ECloc = ECloc_residuals_by_group;
residuals_by_group.Cloc = Cloc_residuals_by_group;
residuals_by_group.Eloc = Eloc_residuals_by_group;

% Perform Bayesian comparisons for each measure
measures = fieldnames(residuals_by_group);

% Define the log10 threshold for significance
log_bf_threshold = 1;  % This corresponds to a Bayes factor threshold of 10

for i = 1:numel(measures)
    measure = measures{i};
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
    
    % Pivot the results to have a column for each comparison
    posthoc_pivot = unstack(posthoc_results, 'log_bayes_factor', 'Comparison');
    
    % Modify variable names to be valid MATLAB identifiers
    valid_names = matlab.lang.makeValidName(posthoc_pivot.Properties.VariableNames);
    posthoc_pivot.Properties.VariableNames = valid_names;
    
    % Save the post hoc test results to an XLSX file
    writetable(posthoc_pivot, [measure '_bayes_posthoc_results.xlsx']);
    disp(['Post hoc test results for ' measure ' saved to XLSX successfully.']);
end

% Display significant regions
disp('Significant regions based on log-transformed Bayes Factors:');
disp(significant_regions);

% Save the significant regions to a text file
output_file = 'significant_regions_bayes.txt';
fid = fopen(output_file, 'w');
fields = fieldnames(significant_regions);
for i = 1:numel(fields)
    fprintf(fid, '%s:\n', fields{i});
    sig_regions = significant_regions.(fields{i});
    for j = 1:size(sig_regions, 1)
        fprintf(fid, 'Region: %d, Comparison: %s, log(Bayes Factor): %.2f\n', sig_regions{j, 1}, sig_regions{j, 2}, sig_regions{j, 3});
    end
    fprintf(fid, '\n');
end
fclose(fid);
disp(['Significant regions saved to ' output_file ' successfully.']);

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

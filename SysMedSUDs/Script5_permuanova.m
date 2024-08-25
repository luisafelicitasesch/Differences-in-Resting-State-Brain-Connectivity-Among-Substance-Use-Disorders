% Setup paths
path_BCT = 'C:\Users\luisa\MATLAB\packages\BCT';
path_GraphVar = 'C:\Users\luisa\MATLAB\packages\GraphVar_2.03a';
path_ComAlg = 'C:\Users\luisa\MATLAB\packages\communityalg-master';
addpath(genpath(path_BCT));
addpath(genpath(path_GraphVar));
addpath(genpath(path_ComAlg));

% Define unique groups from participantGroups
uniqueGroups = {'Controls', 'Cannabis', 'Opiates', 'Nicotine'};

% Define global measures and their corresponding data variables (residuals)
global_measures = {'Tglob', 'Eglob', 'Pglob', 'Cglob'};
residuals_variables = {'residuals_Tglob_by_group', 'residuals_Eglob_by_group', 'residuals_Pglob_by_group', 'residuals_Cglob_by_group'};

% Initialize structures to store results
f_values = zeros(length(global_measures), 1);
p_values = zeros(length(global_measures), 1);
ci_values = zeros(length(global_measures), 1);
group_means = cell(length(global_measures), 1);
results_table = table();

for i = 1:length(global_measures)
    measure = global_measures{i};
    residuals_var = residuals_variables{i};
    
    % Get the residual data for the current measure
    residuals_data = eval(residuals_var);

    % Concatenate all residuals data and create a grouping vector
    residuals_all = [residuals_data.Controls; residuals_data.Cannabis; residuals_data.Opiates; residuals_data.Nicotine];

    % Verify the sizes of each group
    size_controls = size(residuals_data.Controls, 1);
    size_cannabis = size(residuals_data.Cannabis, 1);
    size_opiates = size(residuals_data.Opiates, 1);
    size_nicotine = size(residuals_data.Nicotine, 1);

    % Create the group vector with consistent dimensions
    group = [repelem(1, size_controls)'; repelem(2, size_cannabis)'; ...
             repelem(3, size_opiates)'; repelem(4, size_nicotine)'];

    % Perform permutation-based ANOVA
    [F, p, ci, stats, tbl, dist] = permuanova1(residuals_all, group, 'nperm', 10000);

    % Store results
    f_values(i) = F;
    p_values(i) = p;
    ci_values(i) = ci(1);

    % Store group means from the residuals data
    groupMeans = [mean(residuals_data.Controls), mean(residuals_data.Cannabis), mean(residuals_data.Opiates), mean(residuals_data.Nicotine)];
    group_means{i} = groupMeans;

    % Display results
    disp(['Results for global ' measure ':']);
    disp(['F-value: ' num2str(F)]);
    disp(['p-value: ' num2str(p)]);
    disp(['95% Confidence Interval for F: [' num2str(ci(1)) ', Inf]']);
    disp(['Group Means: ' num2str(groupMeans)]);
    disp(tbl); % Display ANOVA table
    
    % Append to results table
    results_table = [results_table; table({measure}, F, p, ci(1), group_means{i}(1), group_means{i}(2), group_means{i}(3), group_means{i}(4), ...
                    'VariableNames', {'Measure', 'F_value', 'p_value', 'CI_lower', 'Controls_mean', 'Cannabis_mean', 'Opiates_mean', 'Nicotine_mean'})];
end

% FDR correction using mafdr
adj_p_values = mafdr(p_values, 'BHFDR', true);

% Add FDR-corrected p-values to results table
results_table.FDR_p_value = adj_p_values;

% Display the results table with FDR-corrected p-values
disp(results_table);

% Save the results table to a CSV file
writetable(results_table, 'global_measures_results_with_fdr.csv');
disp('Results with FDR correction saved to CSV successfully.');

% Save the results to a .mat file
save('global_measures_permutation_anova_results.mat', 'results_table');
disp('Permutation-based ANOVA results for global measures saved successfully.');

%% Local Measures

% Define local measures and their corresponding data variables (residuals)
local_measures = {'Dloc', 'BCloc', 'ECloc', 'Cloc', 'Eloc'};
residuals_local_variables = {'Dloc_residuals_by_group', 'BCloc_residuals_by_group', 'ECloc_residuals_by_group', 'Cloc_residuals_by_group', 'Eloc_residuals_by_group'};

for k = 1:length(local_measures)
    measure = local_measures{k};
    residuals_var = residuals_local_variables{k};

    % Get the residual data for the current local measure
    residuals_data = eval(residuals_var);
    
    % Extract the data for each group
    Controls = residuals_data.Controls; 
    Cannabis = residuals_data.Cannabis; 
    Opiates = residuals_data.Opiates; 
    Nicotine = residuals_data.Nicotine; 

    % Verify the number of regions in each group
    numRegions = size(Controls, 1); % Regions are rows

    % Initialize vectors to store p-values for each region
    resi_p_values = zeros(numRegions, 1);
    resi_means = zeros(numRegions, 4);  % For storing group means

    % Loop through each region
    for region = 1:numRegions
        % Ensure that the region index does not exceed the number of rows for any group
        if region <= size(Cannabis, 1) && region <= size(Opiates, 1) && region <= size(Nicotine, 1) && region <= size(Controls, 1)
            % Extract data for the current region
            data = [
                Controls(region, :), ...
                Cannabis(region, :), ...
                Opiates(region, :), ...
                Nicotine(region, :)
            ];

            % Create the group vector
            group = [
                repmat(1, size(Controls, 2), 1);
                repmat(2, size(Cannabis, 2), 1);
                repmat(3, size(Opiates, 2), 1);
                repmat(4, size(Nicotine, 2), 1)
            ];

            % Perform permutation-based ANOVA for the current region
            [F, p, ci, stats, tbl, dist] = permuanova1(data(:), group, 'nperm', 10000);

            % Store the p-value
            resi_p_values(region) = p;

            % Store group means
            resi_means(region, :) = stats.means;
        else
            % If the region index exceeds the number of rows for any group, skip this region
            resi_p_values(region) = NaN;
            resi_means(region, :) = NaN;
        end
    end

    % FDR correction
    resi_adj_p_values = mafdr(resi_p_values, 'BHFDR', true);

    % Store results in a table
    resi_results = table((1:numRegions)', resi_p_values, resi_adj_p_values, ...
                         resi_means(:, 1), resi_means(:, 2), resi_means(:, 3), resi_means(:, 4), ...
                         'VariableNames', {'Region', 'p_value', 'adj_p_value', 'Controls_mean', 'Cannabis_mean', 'Opiates_mean', 'Nicotine_mean'});

    % Save the results to a CSV file
    writetable(resi_results, [measure '_permutation_anova_results.csv']);
    disp(['Permutation-based ANOVA results for ' measure ' saved to CSV successfully.']);

    writetable(resi_results, [measure '_permutation_anova_results.xlsx']);
    disp(['Permutation-based ANOVA results for ' measure ' saved to xlsx successfully.']);

    % Save the local results table as a .mat file
    save([measure '_permutation_anova_results.mat'], 'resi_results');
    disp(['Permutation-based ANOVA results for ' measure ' saved successfully.']);
end

%% Collect significant findings into tables
% Collect significant findings into tables

% Initialize a table for significant findings
significant_findings = table();

for k = 1:length(local_measures)
    measure = local_measures{k};
    results_var_name = ['resi_' measure '_results'];
    
    % Retrieve the results variable for the current measure
    if exist(results_var_name, 'var')
        results = eval(results_var_name);

        % Find significant regions before FDR correction
        significant_regions_p = find(results.p_value < 0.05);
        disp(['Significant regions before FDR correction (' measure '):']);
        disp(significant_regions_p);

        % Find significant regions after FDR correction
        significant_regions_adj_p = find(results.adj_p_value < 0.05);
        disp(['Significant regions after FDR correction (' measure '):']);
        disp(significant_regions_adj_p);

        % Store in the table
        significant_findings = [significant_findings; table({measure}, {significant_regions_p}, {significant_regions_adj_p}, ...
            'VariableNames', {'Measure', 'Before_FDR', 'After_FDR'})];
    else
        disp(['No results found for measure: ' measure]);
    end
end

% Save the significant findings table as an .xlsx file
writetable(significant_findings, 'significant_findings.xlsx', 'WriteVariableNames', true);
disp('Significant findings saved to Excel file successfully.');

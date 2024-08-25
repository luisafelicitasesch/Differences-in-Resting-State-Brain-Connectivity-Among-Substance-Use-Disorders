% Setup paths
path_BCT = 'C:\Users\luisa\MATLAB\packages\BCT';
path_GraphVar = 'C:\Users\luisa\MATLAB\packages\GraphVar_2.03a';
path_ComAlg = 'C:\Users\luisa\MATLAB\packages\communityalg-master';
addpath(genpath(path_BCT));
addpath(genpath(path_GraphVar));
addpath(genpath(path_ComAlg));

%% Define unique groups
uniqueGroups = {'Controls', 'Nicotine', 'Cocaine'};

% Define global measures and their corresponding data variables in the desired order
global_measures = {'Cglob', 'Tglob', 'Eglob', 'Pglob'};
data_variables = {'residuals_Cglob_by_group', 'residuals_Tglob_by_group', 'residuals_Eglob_by_group', 'residuals_Pglob_by_group'};
auc_data_variables = {'Cglob_AUC', 'Tglob_AUC', 'Eglob_AUC', 'Pglob_AUC'};

% Initialize structures to store results
f_values = zeros(length(global_measures), 1);
p_values = zeros(length(global_measures), 1);
ci_values = zeros(length(global_measures), 1);
group_means_residuals = cell(length(global_measures), 1);
group_means_auc = cell(length(global_measures), 1);
results_table = table();

for i = 1:length(global_measures)
    measure = global_measures{i};
    data_var = data_variables{i};
    auc_data_var = auc_data_variables{i};
    
    % Get the residual data for the current measure
    residualData = eval(data_var);
    
    % Get the AUC data for the current measure
    aucAllData = eval(auc_data_var);

    % Concatenate all residual data and create a grouping vector
    residual_all = [residualData.Controls; residualData.Nicotine; residualData.Cocaine];

    % Verify the sizes of each group
    size_controls = size(residualData.Controls, 1);
    size_nicotine = size(residualData.Nicotine, 1);
    size_cocaine = size(residualData.Cocaine, 1);

    % Create the group vector with consistent dimensions
    group = [repelem(1, size_controls)'; repelem(2, size_nicotine)'; ...
             repelem(3, size_cocaine)'];

    % Perform permutation-based ANOVA on residual data
    [F, p, ci, stats, tbl, dist] = permuanova1(residual_all, group, 'nperm', 10000);

    % Store results
    f_values(i) = F;
    p_values(i) = p;
    ci_values(i) = ci(1);
    
    % Calculate and store group means from the residuals
    groupMeans_residuals = [mean(residualData.Controls, 'omitnan'), ...
                            mean(residualData.Nicotine, 'omitnan'), ...
                            mean(residualData.Cocaine, 'omitnan')];
    group_means_residuals{i} = groupMeans_residuals;

    % Calculate and store group means from the original AUC data
    groupMeans_auc = [mean(aucAllData.Controls, 'omitnan'), ...
                      mean(aucAllData.Nicotine, 'omitnan'), ...
                      mean(aucAllData.Cocaine, 'omitnan')];
    group_means_auc{i} = groupMeans_auc;

    % Display results
    disp(['Results for global ' measure ':']);
    disp(['F-value: ' num2str(F)]);
    disp(['p-value: ' num2str(p)]);
    disp(['95% Confidence Interval for F: [' num2str(ci(1)) ', Inf]']);
    disp(['Group Means (Residuals): ' num2str(groupMeans_residuals)]);
    disp(['Group Means (AUC): ' num2str(groupMeans_auc)]);
    disp(tbl); % Display ANOVA table
    
    % Append to results table
    results_table = [results_table; table({measure}, F, p, ci(1), ...
                    group_means_auc{i}(1), group_means_auc{i}(2), group_means_auc{i}(3), ...
                    group_means_residuals{i}(1), group_means_residuals{i}(2), group_means_residuals{i}(3), ...
                    'VariableNames', {'Measure', 'F_value', 'p_value', 'CI_lower', ...
                                      'Controls_mean_auc', 'Nicotine_mean_auc', 'Cocaine_mean_auc', ...
                                      'Controls_mean_residual', 'Nicotine_mean_residual', 'Cocaine_mean_residual'})];
end

% FDR correction using mafdr
adj_p_values = mafdr(p_values, 'BHFDR', true);

% Add FDR-corrected p-values to results table
results_table.FDR_p_value = adj_p_values;

% Display the results table with FDR-corrected p-values
disp(results_table);

% Save the results table to a CSV file
writetable(results_table, 'global_measures_results_with_fdr.xlsx');
disp('Results with FDR correction saved to CSV successfully.');

% Store individual results for each measure with FDR correction
Cglob_residual_results = results_table(strcmp(results_table.Measure, 'Cglob'), :);
Tglob_residual_results = results_table(strcmp(results_table.Measure, 'Tglob'), :);
Eglob_residual_results = results_table(strcmp(results_table.Measure, 'Eglob'), :);
Pglob_residual_results = results_table(strcmp(results_table.Measure, 'Pglob'), :);

% Save the results to a .mat file
save('global_measures_permutation_anova_results.mat', 'results_table', 'Cglob_residual_results', 'Tglob_residual_results', 'Eglob_residual_results', 'Pglob_residual_results');
disp('Permutation-based ANOVA results for global measures saved successfully.');


%% Initialize an empty table for all results
all_local_results_table = table('Size', [0 10], ...
    'VariableTypes', {'double', 'cell', 'double', 'double', 'double', 'double', 'double', 'double', 'double', 'double'}, ...
    'VariableNames', {'Region', 'Measure', 'p_value', 'adj_p_value', 'Controls_mean', 'Nicotine_mean', 'Cocaine_mean', 'Controls_mean_auc', 'Nicotine_mean_auc', 'Cocaine_mean_auc'});

%% Local Measures
local_measures = {'Dloc', 'Eloc', 'ECloc', 'Cloc', 'BCloc'};
local_data_variables = {'Dloc_residuals_by_group', 'Eloc_residuals_by_group', 'ECloc_residuals_by_group', 'Cloc_residuals_by_group', 'BCloc_residuals_by_group'};
auc_data_variables = {'Dloc_AUC', 'Eloc_AUC', 'ECloc_AUC', 'Cloc_AUC', 'BCloc_AUC'};

for j = 1:length(local_measures)
    measure = local_measures{j};
    data_var = local_data_variables{j};
    auc_var = auc_data_variables{j};
    
    % Extract the data for each group (residuals)
    Controls = eval([data_var '.Controls']); % 246xN matrix
    Nicotine = eval([data_var '.Nicotine']); % 246xN matrix
    Cocaine = eval([data_var '.Cocaine']); % 246xN matrix
    
    % Extract the data for each group (AUCs)
    Controls_auc = eval([auc_var '.Controls']); % N (participants) x 246 (regions)
    Nicotine_auc = eval([auc_var '.Nicotine']); % N (participants) x 246 (regions)
    Cocaine_auc = eval([auc_var '.Cocaine']); % N (participants) x 246 (regions)

    % Verify the number of regions in each group
    numRegions = size(Controls, 1); % Regions are rows

    % Initialize vectors to store p-values and means for each region
    p_values = zeros(numRegions, 1);
    means = zeros(numRegions, 3);  % For storing group means (residuals)
    means_auc = zeros(numRegions, 3); % For storing group means (AUCs)

    % Loop through each region
    for region = 1:numRegions
        % Extract data for the current region (residuals)
        data = [
            Controls(region, :), ...
            Nicotine(region, :), ...
            Cocaine(region, :)
        ];

        % Extract data for the current region (AUCs)
        data_auc = [
            Controls_auc(:, region)', ...
            Nicotine_auc(:, region)', ...
            Cocaine_auc(:, region)'
        ];

        % Create the group vector
        group = [
            repmat(1, size(Controls, 2), 1);
            repmat(2, size(Nicotine, 2), 1);
            repmat(3, size(Cocaine, 2), 1)
        ];

        % Perform permutation-based ANOVA for the current region
        [F, p, ci, stats, tbl, dist] = permuanova1(data(:), group, 'nperm', 10000);

        % Store the p-value
        p_values(region) = p;

        % Store group means (residuals)
        means(region, :) = stats.means;

        % Store group means (AUCs)
        means_auc(region, :) = [
            mean(Controls_auc(:, region), 'omitnan'), ...
            mean(Nicotine_auc(:, region), 'omitnan'), ...
            mean(Cocaine_auc(:, region), 'omitnan')
        ];
    end

    % FDR correction
    adj_p_values = mafdr(p_values, 'BHFDR', true);

    % Store results in a table
    results = table((1:numRegions)', repmat({measure}, numRegions, 1), p_values, adj_p_values, ...
                         means(:, 1), means(:, 2), means(:, 3), ...
                         means_auc(:, 1), means_auc(:, 2), means_auc(:, 3), ...
                         'VariableNames', {'Region', 'Measure', 'p_value', 'adj_p_value', 'Controls_mean', 'Nicotine_mean', 'Cocaine_mean', 'Controls_mean_auc', 'Nicotine_mean_auc', 'Cocaine_mean_auc'});

    % Append local measure results to all_local_results_table
    all_local_results_table = [all_local_results_table; results];

    % Save the results to a CSV file
    writetable(results, [measure '_permutation_anova_results.xlsx']);
    disp(['Permutation-based ANOVA results for ' measure ' saved to CSV successfully.']);

    % Save the results to a .mat file
    save([measure '_permutation_anova_results.mat'], 'results');
    disp(['Permutation-based ANOVA results for ' measure ' saved successfully.']);
end

% Save all local measures results to a single xlsx file
writetable(all_local_results_table, 'all_local_measures_permutation_anova_results.xlsx');
disp('All local measures permutation ANOVA results saved to Excel file successfully.');

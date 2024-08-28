% Define global measures and corresponding data variables
global_measures = {'Cglob', 'Tglob', 'Eglob', 'Pglob'};
data_variables = {'residuals_Cglob_by_group', 'residuals_Tglob_by_group', 'residuals_Eglob_by_group', 'residuals_Pglob_by_group'};

% Define comparisons
comparisons = {'Controls vs Nicotine', 'Controls vs Opiates', 'Controls vs Cannabis', 'Nicotine vs Opiates', 'Nicotine vs Cannabis', 'Opiates vs Cannabis'};

% Initialize a table to store all results
VariableNames = [{'Measure'}, comparisons]; % Create a cell array for VariableNames
bayes_results_table = table('Size', [numel(global_measures), numel(comparisons) + 1], ...
                            'VariableTypes', [{'string'}, repmat({'double'}, 1, numel(comparisons))], ...
                            'VariableNames', VariableNames);

for i = 1:numel(global_measures)
    measure = global_measures{i};
    data_var = data_variables{i};
    data = eval(data_var);  % Load the residuals data for the current measure
    
    % Extract data for the four groups
    data_controls = data.Controls(~isnan(data.Controls));
    data_nicotine = data.Nicotine(~isnan(data.Nicotine));
    data_opiates = data.Opiates(~isnan(data_opiates));
    data_cannabis = data.Cannabis(~isnan(data.Cannabis));
    
    % Perform pairwise Bayesian t-tests
    bf_CvsNic = bf.ttest2(data_controls, data_nicotine);
    bf_CvsOpi = bf.ttest2(data_controls, data_opiates);
    bf_CvsCan = bf.ttest2(data_controls, data_cannabis);
    bf_NicvsOpi = bf.ttest2(data_nicotine, data_opiates);
    bf_NicvsCan = bf.ttest2(data_nicotine, data_cannabis);
    bf_OpivsCan = bf.ttest2(data_opiates, data_cannabis);
    
    % Log10 transform the Bayes factors
    log_bf_values = log10([bf_CvsNic; bf_CvsOpi; bf_CvsCan; bf_NicvsOpi; bf_NicvsCan; bf_OpivsCan]);
    
    % Store the results in the results_table
    bayes_results_table.Measure{i} = measure;
    bayes_results_table{i, 2:end} = log_bf_values';
end

% Name of the Excel file to save all results
output_xlsx_file = 'global_measures_bayes_factor_results.xlsx';

% Write the results table to the Excel file
writetable(bayes_results_table, output_xlsx_file, 'Sheet', 'Bayesian_Results');
disp(['Bayesian comparison results saved to ' output_xlsx_file ' successfully.']);

% Display the results
disp('Log-transformed Bayes Factors for each comparison:');
disp(bayes_results_table);

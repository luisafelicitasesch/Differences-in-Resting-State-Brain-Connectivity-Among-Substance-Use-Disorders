% Setup paths
path_BCT = 'C:\Users\luisa\MATLAB\packages\BCT';
path_ComAlg = 'C:\Users\luisa\MATLAB\packages\communityalg-master';
addpath(genpath(path_BCT));
addpath(genpath(path_ComAlg));
path_bayes = 'C:\Users\luisa\MATLAB\packages\bayesFactor-master';
addpath(genpath(path_bayes))

% Define global measures and corresponding data variables
global_measures = {'Cglob', 'Tglob', 'Eglob', 'Pglob'};
data_variables = {'residuals_Cglob_by_group', 'residuals_Tglob_by_group', 'residuals_Eglob_by_group', 'residuals_Pglob_by_group'};

% Define comparisons
comparisons = {'Controls vs Nicotine', 'Controls vs Cocaine', 'Nicotine vs Cocaine'};

% Initialize a table to store all results
results_table = table('Size', [numel(global_measures), numel(comparisons) + 1], ...
                      'VariableTypes', [{'string'}, repmat({'double'}, 1, numel(comparisons))], ...
                      'VariableNames', ['Measure', comparisons]);

% Define the log10 threshold for significance
log_bf_threshold = 1;  % This corresponds to a Bayes factor threshold of 10

for i = 1:numel(global_measures)
    measure = global_measures{i};
    data_var = data_variables{i};
    data = eval(data_var);  % Load the residuals data for the current measure
    
    % Extract data for the three groups
    data_controls = data.Controls(~isnan(data.Controls));
    data_nicotine = data.Nicotine(~isnan(data.Nicotine));
    data_cocaine = data.Cocaine(~isnan(data.Cocaine));
    
    % Perform pairwise Bayesian t-tests
    bf_CvsNic = bf.ttest2(data_controls, data_nicotine);
    bf_CvsCoc = bf.ttest2(data_controls, data_cocaine);
    bf_NicvsCoc = bf.ttest2(data_nicotine, data_cocaine);
    
    % Log10 transform the Bayes factors
    log_bf_values = log10([bf_CvsNic; bf_CvsCoc; bf_NicvsCoc]);
    
    % Store the results in the results_table
    results_table.Measure{i} = measure;
    results_table{i, 2:end} = log_bf_values';
end

% Name of the Excel file to save all results
output_xlsx_file = 'global_measures_bayes_factor_results.xlsx';

% Write the results table to the Excel file
writetable(results_table, output_xlsx_file, 'Sheet', 'BayesFactors');
disp(['All results saved to ' output_xlsx_file ' successfully.']);

% Display the results
disp('Log-transformed Bayes Factors for each comparison:');
disp(results_table);

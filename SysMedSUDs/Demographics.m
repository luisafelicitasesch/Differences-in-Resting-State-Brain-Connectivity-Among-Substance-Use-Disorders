% Extract the IDs from participantNames
ids_to_keep = participantNames;

% Assuming demodata contains a column 'ID' with participant IDs
% If the column name is different, adjust the variable name accordingly
demodata_ids = demodata.id;

% Find the indices of rows in demodata where the IDs are in ids_to_keep
rows_to_keep = ismember(demodata_ids, ids_to_keep);

% Filter the rows of demodata
filtered_demodata = demodata(rows_to_keep, :);

%%
% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

% Initialize tables to store the statistics
age_stats = table();
sex_stats = table();
education_stats = table();
macs_stats = table();
group_counts = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
    % Count the number of participants
    num_participants = sum(group_idx);
    
    % Remove rows with NaN values for each variable
    age_data = group_data.Alter;
    age_data = age_data(~isnan(age_data));
    
    sex_data = group_data.Geschlecht;
    sex_data = sex_data(~isnan(sex_data));
    
    education_data = group_data.schule; % Assuming 'schule' is the column name for education
    education_data = education_data(~isnan(education_data));
    
    macs_data = group_data.MACS_sum;
    macs_data = macs_data(~isnan(macs_data));
    
    % Calculate Age statistics
    age_mean = mean(age_data);
    age_median = median(age_data);
    age_sd = std(age_data);
    age_min = min(age_data);
    age_max = max(age_data);
    age_row = table({group}, num_participants, age_mean, age_median, age_sd, age_min, age_max, ...
                    'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    age_stats = [age_stats; age_row];
    
    % Calculate Sex statistics
    sex_male = sum(sex_data == 1);
    sex_female = sum(sex_data == 2);
    sex_row = table({group}, num_participants, sex_male, sex_female, ...
                    'VariableNames', {'Group', 'N', 'MaleCount', 'FemaleCount'});
    sex_stats = [sex_stats; sex_row];
    
    % Calculate Education statistics
    education_mean = mean(education_data);
    education_median = median(education_data);
    education_sd = std(education_data);
    education_min = min(education_data);
    education_max = max(education_data);
    education_row = table({group}, num_participants, education_mean, education_median, education_sd, education_min, education_max, ...
                          'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    education_stats = [education_stats; education_row];
    
    % Calculate MACS statistics
    macs_mean = mean(macs_data);
    macs_median = median(macs_data);
    macs_sd = std(macs_data);
    macs_min = min(macs_data);
    macs_max = max(macs_data);
    macs_row = table({group}, num_participants, macs_mean, macs_median, macs_sd, macs_min, macs_max, ...
                     'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    macs_stats = [macs_stats; macs_row];
end

% Calculate overall statistics for all participants
% Remove rows with NaN values for each variable
age_data_all = filtered_demodata.Alter;
age_data_all = age_data_all(~isnan(age_data_all));

sex_data_all = filtered_demodata.Geschlecht;
sex_data_all = sex_data_all(~isnan(sex_data_all));

education_data_all = filtered_demodata.schule;
education_data_all = education_data_all(~isnan(education_data_all));

macs_data_all = filtered_demodata.MACS_sum;
macs_data_all = macs_data_all(~isnan(macs_data_all));

% Age statistics
age_mean_all = mean(age_data_all);
age_median_all = median(age_data_all);
age_sd_all = std(age_data_all);
age_min_all = min(age_data_all);
age_max_all = max(age_data_all);
age_row_all = table({'All'}, length(age_data_all), age_mean_all, age_median_all, age_sd_all, age_min_all, age_max_all, ...
                    'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
age_stats = [age_stats; age_row_all];

% Sex statistics
sex_male_all = sum(sex_data_all == 1);
sex_female_all = sum(sex_data_all == 2);
sex_row_all = table({'All'}, length(sex_data_all), sex_male_all, sex_female_all, ...
                    'VariableNames', {'Group', 'N', 'MaleCount', 'FemaleCount'});
sex_stats = [sex_stats; sex_row_all];

% Education statistics
education_mean_all = mean(education_data_all);
education_median_all = median(education_data_all);
education_sd_all = std(education_data_all);
education_min_all = min(education_data_all);
education_max_all = max(education_data_all);
education_row_all = table({'All'}, length(education_data_all), education_mean_all, education_median_all, education_sd_all, education_min_all, education_max_all, ...
                          'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
education_stats = [education_stats; education_row_all];

% MACS statistics
macs_mean_all = mean(macs_data_all);
macs_median_all = median(macs_data_all);
macs_sd_all = std(macs_data_all);
macs_min_all = min(macs_data_all);
macs_max_all = max(macs_data_all);
macs_row_all = table({'All'}, length(macs_data_all), macs_mean_all, macs_median_all, macs_sd_all, macs_min_all, macs_max_all, ...
                     'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
macs_stats = [macs_stats; macs_row_all];

% Display the tables
disp('Age Statistics:');
disp(age_stats);

disp('Sex Statistics:');
disp(sex_stats);

disp('Education Statistics:');
disp(education_stats);

disp('MACS Statistics:');
disp(macs_stats);

%% Group comparisons

% Age comparison using ANOVA
[p_age, tbl_age, stats_age] = anova1(filtered_demodata.Alter, filtered_demodata.Gruppe);
title('ANOVA for Age by Group');
fprintf('ANOVA for Age by Group:\n');
disp(tbl_age);
fprintf('p-value = %.20f\n\n', p_age);

% Post-hoc test for age if ANOVA is significant
if p_age < 0.05
    figure;
    multcompare(stats_age, 'CType', 'bonferroni');
    title('Post-hoc Test for Age');
end

% Education comparison using ANOVA
[p_education, tbl_education, stats_education] = anova1(filtered_demodata.schule, filtered_demodata.Gruppe);
title('ANOVA for Education by Group');
fprintf('ANOVA for Education by Group:\n');
disp(tbl_education);
fprintf('p-value = %.20f\n\n', p_education);

% Post-hoc test for education if ANOVA is significant
if p_education < 0.05
    figure;
    multcompare(stats_education, 'CType', 'bonferroni');
    title('Post-hoc Test for Education');
end

% MACS comparison using ANOVA
[p_macs, tbl_macs, stats_macs] = anova1(filtered_demodata.MACS_sum, filtered_demodata.Gruppe);
title('ANOVA for MACS_sum by Group');
fprintf('ANOVA for MACS_sum by Group:\n');
disp(tbl_macs);
fprintf('p-value = %.20f\n\n', p_macs);

% Post-hoc test for MACS if ANOVA is significant
if p_macs < 0.05
    figure;
    multcompare(stats_macs, 'CType', 'bonferroni');
    title('Post-hoc Test for MACS_sum');
end

% Sex comparison using Chi-squared test
sex_data = filtered_demodata.Geschlecht;
sex_groups = filtered_demodata.Gruppe;
[crosstab, chi2, p_sex, labels] = crosstab(sex_data, sex_groups);

% Display results
fprintf('Chi-squared test for Sex by Group:\n');
fprintf('Chi2 = %.2f\n', chi2);

%%
% Assuming demodata and participantGroups are already defined
% demodata.Alter contains the age data
% participantGroups contains the group IDs as a cell array of strings
makeFigurePretty
violinplot(age_reordered, participantGroups)


%%
% Load the data

% Define the group labels and corresponding values in 'demodata.Gruppe'
group_labels = {'Control', 'Nicotine', 'Cannabis', 'Opiates'};
group_values = [4, 3, 2, 1]; % Adjust these values if needed to match the actual group IDs in demodata.Gruppe

% List of _diagnose variables to analyze
diagnose_variables = {
    'Nikton', 'Cannabis', 'Alkohol', 'Opioide', ...
    'Amph', 'XTC', 'Kokain', 'Benzo', 'Prega'
};

% Preallocate data for the bar charts
none_counts = zeros(length(group_labels), length(diagnose_variables));
present_counts = zeros(length(group_labels), length(diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(diagnose_variables));

for var = 1:length(diagnose_variables)
    variable = [diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group = group_labels{i};
        group_value = group_values(i);
        
        % Extract the rows for the current group
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);
        
        % Count occurrences of each diagnosis status
        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Remove the '_diagnose' suffix from diagnose variable names for labeling
diagnose_labels = strrep(diagnose_variables, '_diagnose', '');

% Plot grouped bar charts for each group
figure;
for i = 1:length(group_labels)
    group = group_labels{i};
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar(bar_data', 'grouped');
    title(group, 'Interpreter', 'none');
    xlabel('Diagnosis Variables');
    ylabel('Count');
    set(gca, 'XTickLabel', diagnose_labels, 'XTickLabelRotation', 45);
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'best');
    
    % Customize the axes
    ax = gca;
    ax.Box = 'off'; % Remove the box
    ax.FontSize = 10; % Make the font smaller
    ax.TickDir = 'out'; % Move the ticks outside
    ax.LineWidth = 1; % Line width of the axes
end

% Adjust the figure layout
sgtitle('Diagnosis Status Counts by Group');

%%

% Load the data

% Define the group value for Cannabis in 'demodata.Gruppe'
cannabis_group_value = 2; % Adjust this value if needed to match the actual group ID in demodata.Gruppe

% Find the index for the Cannabis group
cannabis_group_idx = filtered_demodata.Gruppe == cannabis_group_value;

% Find the index for individuals with a 1 on Alkohol_diagnose in the Cannabis group
alkohol_diagnose_idx = filtered_demodata.Alkohol_diagnose == 1;

% Combine the indices to find the specific individuals for Alkohol_diagnose
cannabis_alkohol_diagnose_idx = cannabis_group_idx & alkohol_diagnose_idx;

% Find the row numbers of the individuals
rows_alkohol = find(cannabis_alkohol_diagnose_idx);

% Display the result for Alkohol_diagnose
disp('Rows of individuals in Cannabis group with a 1 on Alkohol_diagnose:');
disp(rows_alkohol);

%%


% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

% Initialize tables to store the statistics
occupation_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
    % Count the number of participants
    num_participants = sum(group_idx);
    
    % Remove rows with NaN values for each variable
    occupation_data = group_data.beruf; % Assuming 'beruf' is the column name for occupation
    occupation_data = occupation_data(~isnan(occupation_data));
    
    % Calculate Occupation statistics
    occupation_mean = mean(occupation_data);
    occupation_median = median(occupation_data);
    occupation_sd = std(occupation_data);
    occupation_min = min(occupation_data);
    occupation_max = max(occupation_data);
    occupation_row = table({group}, num_participants, occupation_mean, occupation_median, occupation_sd, occupation_min, occupation_max, ...
                           'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    occupation_stats = [occupation_stats; occupation_row];
end

% Calculate overall statistics for all participants
% Remove rows with NaN values for each variable
occupation_data_all = filtered_demodata.beruf;
occupation_data_all = occupation_data_all(~isnan(occupation_data_all));

% Occupation statistics
occupation_mean_all = mean(occupation_data_all);
occupation_median_all = median(occupation_data_all);
occupation_sd_all = std(occupation_data_all);
occupation_min_all = min(occupation_data_all);
occupation_max_all = max(occupation_data_all);
occupation_row_all = table({'All'}, length(occupation_data_all), occupation_mean_all, occupation_median_all, occupation_sd_all, occupation_min_all, occupation_max_all, ...
                           'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
occupation_stats = [occupation_stats; occupation_row_all];

% Display the table
disp('Occupation Statistics:');
disp(occupation_stats);

%% Group comparisons for Occupation

% Occupation comparison using ANOVA
[p_occupation, tbl_occupation, stats_occupation] = anova1(filtered_demodata.beruf, filtered_demodata.Gruppe);
title('ANOVA for Occupation by Group');
fprintf('ANOVA for Occupation by Group:\n');
disp(tbl_occupation);
fprintf('p-value = %.20f\n\n', p_occupation);

% Post-hoc test for occupation if ANOVA is significant
if p_occupation < 0.05
    figure;
    multcompare(stats_occupation, 'CType', 'bonferroni');
    title('Post-hoc Test for Occupation');
end

%%
% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

% Initialize tables to store the statistics
sex_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
    % Count the number of participants
    num_participants = sum(group_idx);
    
    % Remove rows with NaN values for sex
    sex_data = group_data.Geschlecht;
    sex_data = sex_data(~isnan(sex_data));
    
    % Calculate Sex statistics
    sex_male = sum(sex_data == 1);
    sex_female = sum(sex_data == 2);
    sex_row = table({group}, num_participants, sex_male, sex_female, ...
                    'VariableNames', {'Group', 'N', 'MaleCount', 'FemaleCount'});
    sex_stats = [sex_stats; sex_row];
end

% Calculate overall statistics for all participants
% Remove rows with NaN values for sex
sex_data_all = filtered_demodata.Geschlecht;
sex_data_all = sex_data_all(~isnan(sex_data_all));

% Sex statistics
sex_male_all = sum(sex_data_all == 1);
sex_female_all = sum(sex_data_all == 2);
sex_row_all = table({'All'}, length(sex_data_all), sex_male_all, sex_female_all, ...
                    'VariableNames', {'Group', 'N', 'MaleCount', 'FemaleCount'});
sex_stats = [sex_stats; sex_row_all];

% Display the table
disp('Sex Statistics:');
disp(sex_stats);

%% Chi-square test for Sex


%% Chi-square test for Sex

% Create contingency table
contingency_table = zeros(2, length(group_labels));
for i = 1:length(group_values)
    group_value = group_values(i);
    group_idx = filtered_demodata.Gruppe == group_value;
    sex_data = filtered_demodata.Geschlecht(group_idx);
    sex_data = sex_data(~isnan(sex_data));
    contingency_table(1, i) = sum(sex_data == 1);
    contingency_table(2, i) = sum(sex_data == 2);
end

% Perform chi-square test
[chi2, p_value, ~] = chi2cont(contingency_table);

% Display results
fprintf('Chi-squared test for Sex by Group:\n');
fprintf('Chi2 = %.2f\n', chi2);
fprintf('p-value = %.20f\n', p_value);

function [chi2, p, df] = chi2cont(tab)
    % Function to perform chi-square test on a contingency table
    observed = tab;
    row_totals = sum(observed, 2);
    col_totals = sum(observed, 1);
    total = sum(row_totals);
    expected = row_totals * col_totals / total;
    chi2 = sum((observed - expected).^2 ./ expected, 'all');
    df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
    p = 1 - chi2cdf(chi2, df);
end

%%
% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Opioids', 'Cannabis', 'Nicotine', 'Controls'};
group_values = [1, 2, 3, 4];

% ANOVA for Age by Group
[p_age, tbl_age, stats_age] = anova1(filtered_demodata.Alter, filtered_demodata.Gruppe);
title('ANOVA for Age by Group');
fprintf('ANOVA for Age by Group:\n');
disp(tbl_age);
fprintf('p-value = %.20f\n\n', p_age);

% Post-hoc test for age if ANOVA is significant
if p_age < 0.05
    figure;
    c_age = multcompare(stats_age, 'CType', 'bonferroni');
    title('Post-hoc Test for Age');
    posthoc_age_table = array2table(c_age, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    posthoc_age_table.Group1 = group_labels(posthoc_age_table.Group1)';
    posthoc_age_table.Group2 = group_labels(posthoc_age_table.Group2)';
    disp('Post-hoc Test Results for Age:');
    disp(posthoc_age_table);
end

% ANOVA for Education by Group
[p_education, tbl_education, stats_education] = anova1(filtered_demodata.schule, filtered_demodata.Gruppe);
title('ANOVA for Education by Group');
fprintf('ANOVA for Education by Group:\n');
disp(tbl_education);
fprintf('p-value = %.20f\n\n', p_education);

% Post-hoc test for education if ANOVA is significant
if p_education < 0.05
    figure;
    c_education = multcompare(stats_education, 'CType', 'bonferroni');
    title('Post-hoc Test for Education');
    posthoc_education_table = array2table(c_education, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    posthoc_education_table.Group1 = group_labels(posthoc_education_table.Group1)';
    posthoc_education_table.Group2 = group_labels(posthoc_education_table.Group2)';
    disp('Post-hoc Test Results for Education:');
    disp(posthoc_education_table);
end

% ANOVA for MACS_sum by Group
[p_macs, tbl_macs, stats_macs] = anova1(filtered_demodata.MACS_sum, filtered_demodata.Gruppe);
title('ANOVA for MACS_sum by Group');
fprintf('ANOVA for MACS_sum by Group:\n');
disp(tbl_macs);
fprintf('p-value = %.20f\n\n', p_macs);

% Post-hoc test for MACS if ANOVA is significant
if p_macs < 0.05
    figure;
    c_macs = multcompare(stats_macs, 'CType', 'bonferroni');
    title('Post-hoc Test for MACS_sum');
    posthoc_macs_table = array2table(c_macs, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    posthoc_macs_table.Group1 = group_labels(posthoc_macs_table.Group1)';
    posthoc_macs_table.Group2 = group_labels(posthoc_macs_table.Group2)';
    disp('Post-hoc Test Results for MACS_sum:');
    disp(posthoc_macs_table);
end

% ANOVA for Occupation by Group
[p_occupation, tbl_occupation, stats_occupation] = anova1(filtered_demodata.beruf, filtered_demodata.Gruppe);
title('ANOVA for Occupation by Group');
fprintf('ANOVA for Occupation by Group:\n');
disp(tbl_occupation);
fprintf('p-value = %.20f\n\n', p_occupation);

% Post-hoc test for occupation if ANOVA is significant
if p_occupation < 0.05
    figure;
    c_occupation = multcompare(stats_occupation, 'CType', 'bonferroni');
    title('Post-hoc Test for Occupation');
    posthoc_occupation_table = array2table(c_occupation, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    posthoc_occupation_table.Group1 = group_labels(posthoc_occupation_table.Group1)';
    posthoc_occupation_table.Group2 = group_labels(posthoc_occupation_table.Group2)';
    disp('Post-hoc Test Results for Occupation:');
    disp(posthoc_occupation_table);
end

%%
% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

% Initialize table to store the statistics for Substitut
substitut_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
    % Count the number of participants
    num_participants = sum(group_idx);
    
    % Remove rows with NaN values for Substitut variable
    substitut_data = group_data.Substitut;
    substitut_data = substitut_data(~isnan(substitut_data));
    
    % Calculate frequency of each unique value in Substitut
    [substitut_values, ~, ic] = unique(substitut_data);
    substitut_counts = accumarray(ic, 1);
    
    % Store the results in a table
    substitut_row = table(repmat({group}, length(substitut_values), 1), substitut_values, substitut_counts, ...
                          'VariableNames', {'Group', 'SubstitutValue', 'Count'});
    substitut_stats = [substitut_stats; substitut_row];
end

% Display the Substitut statistics
disp('Substitut Statistics:');
disp(substitut_stats);

%%

% Define the group labels and corresponding values in 'Gruppe' variable
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

% List of substance variables you want to analyze
substance_variables = {'Heroin', 'Moprhin', 'Phentanyl', 'Tramal', 'Opium', 'Oxycodon', 'Opioide_andere'};

% Initialize a table to store the frequencies for each group
substance_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
    for j = 1:length(substance_variables)
        substance = substance_variables{j};
        
        % Get the data for the current substance within the group
        substance_data = group_data.(substance);
        
        % Remove NaN values from the data
        substance_data = substance_data(~isnan(substance_data));
        
        % Calculate frequency of each unique value
        [substance_values, ~, ic] = unique(substance_data);
        substance_counts = accumarray(ic, 1);
        
        % Store the results in a table
        substance_row = table(repmat({group}, length(substance_values), 1), ...
                              repmat({substance}, length(substance_values), 1), ...
                              substance_values, substance_counts, ...
                              'VariableNames', {'Group', 'Substance', 'Value', 'Count'});
        substance_stats = [substance_stats; substance_row];
    end
end

% Display the substance statistics separated by groups
disp('Substance Statistics by Group:');
disp(substance_stats);

%%
% Load the data (Assume 'filtered_demodata' is already loaded)

% Define the group labels and corresponding values in 'demodata.Gruppe'
group_labels = {'Control', 'Nicotine', 'Cannabis', 'Opiates'};
group_values = [4, 3, 2, 1]; % Adjust these values if needed to match the actual group IDs in demodata.Gruppe

% List of original diagnosis variables in 'filtered_demodata'
original_diagnose_variables = {
    'Nikton', 'Cannabis', 'Alkohol', 'Opioide', ...
    'Amph', 'XTC', 'Kokain', 'Benzo', 'Prega'
};

% Corrected list of diagnosis variable names for labeling
corrected_diagnose_labels = {
    'Nicotine', 'Cannabis', 'Alcohol', 'Opioids', ...
    'Amphetamines', 'Ecstasy', 'Cocaine', 'Benzodiazepines', 'Pregabalin'
};

% Preallocate data for the bar charts
none_counts = zeros(length(group_labels), length(original_diagnose_variables));
present_counts = zeros(length(group_labels), length(original_diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(original_diagnose_variables));

for var = 1:length(original_diagnose_variables)
    variable = [original_diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group = group_labels{i};
        group_value = group_values(i);
        
        % Extract the rows for the current group
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);
        
        % Count occurrences of each diagnosis status
        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Plot grouped bar charts for each group
figure;
for i = 1:length(group_labels)
    group = group_labels{i};
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar_handle = bar(bar_data', 'grouped');
    
    % Customizing the bar appearance
    for j = 1:length(bar_handle)
        bar_handle(j).LineWidth = 1.5;
    end
    
    % Titles and labels
    title(group, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Diagnosis Variables', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    set(gca, 'XTickLabel', corrected_diagnose_labels, 'XTickLabelRotation', 45, 'FontSize', 10);
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'northeastoutside', 'FontSize', 10);
    
    % Customize the axes
    ax = gca;
    ax.Box = 'off'; % Remove the box
    ax.FontSize = 10; % Font size for axes labels
    ax.TickDir = 'out'; % Move the ticks outside
    ax.LineWidth = 1; % Line width of the axes
end

% Adjust the figure layout and add a super title
sgtitle('Diagnosis Status Counts by Group', 'FontSize', 16, 'FontWeight', 'bold');

%% same as above but prettier

% Load the data (Assume 'filtered_demodata' is already loaded)

% Define the group labels and corresponding values in 'demodata.Gruppe'
group_labels = {'Control', 'Nicotine', 'Cannabis', 'Opiates'};
group_values = [4, 3, 2, 1]; % Adjust these values if needed to match the actual group IDs in demodata.Gruppe

% List of original diagnosis variables and corresponding corrected labels
diagnose_variables = {
    'Nikton', 'Cannabis', 'Alkohol', 'Opioide', ...
    'Amph', 'XTC', 'Kokain', 'Benzo', 'Prega'
};
corrected_labels = {
    'Nicotine', 'Cannabis', 'Alcohol', 'Opioids', ...
    'Amphetamines', 'Ecstasy', 'Cocaine', 'Benzodiazepines', 'Pregabalin'
};

% Preallocate data for the bar charts
none_counts = zeros(length(group_labels), length(diagnose_variables));
present_counts = zeros(length(group_labels), length(diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(diagnose_variables));

for var = 1:length(diagnose_variables)
    variable = [diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group_value = group_values(i);
        
        % Extract the rows for the current group
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);
        
        % Count occurrences of each diagnosis status
        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Colors for the bars (limited to three)
bar_colors = [0.2, 0.2, 0.8; 0.8, 0.2, 0.2; 0.9, 0.8, 0.2];

% Plot grouped bar charts for each group
figure('Color', 'w'); % Set background color to white
for i = 1:length(group_labels)
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar_handle = bar(bar_data', 'grouped');
    
    % Customizing the bar appearance
    for j = 1:length(bar_handle)
        bar_handle(j).EdgeColor = 'none'; % Remove the bar outline
        bar_handle(j).FaceColor = bar_colors(j, :); % Apply the limited colors
    end
    
    % Titles and labels
    title(group_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Diagnosis Variables', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    set(gca, 'XTickLabel', corrected_labels, 'XTickLabelRotation', 45, 'FontSize', 10);
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'northeastoutside', 'FontSize', 10);
    
    % Customize the axes
    ax = gca;
    ax.Box = 'off'; % Remove the box
    ax.FontSize = 10; % Font size for axes labels
    ax.TickDir = 'out'; % Move the ticks outside
    ax.LineWidth = 1; % Line width of the axes
    ax.Color = 'w'; % Set the axes background color to white
end

% Adjust the figure layout and add a super title
sgtitle('Diagnosis Status Counts by Group', 'FontSize', 16, 'FontWeight', 'bold');

%% Same but horizontal

% Load the data (Assume 'filtered_demodata' is already loaded)

% Define the group labels and corresponding values in 'demodata.Gruppe'
group_labels = {'Control', 'Nicotine', 'Cannabis', 'Opiates'};
group_values = [4, 3, 2, 1]; % Adjust these values if needed to match the actual group IDs in demodata.Gruppe

% List of original diagnosis variables and corresponding corrected labels
diagnose_variables = {
    'Nikton', 'Cannabis', 'Alkohol', 'Opioide', ...
    'Amph', 'XTC', 'Kokain', 'Benzo', 'Prega'
};
corrected_labels = {
    'Nicotine', 'Cannabis', 'Alcohol', 'Opioids', ...
    'Amphetamines', 'Ecstasy', 'Cocaine', 'Benzodiazepines', 'Pregabalin'
};

% Preallocate data for the bar charts
none_counts = zeros(length(group_labels), length(diagnose_variables));
present_counts = zeros(length(group_labels), length(diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(diagnose_variables));

for var = 1:length(diagnose_variables)
    variable = [diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group_value = group_values(i);
        
        % Extract the rows for the current group
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);
        
        % Count occurrences of each diagnosis status
        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Colors for the bars (limited to three)
bar_colors = [0.2, 0.2, 0.8; 0.8, 0.2, 0.2; 0.9, 0.8, 0.2];

% Plot horizontal bar charts for each group
figure('Color', 'w'); % Set background color to white
for i = 1:length(group_labels)
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar_handle = barh(bar_data', 'grouped'); % Use barh for horizontal bars
    
    % Customizing the bar appearance
    for j = 1:length(bar_handle)
        bar_handle(j).EdgeColor = 'none'; % Remove the bar outline
        bar_handle(j).FaceColor = bar_colors(j, :); % Apply the limited colors
    end
    
    % Titles and labels
    title(group_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Diagnosis Variables', 'FontSize', 12); % Y-label is now for the variables
    xlabel('Count', 'FontSize', 12); % X-label is now for the count
    set(gca, 'YTickLabel', corrected_labels, 'FontSize', 10); % Corrected labels on the Y-axis
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'northeastoutside', 'FontSize', 10);
    
    % Customize the axes
    ax = gca;
    ax.Box = 'off'; % Remove the box
    ax.FontSize = 10; % Font size for axes labels
    ax.TickDir = 'out'; % Move the ticks outside
    ax.LineWidth = 1; % Line width of the axes
    ax.Color = 'w'; % Set the axes background color to white
end

% Adjust the figure layout and add a super title
sgtitle('Diagnosis Status Counts by Group', 'FontSize', 16, 'FontWeight', 'bold');

%% thicker bars
% Load the data (Assume 'filtered_demodata' is already loaded)

% Define the group labels and corresponding values in 'demodata.Gruppe'
group_labels = {'Control', 'Nicotine', 'Cannabis', 'Opiates'};
group_values = [4, 3, 2, 1]; % Adjust these values if needed to match the actual group IDs in demodata.Gruppe

% List of original diagnosis variables and corresponding corrected labels (reversed order)
diagnose_variables = {
    'Prega', 'Benzo', 'Kokain', 'XTC', 'Amph', 'Opioide', 'Alkohol', 'Cannabis', 'Nikton'
};
corrected_labels = {
    'Pregabalin', 'Benzodiazepines', 'Cocaine', 'Ecstasy', 'Amphetamines', 'Opioids', 'Alcohol', 'Cannabis', 'Nicotine'
};

% Preallocate data for the bar charts
none_counts = zeros(length(group_labels), length(diagnose_variables));
present_counts = zeros(length(group_labels), length(diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(diagnose_variables));

for var = 1:length(diagnose_variables)
    variable = [diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group_value = group_values(i);
        
        % Extract the rows for the current group
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);
        
        % Count occurrences of each diagnosis status
        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Colors for the bars (limited to three)
bar_colors = [0.2, 0.2, 0.8; 0.8, 0.2, 0.2; 0.9, 0.8, 0.2];

% Plot horizontal bar charts for each group
figure('Color', 'w'); % Set background color to white
for i = 1:length(group_labels)
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar_handle = barh(bar_data', 'grouped'); % Use default BarWidth to avoid overlapping
    
    % Customizing the bar appearance
    for j = 1:length(bar_handle)
        bar_handle(j).EdgeColor = 'none'; % Remove the bar outline
        bar_handle(j).FaceColor = bar_colors(j, :); % Apply the limited colors
    end
    
    % Titles and labels
    title(group_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Count', 'FontSize', 12); % X-label is now for the count
    set(gca, 'YTickLabel', corrected_labels, 'FontSize', 10); % Corrected labels on the Y-axis
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'northeastoutside', 'FontSize', 10);
    
    % Customize the axes
    ax = gca;
    ax.Box = 'off'; % Remove the box
    ax.FontSize = 10; % Font size for axes labels
    ax.TickDir = 'out'; % Move the ticks outside
    ax.LineWidth = 1; % Line width of the axes
    ax.Color = 'w'; % Set the axes background color to white
end

% Adjust the figure layout and add a super title
sgtitle('Diagnosis Status Counts by Group', 'FontSize', 16, 'FontWeight', 'bold');



%%
% Define the group labels and corresponding values in 'groupCodes' variable
group_labels = {'Cocaine', 'Controls', 'Nicotine'};
group_values = [3, 1, 2];

% Initialize tables to store the statistics
age_stats = table();
sex_stats = table();
education_stats = table();
cocageonset_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = groupCodes == group_value;
    group_data = participants(group_idx, :);
    
    % Count the number of participants
    num_participants = sum(group_idx);
    
    % Remove rows with NaN values for each variable
    age_data = group_data.age;
    age_data = age_data(~isnan(age_data));
    
    sex_data = group_data.sex;
    sex_data = sex_data(~isnan(sex_data));
    
    education_data = group_data.educ_yr; % Assuming 'educ_yr' is the column name for education
    education_data = education_data(~isnan(education_data));
    
    cocageonset_data = group_data.cocageonset;
    cocageonset_data = cocageonset_data(~isnan(cocageonset_data));
    
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
    
    % Calculate cocageonset statistics
    cocageonset_mean = mean(cocageonset_data);
    cocageonset_median = median(cocageonset_data);
    cocageonset_sd = std(cocageonset_data);
    cocageonset_min = min(cocageonset_data);
    cocageonset_max = max(cocageonset_data);
    cocageonset_row = table({group}, num_participants, cocageonset_mean, cocageonset_median, cocageonset_sd, cocageonset_min, cocageonset_max, ...
                     'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    cocageonset_stats = [cocageonset_stats; cocageonset_row];
end

% Calculate overall statistics for all participants
% Remove rows with NaN values for each variable
age_data_all = participants.age;
age_data_all = age_data_all(~isnan(age_data_all));

sex_data_all = participants.sex;
sex_data_all = sex_data_all(~isnan(sex_data_all));

education_data_all = participants.educ_yr;
education_data_all = education_data_all(~isnan(education_data_all));

cocageonset_data_all = participants.cocageonset;
cocageonset_data_all = cocageonset_data_all(~isnan(cocageonset_data_all));

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

% cocageonset statistics
cocageonset_mean_all = mean(cocageonset_data_all);
cocageonset_median_all = median(cocageonset_data_all);
cocageonset_sd_all = std(cocageonset_data_all);
cocageonset_min_all = min(cocageonset_data_all);
cocageonset_max_all = max(cocageonset_data_all);
cocageonset_row_all = table({'All'}, length(cocageonset_data_all), cocageonset_mean_all, cocageonset_median_all, cocageonset_sd_all, cocageonset_min_all, cocageonset_max_all, ...
                     'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
cocageonset_stats = [cocageonset_stats; cocageonset_row_all];

% Display the tables
disp('Age Statistics:');
disp(age_stats);

disp('Sex Statistics:');
disp(sex_stats);

disp('Education Statistics:');
disp(education_stats);

disp('cocageonset Statistics:');
disp(cocageonset_stats);

%% Group comparisons


% Age comparison using ANOVA
[p_age, tbl_age, stats_age] = anova1(participants.age, groupCodes);
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
[p_education, tbl_education, stats_education] = anova1(participants.educ_yr, groupCodes);
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

% cocageonset comparison using ANOVA
[p_cocageonset, tbl_cocageonset, stats_cocageonset] = anova1(participants.cocageonset, groupCodes);
title('ANOVA for cocageonset by Group');
fprintf('ANOVA for cocageonset by Group:\n');
disp(tbl_cocageonset);
fprintf('p-value = %.20f\n\n', p_cocageonset);

% Post-hoc test for cocageonset if ANOVA is significant
if p_cocageonset < 0.05
    figure;
    multcompare(stats_cocageonset, 'CType', 'bonferroni');
    title('Post-hoc Test for cocageonset');
end

% Sex comparison using Chi-squared test
sex_data = participants.sex;
sex_groups = groupCodes;
[crosstab, chi2, p_sex, labels] = crosstab(sex_data, sex_groups);

% Display results
fprintf('Chi-squared test for Sex by Group:\n');
fprintf('Chi2 = %.2f\n', chi2);
fprintf('p-value = %.20f\n', p_sex);

%% Posthoc education
% Ensure that the group_labels are correctly ordered
% Check the output of 'c_education' to verify group indexing

% Define the group labels and corresponding values in 'groupCodes' variable
group_labels = {'Cocaine', 'Controls', 'Nicotine'};
group_values = [3, 1, 2];


if p_education < 0.05
    % Perform post-hoc comparisons using Bonferroni correction
    figure;
    c_education = multcompare(stats_education, 'CType', 'bonferroni');
    title('Post-hoc Test for Education by Group');

    % Display the post-hoc test results
    posthoc_education_table = array2table(c_education, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
    
    % Verify that group_labels are in the correct order by printing c_education
    disp('c_education Output:');
    disp(c_education);

    % Replace group indices with group names using correct ordering
    posthoc_education_table.Group1 = group_labels(posthoc_education_table.Group1)';
    posthoc_education_table.Group2 = group_labels(posthoc_education_table.Group2)';
    
    disp('Post-hoc Test Results for Education:');
    disp(posthoc_education_table);
end

%% Filter ccqn
% Step 1: Convert participant_id to a string array and remove the '-'
% Assume participant_id is numeric and needs to be converted to string
if isnumeric(participants.participant_id)
    participant_ids = string(participants.participant_id); % Convert numeric IDs to string
else
    participant_ids = participants.participant_id; % If already strings, no conversion needed
end

% Remove the '-' from the IDs (this step is redundant if participant_id was numeric)
participant_ids = strrep(participant_ids, '-', ''); 

% Step 2: Convert CCQN.rid to a string array if it isn't already
if isnumeric(CCQN.rid)
    CCQN.rid = string(CCQN.rid); % Convert numeric IDs to string
elseif iscell(CCQN.rid)
    CCQN.rid = string(CCQN.rid); % Convert cell array of strings to string array
end

% Step 3: Find rows in CCQN that have a matching ID in participant_ids
matching_rows = ismember(CCQN.rid, participant_ids);

% Step 4: Keep only the rows in CCQN that have matching IDs
CCQN_filtered = CCQN(matching_rows, :);

% Step 5: (Optional) Display the filtered CCQN
disp('Filtered CCQN:');
disp(CCQN_filtered);

%% Group comparison CCQN
% Define the group labels and corresponding values
group_labels = {'Cocaine', 'Controls', 'Nicotine'};
group_values = [3, 1, 2];

% Initialize tables to store the statistics
ccqn_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group using the groupCodes variable
    group_idx = groupCodes == group_value;
    group_data = CCQN_filtered(group_idx, :);
    
    % Count the number of participants in this group
    num_participants = sum(group_idx);
    
    % Extract ccqnscore and replace NaNs with 0
    variable_data = group_data.ccqnscore;
    variable_data(isnan(variable_data)) = 0; 
    
    % Calculate statistics
    variable_mean = mean(variable_data);
    variable_median = median(variable_data);
    variable_sd = std(variable_data);
    variable_min = min(variable_data);
    variable_max = max(variable_data);
    stats_row = table({group}, num_participants, variable_mean, variable_median, variable_sd, variable_min, variable_max, ...
                      'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
    ccqn_stats = [ccqn_stats; stats_row];
end

% Calculate overall statistics for all participants
variable_data_all = CCQN_filtered.ccqnscore;
variable_data_all(isnan(variable_data_all)) = 0;

% Overall statistics
variable_mean_all = mean(variable_data_all);
variable_median_all = median(variable_data_all);
variable_sd_all = std(variable_data_all);
variable_min_all = min(variable_data_all);
variable_max_all = max(variable_data_all);
overall_row = table({'All'}, length(variable_data_all), variable_mean_all, variable_median_all, variable_sd_all, variable_min_all, variable_max_all, ...
                    'VariableNames', {'Group', 'N', 'Mean', 'Median', 'SD', 'Min', 'Max'});
ccqn_stats = [ccqn_stats; overall_row];

% Display the statistics table
disp('Group Statistics for CCQN_filtered:');
disp(ccqn_stats);

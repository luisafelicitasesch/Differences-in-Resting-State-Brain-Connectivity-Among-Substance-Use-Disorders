%% Define Group Labels and Values
group_labels = {'Cannabis', 'Controls', 'Opiates', 'Nicotine'};
group_values = [2, 4, 1, 3];

%% Function: Calculate Group Statistics
function stats = calculate_group_stats(group_data, variable)
    data = group_data.(variable);
    data = data(~isnan(data));
    stats.mean = mean(data);
    stats.median = median(data);
    stats.sd = std(data);
    stats.min = min(data);
    stats.max = max(data);
    stats.n = length(data);
end

%% Function: Display Group Statistics
function display_group_stats(stats, stat_name)
    disp([stat_name ' Statistics:']);
    disp(struct2table(stats));
end

%% Function: Run ANOVA and Post-Hoc Tests
function run_anova(variable, group, labels, var_name)
    [p, tbl, stats] = anova1(variable, group);
    title(['ANOVA for ' var_name ' by Group']);
    fprintf('ANOVA for %s by Group:\n', var_name);
    disp(tbl);
    fprintf('p-value = %.20f\n\n', p);

    % Post-hoc test if ANOVA is significant
    if p < 0.05
        figure;
        c = multcompare(stats, 'CType', 'bonferroni');
        title(['Post-hoc Test for ' var_name]);
        posthoc_table = array2table(c, 'VariableNames', {'Group1', 'Group2', 'LowerCI', 'MeanDiff', 'UpperCI', 'pValue'});
        posthoc_table.Group1 = labels(posthoc_table.Group1)';
        posthoc_table.Group2 = labels(posthoc_table.Group2)';
        disp(['Post-hoc Test Results for ' var_name ':']);
        disp(posthoc_table);
    end
end

%% Function: Chi-Square Test for Sex
function run_chi_square_test(sex_data, group_data, group_labels)
    contingency_table = zeros(2, length(group_labels));
    for i = 1:length(group_labels)
        group_idx = group_data == i;
        group_sex_data = sex_data(group_idx);
        group_sex_data = group_sex_data(~isnan(group_sex_data));
        contingency_table(1, i) = sum(group_sex_data == 1);
        contingency_table(2, i) = sum(group_sex_data == 2);
    end

    [chi2, p, ~] = chi2cont(contingency_table);
    fprintf('Chi-squared test for Sex by Group:\n');
    fprintf('Chi2 = %.2f\n', chi2);
    fprintf('p-value = %.20f\n', p);
end

%% Filter Demographic Data by Participants with CorrMatrices
ids_to_keep = participantNames;
rows_to_keep = ismember(demodata.id, ids_to_keep);
filtered_demodata = demodata(rows_to_keep, :);

%% Calculate and Display Demographic Statistics
stats_struct = struct();
for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    group_data = filtered_demodata(filtered_demodata.Gruppe == group_value, :);

    % Calculate and store statistics
    stats_struct.(group).age = calculate_group_stats(group_data, 'Alter');
    stats_struct.(group).sex = calculate_group_stats(group_data, 'Geschlecht');
    stats_struct.(group).education = calculate_group_stats(group_data, 'schule');
    stats_struct.(group).macs = calculate_group_stats(group_data, 'MACS_sum');
    stats_struct.(group).occupation = calculate_group_stats(group_data, 'beruf');
end

% Calculate overall statistics for all participants
stats_struct.All.age = calculate_group_stats(filtered_demodata, 'Alter');
stats_struct.All.sex = calculate_group_stats(filtered_demodata, 'Geschlecht');
stats_struct.All.education = calculate_group_stats(filtered_demodata, 'schule');
stats_struct.All.macs = calculate_group_stats(filtered_demodata, 'MACS_sum');
stats_struct.All.occupation = calculate_group_stats(filtered_demodata, 'beruf');

% Display the statistics
fields = fieldnames(stats_struct.All);
for i = 1:length(fields)
    display_group_stats(stats_struct.All.(fields{i}), fields{i});
end

%% Group Comparisons Using ANOVA
run_anova(filtered_demodata.Alter, filtered_demodata.Gruppe, group_labels, 'Age');
run_anova(filtered_demodata.schule, filtered_demodata.Gruppe, group_labels, 'Education');
run_anova(filtered_demodata.MACS_sum, filtered_demodata.Gruppe, group_labels, 'MACS_sum');
run_anova(filtered_demodata.beruf, filtered_demodata.Gruppe, group_labels, 'Occupation');

%% Chi-Square Test for Sex
run_chi_square_test(filtered_demodata.Geschlecht, filtered_demodata.Gruppe, group_labels);

%% Group Comparisons for Specific Substance Frequencies
substance_variables = {'Heroin', 'Moprhin', 'Phentanyl', 'Tramal', 'Opium', 'Oxycodon', 'Opioide_andere'};
substance_stats = struct();

for i = 1:length(substance_variables)
    substance = substance_variables{i};
    for j = 1:length(group_labels)
        group = group_labels{j};
        group_value = group_values(j);
        group_data = filtered_demodata(filtered_demodata.Gruppe == group_value, :);
        substance_stats.(substance).(group) = calculate_group_stats(group_data, substance);
    end
end

% Display substance statistics (if necessary)
disp('Substance Statistics by Group:');
disp(substance_stats);

%% Substitution Substance Frequency Analysis
substitut_stats = table();

for i = 1:length(group_labels)
    group = group_labels{i};
    group_value = group_values(i);
    
    % Extract the rows for the current group
    group_idx = filtered_demodata.Gruppe == group_value;
    group_data = filtered_demodata(group_idx, :);
    
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

%% Charts for SUD Criteria by Group

diagnose_variables = {'Nikton', 'Cannabis', 'Alkohol', 'Opioide', 'Amph', 'XTC', 'Kokain', 'Benzo', 'Prega'};
corrected_labels = {'Nicotine', 'Cannabis', 'Alcohol', 'Opioids', 'Amphetamines', 'Ecstasy', 'Cocaine', 'Benzodiazepines', 'Pregabalin'};
none_counts = zeros(length(group_labels), length(diagnose_variables));
present_counts = zeros(length(group_labels), length(diagnose_variables));
lifetime_counts = zeros(length(group_labels), length(diagnose_variables));

for var = 1:length(diagnose_variables)
    variable = [diagnose_variables{var} '_diagnose'];
    for i = 1:length(group_labels)
        group_value = group_values(i);
        group_idx = filtered_demodata.Gruppe == group_value;
        group_data = filtered_demodata.(variable)(group_idx);

        none_counts(i, var) = sum(group_data == 0);
        present_counts(i, var) = sum(group_data == 1);
        lifetime_counts(i, var) = sum(group_data == 2);
    end
end

% Plot grouped bar charts for each group
figure('Color', 'w');
bar_colors = [0.2, 0.2, 0.8; 0.8, 0.2, 0.2; 0.9, 0.8, 0.2];

for i = 1:length(group_labels)
    subplot(2, 2, i);
    bar_data = [none_counts(i, :); present_counts(i, :); lifetime_counts(i, :)];
    bar_handle = bar(bar_data', 'grouped');

    for j = 1:length(bar_handle)
        bar_handle(j).EdgeColor = 'none';
        bar_handle(j).FaceColor = bar_colors(j, :);
    end

    title(group_labels{i}, 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Diagnosis Variables', 'FontSize', 12);
    ylabel('Count', 'FontSize', 12);
    set(gca, 'XTickLabel', corrected_labels, 'XTickLabelRotation', 45, 'FontSize', 10);
    legend({'None', 'Present', 'Lifetime'}, 'Location', 'northeastoutside', 'FontSize', 10);

    ax = gca;
    ax.Box = 'off';
    ax.FontSize = 10;
    ax.TickDir = 'out';
    ax.LineWidth = 1;
    ax.Color = 'w';
end

sgtitle('Diagnosis Status Counts by Group', 'FontSize', 16, 'FontWeight', 'bold');

%% Helper Function for Chi-Square Test
function [chi2, p, df] = chi2cont(tab)
    observed = tab;
    row_totals = sum(observed, 2);
    col_totals = sum(observed, 1);
    total = sum(row_totals);
    expected = row_totals * col_totals / total;
    chi2 = sum((observed - expected).^2 ./ expected, 'all');
    df = (size(observed, 1) - 1) * (size(observed, 2) - 1);
    p = 1 - chi2cdf(chi2, df);
end

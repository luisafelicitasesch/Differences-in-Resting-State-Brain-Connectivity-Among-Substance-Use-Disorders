%% Imputation education NaNs

% Fit a linear regression model using other variables to predict education years
model = fitlm(participants, 'ResponseVar', 'educ_yr', 'PredictorVars', {'age', 'sex', 'group'});

% Predict the missing values
predicted_values = predict(model, participants(isnan(participants.educ_yr), :));

% Fill the NaNs with the predicted values
participants.educ_yr(isnan(participants.educ_yr)) = predicted_values;


%%
% Extract demographic data directly
sex_reordered = participants.sex;
age_reordered = participants.age;
education_reordered = participants.educ_yr;

% Check if the lengths match the number of participants
if length(Cglob_all) ~= length(participantNames) || ...
   length(Tglob_all) ~= length(participantNames) || ...
   length(Eglob_all) ~= length(participantNames) || ...
   length(Pglob_all) ~= length(participantNames)
    error('Lengths of global measures do not match the number of participants.');
end

% Create separate tables for each global measure

% Table for Cglob
demo_Cglob = table(groupCodes, participantNames, sex_reordered, age_reordered, education_reordered, Cglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Cglob'});
% Display the Cglob table
disp('Table for Cglob:');
disp(demo_Cglob);

% Table for Tglob
demo_Tglob = table(groupCodes, participantNames, sex_reordered, age_reordered, education_reordered, Tglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Tglob'});
% Display the Tglob table
disp('Table for Tglob:');
disp(demo_Tglob);

% Table for Eglob
demo_Eglob = table(groupCodes, participantNames, sex_reordered, age_reordered, education_reordered, Eglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Eglob'});
% Display the Eglob table
disp('Table for Eglob:');
disp(demo_Eglob);

% Table for Pglob
demo_Pglob = table(groupCodes, participantNames, sex_reordered, age_reordered, education_reordered, Pglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Pglob'});
% Display the Pglob table
disp('Table for Pglob:');
disp(demo_Pglob);


%% Fit Linear Models and Calculate Residuals for Pglob, Cglob, Tglob, and Eglob by Group

% Fit the linear model for Pglob and calculate residuals
mdl_Pglob = fitlm(demo_Pglob(:, 3:5), demo_Pglob.Pglob);
residuals_Pglob = mdl_Pglob.Residuals.Raw;

% Separate residuals by group for Pglob
group_names = {'Controls', 'Nicotine', 'Cocaine'};
residuals_Pglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = demo_Pglob.Group == i;
    residuals_Pglob_by_group.(group_name) = residuals_Pglob(group_idx);
end

% Display residuals by group for Pglob
disp('Residuals by group for Pglob:');
disp(residuals_Pglob_by_group);

% Fit the linear model for Cglob and calculate residuals
mdl_Cglob = fitlm(demo_Cglob(:, 3:5), demo_Cglob.Cglob);
residuals_Cglob = mdl_Cglob.Residuals.Raw;

% Separate residuals by group for Cglob
residuals_Cglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = demo_Cglob.Group == i;
    residuals_Cglob_by_group.(group_name) = residuals_Cglob(group_idx);
end

% Display residuals by group for Cglob
disp('Residuals by group for Cglob:');
disp(residuals_Cglob_by_group);

% Fit the linear model for Tglob and calculate residuals
mdl_Tglob = fitlm(demo_Tglob(:, 3:5), demo_Tglob.Tglob);
residuals_Tglob = mdl_Tglob.Residuals.Raw;

% Separate residuals by group for Tglob
residuals_Tglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = demo_Tglob.Group == i;
    residuals_Tglob_by_group.(group_name) = residuals_Tglob(group_idx);
end

% Display residuals by group for Tglob
disp('Residuals by group for Tglob:');
disp(residuals_Tglob_by_group);

% Fit the linear model for Eglob and calculate residuals
mdl_Eglob = fitlm(demo_Eglob(:, 3:5), demo_Eglob.Eglob);
residuals_Eglob = mdl_Eglob.Residuals.Raw;

% Separate residuals by group for Eglob
residuals_Eglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = demo_Eglob.Group == i;
    residuals_Eglob_by_group.(group_name) = residuals_Eglob(group_idx);
end

% Display residuals by group for Eglob
disp('Residuals by group for Eglob:');
disp(residuals_Eglob_by_group);

%% Calculate and Separate Residuals by Group for Various Local Network Measures Across Regions

% Extract demographic data directly from participants
sex_reordered = participants.sex;
age_reordered = participants.age;
education_reordered = participants.educ_yr;

% Assuming BCloc_all, Cloc_all, Eloc_all, ECloc_all, and Dloc_all are already calculated and available

% Define local measures and their names
local_measures_data = {BCloc_all, Cloc_all, Eloc_all, ECloc_all, Dloc_all};
local_measures_names = {'BCloc', 'Cloc', 'Eloc', 'ECloc', 'Dloc'};

% Loop through each local measure
for k = 1:length(local_measures_names)
    measure_name = local_measures_names{k};
    aucData = local_measures_data{k};
    
    % Initialize a cell array to store residuals for each group and region
    residuals_by_group = struct();
    
    % Loop through each region
    numRegions = size(aucData, 2); % Number of regions
    for region = 1:numRegions
        % Extract data for the current region
        AUC_all = aucData(:, region);
        
        % Create a table for the current region
        demo_region = table(groupCodes, participantNames, sex_reordered, age_reordered, education_reordered, AUC_all, ...
                            'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'AUC'});
        
        % Remove constant or nearly constant predictors
        X = demo_region(:, 3:5);
        X = X{:, :}; % Convert to matrix
        X(:, var(X) == 0) = []; % Remove constant predictors
        
        % Check for empty predictor matrix after removing constants
        if isempty(X)
            residuals_region = AUC_all; % If no predictors left, use raw AUC data as residuals
        else
            % Fit the linear model for the current region and calculate residuals
            mdl_region = fitlm(X, AUC_all);
            residuals_region = mdl_region.Residuals.Raw;
        end
        
        % Separate residuals by group for the current region
        unique_groups = unique(groupCodes);
        group_names = {'Controls', 'Nicotine', 'Cocaine'};
        for i = 1:length(unique_groups)
            group_code = unique_groups(i);
            group_idx = (demo_region.Group == group_code);
            group_name = group_names{group_code};
            if ~isfield(residuals_by_group, group_name)
                residuals_by_group.(group_name) = NaN(numRegions, sum(group_idx)); % Preallocate with NaNs
            end
            residuals_by_group.(group_name)(region, :) = residuals_region(group_idx);
        end
    end
    
    % Save the residuals for the current measure
    eval([measure_name '_residuals_by_group = residuals_by_group;']);
end

% Display residuals by group for each local measure
for k = 1:length(local_measures_names)
    measure_name = local_measures_names{k};
    disp(['Residuals by group for ' measure_name ':']);
    eval(['disp(' measure_name '_residuals_by_group);']);
end

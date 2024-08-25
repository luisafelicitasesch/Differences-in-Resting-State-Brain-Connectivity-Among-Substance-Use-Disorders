%% Reordering and Filtering Demographic Data Based on Participant Information

% Extract participant IDs from demodata
demodata_ids = demodata.id;

% Initialize cell arrays to store the reordered data
sex_reordered = cell(length(participantNames), 1);
age_reordered = NaN(length(participantNames), 1);
education_reordered = cell(length(participantNames), 1);
beruf_reordered = cell(length(participantNames), 1);
site_reordered = cell(length(participantNames), 1); % New variable for site

% Loop through participantNames and reorder the data
for i = 1:length(participantNames)
    % Find the index of the current participant in demodata
    idx = find(strcmp(demodata_ids, participantNames{i}), 1);
    
    % If the participant is found in demodata, extract the corresponding data
    if ~isempty(idx)
        % Extract and convert data as necessary
        sex_value = demodata.Geschlecht(idx);
        if iscell(sex_value)
            sex_reordered{i} = sex_value{1};
        elseif isnumeric(sex_value)
            sex_reordered{i} = num2str(sex_value);
        elseif ischar(sex_value) || isstring(sex_value)
            sex_reordered{i} = char(sex_value);
        else
            sex_reordered{i} = 'Unknown';
        end
        
        age_reordered(i) = demodata.Alter(idx); % Assuming numeric
        
        education_value = demodata.schule(idx);
        if iscell(education_value)
            education_reordered{i} = education_value{1};
        elseif isnumeric(education_value)
            education_reordered{i} = num2str(education_value);
        elseif ischar(education_value) || isstring(education_value)
            education_reordered{i} = char(education_value);
        else
            education_reordered{i} = 'Unknown';
        end
        
        beruf_value = demodata.beruf(idx);
        if iscell(beruf_value)
            beruf_reordered{i} = beruf_value{1};
        elseif isnumeric(beruf_value)
            beruf_reordered{i} = num2str(beruf_value);
        elseif ischar(beruf_value) || isstring(beruf_value)
            beruf_reordered{i} = char(beruf_value);
        else
            beruf_reordered{i} = 'Unknown';
        end
        
        % Extract site information based on specific parts of the ID
        % Grouping BC, BG, BO, and BR as site 1, and MC, MG, MO, and MR as site 2
        if contains(participantNames{i}, {'BC', 'BG', 'BO', 'BR'})
            site_reordered{i} = '1';
        elseif contains(participantNames{i}, {'MC', 'MG', 'MO', 'MR'})
            site_reordered{i} = '2';
        else
            site_reordered{i} = 'Unknown';
        end
    else
        % If the participant is not found, set as empty or some placeholder
        sex_reordered{i} = 'Unknown';
        age_reordered(i) = NaN;
        education_reordered{i} = 'Unknown';
        beruf_reordered{i} = 'Unknown';
        site_reordered{i} = 'Unknown'; % For site
    end
end

% Identify participants with any missing data
any_missing = cellfun(@(x) isequal(x, 'Unknown'), sex_reordered) | ...
              isnan(age_reordered) | ...
              cellfun(@(x) isequal(x, 'Unknown'), education_reordered) | ...
              cellfun(@(x) isequal(x, 'Unknown'), beruf_reordered) | ...
              cellfun(@(x) isequal(x, 'Unknown'), site_reordered);

% Display participants with any missing data
disp('Participants with any missing data:');
missing_participants = participantNames(any_missing);
disp(missing_participants);

% Check if there are any participants left after filtering
if all(any_missing)
    error('All participants have missing data, cannot proceed with analysis.');
end

% Filter out participants with any missing data
sex_reordered = sex_reordered(~any_missing);
age_reordered = age_reordered(~any_missing);
education_reordered = education_reordered(~any_missing);
beruf_reordered = beruf_reordered(~any_missing);
site_reordered = site_reordered(~any_missing);
participantNames = participantNames(~any_missing);
participantGroups = participantGroups(~any_missing);

% Filter the global measures to exclude participants with missing data
Cglob_all = Cglob_all(~any_missing);
Tglob_all = Tglob_all(~any_missing);
Eglob_all = Eglob_all(~any_missing);
Pglob_all = Pglob_all(~any_missing);

% Check if there are any remaining NaN or 'Unknown' values after exclusion
disp('Remaining NaN or unknown values after exclusion:');
disp('Age (NaN indices):');
disp(find(isnan(age_reordered)));

disp('Sex (Unknown/empty indices):');
disp(find(cellfun(@(x) isequal(x, 'Unknown'), sex_reordered)));

disp('Education (Unknown/empty indices):');
disp(find(cellfun(@(x) isequal(x, 'Unknown'), education_reordered)));

disp('Beruf (Unknown/empty indices):');
disp(find(cellfun(@(x) isequal(x, 'Unknown'), beruf_reordered)));

disp('Site (Unknown/empty indices):');
disp(find(cellfun(@(x) isequal(x, 'Unknown'), site_reordered)));

% Check if the lengths match the number of participants
if length(Cglob_all) ~= length(participantNames) || ...
   length(Tglob_all) ~= length(participantNames) || ...
   length(Eglob_all) ~= length(participantNames) || ...
   length(Pglob_all) ~= length(participantNames)
    error('Lengths of global measures do not match the number of participants.');
end

%% Create separate tables for each global measure

% Table for Cglob
demo_Cglob = table(participantGroups, participantNames, sex_reordered, age_reordered, education_reordered, beruf_reordered, site_reordered, Cglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Beruf', 'Site', 'Cglob'});
% Display the Cglob table
disp('Table for Cglob:');
disp(demo_Cglob);

% Table for Tglob
demo_Tglob = table(participantGroups, participantNames, sex_reordered, age_reordered, education_reordered, beruf_reordered, site_reordered, Tglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Beruf', 'Site', 'Tglob'});
% Display the Tglob table
disp('Table for Tglob:');
disp(demo_Tglob);

% Table for Eglob
demo_Eglob = table(participantGroups, participantNames, sex_reordered, age_reordered, education_reordered, beruf_reordered, site_reordered, Eglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Beruf', 'Site', 'Eglob'});
% Display the Eglob table
disp('Table for Eglob:');
disp(demo_Eglob);

% Table for Pglob
demo_Pglob = table(participantGroups, participantNames, sex_reordered, age_reordered, education_reordered, beruf_reordered, site_reordered, Pglob_all, ...
                  'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Beruf', 'Site', 'Pglob'});
% Display the Pglob table
disp('Table for Pglob:');
disp(demo_Pglob);


%% Linear Model Fitting and Residual Analysis for Global Network Measures by Group
% Fit the linear model for Pglob and calculate residuals
mdl_Pglob = fitlm(demo_Pglob(:, 3:6), demo_Pglob.Pglob);
residuals_Pglob = mdl_Pglob.Residuals.Raw;

% Separate residuals by group for Pglob
group_names = unique(demo_Pglob.Group);
residuals_Pglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = strcmp(demo_Pglob.Group, group_name);
    residuals_Pglob_by_group.(group_name) = residuals_Pglob(group_idx);
end

% Display residuals by group for Pglob
disp('Residuals by group for Pglob:');
disp(residuals_Pglob_by_group);

% Fit the linear model for Cglob and calculate residuals
mdl_Cglob = fitlm(demo_Cglob(:, 3:6), demo_Cglob.Cglob);
residuals_Cglob = mdl_Cglob.Residuals.Raw;

% Separate residuals by group for Cglob
residuals_Cglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = strcmp(demo_Cglob.Group, group_name);
    residuals_Cglob_by_group.(group_name) = residuals_Cglob(group_idx);
end

% Display residuals by group for Cglob
disp('Residuals by group for Cglob:');
disp(residuals_Cglob_by_group);

% Fit the linear model for Tglob and calculate residuals
mdl_Tglob = fitlm(demo_Tglob(:, 3:6), demo_Tglob.Tglob);
residuals_Tglob = mdl_Tglob.Residuals.Raw;

% Separate residuals by group for Tglob
residuals_Tglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = strcmp(demo_Tglob.Group, group_name);
    residuals_Tglob_by_group.(group_name) = residuals_Tglob(group_idx);
end

% Display residuals by group for Tglob
disp('Residuals by group for Tglob:');
disp(residuals_Tglob_by_group);

% Fit the linear model for Eglob and calculate residuals
mdl_Eglob = fitlm(demo_Eglob(:, 3:6), demo_Eglob.Eglob);
residuals_Eglob = mdl_Eglob.Residuals.Raw;

% Separate residuals by group for Eglob
residuals_Eglob_by_group = struct();

for i = 1:length(group_names)
    group_name = group_names{i};
    group_idx = strcmp(demo_Eglob.Group, group_name);
    residuals_Eglob_by_group.(group_name) = residuals_Eglob(group_idx);
end

% Display residuals by group for Eglob
disp('Residuals by group for Eglob:');
disp(residuals_Eglob_by_group);

%% Calculation and Group Analysis of Residuals for Local Network Measures Across Regions

% Assuming BCloc_all, Cloc_all, Eloc_all, ECloc_all, and Dloc_all are already calculated and available

% Define local measures and their names
local_measures_data = {BCloc_all, Cloc_all, Eloc_all, ECloc_all, Dloc_all};
local_measures_names = {'BCloc', 'Cloc', 'Eloc', 'ECloc', 'Dloc'};

% Identify participants with complete covariate data (sex, age, education, beruf, site)
valid_indices = ~cellfun(@(x) isequal(x, 'Unknown'), sex_reordered) & ...
                ~isnan(age_reordered) & ...
                ~cellfun(@(x) isequal(x, 'Unknown'), education_reordered) & ...
                ~cellfun(@(x) isequal(x, 'Unknown'), beruf_reordered) & ...
                ~cellfun(@(x) isequal(x, 'Unknown'), site_reordered);

% Loop through each local measure
for k = 1:length(local_measures_names)
    measure_name = local_measures_names{k};
    aucData = local_measures_data{k};

    % Initialize a structure to store residuals for each group and region
    residuals_by_group = struct();
    
    % Loop through each region
    numRegions = size(aucData, 2); % Number of regions
    for region = 1:numRegions
        % Extract data for the current region
        AUC_all = aucData(:, region);
        
        % Create a table for the current region, only including valid indices
        demo_region = table(participantGroups(valid_indices), participantNames(valid_indices), sex_reordered(valid_indices), age_reordered(valid_indices), ...
                            education_reordered(valid_indices), beruf_reordered(valid_indices), site_reordered(valid_indices), AUC_all(valid_indices), ...
                            'VariableNames', {'Group', 'ParticipantName', 'Sex', 'Age', 'Education', 'Beruf', 'Site', 'AUC'});
        
        % Fit the linear model for the current region and calculate residuals
        mdl_region = fitlm(demo_region(:, 3:7), demo_region.AUC); % Include Site
        residuals_region = mdl_region.Residuals.Raw;
        
        % Separate residuals by group for the current region
        for i = 1:length(group_names)
            group_name = group_names{i};
            group_idx = strcmp(demo_region.Group, group_name);
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

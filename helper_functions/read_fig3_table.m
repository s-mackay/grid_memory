% [results] = read_fig3_table(table_dir, sheet_name)
% read figure 3 data from .xlsx file.
% returns a struct with data, specifically 
% results.means    = means;
% results.sems     = sems;
% results.x        = x;
% results.pos_true = pos_true;
% results.neg_true = neg_true;
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function [results] = read_fig3_table(table_dir, sheet_name)

if nargin < 1
    table_dir = '.';
end

% % locate the spreadsheet
% fname = fullfile(table_dir, 'source_data_all_figures.xlsx');

% read content of spreadsheet as cell
data = readcell(table_dir, 'Sheet', sheet_name);

% Initialize an empty struct
results = struct();

% Get the number of rows in the data
nrows = size(data, 1);

% Initialize empty cell arrays and structs
means = cell(5, 3);
sems = cell(5, 3);
x = cell(5, 3);
pos_true = cell(5, 3);
neg_true = cell(5, 3);

% Loop over each row
for i = 1:nrows
    % Get the variable name from the first column
    varName = data{i, 1};
    if isempty(varName) || any(ismissing(varName))
        continue; % Skip empty variable names
    end
    
%      % Extract the data for the variable
    varData = data(i, 2:end);
    varData = varData(~cellfun('isempty', varData));
    varData = varData(~cellfun(@ismissing, varData));
%     varData = varData(~cellfun('isempty', varData));

    % Convert to numeric array if all elements are numeric
%     if iscell(varData) && all(cellfun(@isnumeric, varData))
%         varData = varData(~cellfun('isempty', varData));
%         varData = varData(~cellfun(@ismissing, varData));
%     end

    % Parse the variable name to determine where to place the data
    tokens = regexp(varName, '_(\d)_(\d)', 'tokens');
    if isempty(tokens)
        % Check for x variable which has a different format
        tokens_x = regexp(varName, 'x_(\d)', 'tokens');
        if ~isempty(tokens_x)
            reg = str2double(tokens_x{1}{1});
            % Assign the data to the x cell array for all columns
            for col = 1:3
                x{reg, col} = cell2mat(varData);
            end
        end
        continue;
    end

    reg = str2double(tokens{1}{1});
    col = str2double(tokens{1}{2});

    if contains(varName, 'means_rem')
        means{reg, col}(1, :) = cell2mat(varData);
    elseif contains(varName, 'means_forg')
        means{reg, col}(2, :) = cell2mat(varData);
    elseif contains(varName, 'SEMs_rem')
        sems{reg, col}(1, :) = cell2mat(varData);
    elseif contains(varName, 'SEMs_forg')
        sems{reg, col}(2, :) = cell2mat(varData);
    elseif contains(varName, 'pos_true_froms')
        pos_true{reg, col}.froms = cell2mat(varData);
    elseif contains(varName, 'pos_true_tos')
        pos_true{reg, col}.tos = cell2mat(varData);
    elseif contains(varName, 'pos_true_p')
        pos_true{reg, col}.cluster_p = cell2mat(varData);
    elseif contains(varName, 'neg_true_froms')
        neg_true{reg, col}.froms = cell2mat(varData);
    elseif contains(varName, 'neg_true_tos')
        neg_true{reg, col}.tos = cell2mat(varData);
    elseif contains(varName, 'neg_true_p')
        neg_true{reg, col}.cluster_p = cell2mat(varData);
    end
end

results = struct;
results.means    = means;
results.sems     = sems;
results.x        = x;
results.pos_true = pos_true;
results.neg_true = neg_true;


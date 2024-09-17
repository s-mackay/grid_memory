% function results_cell = read_fig2_results_from_xlsx(fname, sheetName)
% reads the data from the specified Excel file and worksheet.
%
%   Inputs:
%       fname (str): Name including path of the .xlsx document.
%       sheet_name (str): Name of the worksheet within fname.
%
%   Outputs:
%       results_cell: Cell array containing the data from the specified 
%       worksheet.
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function results_cell = read_fig2_results_from_xlsx(fname, sheet_name)

    % Read the content of the spreadsheet as a cell array
    data = readcell(fname, 'Sheet', sheet_name);

    % Initialize an empty cell array to store the results
    results_cell = cell(2,2); % Adjust the size based on your data
    
    % Initialize a temporary struct to store the current data
    temp_struct = struct();
    previous_field_name = 'placeholder';
    previous_rowIdx = 0;
    previous_colIdx = 0;
    temp_cell_array = {};
    
    % Loop over each row in the data
    for i = 1:size(data, 1)
        % Get the variable name and data
        varName = data{i, 1};
        varData = data(i, 2:end);
        varData = varData(~cellfun(@(x) any(ismissing(x)), varData));
        % Parse the variable name to get the field name and indices
        tokens = regexp(varName, '(\w+\d?)_(\d+)_(\d+)_(\d+)', 'tokens');
        tokens = tokens{1};
        fieldName = tokens{1};
        rowIdx = str2double(tokens{2});
        colIdx = str2double(tokens{3});
        % here we have several rows per variable (timestamps)
        if startsWith(fieldName, 'rel_ts')
            cell_idx = str2double(tokens{4});
            temp_cell_array{cell_idx} = cell2mat(varData);
            % if this is the last row or the next row starts with a new var
            if i < size(data, 1)
                next_varName = data{i+1, 1};
                next_tokens = regexp(next_varName,...
                    '(\w+\d?)_(\d+)_(\d+)_(\d+)', 'tokens');
                next_cell_idx = str2num(next_tokens{1}{4});
            end
            if i == size(data,1) || next_cell_idx == 1
                %add the cell array we've been collecting to temp_struct
                temp_struct.(fieldName) = temp_cell_array;
                temp_cell_array = cell(0);
            end
        % again several rows (paths to jpg files)
        elseif startsWith(fieldName, 'all_jpg')
            cell_idx = str2double(tokens{4});
            temp_cell_array{cell_idx} = varData{1};
            % if this is the last row or the next row starts with a new var
            if i < size(data, 1)
                next_varName = data{i+1, 1};
                next_tokens = regexp(next_varName,...
                    '(\w+)_(\d+)_(\d+)_(\d+)', 'tokens');
                next_cell_idx = str2num(next_tokens{1}{4});
            end
            if i == size(data,1) || next_cell_idx == 1
                %add the cell array we've been collecting to temp_struct
                temp_struct.(fieldName) = temp_cell_array;
                temp_cell_array = cell(0);
            end
        else
            % everything here just has 1 row of data
            % If the varData is numeric, assign it directly to the field
            if isempty(varData)
                temp_struct.(fieldName) = [];
            elseif ischar(varData{1})
                temp_struct.(fieldName) = varData{1};
            elseif isnumeric(varData{1})
                temp_struct.(fieldName) = cell2mat(varData);
            end
        end
        
        if i == size(data, 1)
            if iscell(temp_cell_array) && ~isempty(temp_cell_array)
                temp_struct.(fieldName) = temp_cell_array;
            end
%             results_cell{previous_rowIdx, previous_colIdx} = temp_struct;
        end
        
        % If this is the last field for the current struct, add it to the cell array
        if strcmp(fieldName, 'spikefile')% && i == size(data, 1)
            results_cell{str2double(tokens{2}), str2double(tokens{3})}...
                = temp_struct;
        end
        
    end
end
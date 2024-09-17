    
% function out_struct = read_figS6_results_from_xlsx(fname, sheetName)
% reads data necessary to create figure S6 from .xlsx. All required
% variables are fields of out_struct
%   Inputs:
%       fname (str): Name including path of the .xlsx document.
%       sheetName (str): Name of the worksheet within fname.
%
%   Output:
%       out_struct (struct): Structure containing the read data.
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function out_struct = read_figS6_results_from_xlsx(fname, sheetName)

    % Read the content of the spreadsheet as a cell array
    data = readcell(fname, 'Sheet', sheetName);

    % Initialize an empty cell array to store the results
    brainregs = {'Amy', 'Hipp', 'EC', 'PHC'};
    types = {'item', 'loc'};
    out_struct = struct;
    out_struct.hists = cell(2, 4); % 2(item / loc) x 4 ( brain regions)
    out_struct.binCenters = cell(2, 4); 
    out_struct.measuredValue = nan(2, 4);
    out_struct.empSize = nan(2, 4);
    
    % Loop over each row in the data
    for i = 1:size(data, 1)
        % Get the variable name and data
        varName   = data{i, 1};
        varData   = data(i, 2:end);
        varData   = varData(~cellfun(@(x) any(ismissing(x)), varData));
        % Parse the variable name to get the field name and indices
        tokens    = regexp(varName, '(\w+)_(\w+)_(\w+)', 'tokens');
        tokens    = tokens{1};
        reg_ind   = strcmp(tokens{3}, brainregs);
        type_ind  = strcmp(tokens{2}, types);
        fieldName = tokens{1};

        % we either get arrays...
        if strcmp(fieldName, 'hists') || strcmp(fieldName, 'binCenters')
            out_struct.(fieldName){type_ind, reg_ind} = cell2mat(varData);
        % ... or doubles
        else % this leaves measuredValue and empSize
            out_struct.(fieldName)(type_ind, reg_ind) = cell2mat(varData);

        end
    end
end
function [spiketable_item, spiketable_loc] = read_Fig2_excel(data_dir)

% since data has to be provided in spreadsheet format, this is a function
% to convert it back to a matlab table.

% initialize more columns to be added to final tables
more_cols = cell(2,6);

%% item cells
for unit = 1:2
    fname = fullfile(data_dir, sprintf('Fig2_item_unit_%i.xlsx', unit));
    sheets = sheetnames(fname);
    % first sheet is called 'info'
    info = readtable(fname, 'Sheet', 'info');
    spiketable_item(unit, 1:5) = info;
    % add following sheets as columns to spiketablevariable
    for i = 2:length(sheets)
        sheet = sheets{i};
        % filenames have to be read as strings, not numbers, into a cell
        if strcmp(sheet, 'fname')
            data = readcell(fname, 'Sheet', sheet, 'NumHeaderLines', 0);
        else
            data = readmatrix(fname, 'Sheet', sheet, 'NumHeaderLines', 0);
        end
        % for rel_ts, i.e. spike times relative to stim onset, the data is
        % converted to a cell array
        if size(data, 2) > 1 && ~strcmp(sheet, 'spikeshapes')
            cellArray = cellfun(@(row) row(~isnan(row)),...
                num2cell(data, 2), 'UniformOutput', false)';
            more_cols{unit, i-1} = cellArray;
        else
            data_array = data';
%             spiektable_item.(sheet) = data_array;
            more_cols{unit, i-1} = data_array;
        end
    end
    %if the last trials had no time stamps, we have to add empty brackets
    %to fill the trials
    if size(more_cols{unit, 1}, 1) < size(more_cols{unit, 2}, 2)
        more_cols{unit, 1}{1, size(more_cols{unit, 2}, 2)} = [];
    end
end
for col = 1:size(more_cols, 2)
    new_var = more_cols(:, col);
    new_name = sheets(col+1);
    spiketable_item = addvars(spiketable_item, new_var, 'NewVariableNames', new_name);
    pause(.2)
end

%% location cells
% here we have less columns because there are no filenames to images
more_cols = cell(2, 5); 
for unit = 1:2
    fname = fullfile(data_dir, sprintf('Fig2_loc_unit_%i.xlsx', unit));
    sheets = sheetnames(fname);
    % first sheet is called 'info'
    info = readtable(fname, 'Sheet', 'info');
    spiketable_loc(unit, 1:5) = info;
    % add following sheets as columns to spiketablevariable
    for i = 2:length(sheets)
        sheet = sheets{i};
        data = readmatrix(fname, 'Sheet', sheet, 'NumHeaderLines', 0);
        % for rel_ts, i.e. spike times relative to stim onset, the data is
        % converted to a cell array
        if size(data, 2) > 1 && ~strcmp(sheet, 'spikeshapes')
            cellArray = cellfun(@(row) row(~isnan(row)),...
                num2cell(data, 2), 'UniformOutput', false)';
            more_cols{unit, i-1} = cellArray;
        else
            data_array = data';
%             spiektable_item.(sheet) = data_array;
            more_cols{unit, i-1} = data_array;
        end
    end
    %if the last trials had no time stamps, we have to add empty brackets
    %to fill the trials
    if size(more_cols{unit, 1}, 1) < size(more_cols{unit, 2}, 2)
        more_cols{unit, 1}{1, size(more_cols{unit, 2}, 2)} = [];
    end
end
for col = 1:size(more_cols, 2)
    new_var = more_cols(:, col);
    new_name = sheets(col+1);
    spiketable_loc = addvars(spiketable_loc, new_var, 'NewVariableNames', new_name);
    pause(.2)
end

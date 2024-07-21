function write_fig2_results_to_xlsx(results_cell, file_name, sheet_name)

% function write_fig2_results_to_xlsx(results_cell, file_name, sheet_name)
% writes the data necessary to plot figure 2 to .xlsx
%
%   Inputs:
%       results_cell: Cell array containing the data to be written.
%       file_name (str): Name including path of the .xlsx document.
%       sheet_name (str): Name of the worksheet within file_name.

    data = {};
    [rows, cols] = size(results_cell);
    
    for r = 1:rows
        for c = 1:cols
            result = results_cell{r,c};
            fields = fieldnames(result);
            
            for i = 1:numel(fields)
                var_name = sprintf('%s_%d_%d', fields{i}, r, c);
                var_data = result.(fields{i});
                
                if iscell(var_data)
                    % Serialize each cell element
                    for j = 1:numel(var_data)
                        row = cell(1,2);
                        cell_var_name = sprintf('%s_%d', var_name, j);
                        cell_var_data = var_data{j};
                        row{1} = cell_var_name;
                        if ~isempty(cell_var_data)
                            row{2} = cell_var_data;
                        end
                        data = [data; row];
                    end
                elseif isnumeric(var_data)
                    for j = 1:size(var_data, 1)
                        row_var_name = sprintf('%s_%d', var_name, j);
                        row_var_data = var_data(j,:);
                        row = [{row_var_name}, row_var_data];
                        data = [data; row];
                    end
                elseif isstr(var_data)
                    var_name = sprintf('%s_%d', var_name, 1);
                    row = [{var_name}, {var_data}];
                    data = [data; row];
                else
                    var_name = sprintf('%s_%d', var_name, 1);
                    row = [{var_name}, num2cell(var_data)];
                    data = [data; row];
                end
            end
        end
    end
    %delete the old sheet so that there is no data remaining. Or rather, since
    %deleting is turning out to be tricky (if you dont have Excel installed), 
    %we're doing the ugly workaround of writing 1000 x 500 empty cells into
    %the sheet
    empty_data = cell(3000, 1000);
    % Write the empty cell array to the specified sheet
    writecell(empty_data, file_name, 'Sheet', sheet_name, 'Range', 'A1');
    writecell(data, file_name, 'Sheet', sheet_name, 'Range', 'A1');
end

% % Example usage
% resultsCell = {
%     struct('stim_type', 1, 'pref_jpg', 'image.jpg', 'x', [1, 2, 3], 'rel_ts_rem', {[1, 2], [3, 4]}, 'rel_ts_forg', [1, 2, 3], 'pref_rem', 0.5, 'pref_forg', 0.6, 'non_pref_rem', 0.7, 'non_pref_forg', 0.8, 'spikeshapes', [1, 2; 3, 4]), 
%     struct('stim_type', 2, 'pref_jpg', 'image2.jpg', 'x', [4, 5, 6], 'rel_ts_rem', {[5, 6], [7, 8]}, 'rel_ts_forg', [4, 5, 6], 'pref_rem', 0.55, 'pref_forg', 0.65, 'non_pref_rem', 0.75, 'non_pref_forg', 0.85, 'spikeshapes', [5, 6; 7, 8])
% };
% writeResultsToExcel(resultsCell, 'source_data_all_figures.xlsx', 'figure_2');
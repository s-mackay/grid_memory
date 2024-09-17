% this function will call figure2_subplanels to create all 4 subplots of
% figure2.
%
% figure2(data_dir, read_from_excel)
% INPUTS:
%   data_dir (str):         Path to source data (default: 'source_data').
%
%   read_from_excel (bool): Flag to read from Excel (true) 
%                           or .mat files (false) (default: false).
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function figureS1(data_directory, read_from_excel)

if nargin < 1
    data_directory  = 'source_data';
end
if nargin < 2
    read_from_excel = false;
end

% location of .xlsx file
all_figure_xlsx_fname = sprintf('%s/source_data_all_figures.xlsx',...
    data_directory);

% read or prepare to write .xlsx data
if read_from_excel
    results_cell = read_fig2_results_from_xlsx(all_figure_xlsx_fname,...
        'figure_S1');
else
    results_cell = cell(2, 2); %will be completed along the way
    load(sprintf('%s/Fig2_data.mat', data_directory), 'spiketable_item',...
        'spiketable_loc');
end

%item neurons
for id = [1, 2]
    if read_from_excel
        figureS1_item(data_directory, id, results_cell{1, id}, true);
    else
        results = figureS1_item(data_directory, id, spiketable_item);
        results_cell{1, id} = results;
    end
end
%location neurons
for id = [1, 2]
    if read_from_excel
        figureS1_loc(data_directory, id, results_cell{2, id}, true);
    else
        results = figureS1_loc(data_directory, id, spiketable_loc);
        results_cell{2, id} = results;
    end
end
% this is how data is written to source_data_all_figures.xlsx
% if ~read_from_excel
%     write_fig2_results_to_xlsx(results_cell, all_figure_xlsx_fname, 'figure_S1');
% end
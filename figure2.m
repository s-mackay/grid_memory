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

function figure2(data_dir, read_from_excel)

if nargin < 1
    data_dir  = 'source_data';
end
if nargin < 2
    read_from_excel = false;
end

all_figure_xlsx_fname = sprintf('%s/source_data_all_figures.xlsx',...
        data_dir);

if read_from_excel
    % load excel data
    results_cell = read_fig2_results_from_xlsx(all_figure_xlsx_fname,...
        'figure_2');
    for item_or_loc = [1, 2]
        for id = [1, 2]
            figure2_from_xlsx_subpanels(results_cell{item_or_loc, id},...
                item_or_loc, id);
        end
    end
else
    load(sprintf('%s/Fig2_data.mat', data_dir), 'spiketable_item',...
        'spiketable_loc');
    spiketables = {spiketable_item, spiketable_loc};
    results_cell = cell(2, 2);
    for item_or_loc = 1:2
        for id = 1:2
            results_cell{item_or_loc, id} = ...
                figure2_subpanels(data_dir, id, item_or_loc, ...
                spiketables{item_or_loc});
        end
    end
%     % this is how data is written to source_data_all_figures.xlsx
%     write_fig2_results_to_xlsx(results_cell, all_figure_xlsx_fname, 'figure_2');
end
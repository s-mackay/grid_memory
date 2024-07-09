
function all_figures_and_stats(data_directory, read_from_excel)

%% function all_figures_and_stats(data_directory, read_from_excel)
% generates all figures and stats
% 
% parameters:
% data_directory (str): e.g. '/path/to/git_dir/grid_memory'
% read_from_excel (bool): false (read .mat files, recommended)
%                         true (read .xlsx files, not recommended)

if nargin < 1
    git_dir = input(['Please enter the path to the git repo, in which the',...
        'source_data directory is located. E.g. \n /path/to/dir', 's']);
    %remove any single or double quotes the user may have added 
    git_dir = strrep(git_dir, '''', '');
    git_dir = strrep(git_dir, '"', '');
    data_directory = [git_dir, '/source_data'];
end
if nargin < 2
    read_from_excel = false;
end

% add dir 'source_data' to data path if not already included
if ~contains(data_directory, 'source_data')
    data_directory = fullfile(data_directory, 'source_data');
end

%% Figure 2
% load excel data
all_figure_xlsx_fname = sprintf('%s/source_data_all_figures.xlsx',...
    data_directory);
if read_from_excel
    results_cell = read_fig2_results_from_xlsx(all_figure_xlsx_fname,...
        'figure_2');
else
    load(sprintf('%s/Fig2_data.mat', data_directory), 'spiketable_item',...
        'spiketable_loc');
end

% generate all sub-panels

if read_from_excel
    for item_or_loc = [1, 2]
        for id = [1, 2]
            figure2_from_xlsx(results_cell{item_or_loc, id},...
                item_or_loc, id);
        end
    end
else
    spiketables = {spiketable_item, spiketable_loc};
    results_cell = cell(2, 2);
    for item_or_loc = [1, 2]
        for id = [1, 2]
            results_cell{item_or_loc, id} = ...
                figure2(data_directory, id, item_or_loc, ...
                spiketables{item_or_loc});
        end
    end
%     % this is how data is written to source_data_all_figures.xlsx
%     write_fig2_results_to_xlsx(results_cell, all_figure_xlsx_fname, 'figure_2');
end

%% Figure 3
alpha_resp = .001;  plot_visibility = 'on'; narrow_ylims = false;
seed = 12;          no_overlaps = 0;        which_half = 0;
binned_test = true; 

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, read_from_excel, data_directory);


%% Linear Mixed Effects Models
do_plots = false;
gm_trialLength_RT_corrFalse_mixedEffectsModel(data_directory, ...
    read_from_excel, do_plots)

%% various stats stated in the text
% binomial tests, fractions of stimuli responded to, anterior vs. posterior
% hippocampal recording sites
show_plots = false;
stats_grid_memory(data_directory, read_from_excel, show_plots)

%% counting/distractor task
gm_counting_consolidation(data_directory, read_from_excel);

%% Figure S1 (complete raster plots)
% generate all sub-panels
if read_from_excel
    results_cell = read_fig2_results_from_xlsx(all_figure_xlsx_fname,...
        'figure_S1');
else
    results_cell = cell(2, 2); %will be completed along the way
end
all_figure_xlsx_fname = sprintf('%s/source_data_all_figures.xlsx',...
    data_directory);
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

%% Figure S2 (selectivity)
figureS2(data_directory, read_from_excel)

%% Figure S3 (Figure 3 with adjusted ylims for 3rd column)
alpha_resp = .001;  plot_visibility = 'on'; narrow_ylims = true;
seed = 12;          no_overlaps = false;    which_half = 0;
binned_test = true; %read_from_excel = false;

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, false, data_directory);

%% Figure S4 (Figure 3 but excluding neurons with both response types)
alpha_resp = .001;  plot_visibility = 'on'; narrow_ylims = false;
seed = 12;          no_overlaps = 2;    which_half = 0;
binned_test = true; 

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, read_from_excel, data_directory);

%% Figure S5 (retrieval)
figureS5(data_directory, read_from_excel)

%% Figure S6 (empirical size)
figureS6(data_directory, read_from_excel)

end

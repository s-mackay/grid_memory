
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
    %% do you want to read .mat files or .xlsx?
    % This (read_from_excel=true) will not work for figure 3, because the
    % source files were too large in .xlsx format
    read_from_excel = ...
        logical(input(['Would you like to read .xlsx files instead of .mat?',...
            '\nThis will be ignored for Figure 3 since file sizes were',...
            '\ntoo large here, and is generally very slow.\n',...
            '1 = read .xlsx (not recommended)\n',...
            '0 = read .mat (recommended)\n>']));
    fprintf('read_from_excel = %s\n', string(read_from_excel));
end

% add dir 'source_data' to data path if not already included
if ~contains(data_directory, 'source_data')
    data_directory = fullfile(data_directory, 'source_data');
end

%% Figure 2
% load excel data
if read_from_excel
    [spiketable_item, spiketable_loc] = read_Fig2_excel(data_directory);
else
    load(sprintf('%s/Fig2_data.mat', data_directory), 'spiketable_item',...
        'spiketable_loc');
end

% generate all sub-panels
spiketables = {spiketable_item, spiketable_loc};
for item_or_loc = [1, 2]
    for id = [1, 2]
        figure2(data_directory, id, item_or_loc, spiketables{item_or_loc});
    end
end

%% Figure 3
alpha_resp = .001;  plot_visibility = 'on'; narrow_ylims = false;
seed = 12;          no_overlaps = 0;        which_half = 0;
binned_test = true; %read_from_excel = false;

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, read_from_excel, data_directory);


%% Linear Mixed Effects Models
gm_trialLength_RT_corrFalse_mixedEffectsModel(data_directory, read_from_excel)

%% various stats stated in the text
% binomial tests, fractions of stimuli responded to, anterior vs. posterior
% hippocampal recording sites
stats_grid_memory(data_directory, read_from_excel)

%% counting/distractor task
gm_counting_consolidation(data_directory, read_from_excel);

%% Figure S1 (complete raster plots)
% generate all sub-panels
%item neurons
for id = [1, 2]
    figureS1_item(data_directory, id, spiketable_item);
end
%location neurons
for id = [1, 2]
    figureS1_loc(data_directory, id, spiketable_loc);
end

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
binned_test = true; %read_from_excel = false;

figure3(alpha_resp, plot_visibility, narrow_ylims, seed, no_overlaps,...
    which_half, binned_test, read_from_excel, data_directory);

%% Figure S5 (retrieval)
figureS5(data_directory, read_from_excel)

%% Figure S6 (empirical size)
% have to set the 2nd parameter, read_from_excel, to false because the
% .xlsx source file was too large to be included in the git repo
response_permutation(data_directory, false)

end

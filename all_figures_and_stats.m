%% function all_figures_and_stats(git_directory, read_from_excel)
% generates all figures and stats
% 
% parameters:
% git_directory    (str): e.g. '/path/to/git_dir/grid_memory'
% read_from_excel (bool): false (read .mat files, recommended)
%                         true (read .xlsx files, not recommended)
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function all_figures_and_stats(git_directory, read_from_excel)

if nargin < 1
    git_directory = input(['Please enter the path to the git repo, in which',...
        ' the source_data directory is located. E.g. \n /path/to/dir\n',...
        'If you are already in the git directory, hit enter\n'], 's');
    %remove any single or double quotes the user may have added 
    git_directory = strrep(git_directory, '''', '');
    git_directory = strrep(git_directory, '"', '');
    
    if isempty(git_directory)
        git_directory = pwd;
    end
end

% add dir 'source_data' to git_directory to obtain data_directory
if contains(git_directory, 'source_data')
    git_directory = strrep(git_directory, '/source_data', '');    
end
data_directory = fullfile(git_directory, 'source_data');

addpath(genpath(git_directory));

if nargin < 2
    read_from_excel = false;
end

%% Figure 2
figure2(data_directory, read_from_excel);

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
figureS1(data_directory, true)

%% Figure S2 (selectivity)
figureS2(data_directory, read_from_excel)

%% Figure S3 (Figure 3 with adjusted ylims for 3rd column)
figureS3(data_directory, read_from_excel)

%% Figure S4 (Figure 3 but excluding neurons with both response types)
figureS4(data_directory, read_from_excel)

%% Figure S5 (retrieval)
figureS5(data_directory, read_from_excel)

%% Figure S6 (empirical size)
figureS6(data_directory, read_from_excel)

end

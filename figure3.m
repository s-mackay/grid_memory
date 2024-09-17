% This function plots figure 3 of Mackay et al. 2024 
% (DOI:10.1038/s41467-024-52295-5)
%
% function figure3(alpha_resp, plot_visibility, ...
%     narrow_ylims, seed, no_overlaps, which_half, binned_test,...
%     read_from_excel, data_dir)
%
% INPUTS:
%   alpha_resp (numeric):   Threshold for detecting neurons (default: 0.001).
%   plot_visibility (char): Visibility of plot ('on' or 'off') (default:
%                           'on')
%   narrow_ylims (bool):    zoom in on column 3? (default: false).
%   seed (numeric):         Random seed for reproducibility (default: 12).
%   no_overlaps (numeric):  Handling of units showing both response types,
%                           item and location (0: allow overlaps, 
%                           1: keep stronger response, 
%                           2: exclude dual responses) (default: 0).
%   which_half (numeric):   Subset of trials to include (0: all trials, 
%                           1: first half, 2: second half) (default: 0).
%   binned_test (numeric):  Type of statistical test (1: bw ranksum, 
%                           2: ranksum) (default: 1).
%   read_from_excel (bool): Flag to read from Excel (true) 
%                           or .mat files (false) (default: false).
%   data_dir (str):         Path to source data (default: 'source_data').
%
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function figure3(alpha_resp, plot_visibility, ...
    narrow_ylims, seed, no_overlaps, which_half, binned_test,...
    read_from_excel, data_dir)

if nargin < 9; data_dir ='source_data';end % assuming your working 
                                           %  directory is the git repo
if nargin < 8; read_from_excel = false;end % if set to false, the .mat
                                           %  files will be read.
                                           %  Recommended (much faster).
if nargin < 7; binned_test = 1;        end % 1= bw ranksum, 2 = ranksum
if nargin < 6; plot_visibility = 'on'; end
if nargin < 5; narrow_ylims = false;   end % tailored to last column
if nargin < 4; seed = 12;              end % selected arbitrarily, fixed to
                                           %  allow for exact reproduction
if nargin < 3; no_overlaps = 0;        end % 0: overlaps are allowed
                                           % 1: if a unit shows both resps,
                                           %    only the stronger one will 
                                           %    be included. 
                                           % 2: exclude units that respond
                                           %    to items and locations
if nargin < 2; which_half = 0;         end % 0: all trials
                                           % 1: first half
                                           % 2: second half
if nargin < 1; alpha_resp = 0.001;     end % separates responsive from 
                                           %  non-responsive  
                                                                                 
% set seed                                           
randn('state', seed);
rand ('state', seed);

% response window 0-1000ms from stimulus onset
resp_win = [0, 1000];
% resposne window is the same for all 4 brain regions
resp_wins = repmat(resp_win, 4, 1);

if binned_test % binwise ranksum test
    do_binned = true;
    testname = 'binwiseRS'; % binwise ranksum
else
    do_binned = false;
    testname = 'NoBinsRS'; % ranksum test with no bins (just 1 window)
end

if no_overlaps == 1 % if both response types, count only the stronger one
    overl_str = '_ovrlpStronger';
elseif no_overlaps == 2 % if both resp types, exclude the unit
    overl_str = '_ovrlpExcl';
elseif no_overlaps == 0 % allow overlaps, count both response types
   overl_str = '';
end  

if narrow_ylims
    ylims = [-.05, .1];
    yZoomStr = '_yZoom'; % will be added to figure title
else
    ylims = [-.5 4];
    yZoomStr = '';
end
% when looking at only half of the trials, the figures will have
% appropriate filenames
half_strs = {'', '_1stHalf', '_2ndHalf'};
half_str = half_strs{which_half+1};

% what kind of figure is it
if no_overlaps==0 && narrow_ylims==0 && which_half==0 && binned_test==1
    load_figure = '3'; 
    plotname = '3';
elseif no_overlaps==2 && narrow_ylims==0 && which_half==0 && binned_test==1
    load_figure = 'S4'; 
    plotname  = 'S4';
elseif no_overlaps==0 && narrow_ylims==1 && which_half==0 && binned_test==1
    % here we use the same data as in figure 3, but save as figure S3
    % because scaling is different
    load_figure = '3';
    plotname = 'S3';
else
    load_figure = sprintf('3_variation_with%s%s%s_%s', overl_str,...
        yZoomStr, half_str, testname);
    plotname = sprintf('figure%s', load_figure);
end

% for loading or saving source data from/to excel
% if ~narrow_ylims %because then we can just use figure3 data 
    sheet_name = sprintf('figure_%s', load_figure);
% end

% we define only one alpha level to separate all cells into responsive or
% non-responsive
alpha_nonresp = alpha_resp;


% if you want little labels on the effects, stating effect size etc, set
% this to true
add_labels = false;

ttype_resps = 'enc';     % trial type that responses are based on
n_stims = 'allRespElic'; % all resposne eliciting
fs = 14;                 % font size for figure

source_data_fname = fullfile(data_dir, 'source_data_all_figures.xlsx');
if read_from_excel
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');
    results = read_fig3_table(source_data_fname, sheet_name);
else
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
    % collect summary data for excel source files, such as
    means = cell(5, 3);
    sems  = cell(5, 3); % standard errors of the means
    x     = cell(5, 3); % x axis in ms
    pos_true = cell(5, 3); % significant clusters where rememb.>forgotten
    neg_true  = cell(5, 3); % significant clsuters wehre forgotten>rememb.
end

% cd(data_dir)
nperms = 10000;
%nperms = 100;% for testing
% in permutation test: alpha level for detecting a significant
% difference at one single time point
alpha_diff = .05; 
% in permutation test: alpha level for detecting a significantly large
% cluster
alpha_cluster = .05;
% stepsize_ms = 5;
% kern_length = 5*kern_width; 
% in order to avoid periods in filenames, the strings i use in filenames
% are, e.g. 0010 for .001

if do_binned
    lowest_p_item = response_info.lowest_p_item;
    lowest_p_loc  = response_info.lowest_p_loc;
    % convolved firing rates from the following 3 files:
    % 1. convolved firing rates for item responses
    fname1 = 'convolved_spiking_zvals_enc_allRespElic_image_p001_encResps_seed12';
    % 2. convolved firing rates for location responses
    fname2 = 'convolved_spiking_zvals_enc_allRespElic_position_p001_encResps_seed12';
else
    respWinStr = sprintf('A%i-%i_H%i-%i_E%i-%i_P%i-%i', resp_wins(1,1),...
        resp_wins(1,2),resp_wins(2,1),resp_wins(2,2),resp_wins(3,1),...
        resp_wins(3,2),resp_wins(4,1),resp_wins(4,2));
    lowest_p_item = response_info.lowest_p_item_1bin;
    lowest_p_loc  = response_info.lowest_p_loc_1bin;
    % convolved firing rates from the following files:
    % 1. convolved firing rates for item responses
    fname1 = sprintf(...
        'convolved_spiking_zvals_enc_allRespElic_image_p001_encResps_ranksumNoBins_%s_seed12',...
        respWinStr);
    % 2. convolved firing rates for location responses
    fname2 = sprintf(...
        'convolved_spiking_zvals_enc_allRespElic_position_p001_encResps_ranksumNoBins_%s_seed12',...
        respWinStr);
end

% 1. convolved firing rates for non-responsive cells
fname3 = 'convolved_spiking_zvals_enc_all_allUnitsAllTrials_seed12';

% load .mat files:
conv_file1 = sprintf('%s/%s%s.mat', data_dir, fname1, half_str);
conv_file2 = sprintf('%s/%s%s.mat', data_dir, fname2, half_str);
conv_file3_c = sprintf('%s/%s%s_c.mat', data_dir, fname3, half_str);
conv_file3_f = sprintf('%s/%s%s_f.mat', data_dir, fname3, half_str);
conv1 = load(conv_file1,'y_unit_time_cf', 'x_all');
conv2 = load(conv_file2,'y_unit_time_cf', 'x_all');
conv3_c = load(conv_file3_c,'y_unit_time_c', 'x_all');
conv3_f = load(conv_file3_f,'y_unit_time_f', 'x_all');
y_unit_time_cf = cat(3, conv3_c.y_unit_time_c, conv3_f.y_unit_time_f);
conv3.y_unit_time_cf = y_unit_time_cf;
conv3.x_all = conv3_f.x_all;    

% these are the axis lims for the standard plots. ylims are overwritten
% later if narrow_ylims is true
xlims = [-500 1750];
% ylims = [-.5 4];
clrs = [.2, .2, .8; .8 0 .2]; % colors


rgns  = {'Amy', 'Hipp', 'EC', 'PHC'};
if alpha_resp == .001 && alpha_nonresp == alpha_resp
    p_str = sprintf('%.3f', alpha_resp);
    p_str_nonresp = sprintf('%.3f', alpha_nonresp);
else
    p_str = sprintf('%.5f', alpha_resp);
    p_str_nonresp = sprintf('%.5f', alpha_nonresp);    
end

reg_titles = {{'item neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {'location neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {sprintf('both %s \\geq %s', '{\it P}', p_str_nonresp)},...
    {'item neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {'location neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {sprintf('unresponsive units (%s \\geq %s)', '{\it P}', p_str_nonresp)},...
    {'item neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {'location neurons', sprintf('%s < %s', '{\it P}', p_str)},...
    {sprintf('unresponsive units (%s \\geq %s)', '{\it P}', p_str_nonresp)}};

set(groot,'defaultAxesColorOrder',[0,0,200; 200,0,0; 68,1,84]/255)

h1 = figure('visible', plot_visibility);
set(h1, 'Units', 'pixels', 'Position', [0 1920 1280 1024]);

x_final = [find(conv1.x_all >= xlims(1), 1), ...
           find(conv1.x_all >= xlims(2), 1)];
x_inds = x_final(1):x_final(2);
ncols = 3;
nrows = 5;

alpha_r_selection = alpha_resp;

% iterate over the 4 brain regions, Amygdala, Hippocampus, EC, PHC
for reg = 1:4
    
    %%%%%%%%%%%%%%%%% FIRST PLOT, LEFT COLUMN - ITEM RESPONSES %%%%%%%%%%%%
    curr_ax = subplot(nrows,ncols,(reg-1)*ncols+1); 
    % Responsive Units
    if no_overlaps == 1
        % each unit can only either be an item unit OR location unit 
        % (depending on which response is stronger), but not both. Units 
        % with bothresponse types are included, but only in one category
        incl_these = lowest_p_item < alpha_r_selection &...
                     lowest_p_item < lowest_p_loc & ...
                     response_info.brain_regs == reg;
    elseif no_overlaps == 2
        % only units that show a response in one of the categories (item or
        % location), but not both
        incl_these = lowest_p_item < alpha_r_selection &...
                     lowest_p_loc >= alpha_r_selection &...
                     response_info.brain_regs == reg;
    elseif no_overlaps == 0 % default
        % units can be included in both categories, item and location
        % responsive
        incl_these = lowest_p_item < alpha_r_selection &...
                     response_info.brain_regs == reg &...
                     ~isnan(conv1.y_unit_time_cf(:,1,1));
    end

    fprintf('\n%i item responsive neurons in %s\n', sum(incl_these), rgns{reg})
    set(gca, 'fontsize', fs)
    plot_info.ylims = ylims; plot_info.reg = reg; plot_info.stimType = 1; %1=img, 2=pos
    plot_info.add_labels = add_labels; plot_info.plot_colors = clrs;
    if sum(incl_these) >0 && ~narrow_ylims
        if read_from_excel
            handles = plot_convolved(results, [reg, 1], plot_info, curr_ax);
        else
            [handles, results] = perm_multcompare_gm_publ(...
                conv1.y_unit_time_cf(incl_these,x_inds,1),...
                conv1.y_unit_time_cf(incl_these,x_inds,2),...
                nperms, alpha_cluster, alpha_diff,...
                conv1.x_all(x_inds), [], curr_ax,7, fs, plot_info);
            means{reg, 1}    = results.means; %mean signals rem., forg.
            sems{reg, 1}     = results.sem; %standard errors of the means
            x{reg, 1}        = conv1.x_all(x_inds);
            pos_true{reg, 1} = results.pos_true; % significant clusters 
            neg_true{reg, 1}  = results.neg_true; % significant clsuters 
        end
    else
        ylim(ylims);xlim([-300 1300])
        handles = gca;
        % fill summary stats with empty brackets
        [means{reg, 2}, sems{reg, 2}, x{reg, 2}] = deal([]);
        [pos_true{reg, 2}.froms, pos_true{reg, 2}.tos,...
            neg_true{reg, 2}.froms, neg_true{reg, 2}.tos] = deal([]);
    end

    if reg==1
        title(handles(1), reg_titles{(reg-1)*3+1},...
            'interpreter', 'tex')
    end
    text(handles(1), xlims(2), mean(ylims)+.35*diff(ylims),...
        sprintf('%s = %i / %i', '{\it n}', sum(incl_these),...
        sum(response_info.brain_regs==reg)), 'horizontalalignment', ...
        'right','interpreter', 'tex', 'fontsize', fs+2)
    ylim(handles(1), ylims)
    set(handles(:), 'fontsize', fs)

    if length(handles) > 3; set(handles(4), 'fontsize', 18); end
    set(h1, 'CurrentAxes', handles(1))
    xl = xlim;
    if narrow_ylims
        textXPos = xl(1)-diff(xl)/2.5;
    else
        textXPos = xl(1)-diff(xl)/4;
    end
    text(textXPos, mean(ylims), rgns{reg}, 'fontsize', fs+2,...
        'horizontalalignment', 'right', 'verticalalignment', 'middle');
    set(gca, 'fontsize', fs)
    if reg<4; set(gca,'xticklabel',{[]}); end
    xticks(-1000:500:2000)
    
    %%%%%%%%%%%%%%%% NEXT PLOT, CENTER COLUMN - LOC RESPS %%%%%%%%%%%%%%%%%

    subplot(nrows,ncols,(reg-1)*ncols+2); 
    % responsive units
    if no_overlaps == 1
        incl_these = lowest_p_loc < alpha_r_selection &...
                     lowest_p_loc <= lowest_p_item & ...
                     response_info.brain_regs == reg;
    elseif no_overlaps == 2
        incl_these = lowest_p_loc < alpha_r_selection &...
                     lowest_p_item >= alpha_r_selection &...
                     response_info.brain_regs == reg;
    elseif no_overlaps == 0
        incl_these_ = lowest_p_loc < alpha_r_selection &...
                      response_info.brain_regs == reg & ...
                      ~isnan(conv2.y_unit_time_cf(:,1,1));
        incl_these = incl_these_ & ~isnan(conv2.y_unit_time_cf(:,1,1));
    end  
    
    set(gca, 'fontsize', fs)
    plot_info.ylims = ylims; plot_info.reg = reg; 
    plot_info.stimType = 2; %1=img, 2=pos
    plot_info.add_labels = add_labels;
    plot_info.plot_colors = clrs;
    if sum(incl_these) >0 && ~narrow_ylims
        if read_from_excel
            handles = plot_convolved(results, [reg, 2], plot_info);
        else
            [handles, results] = perm_multcompare_gm_publ(...
            conv2.y_unit_time_cf(incl_these,x_inds,1),...
            conv2.y_unit_time_cf(incl_these,x_inds,2), ...
            nperms, alpha_cluster, alpha_diff, conv2.x_all(x_inds),...
            [], gca,7, fs, plot_info);
            means{reg, 2}    = results.means; %mean signals rem., forg.
            sems{reg, 2}     = results.sem; %standard errors of the means
            x{reg, 2}        = conv1.x_all(x_inds);
            pos_true{reg, 2} = results.pos_true; % significant clusters 
            neg_true{reg, 2} = results.neg_true; % significant clsuters 
        end
    else
        ylim(ylims);xlim([-300 1700])
        handles = gca;
        % fill summary stats with empty brackets
        [means{reg, 2}, sems{reg, 2}, x{reg, 2}] = deal([]);
        [pos_true{reg, 2}.froms, pos_true{reg, 2}.tos,...
            pos_true{reg, 2}.cluster_p, neg_true{reg, 2}.froms,...
            neg_true{reg, 2}.tos, neg_true{reg, 2}.cluster_p] = deal([]);
    end
    
    fprintf('%i loc responsive neurons in %s\n', sum(incl_these), rgns{reg})
    if reg==1
        title(handles(1), reg_titles{(reg-1)*3+2},...
        'interpreter', 'tex')
    end
    if exist('incl_these_', 'var')
        incl_these = incl_these_;
    end
    text(handles(1),xlims(2), mean(ylims)+.35*diff(ylims), ...
        sprintf('%s = %i / %i', '{\it n}', sum(incl_these),...
        sum(response_info.brain_regs==reg)),...
        'horizontalalignment', 'right','interpreter', 'tex', 'fontsize', fs+2)
    ylim(handles(1), ylims)
    set(handles(:), 'fontsize', fs)
    %handles1-5:?, colorbar, main legend
    if length(handles) > 3; set(handles(4), 'fontsize', 18); end
    set(h1, 'CurrentAxes', handles(1))
    set(gca, 'fontsize', fs)
    if reg<4; set(gca,'xticklabel',{[]}); end
    set(gca,'yticklabel',{[]})
    xticks(-1000:500:2000)
    ylabel('')
    
    %%%%%%%%%%%%%%%% NEXT PLOT, RIGHT COLUMN> UNRESPONSIVE %%%%%%%%%%%%%%%

    subplot(nrows,ncols,(reg-1)*ncols+3);
    incl_these = lowest_p_item >= alpha_nonresp & ...
                 lowest_p_loc  >= alpha_nonresp & ...
                 response_info.brain_regs == reg;
    data1 = conv3.y_unit_time_cf(incl_these, x_inds, 1);
    data2 = conv3.y_unit_time_cf(incl_these, x_inds, 2);
    fprintf('%i non-responsive neurons in %s\n', sum(incl_these), rgns{reg})

    %UNresponsive units
    if reg==4
        leg_strs = {'remembered', 'forgotten'};
    else
        leg_strs = [];
    end
    plot_info.ylims = ylims; plot_info.reg = reg; plot_info.stimType = 3; %1=item, 2=pos 3=unresp
    plot_info.add_labels = add_labels; plot_info.plot_colors = clrs;
    if read_from_excel
            handles = plot_convolved(results, [reg, 3], plot_info);
    else
        [handles, results] = perm_multcompare_gm_publ(data1, data2, nperms,...
            alpha_cluster, alpha_diff, conv3.x_all(x_inds),...
            leg_strs, gca, 7, fs, plot_info);
            means{reg, 3}    = results.means; %mean signals rem., forg.
            sems{reg, 3}     = results.sem; %standard errors of the means
            x{reg, 3}        = conv1.x_all(x_inds);
            pos_true{reg, 3} = results.pos_true; % significant clusters 
            neg_true{reg, 3} = results.neg_true; % significant clsuters 
    end

    if reg==1
        title(handles(1), reg_titles{(reg-1)*3+3},...
        'interpreter', 'tex')
    end
    text(handles(1), xlims(2), mean(ylims)+.35*diff(ylims),...
        sprintf('%s = %i / %i', '{\it n}', sum(incl_these),...
        sum(response_info.brain_regs==reg)),...
        'horizontalalignment', 'right','interpreter', 'tex', 'fontsize', fs+2)
    ylim(handles(1), ylims)
    set(handles(:), 'fontsize', fs)
    xticks(-1000:500:2000)
    if reg<4; set(gca,'xticklabel',{[]}); end
%     ylabel('')
    if ~narrow_ylims
        set(gca,'yticklabel',{[]})
        ylabel('')
    end
    
end

%last plot center bottom, aligned to confirmation tap
subplot(nrows, ncols, 14);
xlims = [-1000 1000];
cond_  = 'encTap';
if ~read_from_excel 
    conv_file2 = sprintf('%s/convolved_spiking_zvals_%s_%s_position_p%.3i_%sResps_seed%i%s.mat',...
        data_dir, cond_, n_stims, alpha_r_selection*1000, ttype_resps, seed, half_str);
    conv2 = load(conv_file2,'y_unit_time_cf', 'x_all');
end

x_final = [find(conv2.x_all >= xlims(1), 1), ...
       find(conv2.x_all >= xlims(2), 1)];
x_inds = x_final(1):x_final(2);
if no_overlaps ==1  %reg is still 4
    incl_these = lowest_p_loc < alpha_r_selection & ...
                 lowest_p_loc < lowest_p_item & ...
                 response_info.brain_regs == reg;
elseif no_overlaps ==2
    incl_these = lowest_p_loc < alpha_r_selection & ...
                 lowest_p_item >= alpha_r_selection &...
                 response_info.brain_regs == reg;
elseif no_overlaps ==0
    incl_these = lowest_p_loc < alpha_r_selection &...
                 response_info.brain_regs == reg;
end
fprintf('%i loc responsive neurons in %s (aligned to tap)\n',...
    sum(incl_these), rgns{reg})
set(gca, 'fontsize', fs)
plot_info.ylims = ylims; plot_info.reg = reg; plot_info.stimType = 2; %1=img, 2=pos
plot_info.add_labels = add_labels; plot_info.plot_colors = clrs;
if sum(incl_these) >0 && ~narrow_ylims
    if read_from_excel
            handles = plot_convolved(results, [5, 2], plot_info);
    else
        [handles, results] = perm_multcompare_gm_publ(...
            conv2.y_unit_time_cf(incl_these,x_inds,1),...
            conv2.y_unit_time_cf(incl_these,x_inds,2),...
            nperms, alpha_cluster, alpha_diff, conv2.x_all(x_inds), ...
            [], gca,7, fs, plot_info);
        means{5, 2}    = results.means; %mean signals rem., forg.
        sems{5, 2}     = results.sem; %standard errors of the means
        x{5, 2}        = conv2.x_all(x_inds);
        pos_true{5, 2} = results.pos_true; % significant clusters 
        neg_true{5, 2}  = results.neg_true; % significant clsuters
    end
else
    ylim(ylims);xlim([-300 1700])
    handles = gca;
end
text(handles(1),xlims(2), mean(ylims)+.35*diff(ylims), ...
    sprintf('%s = %i / %i', '{\it n}', sum(incl_these), ...
    sum(response_info.brain_regs==reg)),'horizontalalignment',  'right',...
        'interpreter', 'tex', 'fontsize', fs+2)
ylim(handles(1), ylims)
set(handles(:), 'fontsize', fs)
%handles1-5:, colorbar, main legend
if length(handles) > 3; set(handles(4), 'fontsize', 18); end
set(h1, 'CurrentAxes', handles(1))
%     tx = text(-400, 3.5, rgns{reg}, 'fontsize', fs+2);
set(gca, 'fontsize', fs)
if reg<4; set(gca,'xticklabel',{[]}); end
xticks(-1000:500:1200)
text(-3500, 2.25, 'PHC location neurons', 'fontsize', fs)
text(-3500, 1.25, 'aligned to confirmation tap',...
    'fontsize', fs)
if ~do_binned
    text(1200, 2.25, 'Ranksum Test with just ONE BIN:', 'fontsize', fs)
    if ischar(resp_win)
        txt = 'entire presentation duration';
    elseif isnumeric(resp_win) || isempty(resp_win)
        txt = sprintf('Amy: %i-%i, H: %i-%i, E: %i-%i, P: %i-%ims',...
        resp_wins(1,1),resp_wins(1,2),resp_wins(2,1),resp_wins(2,2),...
        resp_wins(3,1),resp_wins(3,2),resp_wins(4,1),resp_wins(4,2));
    elseif which_half == 1
        txt = '1st half';
    elseif which_half == 2
        txt = '2nd half';
    end

    text(1200, 1.25, txt, 'fontsize', fs-2)
end

if ~read_from_excel && ~narrow_ylims && no_overlaps ~= 1
    fig3_source_to_xlsx(means, sems, x, pos_true, neg_true,...
        source_data_fname, sheet_name);
end


fname = sprintf('%s/figures/Figure%s.png', data_dir, plotname);
print(gcf, '-dpng', fname, '-r600');
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);

end
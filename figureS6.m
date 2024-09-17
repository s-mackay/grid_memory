% figureS6(data_dir, read_from_excel)
%
% INPUTS:
%   data_dir (str):         Path to source data (default: 'source_data').
%
%   read_from_excel (bool): Flag to read from Excel (true) 
%                           or .mat files (false) (default: false).
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function figureS6(data_dir, read_from_excel)

%formerly called response_permutation
file = 'labelshuf_10k_perms_nResponses';
file_xlsx = 'source_data_all_figures.xlsx';
sheet_name = 'figure_S6';
if nargin <1;   data_dir = 'source_data';       end
if nargin <2;   read_from_excel = false;        end
tic;
if read_from_excel
    results_struct = read_figS6_results_from_xlsx(file_xlsx, sheet_name);
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');
else
    load(sprintf('%s/%s.mat', data_dir, file), 'has_item_response_10k',...
        'has_loc_response_10k');
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
    results_cell = cell(32, 2);
end

alpha_resp = .001;
plot_visibility = 'on';
brainregs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
b_r = {'Amy', 'Hipp', 'EC', 'PHC'};
n_shuffles = 10000;

h1 = figure('visible', plot_visibility);
set(h1, 'Units', 'pixels', 'Position', [0 1920 1000 600]);

for i =1:4
    if read_from_excel
        bincenters_item = results_struct.binCenters{1, i};
        bincenters_loc  = results_struct.binCenters{2, i};
        hists_item      = results_struct.hists{1, i};
        hists_loc       = results_struct.hists{2, i};
        measured_item   = results_struct.measuredValue(1, i);
        measured_loc    = results_struct.measuredValue(2, i);
        empSize_item    = results_struct.empSize(1, i);
        empSize_loc     = results_struct.empSize(2, i);
    else
        units_reg = response_info.brain_regs == i;
        % true measured counts
        item_n = sum(units_reg & response_info.lowest_p_item < alpha_resp);
        loc_n  = sum(units_reg & response_info.lowest_p_loc  < alpha_resp);
        % get distribution of counts from labelshuffled data
        distr = sum(has_item_response_10k(units_reg, :), 1);
        % how many of the shuffled versions have more responsive units than we
        % find in true measured data?
        n_larger = sum(distr > item_n);
        % how many have equal counts?
        n_equal = sum(distr == item_n);
        % calculate fraction of more extreme outcomes in shuffled data. We
        % count half of equal counts aswell. "full" converts a sparse matrix to
        % array
        empSize_item = full((n_larger + .5*n_equal) / length(distr));
        
        edges_i = 0:max(distr)+1;
        hists_item = histcounts(distr, edges_i);
        bincenters_item = (edges_i(1:end-1) + .5)/sum(units_reg);
        measured_item = item_n/sum(units_reg);
        
        %same for location responses
        distr = sum(has_loc_response_10k(units_reg, :), 1);
        n_larger = sum(distr > loc_n);
        n_equal = sum(distr== loc_n);
        empSize_loc = full((n_larger + .5*n_equal) / length(distr));
        edges_l = 0:max(distr)+1;
        hists_loc = histcounts(distr, edges_l);
        bincenters_loc = (edges_l(1:end-1)+ .5) / sum(units_reg);
        measured_loc = loc_n/sum(units_reg);
        
        % create a cell array that will be saved to .xlsx. the first column 
        % will be var names, the 2nd column are numeric arrays or single 
        % values. There will only be 1d arrays
        results_cell{   i, 1} = sprintf('hists_item_%s', b_r{i});
        results_cell{   i, 2} = hists_item;
        results_cell{ 4+i, 1} = sprintf('binCenters_item_%s', b_r{i});
        results_cell{ 4+i, 2} = bincenters_item;
        results_cell{ 8+i, 1} = sprintf('hists_loc_%s', b_r{i});
        results_cell{ 8+i, 2} = hists_loc;
        results_cell{12+i, 1} = sprintf('binCenters_loc_%s', b_r{i});
        results_cell{12+i, 2} = bincenters_loc;
        results_cell{16+i, 1} = sprintf('measuredValue_item_%s', b_r{i});
        results_cell{16+i, 2} = measured_item;
        results_cell{20+i, 1} = sprintf('measuredValue_loc_%s', b_r{i});
        results_cell{20+i, 2} = measured_loc;
        results_cell{24+i, 1} = sprintf('empSize_item_%s', b_r{i});
        results_cell{24+i, 2} = empSize_item;
        results_cell{28+i, 1} = sprintf('empSize_loc_%s', b_r{i});
        results_cell{28+i, 2} = empSize_loc;
    end
    
    % plot item responses
    subplot(2,4,i)
    bar(bincenters_item, hists_item, 'EdgeColor', 'none', 'BarWidth', 1);
    hold on
    xline(measured_item, 'r');
    fprintf('in %s, empirical size for item responses is %.3g.\n',...
        brainregs{i}, empSize_item);
    title(sprintf('%s', brainregs{i}));
    
    xl = xlim;
    yl = ylim;
    if i==1
        legend({sprintf('%i shuffles', n_shuffles),...
            'measured value'}, 'location', 'east')
        ylabel('iteration count')
        text(-.045, mean(yl), sprintf('Item\nCells'), 'FontSize', 11,...
            'HorizontalAlignment', 'right');
    end
    if empSize_item == 0
        text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n< 0.0001'), ...
            'horizontalalignment', 'right');
    else
        text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n= %.2g',...
            empSize_item), 'horizontalalignment', 'right');
    end
    xlim(xl); ylim(yl);
    
    % plot location responses
    subplot(2,4,i+4)
    bar(bincenters_loc, hists_loc, 'EdgeColor', 'none', 'BarWidth', 1);
    hold on
    xline(measured_loc, 'r');
    fprintf('in %s, empirical size for location responses is %.3g.\n',...
        brainregs{i}, empSize_loc);
    title(sprintf('%s', brainregs{i}));
    xlim(xl);
    ylim(yl);
    text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n= %.2g',...
        empSize_loc), 'horizontalalignment', 'right');
    if i == 1
        ylabel('iteration count')
        text(-.045, mean(yl), sprintf('Location\nCells'), 'FontSize', 11,...
            'HorizontalAlignment', 'right');
    end

end

%% optional: save results to .xlsx
% if ~read_from_excel
%     %delete the old sheet so that there is no data remaining. Or rather, since
%     %deleting is turning out to be tricky (if you dont have Excel installed), 
%     %we're doing the ugly workaround of writing 1000 x 500 empty cells into
%     %the sheet
%     sheet_name = 'figure_S6';
%     emptyData = cell(3000, 1000);
%     % Write the empty cell array to the specified sheet
%     writecell(emptyData, file_xlsx, 'Sheet', sheet_name, 'Range', 'A1');
%     writecell(results_cell,  file_xlsx, 'Sheet', sheet_name, 'Range', 'A1');
% end

%% print figures
fname_png = sprintf('%s/figures/FigureS6.png', data_dir);
fname_eps = sprintf('%s/figures/FigureS6.eps', data_dir);
print(gcf, '-dpng', fname_png, '-r600');
print('-painters', gcf, '-depsc', fname_eps);
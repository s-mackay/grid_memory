% figureS2 (data_dir, read_from_excel)
%   creates figure S2 (line pltos, "Selectivity of responses to items and 
%                                                       spatial locations"
%
%   Inputs:
%       data_dir (char): Directory path where data files are located.
%       read_from_excel (logical): Indicates whether to read data from an Excel file (true) or from files in data_dir (false).
%
%   Example usage:
%       figureS2('/path/to/data', true)  % Read data from Excel file
%       figureS2('/path/to/data', false) % Read data from .mat files
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function figureS2(data_dir, read_from_excel)

if read_from_excel
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');
else
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
end

% the edges are irregular. Each session can have a different sized item
% pool, such that a response to a single item may correspond to 1/5 of
% stimuli in one session, but 1/6 or 1/7 in another. We include all
% possible steps for the range of item pools from n=4 to n=8.
% exact fractions don't work for data saved in .xlsx format, therefore
% switching to unique values
% edges_i = [0, unique([(1:4)/4, (1:5)/5, (1:6)/6, (1:7)/7, (1:8)/8])];
edges_i = unique(response_info.frac_pref_items);
% for locations, it's more straight forward. steps of 1/9.
% edges_l = (0:9)/9;
edges_l = unique(response_info.frac_pref_locs);
item_data = zeros(length(edges_i), 4);
pos_data = zeros(length(edges_l), 4);
pos_data_hist = zeros(8, 4);
l_mat = lines(4); % color scheme
fs = 12; % font size
reglabels = {'Amy', 'Hipp', 'EC', 'PHC'};


for reg = 1:4
    % include only units in one brain region
    which_units = response_info.brain_regs==reg;
    % ITEMS
    % plot the fractions of items responded to
    data = response_info.frac_pref_items(which_units);
    % sort fractions ascending
    frac_img_resp_to = sort(data(data>0));
    % how many units responded to more than 50% of items?
    % find works because fractions are sorted
    overHalf = find(frac_img_resp_to > .5, 1);
    frac_underHalf = (overHalf-1)/length(frac_img_resp_to);
    if isempty(frac_underHalf) && ~isempty(frac_img_resp_to)
        frac_underHalf=1;
    end
    fprintf('In %s, %.2f%% of ',reglabels{reg}, frac_underHalf*100);
    fprintf('responsive units responded to at most half of the items\n');
    item_data(:, reg) = [cumsum(histcounts(data, edges_i))/length(data), 1];
    % LOCATIONS
    data = response_info.frac_pref_locs(which_units);
    frac_pos_resp_to = sort(data(data>0));
    overHalf = find(frac_pos_resp_to > .5, 1);
    frac_underHalf = (overHalf-1)/length(frac_pos_resp_to);
    if isempty(frac_underHalf) && ~isempty(frac_pos_resp_to)
        frac_underHalf=1;
    end
    if reg == 4
        fprintf('In %s, %.2f%% of', reglabels{reg}, frac_underHalf*100);
        fprintf('responsive units responded to at most half of the locations\n');
    end
    pos_data(:, reg) = [cumsum(histcounts(data, edges_l))/length(data), 1];
    pos_data_hist(:, reg) = histcounts(data, edges_l(2:end));
end

% create plots
h2 = figure;
set(h2,'PaperUnits','centimeters');
% set(h2, 'PaperSize', [18 8], 'PaperPosition', [0 0 18 8]);
set(h2, 'PaperSize', [12.6 7.5], 'PaperPosition', [0 0 12.6 7.5]);

subplot(121) % item responses
plot(edges_i, item_data(:, 1), 'linewidth', 1.5,'color', l_mat(1,:))
hold on
plot(edges_i, item_data(:, 2), 'linewidth', 1.5,'color', l_mat(2,:))
plot(edges_i, item_data(:, 3), 'linewidth', 1.5,'color', l_mat(3,:))
plot(edges_i, item_data(:, 4), 'linewidth', 1.5,'color', l_mat(4,:))
text(-.4, 1.02, 'A', 'fontweight', 'bold', 'fontsize', 14)
legend({'item Amy', 'item Hipp', 'item EC', 'item PHC'},...
    'location', 'southeast')
title('item responses', 'fontsize', 12, 'fontweight', 'normal')
ylabel('fraction of cells')
xlabel({'fraction of images', 'responded to'})
grid on; %ax = gca;  ax.YGrid = 'on';
set(gca, 'fontsize', fs)
legend('boxoff')
ylims = ylim;

subplot(122) % location responses
plot(edges_l, pos_data(:, 1), 'linewidth', 1.5,'color', l_mat(1,:))
hold on
plot(edges_l, pos_data(:, 2), 'linewidth', 1.5,'color', l_mat(2,:))
plot(edges_l, pos_data(:, 3), 'linewidth', 1.5,'color', l_mat(3,:))
plot(edges_l, pos_data(:, 4), 'linewidth', 1.5,'color', l_mat(4,:))
text(-.4, 1.02, 'B', 'fontweight', 'bold', 'fontsize', 14)
legend({'loc amy', 'loc hipp', 'loc EC', 'loc PHC'},...
    'location', 'southeast')
title('location responses', 'fontsize', 12, 'fontweight', 'normal')
xlabel({'fraction of grid locations', 'responded to'})
grid on; %ax = gca; ax.YGrid = 'on';
set(gca, 'fontsize', fs)
legend('boxoff')
ylim(ylims)

%print figures
fname = sprintf('%s/figures/FigureS2.png', data_dir);
fname_eps = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps);
print(gcf, '-dpng', fname, '-r300')

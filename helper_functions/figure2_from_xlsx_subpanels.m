function figure2_from_xlsx_subpanels(results, item_or_loc, id)

% function figure2_from_xlsx(results, item_or_loc, id)
% Produces a figure from summary data in an Excel sheet
%
% This function generates a figure based on the summary data provided in
% the input struct 'results', which is loaded from an Excel file. The
% figure produced corresponds to figure 2 in a publication, showing either
% item or location cells as specified by the inputs.
%
% INPUTS:
%   results:     struct containing summary data loaded from xlsx
%   item_or_loc: int, 1 for item, 2 for loc
%   id:          int, can be 1 or 2 since figure 2 shows 2 item and 
%                2 location cell examples
%
% USAGE:
%  
%   % To produce a figure for the first item cell example:
%   figure2_from_xlsx(results, 1, 1);
%
%   % To produce a figure for the second location cell example:
%   figure2_from_xlsx(results, 2, 2);
%
% NOTE:
% This is a reduced version of figure2. No statistical analysis is 
% performed, it simply produces a figure from the summary data in the 
% spreadsheet.

%% configure figure

fs = 16; %font size
h1 = figure('visible', 'on');
set(h1, 'PaperUnits','centimeters');
set(h1, 'PaperSize', [10, 13], 'PaperPosition', [0, 0, 10, 13]);
padding     = .11; % white frame around everything
pl_h        = (1 - 4 * padding) / 3; % we want 3 rows of plots
full_w      = 1 - 2.8 * padding; % for raster plot and convolved plot
thumbnail_w = (full_w - 1.5 * padding) * (1/3); %thumbnail width
dens_w      = (full_w - 1.5 * padding) * (2/3); %density plot width

%% start plotting

% stimulus image
axes('Position', [padding, 3 * padding + 2 * pl_h, thumbnail_w, pl_h]);
imshow(results.pref_jpg)

% spike shape
axes('Position', [1.5 * padding + 1.5 * padding + thumbnail_w,...
    3 * padding + 2 * pl_h, dens_w, pl_h]);
load(results.spikefile, 'spikeshapes')
density_plot(spikeshapes)
xticks([-.5, 0, .5, 1, 1.34])
xticklabels({'-.5', '0', '.5', '1', 'ms'})
xlabel('')
set(gca, 'Fontsize', fs)

% raster plot
rel_ts_r = results.rel_ts_rem;
rel_ts_f = results.rel_ts_forg;
axes('Position', [1.5 * padding, 2* padding + pl_h, full_w, pl_h])
n_trials_r_f = [length(results.rel_ts_rem), length(results.rel_ts_forg)];
% plot_resp_byID(id, item_or_loc, spiketable)
% individual spikes
% plot red lines, one trial, i.e. row, at a time
for trF = 1:length(rel_ts_f)
    ypos = length(rel_ts_f)- trF + 1;
    xvals = [rel_ts_f{trF}; rel_ts_f{trF}];
    yvals = repmat([ypos-0.5; ypos+0.5], 1, length(rel_ts_f{trF}));
    line(xvals, yvals, 'LineWidth', 0.8, 'Color', 'r');
end
% plot horizontal grey lines across the whole panel
lines_ypos = [length(rel_ts_f)+length(rel_ts_r)-[cumsum(n_trials_r_f(:,1));...
    sum(n_trials_r_f(:,1))+cumsum(n_trials_r_f(:,2))]]+.5;
xvals = repmat([-500; 1500], [1, numel(lines_ypos)-1]);
line(xvals, repmat(lines_ypos(1:end-1)', 2, 1), 'color', [.5,.5,.5])
% plot blue lines, one trial, i.e. row, at a time
for trC = 1:length(rel_ts_r)
    ypos = trF+length(rel_ts_r)-trC+1;
    xvals = [rel_ts_r{trC}; rel_ts_r{trC}];
    yvals = repmat([ypos-0.5; ypos+0.5], 1, length(rel_ts_r{trC}));
    line(xvals, yvals, 'LineWidth', 0.8, 'Color', 'b');
end
%plot black vertical line at x=0
line([0 0],[-1.7 length(rel_ts_r)+length(rel_ts_r)+1.5], 'Color', 'k');
ylim([0, length(rel_ts_r)+length(rel_ts_f)+1])
xlim([-500 1500])
ylabel('trials'); %xlabel('ms');
set(gca, 'Fontsize', fs)

% next panel, plot lines
axes('Position', [1.5*padding, padding, full_w, pl_h]);
hold on
plot(results.x, results.non_pref_rem,  '--b');
plot(results.x, results.non_pref_forg, '--r');
hclp = plot(results.x, results.pref_rem,  'b', 'linewidth', 1.5);
hflp = plot(results.x, results.pref_forg, 'r', 'linewidth', 1.5);
maxval = max([results.pref_rem, results.pref_forg]);
line([0 0],[-1.7 maxval*1.2], 'Color', 'k'); % stim presentation time
xlim([-500 1500])
if id == 1 && item_or_loc == 1
    yls = ylim; ylim([-.05*yls(2) 55])
elseif id == 2 && item_or_loc == 1
    yls = ylim; ylim([-.05*yls(2) 20])
elseif id == 1 && item_or_loc == 2
    yls = ylim; ylim([-.05*yls(2) 18])
elseif id == 2 && item_or_loc == 2
    yls = ylim; ylim([-.05*yls(2) 20])
else
    yls = ylim; ylim([-.05*yls(2) yls(2)])
end
lgnd = legend([hclp, hflp], {'remembered', 'forgotten'}, 'Location', 'northeast',...
    'numcolumns', 1, 'fontsize', fs-2);
legend('boxoff')
lgnd.Position =[.475 .2 .5 .1];
xlabel('ms')
ylabel('Hz')
set(gca, 'Fontsize', fs)
xlim([-500 1500])
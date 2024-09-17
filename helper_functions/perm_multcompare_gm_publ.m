%% function to run pairwise comparision and one-sample t-test
% against zero and plot all results into the same plot
% function [handles] = perm_multcompare(condA, condB, nperms, ...
%                          clusteralpha, alpha, xaxis, plotcolors,...
%                          conditionnames, ax1, markSize, LfontSiz,...
%                          plot_info)
%
% Portions of this code were adapted from work by Thomas Reber, with 
% permission. Original work can be found at 
% https://github.com/rebrowski/neuralAdapatationInMTL

function [handles, results] = perm_multcompare_gm_publ(...
    condA, condB,nperms, clusteralpha, alpha, xaxis, ...
    conditionnames, ax1, markSize, LfontSiz, plot_info)

if ~exist('nperms', 'var')
    nperms = 10000;
end

if ~exist('clusteralpha', 'var')
    clusteralpha = 0.005;
end

if ~exist('alpha', 'var')
    alpha = 0.05;
end

if ~exist('xaxis', 'var')
    xaxis = 1:size(condA, 2);
end

if ~exist('conditionnames', 'var') || isempty(conditionnames)
    conditionnames = {'A', 'B'};
    legendoff = true;
else
    legendoff = false;
end

if ~exist('markSize', 'var')
    markSize = 10;
end

if ~exist('LfontSiz', 'var')
    LfontSiz = markSize;
end

if ~exist('adjust_y_axis', 'var')
    adjust_y_axis = false;
end

results = struct;

%% run the paired t-test
doplot = false;
if size(condA, 1) > 1
    [pos_true, neg_true, ~, ~, ~, ~, ~] = ...
        perm_ttest(condA, condB, nperms, [],clusteralpha, alpha, doplot);
else
    pos_true = struct;
    pos_true.froms = [];
    pos_true.tos = [];
    pos_true.sizes = [];
    pos_true.sumts = [];
    pos_true.cluster_p = [];
    neg_true = pos_true;
end
    

%% plot curves
if ~exist('ax1', 'var') || isempty(ax1)
    ax1 = gca;
end

plotcolors = plot_info.plot_colors;
lineh(1) = plot_signals(condA, xaxis, plotcolors(1,:));
lineh(2) = plot_signals(condB, xaxis, plotcolors(2,:));
xlim([xaxis(1) xaxis(end)]);
ylabel(sprintf('%s %s boot SEM', '{\itM}', '\pm'), 'interpreter', 'tex');

grid on
pos = get(ax1, 'Position');
set(ax1, 'Position', pos);

if ~legendoff
    [~, objh] = legend(lineh, conditionnames, 'Position', ...
        [.79, .13, .09, .02], 'FontSize', LfontSiz); 
    set(objh,'linewidth',2);
end

% make sure the legend is not written on top of data
if ~exist('plot_info', 'var') || isempty(plot_info.ylims)
    yl = ylim;
    if diff(yl) <.5 && adjust_y_axis
        yl(1) = -.04;
    else
        yl(1) = -.5;
    end
    ylim([yl(1) yl(2)+diff(yl)*0.1]);
elseif any(isnan(plot_info.ylims))
    %dont take action. need this in cross corr script
else
    ylim(plot_info.ylims);
end
yl = ylim;

%% annotate sign clusters paired test
y_star = yl(1)+diff(yl)/15; %-0.2;%yl(1);
if size(condA, 1)>1 && size(condB, 1)>1
for c = 1:length(pos_true.froms)
    
    plot(xaxis([pos_true.froms(c), pos_true.tos(c)]), [yl(1), yl(1)], ...
        'linewidth', 5, 'color', [.4 .4 .4]);
    x_ = xaxis(round(mean([pos_true.tos(c) pos_true.froms(c)])));
    if pos_true.cluster_p(c) < .001
        text(x_,y_star, '***', 'HorizontalAlignment', 'center', ...
            'fontname', 'monospaced');
    elseif pos_true.cluster_p(c) < .01
        text(x_,y_star, '**', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    elseif pos_true.cluster_p(c) < .05
        text(x_,y_star, '*', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    end


    meanA = nanmean(condA(:,pos_true.froms(c):pos_true.tos(c)));
    meanB = nanmean(condB(:,pos_true.froms(c):pos_true.tos(c)));
    stdA = nanstd(condA(:, pos_true.froms(c):pos_true.tos(c)));
    stdB = nanstd(condB(:, pos_true.froms(c):pos_true.tos(c)));
    d = (meanA-meanB)./sqrt((stdA.^2+stdB.^2)./2);
    if plot_info.add_labels
        text(x_, yl(1)+diff(yl)/15*(c+1),...
            sprintf('max(d)=%.3f mean=%.3f, sd=%.3f', max(d), mean(d),...
            std(d)), 'HorizontalAlignment', 'center')
    end
    fprintf(sprintf('    max(d)=%.3f mean=%.3f, sd=%.3f\n', max(d), mean(d), std(d)))
    fprintf('    %i - %ims, pval true %.8g\n', ...
        xaxis(pos_true.froms(c)), xaxis(pos_true.tos(c)), pos_true.cluster_p(c))
end

for c = 1:length(neg_true.froms)
    plot(xaxis([neg_true.froms(c), neg_true.tos(c)]), [yl(1), yl(1)], ...
        'linewidth', 5, 'color', [.5 .3 .3]);
    x_ = xaxis(round(mean([neg_true.tos(c) neg_true.froms(c)])));
    if neg_true.cluster_p(c) < .001
        text(x_, y_star, '***', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    elseif neg_true.cluster_p(c) < .01
        text(x_, y_star, '**', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    elseif neg_true.cluster_p(c) < .05
        text(x_, y_star, '*', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    end
%     results = [results sprintf('%s < %s \\t%d to %d\\t %0.3f\\n', ...
%                    conditionnames{1}, conditionnames{2}, ...
%                    xaxis(neg_true.froms(c)), xaxis(neg_true.tos(c)), ...
%                    neg_true.cluster_p(c))];
    meanA = mean(condA(:,neg_true.froms(c):neg_true.tos(c)));
    meanB = mean(condB(:,neg_true.froms(c):neg_true.tos(c)));
    stdA = std(condA(:, neg_true.froms(c):neg_true.tos(c)));
    stdB = std(condB(:, neg_true.froms(c):neg_true.tos(c)));
    d = (meanA-meanB)./sqrt((stdA.^2+stdB.^2)./2);
    if plot_info.add_labels
        text(x_, yl(1)+diff(yl)/15*(c+length(pos_true.froms)+1),...
            sprintf('max(d)=%.3f, mean=%.3f sd=%.3f', max(d), mean(d), ...
            std(d)), 'HorizontalAlignment', 'center')
    end
    fprintf(sprintf('    max(d)=%.3f, mean=%.3f sd=%.3f\n', max(d), mean(d), std(d)))
     fprintf('    %i - %ims, pval neg true %.8g\n', ...
         xaxis(neg_true.froms(c)), xaxis(neg_true.tos(c)), neg_true.cluster_p(c))
end


end
%calculate some summary data for saving to excel. This is normally done in
%plot_signals.
both_singals = {condA, condB};
[m, sem] = deal(nan(2, size(condA, 2)));
for i_s = 1:length(both_singals)
    signal = both_singals{i_s};
    nsigs = size(signal,1);
    m(i_s, :) = mean(signal, 1);
    sd = std(signal, 1);
    sem(i_s, :) = sd./sqrt(nsigs);
end
results.means     = m;
results.sem       = sem;
results.pos_true  = pos_true;
results.neg_true  = neg_true;
handles = ax1;
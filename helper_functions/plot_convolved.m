% plot_convolved(results, which_panel, plot_info, curr_ax)
%
% this function plots the main figures from summary source data saved in
% .xlsx format, rather than raw data
% it is a combination of elements of perm_multcompare_gm_publ.m and
% plot_signals.m, and was put together in response to final requests to
% combine all plotted data into one excel sheet. 
%
% Portions of this code were adapted from work by Thomas Reber, with 
% permission. Original work can be found at 
% https://github.com/rebrowski/neuralAdapatationInMTL
% -------------------------------------------------------------------------

function handles = plot_convolved(results, which_panel, plot_info, curr_ax)

m = cell2mat(results.means(which_panel(1), which_panel(2)));
sem = cell2mat(results.sems(which_panel(1), which_panel(2)));
x = cell2mat(results.x(which_panel(1), which_panel(2)));
error_type = 'boot SEM';
pos_true = results.pos_true{which_panel(1), which_panel(2)};
neg_true = results.neg_true{which_panel(1), which_panel(2)};
plot_colors = plot_info.plot_colors;
condA = m(1, :);
condB = m(2, :);
yl = plot_info.ylims;

if ~exist('curr_ax', 'var') || isempty(curr_ax)
    curr_ax = gca;
end
pos = get(curr_ax, 'Position');
set(curr_ax, 'Position', pos);
               
hold on
if sum(diff(sem(1,:))) ~= 0 % skip when there's 0 or 1 trace only
    for i = 1:2
        plot_color = plot_colors(i, :);
        signal = m(i, :);
        error_upper = signal + sem(i, :);
        error_lower = signal - sem(i, :);
        X = [x, fliplr(x)];
        Y = [error_upper, fliplr(error_lower)];
        f = fill(X, Y, plot_color);
        set(f,'EdgeColor','none'),
        alpha(0.25);
    end
end

% plot means
plot(x, m(1, :), 'color', plot_colors(1, :));
plot(x, m(2, :), 'color', plot_colors(2, :));

y_star = yl(1)+diff(yl)/15;
ylim(yl)
box off
hold on
%put in bars marking significant clusters where the response during
%remembered trials is stronger than during subsequently forgotten trials
for c = 1:length(pos_true.froms)
    plot(x([pos_true.froms(c), pos_true.tos(c)]), [yl(1), yl(1)], ...
        'linewidth', 5, 'color', [.4 .4 .4]);
    x_ = x(round(mean([pos_true.tos(c) pos_true.froms(c)])));
    if pos_true.cluster_p(c) < .001
        text(x_, y_star, '***', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    elseif pos_true.cluster_p(c) < .01
        text(x_, y_star, '**', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    elseif pos_true.cluster_p(c) < .05
        text(x_, y_star, '*', 'HorizontalAlignment', 'center',...
            'fontname', 'monospaced');
    end
    ylim(yl)
end

% put in bars marking significant clusters where the response during
% subsequently forgotten trials was stronger
for c = 1:length(neg_true.froms)
    plot(x([neg_true.froms(c), neg_true.tos(c)]), [yl(1), yl(1)], ...
        'linewidth', 5, 'color', [.5 .3 .3]);
    x_ = x(round(mean([neg_true.tos(c) neg_true.froms(c)])));
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
    ylim(yl)
end
xlim([x(1), x(end)])

handles = curr_ax;
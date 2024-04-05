function response_permutation(path, read_from_excel)

file = 'labelshuf_10k_perms_nResponses';
tic;
if read_from_excel
    has_item_response_10k = readmatrix(sprintf('%s/%s.xlsx', path, file),...
        'Sheet', 'has_item_response');
    toc;
    has_loc_response_10k = readmatrix(sprintf('%s/%s.xlsx', path, file),...
        'Sheet', 'has_loc_response');
    toc;
    response_info = readtable(sprintf('%s/response_info.xlsx', path));
else
    load(sprintf('%s/%s.mat', path, file), 'has_item_response_10k',...
        'has_loc_response_10k');
    load(sprintf('%s/response_info.mat', path), 'response_info');
end

alpha_resp = .001;
plot_visibility = 'on';
brainregs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};

h1 = figure('visible', plot_visibility);
set(h1, 'Units', 'pixels', 'Position', [0 1920 1000 600]);
for i =1:4
    units_reg = response_info.brain_regs == i;
    % true measured counts
    item_n = sum(units_reg & response_info.lowest_p_item < alpha_resp);
    pos_n  = sum(units_reg & response_info.lowest_p_loc  < alpha_resp);
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
    emp_size = full((n_larger + .5*n_equal) / length(distr));
    subplot(2,4,i)
    edges = 0:max(distr)+1;
    bar((edges(1:end-1) + .5)/sum(units_reg), histcounts(distr, edges),...
        'EdgeColor', 'none', 'BarWidth', 1);
    hold on
    xline(item_n/sum(units_reg), 'r');
    fprintf('in %s, empirical size for item responses is %.3g.\n',...
        brainregs{i}, emp_size);
    units_reg = response_info.brain_regs == i;
    title(sprintf('%s', brainregs{i}));
    
    xl = xlim;
    yl = ylim;
    if i==1
        legend({sprintf('%i shuffles', size(has_item_response_10k, 2)),...
            'measured value'}, 'location', 'east')
        assert(size(has_item_response_10k, 2) == size(has_loc_response_10k, 2));
        ylabel('iteration count')
        text(-.045, mean(yl), sprintf('Item\nCells'), 'FontSize', 11,...
            'HorizontalAlignment', 'right');
    end
    if emp_size == 0
        text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n< 0.0001'), ...
            'horizontalalignment', 'right');
    else
        text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n= %.2g',...
            emp_size), 'horizontalalignment', 'right');
    end
    xlim(xl); ylim(yl);

    % same for location responses
    distr = sum(has_loc_response_10k(units_reg, :), 1);
    n_larger = sum(distr > pos_n);
    n_equal = sum(distr== pos_n);
    emp_size = full((n_larger + .5*n_equal) / length(distr));
    subplot(2,4,i+4)
    edges = 0:max(distr)+1;
    bar((edges(1:end-1)+ .5)/sum(units_reg), histcounts(distr, edges),...
        'EdgeColor', 'none', 'BarWidth', 1);
    hold on
    xline(pos_n/sum(units_reg), 'r');
    fprintf('in %s, empirical size for location responses is %.3g.\n',...
        brainregs{i}, emp_size);
    title(sprintf('%s', brainregs{i}));
    xlim(xl);
    ylim(yl);
    text(xl(2)*.9, yl(2)*.9, sprintf('empirical size\n= %.2g',...
        emp_size), 'horizontalalignment', 'right');
    if i == 1
        ylabel('iteration count')
        text(-.045, mean(yl), sprintf('Location\nCells'), 'FontSize', 11,...
            'HorizontalAlignment', 'right');
    end
end
fname_png = sprintf('%s/figures/FigureS6.png', path);
fname_eps = sprintf('%s/figures/FigureS6.eps', path);
print(gcf, '-dpng', fname_png, '-r600');
print('-painters', gcf, '-depsc', fname_eps);
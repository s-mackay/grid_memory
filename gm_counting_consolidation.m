function gm_counting_consolidation(data_dir, read_from_excel)

% function gm_counting_consolidation(data_dir, read_from_excel)
% here we consider all units responding to exactly one item or location
% we compare their firing rates during the distraction task for
% subsequently correct vs. incorrect responses.

if read_from_excel
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');
    counting_table = readtable(sprintf('%s/fr_counting.xlsx', data_dir));
else
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
    load(sprintf('%s/fr_counting.mat', data_dir), 'counting_table');
end

regs = response_info.brain_regs;

regions = {'Amy ', 'Hipp', 'EC  ', 'PHC ', 'All '};
h1 = figure('visible', 'on');
set(h1, 'Units', 'pixels', 'Position', [0 1920 1100 600]);
asx = {'', '*'}; % asterisk (or not)
yl = [-1 1.4]; %ylims
xl = [.5 2.5]; %xlims
fs = 12; % font size


fprintf('\n---------- COUNTING / DISTRACTOR TASK ----------')
fprintf('\nCOMPARING NEURONAL FIRING WHILE RETENTION VS. FORGETTING\n\n')
fprintf('counting task during retention vs forgetting (pref. item):\n')
for reg = 1:5
    if reg == 5 % the last 'brain region' is all regions pooled
        incl_these_i = find(response_info.n_pref_items==1);
        incl_these_l = find(response_info.n_pref_locs==1);
    else
        incl_these_i = find(regs==reg & response_info.n_pref_items==1);
        incl_these_l = find(regs==reg & response_info.n_pref_locs ==1);
%         disp(sprintf('n item responsive %s: %i', regions{reg}, ...
%             sum(regs==reg & response_info.n_pref_items > 0)));
%         disp(sprintf('n loc  responsive %s: %i', regions{reg}, ...
%             sum(regs==reg & response_info.n_pref_locs > 0)));
    end
    % get fr during counting for all units that have exactly 1 stim
    % response, separated out into counting while remembereing and while
    % forgetting
    fr_i = [counting_table.fr_item_c(incl_these_i),...
            counting_table.fr_item_f(incl_these_i)];
    fr_l = [counting_table.fr_loc_c(incl_these_l),...
            counting_table.fr_loc_f(incl_these_l)];
    % does activity while retaining and forgetting differ?
    [~, p_s, ci, stats] = ttest(fr_i(:, 1), fr_i(:, 2));
    fprintf('%s: p = %.2g, t = %.2g, df = %i, ci = [%.2g, %.2g]\n',...
        regions{reg}, p_s, stats.tstat, stats.df, ci(1), ci(2))
    
    subplot(2,5,reg)
    boxplot(fr_i(:,1), 'positions', 1);% ylim(yl);
    hold on;
    boxplot(fr_i(:,2), 'positions', 2);
    ylim(yl); xlim(xl); 
    title(sprintf('%s%s', asx{(p_s<.05)+1}, regions{reg}))
    if reg == 1
        text(-0.35, mean(yl), sprintf('Item\nNeurons'), 'fontsize',...
            fs, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
        % use text instead of ylabel bc ylabel will squish the subplot
        text(.01, mean(yl), 'mean FR rel. to bsl.', 'Rotation', 90,...
            'FontSize', fs, 'HorizontalAlignment', 'center');
%         ylabel('mean fr rel. to baseline')
    end
    xticks([1,2])
    xticklabels({'rem.', ' forg.'});
    % set manually to skip '0.5' (takes up too much space)
    yticks(-1:.5:1)
    yticklabels({'-1', '', '0', '', '1'})
    xlabel(sprintf('n = %i', length(incl_these_i)));
    plot([1.1, 1.9], [1, 1], 'k')
    set(gca, 'fontsize', fs)
    if p_s <= .05
        text(1.5, 1.2, sprintf('*%s = %.3f', '{\itp}', p_s), ...
            'HorizontalAlignment', 'center', 'FontWeight', 'bold',...
            'fontsize', fs);
    else
        text(1.5, 1.2, sprintf('n.s., %s = %.2f', '{\itp}', p_s),...
            'HorizontalAlignment', 'center', 'fontsize', fs);
    end
    grid on
    
    subplot(2,5,5+reg)
    if reg == 1
        % label the row of subplots from here
        ylim(yl); xlim(xl);
        text(-.9, mean(yl), sprintf('Location\nNeurons'), 'fontsize',...
            fs, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        axis off
        box off
    % we only find a significant population of location neurons in PHC
    elseif reg == 4
        % does activity differ during counting task, depending on
        % subsequent memory performance?
        [~, p_p, ci, stats] = ttest(fr_l(:,1), fr_l(:,2));
        if size(fr_l, 1)==1
            fprintf('%s: not enough location selective neurons', regions{reg})
        else
            fprintf('%s (location): p = %.2g, t = %.2g, df = %i, ci = [%.2g, %.2g]\n',...
                regions{reg}, p_p, stats.tstat, stats.df, ci(1), ci(2))
        end
        boxplot(fr_l(:,1), 'positions', 1); %ylim(yl);
        title(sprintf('%s%s', asx{(p_p<.05)+1}, regions{reg}))
        hold on
        boxplot(fr_l(:,2), 'positions', 2); ylim(yl);xlim(xl);
        text(.01, mean(yl), 'mean FR rel. to bsl.', 'Rotation', 90,...
            'FontSize', fs, 'HorizontalAlignment', 'center');
        xticks([1,2])
        xticklabels({'rem.', 'forg.'})
        yticks(-1:.5:1)
        yticklabels({'-1', '', '0', '', '1'})
        xlabel(sprintf('n = %i', length(incl_these_l)))
        plot([1.1, 1.9], [1, 1], 'k')
        set(gca, 'fontsize', fs);
        if p_p <= .05
            text(1.5, 1.2, sprintf('*%s = %.3f', '{\itp}', p_p), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold',...
                'fontsize', fs);
        else
            text(1.5, 1.2, sprintf('n.s., %s = %.2f', '{\itp}', p_p),...
                'HorizontalAlignment', 'center', 'fontsize', fs);
        end
        grid on
    else
        box off
        axis off
    end
end
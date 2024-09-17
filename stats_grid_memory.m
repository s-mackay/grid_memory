% stats_grid_memory (data_dir, read_from_excel, show_plots)
% this function prints stats into the console and 
% generates some figures
% 
%   INPUTS:
%          data_dir (str): path to source data
%  read_from_excel (bool): load data from .xlsx? (as opposed to .mat)
%       show_plots (bool): visible plots?
%
% Mackay et al. 2024 (DOI:10.1038/s41467-024-52295-5)
% License: MIT License (see LICENSE file for details)
% -------------------------------------------------------------------------

function stats_grid_memory(data_dir, read_from_excel, show_plots)

if nargin < 1
    data_dir = '.'; % or indicate where data is saved
end

if read_from_excel
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');
else
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
end

if show_plots
    plot_visibility = 'on';
else
    plot_visibility = 'off';
end

% get fraction of responsive cells per brain region
brain_regs = {'Amygdala', 'Hippocampus', 'EC', 'PHC'};
alpha_resp = .001;


%% Calculating Binomial Tests and Sparsely Firing Units
% to check for significant numbers of responsive cells

% collect fractions in array (will later convert to table)
binom_item = nan(4, 3);
binom_loc  = nan(4, 3);
resp2half_items = nan(4,3);
resp2half_locs = nan(4,3);

for r = 1:length(brain_regs)
    in_brain_reg = response_info.brain_regs == r;
    n_in_brain_reg  = sum(in_brain_reg);
    responsive_loc  = response_info.n_pref_locs > 0;
    responsive_item = response_info.n_pref_items > 0;
    n_responsive_i = sum(in_brain_reg & responsive_item);
    n_responsive_l = sum(in_brain_reg & responsive_loc);
    % number of units responding to at most half of items / locations
    n_resp2half_i = sum(in_brain_reg & response_info.n_pref_items > 0 &...
        response_info.frac_pref_items <=.5);
    n_resp2half_l = sum(in_brain_reg & response_info.n_pref_locs > 0 & ...
        response_info.frac_pref_locs <=.5);
    prcnt_resp2half_i = n_resp2half_i / n_responsive_i * 100;
    prcnt_resp2half_l = n_resp2half_l / n_responsive_l * 100;
    % using this sum instead of the commented out alternative in the line
    % below because the values get very small and are rounded to zero in
    % the other version (namely "1-(binocdf(n_responsive_i-1,
    % n_in_brain_reg, alpha_resp));")
    binocdf_i = sum(binopdf(n_responsive_i : n_in_brain_reg,...
        n_in_brain_reg, alpha_resp));
%     binocdf_i = 1-(binocdf(n_responsive_i-1, n_in_brain_reg, alpha_resp));
    binocdf_l = sum(binopdf(n_responsive_l : n_in_brain_reg,...
        n_in_brain_reg, alpha_resp));
%     binocdf_l = 1-(binocdf(n_responsive_l-1, n_in_brain_reg, alpha_resp)

    binom_item(r, :) = [n_responsive_i, n_in_brain_reg, binocdf_i];
    binom_loc(r, :)  = [n_responsive_l, n_in_brain_reg, binocdf_l];
    resp2half_items(r,:) = [n_resp2half_i, n_responsive_i, prcnt_resp2half_i];
    resp2half_locs(r,:)  = [n_resp2half_l, n_responsive_l, prcnt_resp2half_l];
end
table_binom_item = array2table(binom_item, 'RowNames', brain_regs,...
    'VariableNames', {'n responsive', 'n total', 'p (binomial test)'});
table_binom_loc  = array2table(binom_loc,  'RowNames', brain_regs,...
    'VariableNames', {'n responsive', 'n total', 'p (binomial test)'});

table_half_item = array2table(resp2half_items, 'RowNames', brain_regs,...
    'VariableNames', {'n responding to <=50% of items', 'n responsive', ...
    '% responding to <=50% of items'});
table_half_loc = array2table(resp2half_locs, 'RowNames', brain_regs,...
    'VariableNames', {'n units responding to <=50% of locs', 'n responsive', ...
    '% responding to <=50% of locs'});

format short g
fprintf('\nNEURONS WITH ITEM RESPONSES:\n')
fprintf('Calculating one binomial test per brain region\n')
disp(table_binom_item)
fprintf('NEURONS WITH LOCATION RESPONSES:\n')
fprintf('Calculating one binomial test per brain region\n')
disp(table_binom_loc)
disp('---------------------------------------------------------------')
fprintf('\n NEURONS RESPONDING TO AT MOST 50%% OF ITMES:\n\n')
disp(table_half_item)
fprintf('\n NEURONS RESPONDING TO AT MOST 50%% OF LOCATIONS:\n\n')
disp(table_half_loc)

%% Difference between anterior and posterior hippocampal region
%we need all inds for each patient out of spikecell
pat_IDs = unique(response_info.subject);
%one column per patient, one row per unit
pat_inds = zeros(height(response_info), size(pat_IDs, 1));

for pt = 1:length(pat_IDs)
    % in column number sess, set all neurons belonging to session to 1
    pt_neurons = response_info.subject == pat_IDs(pt);
    pat_inds(pt_neurons, pt) = 1;
end

AH_inds = ~cellfun('isempty', strfind(response_info.channel_name, 'AH'));
MH_inds = ~cellfun('isempty', strfind(response_info.channel_name, 'MH'));
% these will hold all fractions. first column is AH, 2nd column MH
[frac_itemUnits, frac_locUnits] = deal(zeros(size(pat_inds, 2), 2));
item_lpc = response_info.lowest_p_item;
loc_lpc  = response_info.lowest_p_loc;
% for each column of sess_inds
for sj = 1: length(pat_IDs)
    AH_pt = AH_inds & pat_inds(:, sj);
    MH_pt = MH_inds & pat_inds(:, sj);
    frac_itemUnits(sj, 1) = sum(item_lpc(AH_pt)<=.001) / sum(AH_pt);
    frac_itemUnits(sj, 2) = sum(item_lpc(MH_pt)<=.001) / sum(MH_pt);
    frac_locUnits (sj, 1) = sum(loc_lpc(AH_pt)<=.001) / sum(AH_pt);
    frac_locUnits (sj, 2) = sum(loc_lpc(MH_pt)<=.001) / sum(MH_pt);
end


clrs = jet(length(pat_IDs));
ms = '*+x><^o+x><^o+x><';

legStr = cell(1, length(pat_IDs));

fu = {frac_itemUnits, frac_locUnits};
type_str = {'item', 'location'};

for types = 1:2
    fprintf('\n\n ---- %s NEURONS (anterior vs posterior hippocampus) ----\n', ...
            upper(type_str{types}));    
    figure('visible', plot_visibility);
    bh = bar(mean(fu{types}));
    bh.FaceColor= [.9,.9,.9];
    hold on
    er = errorbar([1,2], mean(fu{types}), std(fu{types}), std(fu{types}));
    er.Color = [0,0,0];
    er.LineStyle = "none";
    hold on
    hl = gobjects(1, length(pat_IDs));
    for pt = 1:length(pat_IDs)
        plot(1, fu{types}(pt, 1), ms(pt), 'color', clrs(pt,:))
        hl(pt) = plot(2, fu{types}(pt, 2), ms(pt), 'color',...
            clrs(pt,:), 'DisplayName', sprintf('P %.3i', pat_IDs(pt)));
        hold on
        plot(1:2, fu{types}(pt, :), 'color', clrs(pt,:),...
            'HandleVisibility','off');
        legStr{pt} = sprintf("Patient %i", pt); % have to anonymize:
        %legStr{i} = num2str(upid(pt));

    end
    legend(hl, legStr);
    test_y = max(max(fu{types}))*1.1;
    p_y = max(max(fu{types}))*1.14;
    yl = max(max(fu{types}))*1.2;
    yl_ = ylim;
    ylim([yl_(1), yl])
    plot([1,2], [test_y, test_y], 'k', 'HandleVisibility', 'off')
    % do test, print stats
    [~, p, ci, stats] = ttest(fu{types}(:,1), fu{types}(:,2));
    fprintf(sprintf('paired ttest %s: p=%.4f, t=%.3f, df=%i, ci=[%.2f %.2f]\n', ...
        type_str{types}, p, stats.tstat, stats.df, ci(1), ci(2)))
    test_str = 'paired ttest';

    if p<=.05
        text(1.5, p_y, sprintf('* p=%.4f', p))
        fprintf('\n(Significant difference between fractions of\n');
        fprintf(' responsive neurons in anterior vs. posterior\n hippocampal');
        fprintf(' target location)')
    else
        text(1.5, p_y, 'n.s.')
        fprintf('(No significant difference n.s. between fractions of');
        fprintf(' responsive \nneurons in anterior vs. posterior hippocampal');
        fprintf(' target location)\n')
    end
    xlabel('Anterior vs Posterior Hippocampus')
    ylabel(sprintf('fraction of %s neurons', type_str{types}))
    xlim([.2 2.8])
    xticks([1 2])
    xticklabels({'AH','PH'})
    legend box off
    title({sprintf('Are there more %s units in', type_str{types}),...
        sprintf('anterior or posterior hippocampus? (%s)', test_str)})
%     print(gcf, '-dpng', sprintf('outputpath/comparison_%sUnits_means_perPatient_%s',...
%         type_str{types}(1:4), strrep(test_str, ' ', '_')), '-r300');
    
    %% chi square test: are there more location neurons in PHC...
    %      than in other brain regions?
    
    %chi-square test 
    % https://de.mathworks.com/matlabcentral/answers/96572-how-can-i-perform-a-chi-square-test-to-determine-how-statistically-different-two-proportions-are-in
    %observed data. resp_per_unit is n_units x stim_pos x .05_.01_.001
    % we want to test PHC against Amy, Hipp and EC.
    [p_chi2_phc_vs_amy_hipp_ec,chi2stat] = deal(zeros(1,3));
    for i = 1:3
        %observed data
        neurons_phc   = response_info.brain_regs == 4;
        neurons_other = response_info.brain_regs == i;
        % n location neurons in PHC
        n1 = sum(response_info.lowest_p_loc <= .001 & neurons_phc); 
        % n neurons in PHC
        N1 = sum(neurons_phc);
        % n location neurons in other brain region
        n2 = sum(response_info.lowest_p_loc <= .001 & neurons_other);
        % n neurons in other brain region
        N2 = sum(neurons_other);
        disp([n1 N1 n2 N2])

        %pooled estimate of proportion
        p0 = (n1+n2)/(N1+N2);
        %expected counts under H0
        n10 = N1*p0; n20 = N2*p0;
        %chi-square test, by hand
        observed = [n1  N1-n1  n2  N2-n2 ];
        expected = [n10 N1-n10 n20 N2-n20];
        chi2stat(i) = sum((observed-expected).^2./expected);
        %next calculate the p value
        p_chi2_phc_vs_amy_hipp_ec(i) = 1-chi2cdf(chi2stat(i), 1);
    end

    fprintf('\nComparing the fraction of position responsive units in the PHC\n')
    fprintf('    to the other brain regions amy, hipp and EC.\n')
    fprintf('    pvals: %.2g (amy), %.2g (hipp), %.2g (EC)\n', ...
        p_chi2_phc_vs_amy_hipp_ec(1),p_chi2_phc_vs_amy_hipp_ec(2),...
         p_chi2_phc_vs_amy_hipp_ec(3))
    fprintf('    chi^2: %.2f (amy), %.2f (hipp), %.2f (EC)\n', chi2stat(1),...
        chi2stat(2), chi2stat(3))
end
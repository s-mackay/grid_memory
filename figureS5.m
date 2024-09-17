% figureS5(data_dir, read_from_excel)
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

function figureS5(data_dir, read_from_excel)

plot_visibility = 'on';
alpha_resp = .001;
alpha_nonresp = .001;
fs = 12;
if nargin <1;   data_dir = 'source_data';       end
if nargin <2;   read_from_excel = false;        end
rgns  = {' Amy', 'Hipp', '  EC', ' PHC'};
types = {'items', 'locs '};

if read_from_excel
    fname_excel = sprintf('%s/source_data_all_figures.xlsx',...
        data_dir);
    % 4 (regions) x 2 (item/loc) x 2 (resp/non-resp) 
    means_r_nr = cell(4, 2, 2);
    t_in = readtable(fname_excel, 'Sheet', 'figure_S5',...
        'PreserveVariableNames', true);
    % amygdala
    means_r_nr{1,1,1} = t_in.amy_item_resp(1:find(~isnan(t_in.amy_item_resp), 1, 'last'));
    means_r_nr{1,1,2} = t_in.amy_item_unresp(1:find(~isnan(t_in.amy_item_unresp), 1, 'last'));
    means_r_nr{1,2,1} = t_in.amy_loc_resp(1:find(~isnan(t_in.amy_loc_resp), 1, 'last'));
    means_r_nr{1,2,2} = t_in.amy_loc_unresp(1:find(~isnan(t_in.amy_loc_unresp), 1, 'last'));
    % hippocampus
    means_r_nr{2,1,1} = t_in.hipp_item_resp(1:find(~isnan(t_in.hipp_item_resp), 1, 'last'));
    means_r_nr{2,1,2} = t_in.hipp_item_unresp(1:find(~isnan(t_in.hipp_item_unresp), 1, 'last'));
    means_r_nr{2,2,1} = t_in.hipp_loc_resp(1:find(~isnan(t_in.hipp_loc_resp), 1, 'last'));
    means_r_nr{2,2,2} = t_in.hipp_loc_unresp(1:find(~isnan(t_in.hipp_loc_unresp), 1, 'last'));
    % ec
    means_r_nr{3,1,1} = t_in.ec_item_resp(1:find(~isnan(t_in.ec_item_resp), 1, 'last'));
    means_r_nr{3,1,2} = t_in.ec_item_unresp(1:find(~isnan(t_in.ec_item_unresp), 1, 'last'));
    means_r_nr{3,2,1} = t_in.ec_loc_resp(1:find(~isnan(t_in.ec_loc_resp), 1, 'last'));
    means_r_nr{3,2,2} = t_in.ec_loc_unresp(1:find(~isnan(t_in.ec_loc_unresp), 1, 'last'));
    % phc
    means_r_nr{4,1,1} = t_in.phc_item_resp(1:find(~isnan(t_in.phc_item_resp), 1, 'last'));
    means_r_nr{4,1,2} = t_in.phc_item_unresp(1:find(~isnan(t_in.phc_item_unresp), 1, 'last'));
    means_r_nr{4,2,1} = t_in.phc_loc_resp(1:find(~isnan(t_in.phc_loc_resp), 1, 'last'));
    means_r_nr{4,2,2} = t_in.phc_loc_unresp(1:find(~isnan(t_in.phc_loc_unresp), 1, 'last'));
    
    % also need info from figure S2
    response_info = readtable(sprintf('%s/source_data_all_figures.xlsx',...
        data_dir), 'Sheet', 'figure_S2');

else
    fname = sprintf('%s/mean_fr_recall_resp_unresp.mat', data_dir);
    load(fname, 'means_r_nr');
    load(sprintf('%s/response_info.mat', data_dir), 'response_info');
end



disp(means_r_nr)
nrows = size(means_r_nr, 1)+1;
ncols = size(means_r_nr, 2);
h4 = figure;
set(h4, 'Units', 'pixels', 'Position', [0 0 750 1000], 'visible', plot_visibility);

xtl = {['non-pref\newline', 'items\newline',...
    sprintf('%s \\geq %.3f','{\itP}', alpha_nonresp)],...
    ['preferred\newline', 'items\newline',...
    sprintf('%s < %.3f','{\itP}', alpha_resp)],...
    ['non-pref\newline', 'locations\newline',...
    sprintf('%s \\geq %.3f','{\itP}', alpha_nonresp)],...
    ['preferred\newline', 'locations\newline',...
    sprintf('%s < %.3f','{\itP}', alpha_resp)]};
subPl = 1;
ylims = [floor(min(reshape(cellfun(@(x)max(x(:)), means_r_nr), 1, []))),...
          ceil(max(reshape(cellfun(@(x)max(x(:)), means_r_nr), 1, [])))+1.5];

for i =1:nrows-1
    for j = 1:ncols
        subplot(nrows, ncols, subPl)
        data = [means_r_nr{i,j,1}; means_r_nr{i,j,2}];
        [p_r, ~, stats] = signrank(means_r_nr{i,j,1},means_r_nr{i,j,2},...
            'method', 'approximate');
        if i == 4 || j == 1
            fprintf('%s (%s): p=%.2g, z=%.2f, n=%i\n', ...
                rgns{i}, types{j}, p_r, stats.zval, ...
                numel(means_r_nr{i,j,1}));  
        end
        grouping = [ones(size(means_r_nr{i,j,1}))*2; ...
                    ones(size(means_r_nr{i,j,2}))];
        left_x = 1;
        right_x = 2;
        if size(data, 2)==1 %for n=1 the subplot size is too wide for some reason
            boxplot(repmat(data, 4,1), repmat(grouping, 4, 1))
            %will be overwritten now
        end
        boxplot(data, grouping,'Widths', .2, ...
            'positions', [left_x(1), right_x(1)], 'color', 'kk');     
        maxval = max(data)+.04*diff(ylims);
        xl = xlim;
        hold on
        plot([left_x, right_x], repmat(maxval, 2, 1), 'k')
        if p_r <=.001
            text(mean([left_x; right_x]), maxval, '***', ...
                'horizontalalignment', 'center', 'fontsize', fs);
        elseif p_r <=.01
            text(mean([left_x; right_x]), maxval, '**', ...
                'horizontalalignment', 'center', 'fontsize', fs);
        elseif p_r <=.05
            text(mean([left_x; right_x]), maxval, '*', ...
                'horizontalalignment', 'center', 'fontsize', fs);
        else
            text(mean([left_x; right_x]), maxval+.1*diff(ylims), 'n.s.', ...
                'horizontalalignment', 'center', 'fontsize', fs);
        end
        ylim(ylims)
        if i==nrows
%             set(gca, 'xticklabel', xtl, 'TickLabelInterpreter', 'default')
        else
            set(gca, 'xticklabel', '');
        end
        if j==1
            textXPos = xl(1)-diff(xl)/5;%8;
            text(textXPos, mean(ylims), rgns{i}, 'fontsize', fs+2,...
            'horizontalalignment', 'right', 'verticalalignment', 'middle');
            ylabel('fr (z-score)')
        end
        if i==1 && j==1
            title('item cells')
        elseif i==1 && j==2
            title('location cells')
        end
        set(gca, 'fontsize', fs)
        text(xl(2), mean(ylims)+.4*diff(ylims),...
        sprintf('%s = %i / %i ','{\itn}', size(data, 1)/2,...
            sum(response_info.brain_regs==i)), 'horizontalalignment',...
            'right','interpreter', 'tex', 'fontsize', fs)
        grid on
        subPl = subPl + 1;
    end
end

%put xticklabels in separate plots so they dont squish bottom subplots
for i = 1:2
    subplot(nrows, ncols, subPl)
    xlim(xl)
    text(1, 1, xtl{(i-1)*2+1}, 'interpreter', 'default',...
        'horizontalalignment', 'center', 'fontsize', fs)
    text(2, 1, xtl{(i-1)*2+2}, 'interpreter', 'default',...
        'horizontalalignment', 'center', 'fontsize', fs)
    ylim([0,1.2])
    axis off
    subPl = subPl + 1;
end

when=' 1s before tap';

try
    % this is a problem with some Matlab versions
    sgtitle({sprintf('Firing at recall,%s', when), ...
        sprintf('inclusion %s < %.3f @enc','{\itP}', alpha_resp)})
catch
    suptitle({sprintf('Firing at recall,%s', when), ...
        sprintf('inclusion %s < %.3f @enc','{\itP}', alpha_resp)})
end


fname=sprintf('%s/figures/FigureS5.png', data_dir);
print(gcf, '-dpng', fname, '-r600');
fname_eps2 = [fname(1:end-3), 'eps'];
print('-painters', gcf, '-depsc', fname_eps2);

end
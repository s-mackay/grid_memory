function [pos_true, neg_true, pos_dist, neg_dist, mask, handles, pvals_true] =...
    perm_ttest(condA, condB, nperms, mode, clusteralpha, alpha, doplot, ...
                                                   xaxis, plotcolors, tail)
% function [pos_true, neg_true, pos_dist, neg_dist, mask, handles, pvals_true] =...
%   perm_ttest(condA, condB, nperms, mode, clusteralpha, alpha, doplot, ...
%                                                xaxis, plotcolors, tail)
%
% Computes paired-ttests of condA vs condB and ranks the resulting
% t-value in a distribution of t-values obtained by
% label-shuffling, condA and condB have to be of the same size, 
% if condB is one value, a one-sample t-test of condA against this
% value is computed
%
% condA and condB are N x M arrays of the same size, N: trials/conditions, M: timepoints of
% e.g. inst firing rates or raw EEG data (mV), within unit of
% observation (subject/unit)
%
% nperms (int): number of random assignments of condition-labels
% mode (str): 'paired', 'onesample' or 'independent-amples'
% clusteralpha (d): threshold p-value to include a timepoint in a cluster
% alpha (d): exclude clusters with p>=alpha from result
% doplot (bool): plot data
% xaxis (array): x values for plot
% plotcolors (array): provide 3 columns with color values for each trace
% tail(str): default is 'both'
%
% Portions of this code were adapted from work by Thomas Reber, with 
% permission. Original work can be found at 
% https://github.com/rebrowski/neuralAdapatationInMTL
% -------------------------------------------------------------------------


    if ~exist('mode', 'var') || isempty(mode)
        mode = 'paired';   
        %mode = 'onesample';
        %mode = 'independent-samples';
    end
    
    if ~exist('alpha', 'var') || isempty(alpha)     
        alpha = 0.05;
    end
    if ~exist('clusteralpha', 'var') || isempty(clusteralpha)
        clusteralpha = 0.05;
    end

    if ~exist('nperms', 'var') || isempty(nperms)
        nperms = 1000;
    end

    if ~exist('doplot', 'var') || isempty(doplot)
        doplot = true;
    end
    if ~doplot 
        handles = [];
    end

    if ~exist('xaxis', 'var') || isempty(xaxis)
        xaxis = 1:size(condA,2);
    end
    
    if ~exist('plotcolors', 'var') || isempty(plotcolors)
        plotcolors =    [0.3639 0.5755 0.7484;
                         0.9153 0.2816 0.2878];
 
    end
    
    if size(condA,1) ~= size(condB,1) | size(condA,2) ~= size(condB,2) ...
        % we could pass just one value as condB in case of
        % onesample test against this value
            if size(condB,1) == 1 & size(condB,2) == 1
                condB = repmat(condB, size(condA,1), size(condA, ...
                                                          2));
                mode = 'onesample'
            end
                            
    end

    N = size(condA,1);    
    N2 = size(condB,1);        
    nSmpls = size(condA,2);
    nSmpls2 = size(condB,2);
    assert(nSmpls == nSmpls2);
    labels_true = [repmat(1,1,N) repmat(0,1,N2)];
    
    data = [condA; condB];
    nanDatA = sum(isnan(condA),2);
    nanDatB = sum(isnan(condB),2);
    if strcmp(mode, 'paired')
        nanDat = repmat(nanDatA | nanDatB, 2, 1);
    else
        nanDat = isnan(data);
    end

    if sum(nanDat) > 0
        if strcmp(mode, 'paired')
            temp = [condA, condB];
            nanDat = repmat(sum(isnan(temp),2) > 0, 2,1);
            clear temp
        end
        data = data(~nanDat,:);
        labels_true = labels_true(~nanDat);
        
        N = sum(labels_true == 1);
        N2 = sum(labels_true == 0);
    end
    
    %% do t-tests with true labels
    if ~strcmp(mode, 'independent-samples')
        diffs = data(labels_true==1,:)-data(labels_true==0,:);
        t = mean(diffs)./(std(diffs)./sqrt(N));
        p = 2 * tcdf(-abs(t), (N-1));
        t(find(p>=clusteralpha)) = 0;    
        tvals_true = t;
        pvals_true  = p;
    
    elseif strcmp(mode, 'independent-samples')
        if ~exist('tail', 'var')
            tail = 'both';
        end
        [h p ci stats] = ttest2(data(labels_true==1, :), ...
                                data(labels_true==0, :), ...
                                'tail', tail);
        tvals_true = stats.tstat;
        tvals_true(find(p>clusteralpha)) = 0;
        pvals_true = p;    
    end

    
    %% do t-tests with false labels nperms times
    tvals_false = zeros(nperms, nSmpls); 
    for i = 1:nperms

        %% shuflle labels
        % the labels of condA and B should be swapped randomly
        % within unit of observation in case of paired ttest:
        if ~strcmp(mode, 'independent-samples')
            labels_false = [repmat(1,1,floor(N/2)) repmat(0,1,ceil(N/2))];
            labels_false = labels_false(randperm(N));
            labels_false = [labels_false ~labels_false];

            diffs = data(labels_false==1,:)-data(labels_false==0,:);
            t = mean(diffs)./(std(diffs)./sqrt(N));
            p = 2 * tcdf(-abs(t), (N-1));
            t(find(p>=clusteralpha)) = 0;
            
            tvals_false(i,:) = t;
        else
            % the labels of condA and B should be swapped randomly
            % while keeping the Ns of the groups constant

            labels_false = labels_true(randperm(N+N2));
            [h p ci stats] = ttest2(data(labels_false==1, :), ...
                                    data(labels_false==0, :), ...
                                    'tail', tail);

            stats.tstat(p>clusteralpha) = 0;
            tvals_false(i,:) = stats.tstat;

        end
        
    end %nperms
    
%     if version < 9
%         matlabpool open
%     end
    
    %% get the clusters
    [pos_true neg_true] = get_clusters(tvals_true, 'all');
    pos_false.froms = [];
    neg_false.froms = [];
    pos_false.tos = [];
    neg_false.tos = [];
    pos_false.sizes = [];
    neg_false.sizes = [];
    pos_false.sumts = [];
    neg_false.sumts = [];

    fdnms = fieldnames(pos_false);
    
    %% concatenate clusters resulting form different permutationsn    
    for i=1:nperms
        [pos_false_ neg_false_] = get_clusters(tvals_false(i,:), ...
                                                    'maxsum');    
        for k = 1:length (fdnms)
            pos_false.(fdnms{k}) = [pos_false.(fdnms{k}) ...
                                pos_false_.(fdnms{k})];
            neg_false.(fdnms{k}) = [neg_false.(fdnms{k}) ...
                                neg_false_.(fdnms{k})];
        end
    end
    
    %% rank true clusters in distribution of false clusters
    pos_true.cluster_p = ones(length(pos_true.sumts),1);
    neg_true.cluster_p = ones(length(neg_true.sumts),1);
        
    % do distributions of false summmed tvals with size nperms
    pos_dist = []; 
    neg_dist = [];
    %% neg
    if length(neg_false.sumts) < nperms
        neg_dist = [neg_false.sumts zeros(1,nperms- ...
                                           length(neg_false.sumts))]; 
    else
        neg_dist = sort(neg_false.sumts, 'ascend'); % large
                                                    % negative sums
        neg_dist = neg_dist(1:nperms);
    end
    [neg_dist nidx] = sort(neg_dist, 'ascend'); 

    for i = 1:length(neg_true.sumts)
        [val rank] = min(abs(neg_dist - neg_true.sumts(i)));                
        neg_true.cluster_p(i) = rank/nperms;
    end

    %% pos
    if length(pos_false.sumts) < nperms
        pos_dist = [pos_false.sumts zeros(1,nperms- ...
                                           length(pos_false.sumts))]; 
    else
        pos_dist = sort(pos_false.sumts, 'descend');
        pos_dist = pos_dist(1:nperms);
    end
    pos_dist = sort(pos_dist, 'descend'); % large positive sums
                                          % should be at the beginnig
    for i = 1:length(pos_true.sumts)
        [val rank] = min(abs(pos_dist - pos_true.sumts(i)));                
        pos_true.cluster_p(i) = rank/nperms;
        pos_true.cluster_p(i) = rank/nperms;
    end
    
    
    %% remove insignificant clusters from result
    fdnms = fieldnames(pos_true);
    sign_p = find(pos_true.cluster_p<alpha);
    sign_n = find(neg_true.cluster_p<alpha);
    for i = 1:length(fdnms)
        pos_true.(fdnms{i}) = pos_true.(fdnms{i})(sign_p);
        neg_true.(fdnms{i}) = neg_true.(fdnms{i})(sign_n);
    end
   
    %% create a mask
    mask = zeros(1,nSmpls);
    for i = 1:length(pos_true.froms)
        mask(pos_true.froms(i):pos_true.tos(i)) = 1;
    end
    for i = 1:length(neg_true.tos)
        mask(neg_true.froms(i):neg_true.tos(i)) = -1;
    end

    %% create a plot
    if doplot
        ax1 = gca;
        lineh(1) = plot_signals(condA, xaxis, plotcolors(1,:));
        lineh(2) = plot_signals(condB, xaxis, plotcolors(2,:));
        yl = ylim;
        xlim([xaxis(1) xaxis(end)]);
        % makes some space for annotations
        yl = ([yl(1)-abs(diff(yl))* 0.1, yl(2)]);
        ylim(yl);
        grid on
        % do color-coded p-values at the bottom of the plot 
        pos = get(ax1, 'Position');
        
        set(ax1, 'Color', 'none');
        cbpos = [pos(1:3) pos(4)*0.05];
        ax2 = axes('Position', cbpos);
        %colormap(ax2, flipud(parula));
        colormap(ax2, flipud(hot));
        psh = imagesc(xaxis, 1, pvals_true); caxis([0 0.5]); axis off; axis xy;
        keyboard
        cbh = colorbar;
        set(cbh, 'Position', [cbpos(1)+cbpos(3)*1.005, ...
                            cbpos(2), ...
                            pos(3)*0.03, ...
                            cbpos(4)*2])
        axis off
        yl2 = ylim;
        hgt = abs(diff(ylim));
        % annotate sign clusters
        anotpos = [pos(1), ...
                   pos(2) + pos(4)*0.05, ...
                   pos(3), ...
                   pos(4)*0.03];
        ax3 = axes('Position', anotpos);
        set(ax3, 'color', 'none')
        xlim([xaxis(1) xaxis(end)]);
        ylim([0 1])
        axis off
        for c = 1:length(pos_true.froms)
            k = patch(xaxis([pos_true.froms(c) pos_true.tos(c)...
                         pos_true.tos(c) pos_true.froms(c)]), ...
                      [0 0 1 1], [0.7 0.7 0.7]);
            set(k, 'EdgeColor', 'none')
            x_ = xaxis(round(mean([pos_true.tos(c) pos_true.froms(c)])));
            text(x_,0.5, '*', 'HorizontalAlignment', 'center');
        end

        for c = 1:length(neg_true.froms)
            k = patch(xaxis([neg_true.froms(c) neg_true.tos(c)...
                         neg_true.tos(c) neg_true.froms(c)]), ...
                      [0 0 1 1], [0.7 0.7 0.7]);
            set(k, 'EdgeColor', 'none');
            x_ = xaxis(round(mean([neg_true.tos(c) neg_true.froms(c)])));

            text(x_, 0.5, '*', 'HorizontalAlignment', 'center'); 
        end

        handles = [ax1 ax2 ax3 cbh lineh];
        
    end    
end
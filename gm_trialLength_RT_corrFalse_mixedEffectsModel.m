function gm_trialLength_RT_corrFalse_mixedEffectsModel(data_dir,...
                                                 read_from_excel, do_plots)

% Fitting a linear mixed effects model. 
% reference:
% https://www.mathworks.com/help/stats/linear-mixed-effects-model-workflow.html
% 
% INPUTS:
%         data_dir (str): path to source data
% read_from_excel (bool): row number within spiketable of unit to plot
%        do_plots (bool): produce some plots
%  
% Example: set size
% 
% setSize ~ corrFalse:                We are modeling the variable 
%                                       setSize as a function of the fixed 
%                                       effect of memory performance, coded
%                                       as corrFalse
% (corrFalse | patientID):            Models the variability in the
%                                       intercept and the effect of
%                                       corrFalse across patients
% (corrFalse | patientID:sessionInd): Allows for variablility in the effect
%                                       of corrFalse across different
%                                       patients and different sessions
%                                       within patients


% if you want to produce plots and save them, define this
% outputDir = '';

%% load data
excel_fname = sprintf('%s/source_data_mixed_effects.xlsx', data_dir);
if read_from_excel
    ds = readtable(excel_fname);
else
    load(strrep(excel_fname, 'xlsx', 'mat'), 'ds');
end

% converting some of the variables to categorical
for var = {'patientID', 'sessionInd', 'corrFalse'}
    ds.(var{1}) = categorical(ds.(var{1}));
end



fprintf('\n<strong>----------- SET SIZE --------</strong>')
ss_model = fitlme(ds,     'setSize ~ corrFalse + (corrFalse | patientID) + (corrFalse | patientID:sessionInd)');
disp(ss_model)
fprintf("\n<strong>------- TRIAL DURATION --------</strong>")
tl_model = fitlme(ds, 'trialLength ~ corrFalse + (corrFalse | patientID) + (corrFalse | patientID:sessionInd)');
disp(tl_model)
fprintf('\n<strong>------ REACTION TIME --------</strong>')
rt_model = fitlme(ds,          'rt ~ corrFalse + (corrFalse | patientID) + (corrFalse | patientID:sessionInd)');
disp(rt_model)

if do_plots
    %% plotting average set sizes per patient as an example
    
    % choosing set size as the variable with the most consistent difference
    % get unique patient IDs as double (conversion to double just so they
    % are not sorted alphabetically)
    upid = double(unique(ds.patientID));
    clrs = jet(length(upid));
    ms = '*+x><^o+x><^o+x><';
    handles = nan(size(upid)); % handles for legend

    figure
    legStr = cell(1, length(upid));

    for i = 1:length(upid)
        sesss = find(double(ds.patientID)==upid(i));
        h = plot(1:2, [ds.setSize(sesss(1:2:end)),...
            ds.setSize(sesss(2:2:end))],...
            ms(i), 'color', clrs(i,:), 'DisplayName', sprintf('P %.3i', upid(i)));
        handles(i) = h(1);
        hold on
        plot(1:2, [mean(ds.setSize(sesss(1:2:end))),...
            mean(ds.trialLength(sesss(2:2:end)))]);
        legStr{i} = num2str(upid(i));

    end
    legend(handles, legStr);
    xlim([.5 2.5])
    xticks([1 2])
    xticklabels({'not remembered','remembered'})
    title('Average set size per session (markers) and patient (lines)')
    
    %% plotting residuals for each of the 3 variables
    %  trial length, set size and reaction time 
    width = 20; height = 9;
    h1 = figure;
    set(h1, 'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    legStr = cell(1, length(upid));
    for i = 1:length(upid)
        subplot(2,length(upid),i);
        sesss = find(double(ds.patientID)==upid(i));
        plot(1:2, [ds.trialLength(sesss(1:2:end)), ds.trialLength(sesss(2:2:end))], '-oc')
        title(sprintf('P%i', ds.patientID(sesss(1))));
        ylim([min(ds.trialLength), max(ds.trialLength)])
        hold on
        plot(1:2, [mean(ds.trialLength(sesss(1:2:end))),...
            mean(ds.trialLength(sesss(2:2:end)))], 'b');
        xlim([.5 2.5])
        xticklabels({'c','f'})
        legStr{i} = num2str(upid(i));
        if i==1; ylabel('mean trial duration [s]'); end
    end   
    subplot(2, length(upid), i+1:2*length(upid))
    plotResiduals(tl_model, 'fitted')
    try
        sgtitle('Trial length')
    catch
        suptitle('Trial length')
    end
%     print(h1, '-dpng', sprintf('%s/lmeModel_trialLength', outputDir));


    h2 = figure;
    set(h2, 'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    for i = 1:length(upid)
        subplot(2,length(upid),i);
        sesss = find(double(ds.patientID)==upid(i));
        plot(1:2, [ds.setSize(sesss(1:2:end)), ds.setSize(sesss(2:2:end))], '-oc')
        title( sprintf('P%i', ds.patientID(sesss(1))));
        ylim([min(ds.setSize), max(ds.setSize)])
        hold on
        plot(1:2, [mean(ds.setSize(sesss(1:2:end))),...
            mean(ds.setSize(sesss(2:2:end)))], 'b');
        xlim([.5 2.5])
        xticklabels({'c','f'})
        if i==1; ylabel('set size'); end
    end   
    subplot(2, length(upid), i+1:2*length(upid))
    plotResiduals(ss_model, 'fitted')
    try
        sgtitle('Set size')
    catch
        suptitle('Set size')
    end
%     print(gcf, '-dpng', sprintf('%s/lmeModel_setSize', outputDir));


    h3 = figure;
    set(h3, 'PaperSize', [width height], 'PaperPosition', [0 0 width height]);
    for i = 1:length(upid)
        subplot(2,length(upid),i);
        sesss = find(double(ds.patientID)==upid(i));
        plot(1:2, [ds.rt(sesss(1:2:end)), ds.rt(sesss(2:2:end))], '-oc')
        title( sprintf('P%i', ds.patientID(sesss(1))));
        ylim([min(ds.rt), max(ds.rt)])
        hold on
        plot(1:2, [mean(ds.rt(sesss(1:2:end))),...
            mean(ds.rt(sesss(2:2:end)))], 'b');
        xlim([.5 2.5])
        xticklabels({'c','f'})
        if i==1; ylabel('reaction time [s]'); end
    end   
    subplot(2, length(upid), i+1:2*length(upid))
    plotResiduals(rt_model, 'fitted')
    try
        sgtitle('Reaction time')
    catch
        suptitle('Reaction time')
    end
%     print(gcf, '-dpng', sprintf('%s/lmeModel_RT', outputDir));
end

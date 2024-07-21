function fig3_source_to_xlsx(means, sems, x, pos_true, neg_true, ...
    source_data_fname, sheet_name)

% wrote this in an effort to fit data for all figures into one excel file
% with one sheet per figure. this is the fig 3 part.

%means is a cell wth a 2-row-array for each panel. same with sems. x just
%has one row because it's shared.
%each row gets its own label because you can't really fit several 2d arrays
%into one excel sheet well.

% rem, remembered
% forg, forgotten

data = {};
for reg = 1:5
    for col = 1:3
        if reg < 5 || col == 2
            % solid traces, means
            var_name = sprintf('means_rem_%i_%i', reg, col);
            if ~isempty(means{reg, col})
                var_data_m1 = means{reg, col}(1, :);
                var_data_m2 = means{reg, col}(2, :);
                var_data_s1 = sems{reg, col}(1, :);
                var_data_s2 = sems{reg, col}(2, :);
            end
            row = [{var_name}, {var_data_m1}];
            data = [data; row];

            var_name = sprintf('means_forg_%i_%i', reg, col);
%             var_data = means{reg, col}(2, :);
            row = [{var_name}, {var_data_m2}];
            data = [data; row];

            %shaded areas, SEMs
            var_name = sprintf('SEMs_rem_%i_%i', reg, col);
%             var_data = sems{reg, col}(1, :);
            row = [{var_name}, {var_data_s1}];
            data = [data; row];

            var_name = sprintf('SEMs_forg_%i_%i', reg, col);
%             var_data = sems{reg, col}(2, :);
            row = [{var_name}, {var_data_s2}];
            data = [data; row];

            % x values
            var_name = sprintf('x_%i', reg);
            var_data = x{reg, col};
            row = [{var_name}, {var_data}];
            data = [data; row];

            var_name = sprintf('pos_true_froms_%i_%i', reg, col);
            var_data = pos_true{reg, col}.froms;
            row = [{var_name}, {var_data}];
            data = [data; row];

            var_name = sprintf('pos_true_tos_%i_%i', reg, col);
            var_data = pos_true{reg, col}.tos;
            row = [{var_name}, {var_data}];
            data = [data; row];
            
            var_name = sprintf('pos_true_p_%i_%i', reg, col);
            var_data = pos_true{reg, col}.cluster_p;
            row = [{var_name}, {var_data}];
            data = [data; row];

            var_name = sprintf('neg_true_froms_%i_%i', reg, col);
            var_data = neg_true{reg, col}.froms;
            row = [{var_name}, {var_data}];
            data = [data; row];

            var_name = sprintf('neg_true_tos_%i_%i', reg, col);
            var_data = neg_true{reg, col}.tos;
            row = [{var_name}, {var_data}];
            data = [data; row];
            
            var_name = sprintf('neg_true_p_%i_%i', reg, col);
            var_data = neg_true{reg, col}.cluster_p;
            row = [{var_name}, {var_data}];
            data = [data; row];
        end
    end
end

%delete the old sheet so that there is no data remaining. Or rather, since
%deleting is turning out to be tricky (if you dont have Excel installed), 
%we're doing the ugly workaround of writing 1000 x 500 empty cells into
%the sheet
emptyData = cell(1000, 3000);
% Write the empty cell array to the specified sheet
writecell(emptyData, source_data_fname, 'Sheet', sheet_name);
    
%write actual data and save to excel sheet
writecell(data, source_data_fname, 'Sheet', sheet_name, 'Range', 'A1');

end
function [struct_out] = read_convolved_fr_excel(path)


x_all = readmatrix(path, 'Sheet', 'x_all', 'NumHeaderLines', 0);

%the shape of the data will be [3681, size(x_all, 2)]
first_cell = 'A1';
cn = columnNumberToExcelColumn(size(x_all, 2));
last_cell =  sprintf('%s%i', cn, 3681);

dataRange = sprintf('%s:%s', first_cell, last_cell);
%the DataRange argument makes sure the any leading empty rows are filled
%with nans.
y_c = readmatrix(path, 'Sheet', 'y_unit_time_correct',...
    'NumHeaderLines', 0, 'DataRange', dataRange);
y_f = readmatrix(path, 'Sheet', 'y_unit_time_false',...
    'NumHeaderLines', 0, 'DataRange', dataRange);

% when writing to excel, if the last rows are all nan, they don't get
% written. This will happen when the last cells are not response eliciting
% and the file contains valid data only for reponse elicitinv cells.
% Add the last NaN rows back in.
% if size(y_c, 1) < 3681
%     y_c(size(y_c, 1)+1 : 3681, :) = NaN;
%     y_f(size(y_f, 1)+1 : 3681, :) = NaN;
% end

y_unit_time_cf = cat(3, y_c, y_f);

struct_out.x_all = x_all;
struct_out.y_unit_time_cf = y_unit_time_cf;


function columnName = columnNumberToExcelColumn(n)
    dividend = n;
    columnName = '';
    while dividend > 0
        modulo = mod(dividend - 1, 26) + 1;
        columnName = [char('A' + modulo - 1), columnName];
        dividend = floor((dividend - modulo) / 26);
    end
end
end
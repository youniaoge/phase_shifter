input_file='4096_s21_deg_test.csv';
data = readtable(input_file);
% Keep the first column and then select only columns with even indices (2, 4, 6, etc.)
data_filtered = data(:, [1, 2:2:end]);
num_cols = width(data_filtered);  % Get the number of columns
new_column_names = cell(1, num_cols);  % Initialize an empty cell array for column names
% First column is "freq"
new_column_names{1} = 'freq';
% Generate the rest of the names in the binary format
for i = 2:num_cols
    binary_str = dec2bin(i-2, 12);  % Convert to binary string with 12 bits
    new_column_names{i} = sprintf('%s %s', binary_str(1:6), binary_str(7:end));
end
% Update the table column names
data_filtered.Properties.VariableNames = new_column_names;
% Create the output file name
[~, name, ext] = fileparts(input_file);
output_file = strcat(name, '_processed', ext);
    
% Write the updated data to a new Excel file with the new column names
writetable(data_filtered, output_file);
disp(['Modified file saved as: ', output_file]);

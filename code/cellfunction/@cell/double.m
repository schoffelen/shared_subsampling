function output = double(input)

output = cellfun(@double, input, 'uniformoutput', false);
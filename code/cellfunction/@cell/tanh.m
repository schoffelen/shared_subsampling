function output = tanh(input)

output = cellfun(@tanh, input, 'uniformoutput', false);
function output = hilbert(input)

% output the hilbert transform (per row)
output = cell(size(input));
for k = 1:numel(input)
  output{k} = hilbert(input{k}.').';
end

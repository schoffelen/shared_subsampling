function update_bar_colors(b, colors)
% updates the colors of a bar plot

% b : handle of the bar plot
% colors : colors to use, has to match data size

for ii = 1:length(b)
    b(ii).FaceColor = colors(ii, :);
%     b(ii).EdgeColor = "none";

    % also change the width
    b(ii).BarWidth = 1.0;

end

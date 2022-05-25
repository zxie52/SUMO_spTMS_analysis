function shadedError(x, y, dy, Color)
% x for the x axis value
% y for the y axis value
% dy for the error for the y value
% Color: the string for the color

x = x';
y = y';
dy = dy';

fill([x;flipud(x)],[y-dy;flipud(y+dy)], Color, 'FaceAlpha', 0.3, 'linestyle','none');
% line(x,y, 'Color', Color, 'LineWidth', 5);
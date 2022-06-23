function ersp_fig_setup()
% Set up title font size
titleFontSize = 25;
axisFontSize = 22;

% Change the y axis to normal
set(gca, 'YDir', 'normal');
% Not Displaying the figure in matlab gui
% set(0,'DefaultFigureVisible','on');
% Opening the figure in the full size
set(gcf, 'Position', get(0, 'Screensize'));
colorbar;
caxis([-20 10]); % Set up the range for power(-10, 20)
xlabel("Time across the full trial",'Fontsize', axisFontSize);
ylabel("frequency",'Fontsize', axisFontSize);
end
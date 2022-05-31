function [titleFontSize, axisFontSize, textFontSize] = erp_plotting_setup()
figure;
hold on;

% set up the gcf and gca before plotting
set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'linewidth',1);
set(gca, 'Fontsize', 12);

% Set up title font size
titleFontSize = 28;
axisFontSize = 20;
textFontSize = 16;
end
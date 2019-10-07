
F1 = mkfig(1, 0.75, 1.5, false);clf

ax4 = tight_subplot(1, 1, 0.05, 0.01, 0.01, false);


h = colorbar(ax4, 'Location', 'EastOutside', 'Position', [0.03194 0.05 0.241 0.9]);
set(h, 'YAxisLocation','right')

colormap gray
% convert the threshold to nanometers

thresh_color = thresh*(7/20)*(1000/1)*1.005; %should give 22
caxis([-thresh_color, thresh_color])
% set(h, 'Ticks', [-20, -10, 0, 10, 20])

text(.75, .3, '$Z$ [nm]', 'rot', 90, 'interpreter', 'latex',...
    'Units', 'normalized', 'FontSize', 12)

ax4.Visible = 'off';

save_fig(F1, 'latex/figures/colorbar', false)
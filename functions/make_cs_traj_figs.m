function [figs, axs] = make_cs_traj_figs(figbase, how_many)
% [figs, axs] = make_cs_traj_figs(figbase, how_many)
% 
% The order is:
%   axs =  ax_ze, how_many >=1
%   axs += ax_ze, how_many >=2
%   axs += ax_y,  how_many >=3
%   axs += ax_y,  how_many >=4

  Fig_uz = figure(20+figbase); clf
  ax1 = gca();
  figs = {Fig_uz};
  axs = {ax1};
  if how_many >1
    Fig_ze = figure(30+figbase); clf
    ax2 = gca();
    figs = {Fig_uz, Fig_ze};
    axs = {ax1, ax2};
  end
  if how_many > 2
    Fig_x = figure(40+figbase); clf
    ax3 = gca();
    figs = {Fig_uz, Fig_ze, Fig_x};
    axs = {ax1, ax2, ax3};
  end
  
  if how_many > 3
    Fig_y = figure(50+figbase); clf
    ax4 = gca();
    figs = {Fig_uz, Fig_ze, Fig_x, Fig_y};
    axs = {ax1, ax2, ax3, ax4};
  end
end
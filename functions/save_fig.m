

function save_fig(fig, fpath_no_ext, for_presentation)
% save a figure as either pdf (for_presentation=true) or svg
% (for_presentation=false). In both cases, the background color will be set to
% the color defined in the function fig_color_presentation() and fig_color()
% respectively. 
%
% for_presentation is optional, default is false.

  if nargin <3
    for_presentation = false;
  end

  if for_presentation
    ext = '.pdf';
    bg = fig_color_presentation();
  else
    ext = '.svg';
    bg = [1, 1, 1]; %fig_color();
  end

  % Check to see if fpath has a file extension:
  % Cant see how to match a period. Docs say \x{N} will match character with 
  % hexedecimal of value N. A period is 2E in Hex...
  expr = '\x{2E}[A-Za-z]{4,4}$';
  idx = regexp(fpath_no_ext, expr);
  if ~isempty(idx)
    ext_old = fpath_no_ext(idx:end);
    fpath = [fpath_no_ext(1:idx-1), ext];
    warning('File extension detected. Replacing %s with %s\n', ext_old, ext)
  else
    fpath = [fpath_no_ext, ext];
  end

  
  
  % set the figure to have the right background and not be insane.  
  set(fig, 'InvertHardcopy', 'off', 'Color', bg);

  % set the axes to have the right background. Not sure if the figure can have
  % children besides axes...
  childs = get(fig, 'Children');
  for k=1:length(childs)
    if isa(childs(k), 'matlab.graphics.axis.Axes') || isa(childs(k), 'matlab.graphics.illustration.Legend')
      set(childs(k), 'Color', bg);
    else
      fprintf('Dont know what to do with child of type %s\n', class(childs(k)));
    end
  end

  if strcmp(ext, '.pdf')
    print(fig, '-dpdf', fpath)
  elseif strcmp(ext, '.svg')
    saveas(fig, fpath)
  else
    warning('Something went wrong. Dont know what to do with extenstion %s\n', ext)
  end



end


function hands = mimo_bode_mag(HH, freq_s, ha, varargin)
    nu = size(HH, 2);
    no = size(HH, 1);
    
    if isempty(ha)
        [ha] = tight_subplot(no, nu, .01, [.061, 0.03], [.065, .02]);
        ha = reshape(ha', [], no)';
    end
    
    hands = gobjects(no, nu);
    ylm = [-95, 30];
    coords = {'X', 'Y', 'Z'};
    for ny = 1:size(HH, 1)
        for nu=1:size(HH, 2)
            try
                hands(ny, nu) = frf_bode_mag(HH(ny,nu), freq_s, ha(ny,nu), 'Hz', varargin{:});
            catch
                keyboard
            end
            set(ha(ny,nu), 'YLim', ylm);
            if nu >1
                set(ha(ny,nu), 'YLabel', [], 'YTickLabel', []);
            else
                ylabel(ha(ny,nu), sprintf('Output %s (Mag [dB])', coords{ny}));
            end
            if ny ==1
                title(ha(ny,nu), sprintf('Input %s', coords{nu}));
            end
            if ny <3
                set(ha(ny,nu),  'XLabel', [], 'XTickLabel', []);
            end
            set(ha(ny,nu), 'XTick', [10,100,1000,10000]);
            
        end
    end
  
end
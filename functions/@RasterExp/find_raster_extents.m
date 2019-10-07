% [kmin, kmax] = find_raster_extents(rdat)
%
% Find the x-axis extent of the raster data. kmin and kmax are the indeces
% where the raster data actually exists. 

function [kmin, kmax] = find_raster_extents(rdat)


    kmin = 1;
    kmax = size(rdat.pixelifsampled, 2);

    kstarts = zeros(kmax,1);
    kends = zeros(kmax,1);
    
    for i_row=1:rdat.npix
       inds_not = find(rdat.pixelifsampled(i_row,:) ~= 0);
       kstarts(i_row) = inds_not(1);
       kends(i_row) = inds_not(end);
    end
    
    kmin = max(kstarts);
    kmax = min(kends);
%     for i=1:rdat.npix
%        if mean(rdat.pixelifsampled(:,i)) < 0.5
%            kmin = i+1;
%        else
%            break
%        end
%     end
% 
%     for i = kmin+1:rdat.npix
%        if mean(rdat.pixelifsampled(:,i)) < 0.5
%            kmax = i-1;
%            break
%        end
%     end




end


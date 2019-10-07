% [pix_mask_mat, pix_idx] = mu_path_mask(mupath_len, n, m, samplingRatio, repeat_sampling)
%
% Return Mu path sampling mask
%
% Arguments
% ----------
%   mupath_len : size (in pixels) of mu path pattern
%   n,m : image size
%   samplingRatio : (Percentage) - total pixels to sample = samplingRatio*n*m
%   repeat_sampling : (true|false) If true, repeat sampling allowed. 
%                     If false: repeat sampling not allowed.
% Returns
% --------
%   pix_idx : a vector of sampled pixel indeces.
%   pix_mask_mat : an n by m mu path pattern mask, which contains 1's in the
%                  mu-path areas, and zeros elsewhere.


function [pix_mask, pix_idx, npaths] = mu_path_mask(mupath_len, n, m, samplingRatio, repeat_sampling)
   
    if nargin < 5
        repeat_sampling = false;
    end

    pix_mask = zeros(n,m);
    npaths = 0;
    while (sum(sum(pix_mask))<samplingRatio*n*m)
        
        if repeat_sampling
            
            rand_i = randi(n);
            rand_j = randi([2-mupath_len m]);            
            pix_mask(rand_i,max(rand_j,1):min(rand_j+mupath_len-1,n)) = 1;            
            
        else            
            rand_i = randi(n);
            rand_j = randi(m-mupath_len+1);
            if sum(pix_mask(rand_i,rand_j:rand_j+mupath_len-1)) < 0.5
                pix_mask(rand_i,rand_j:rand_j+mupath_len-1) = 1;
                npaths = npaths + 1;
            end            
        end
    end
    
    pix_idx = find(CsTools.pixmat2vec(pix_mask)>0.5);
   
end





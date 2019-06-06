function im = interp_missing(self, im)
    mask = self.pix_mask;
    for k=1:size(mask, 2)
        sm = sum(mask(:, k));
        if sm > 0
            break
        end
    end
    idx_start = k;
    for k=size(mask, 2):-1:1
        sm = sum(mask(:, k));
        if sm > 0
            break
        end
    end
    idx_end = k;
    
    mask_slice = mask(:, idx_start:idx_end);
    im_slice = im(:, idx_start:idx_end);
    n = size(im_slice, 1);
    
    
    for k=1:size(im_slice, 1)
        idx_row_have = find(mask_slice(k, :) == 1);
        row_have = im_slice(k, idx_row_have);
        row_slice = interp1(idx_row_have, row_have, idx_start:idx_end,...
            'nearest', 'extrap');
        
        im(k, idx_start:idx_end) = row_slice;
    end
    
    
end
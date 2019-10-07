
function data_writer(data, fname, numeric_fmt)
    if nargin < 3
        numeric_fmt = '%.3g';
    end
   assert(isa(data, 'cell')); 
   fid = fopen(fname, 'w+');
   
   
   for k=1:length(data)
      dat_k = data{k};
      assert(isstruct(dat_k));
      write_struct(dat_k, fid, numeric_fmt)
   end
    
end


function write_struct(data, fid, fmt)
    
   flds = fieldnames(data);
   
   for k=1:length(flds)
      fld = flds{k};
      dat_k =  data.(fld);
      if isnumeric(dat_k) || islogical(dat_k)
          write_number(dat_k, fld, fid, fmt)
      elseif ischar(dat_k)
          write_string(dat_k, fld, fid);
      elseif isstruct(dat_k)
          write_struct(dat_k, fid, fmt)
      else
          error(['I dont know what to do with type %s, in field %s. ',...
              'Expected numeric or struct.'], class(dat_k), fld)
      end
   end
    
end

function write_string(data, name, fid)
    
   fprintf(fid, '%%<*%s>\n', name);
   fprintf(fid, '%s\n', data);
   fprintf(fid, '%%</%s>\n', name);
    
end

function write_number(data, name, fid, fmt)
    
    fmt_str = sprintf('%s\n', fmt);
   fprintf(fid, '%%<*%s>\n', name);
   fprintf(fid, fmt_str, data);
   fprintf(fid, '%%</%s>\n', name);
    
end
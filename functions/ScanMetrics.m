classdef ScanMetrics < handle
  properties
    damage;
    quality;
    psnr;
    ssim;
    time;
    rate;
    coverage;
    type;
    ypix;
  end
  
  methods
    function self = ScanMetrics(varargin)
        self.type = {};
    end
    function append_metrics(self, varargin)
      p = inputParser();
      p.addParameter('damage', Inf);
      p.addParameter('quality', Inf);
      p.addParameter('psnr', 0);
      p.addParameter('ssim', 0);
      p.addParameter('rate', 0);
      p.addParameter('time', 0);
      p.addParameter('type', 'na');
      p.addParameter('coverage', 0);
      p.addParameter('ypix', 0);
      
      p.parse(varargin{:});
      
      self.damage(end+1) = p.Results.damage;
      self.quality(end+1) = p.Results.quality;
      self.psnr(end+1) = p.Results.psnr;
      self.ssim(end+1) = p.Results.ssim;
      self.rate(end+1) = p.Results.rate;
      self.time(end+1) = p.Results.time;
      self.coverage(end+1) = p.Results.coverage;
      self.type{end+1} = p.Results.type;
      self.ypix(end+1) = p.Results.ypix;
      
    end
  end

end
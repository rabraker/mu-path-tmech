classdef frf2ss

    properties
       frf;
       freqs_rad;
       frf_even
       freqs_rad_even;
       impulse_experimental;
       Nd; 
       opts;
       sigmas;
       ss_full;
       ni;
       no;
    end
    
    methods
        function self = frf2ss(frf, freqs_rad, Nd, opts)
        % self = frf2ss(frf, freqs_rad, Nd, opts)
        % This should be a 3d array with No
        %  rows, Ni columns and length(t) pages.
            sz = size(frf);
            if length(sz) == 2
              frf = reshape(frf, 1, 1, []);
              ni = 1;
              no = 1;
            else
              ni = sz(2); % ni cols
              no = sz(1); % no rows
              
            end
            self.ni = ni;
            self.no = no;
            
            self.opts = opts;
            self.freqs_rad = freqs_rad;
            self.frf = frf;
            
            self.Nd = Nd;
            t = [0:opts.Ts:opts.impulse_length_sec]';
            N = length(t);

            % For the fourier transform, N we need evenly spaced frequencies, where N
            % is length of impulse response.
            freqs_even = [0:(1/(N*opts.Ts)):floor((1/opts.Ts)*0.5)]'; %Frequencies up to nyquist.
            ws_even    = freqs_even*2*pi;
            self.freqs_rad_even = ws_even;
            
            self.impulse_experimental = []; %zeros(ni, no, length(ws_even));
            self.frf_even = zeros(ni, no, length(ws_even));
            % Interpolate the FRF we have onto the evenly spaced frequencies.
            K_freq_max = find(ws_even <= freqs_rad(end), 1, 'last');
            K_freq_min = find(ws_even <= freqs_rad(1), 1, 'last');
            gdelay = (exp(j*ws_even*opts.Ts)).^Nd;
            
            for k_in = 1:ni  % no cols
              for j_out = 1:no % ni rows
                frf_iter = squeeze(frf(j_out, k_in, :));
                [frf_iter, freqs_rad] = monotonicFRF(frf_iter, freqs_rad);
                frf_iter = interp1(freqs_rad, frf_iter, ws_even(1:K_freq_max), 'spline');
                frf_iter = frf_iter.*gdelay(1:K_freq_max);
                
                % Window the FRF for frequencies past what we have, because
                % interpolation becomes inaccurate.
                frf_iter(K_freq_max:length(ws_even)) = 0;

                self.frf_even(j_out, k_in, :) = frf_iter;
                % Get the impulse response
                self.impulse_experimental(j_out, k_in,:) = frf2impulse(frf_iter);
                
              end
            end

            % Do the ERA on the impulse response
            [Ad, Bd, Cd, sigmas] = impulse2ss(self.impulse_experimental, opts.r, opts.s);
            self.ss_full.A = Ad;
            self.ss_full.B = Bd;
            self.ss_full.C = Cd;
            self.sigmas = diag(sigmas);
        end
        
        function sys = realize(self, Ns)
            A = self.ss_full.A(1:Ns,1:Ns);
            B = self.ss_full.B(1:Ns, 1:self.ni);
            C = self.ss_full.C(1:self.no, 1:Ns);
            D = zeros(self.no, self.ni);
            sys = ss(A, B, C, D, self.opts.Ts, 'InputDelay', round(self.Nd));
            
        end
        
    end
end














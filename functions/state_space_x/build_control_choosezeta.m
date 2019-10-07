function [Q, R0, S, P_x] = build_control_choosezeta(sys_recyc, can_cntrl, tol)
% [Q, R0, S] = build_control(sys_recyc, can_cntrl)
% Construct the quadratic cost matrices using the root-locus
% technique (fictitious zeros). 
% 
% Inputs
% ------
%  sys_recyc : the deltaUk ss system
%  can_control : an instance of CanonControlParams_XX.m
%
% tol : tolerance for ctrbf, which will be used if pole placment fails.
% Unfortunately, the default tolerance is often not sufficient to promote a
% decomposition for all cases where pole placement fails.
  
  if ~exist('tol', 'var')
      tol = 1e-8;
  end
  
  if sys_recyc.InputDelay ~= 0 || sys_recyc.IOdelay ~=0
    error(['All Delays of sys_obs must be part of the states, not ' ...
           'included as a property']);
  end
  
  
  P_x  = getCharDes(sys_recyc, can_cntrl.gam_s, can_cntrl.pint,...
      can_cntrl.zeta_s, can_cntrl.rho_s, can_cntrl.rad);

  try
      [Chat, Dhat] = place_zeros(sys_recyc, P_x);
  catch
      [Chat, Dhat] = place_zeros_kalman_decomp_local(sys_recyc, P_x, tol);
  end

  
  Q = Chat'*Chat;
  S = Chat'*Dhat;
  R0 = Dhat'*Dhat;
  

end


function [Chat, Dhat] = place_zeros_kalman_decomp_local(sys, P_x, tol)
% The trouble with controllability is that the recycled system (delta Uk
% puts a zero at z=0, but we also have delay poles there, so pole-zero
% cancellation. This function tries to decompose the controllable and
% uncontrollable subspaces.

    ns = size(sys.b, 1);
    for i=1:100
        % We dont care about Cc.
        [ab, bb, cb, T, L] = ctrbf(sys.a, sys.b, sys.c, tol);
        nc = sum(L);
        if nc < ns
            fprintf(['Controllable/uncontrollable decomposition succeeded ',...
                'with tol=%f\n'], tol)
            break
        end
        tol = tol * 10;
    end
    nnc = length(bb) - nc;
    Cnc = cb(1:nnc);
    % In this specialized situation, we expect that 
    % 1) there was one uncollable mode and 
    % 2) The mode was at z==0. Lets check.
    Anc = ab(1:nnc, 1:nnc);
    Bnc = bb(1:nnc);
    
    if nnc ~= 1
        error(['Extracted more than one uncolloable mode. ',...
            'Dont know what to do with that. size(Bnc) = %d'], size(Bnc,1));
    end
    
    if abs(eig(Anc)) > 1e-8
        error(['Eigenvalue of uncolloable mode does not appear to at z==0 ',...
            'Mode is at %f'], eig(Anc));
    end
    
    
    Ac = ab(nnc+1:end, nnc+1:end);
    Bc = bb(nnc+1:end);
    Cc = cb(nnc+1:end);
    P_x_c = P_x(P_x ~= 0);

    if length(P_x_c) ~= size(Bc, 1)
        error('Size of Bc and P_x_c do not match.');
    end
    sys_tmp = ss(Ac, Bc, Cc, sys.D, sys.Ts);
    [Chat_bar, Dhat] = place_zeros(sys_tmp, P_x_c);
%     Kc = place(Ac, Bc, P_x_c);
    Chat = [0, Chat_bar] * T;
    
end
function [phase_ini, amp, err ] = algo3ph_zero_it_par(k, contr, phase_ini, amp, err, x_phase0, use_invA, end_j00, ...
    method, V1, EndOfLoop0, decomp_LU, parallel_comp)
%[phase_ini, amp, err ] = algo3ph_zero_it_par(k, contr, phase_ini, amp, err, x_phase0, use_invA, end_j00, ...
%    method, V1, EndOfLoop0, decomp_LU, parallel_comp)
%   it does parallel_comput case or not
%
% Maxime PINSARD

if (use_invA ||  decomp_LU)
% %     case 1 % 1 for inv. of matrix then multiply
        
        A00 = [length(x_phase0), sum(cos(x_phase0/180*pi)), sum(sin(x_phase0/180*pi)); ...
            sum(cos(x_phase0/180*pi)), sum(cos(x_phase0/180*pi).^2), sum(cos(x_phase0/180*pi).*sin(x_phase0/180*pi)); ...
            sum(sin(x_phase0/180*pi)), sum(cos(x_phase0/180*pi).*sin(x_phase0/180*pi)), sum(sin(x_phase0/180*pi).^2)];
        if decomp_LU
            [L,U] = lu(A00);
            inv_A00 ={L, U};
        elseif use_invA
            inv_A00 = A00^(-1);
            
        end
        
else %  0 for inv. of system
    inv_A00 = 0;
end

switch parallel_comp
    case 1
        hh=waitbar(0,sprintf('it. # %d : fit-I-SHG calculation (parallel computing implies no view or progress)...', k));
        parforArg = Inf; % max. of workers available
%         dividerWaitbar0 = 0;

    case 0
        hh=waitbar(0,sprintf('it. # %d : fit-I-SHG calculation ...', k));
        parforArg = 0;
end

dividerWaitbar0=10^(floor(log10(EndOfLoop0))-1);

if method == 1 % 2D matrices
    end_i = size(contr,1);
    
    parfor(i=1:end_i, parforArg) % dimension Y
        % WARNING : replace by a for loop to debug !
        
        for j = 1:end_j00 % dimension X
            
            [phase_ini(i, j), amp(i, j), err(i, j)] = fit_I_SHG_2(x_phase0, squeeze(contr(i,j, :)), inv_A00, use_invA, decomp_LU, [], 0, 0);
            
        end
        if (~parallel_comp && (round((end_i-i)/dividerWaitbar0)==(end_i-i)/dividerWaitbar0))
            
            waitbar((end_i-i)/EndOfLoop0,hh, sprintf('fit-I-SHG calculation : %.3g %%', (end_i-i)/EndOfLoop0*100));
        end
    end
    
else % 1D matrices
    
    parfor(j = 1:end_j00, parforArg) % dimension XxY
         % WARNING : replace by a for loop to debug !
        
        [phase_ini( j), amp( j), err(j)] = fit_I_SHG_2(x_phase0, V1(j, :)', inv_A00, use_invA, decomp_LU, [], 0, 0);
        
        % phase is between -1 , 1
        if (~parallel_comp && (round((end_j00-j)/dividerWaitbar0)==(end_j00-j)/dividerWaitbar0))
            waitbar((end_j00-j)/EndOfLoop0,hh, sprintf('it. # %d : fit-I-SHG calculation : %.3g %%', k, (end_j00-j)/EndOfLoop0*100));
        end
        
    end
    
end

 delete(hh);
end
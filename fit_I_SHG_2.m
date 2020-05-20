function [phi_opt, amp_opt, err_rel] = fit_I_SHG_2(x_phase,contr, inv_A, use_invA, decomp_LU, B, algo_vib, alpha_temp)
% [phi_opt, amp_opt, err_rel] = fit_I_SHG_2(x_phase,contr, inv_A, meth)
%
% edited 2015 - 10, Maxime Pinsard
%
% Prend 'contr' qui est la pile des 12 images, et la position
% correspondante. On a besoin des pixels adjacents pour moyenner
%
% Va creer 2 vecteurs de cos et sin avec la variation de phase imposee X(1)
% et X(2)
% solve the equation X*beta = pixel, et trouve la phase avec
% atan2(beta(2)/beta(1));
%
% On calcule la moyenne de la région 3x3 autour du pixel considéré en
% fonction de la phase.
% Poids attribué : 5 centre, 2 croix, 1 diago
% % if isempty(B) % % not input
    tmp = squeeze(contr);
    % On mesure l'amplitude et le décalage (en y) du cos obtenu
    amp = (max(tmp)-min(tmp))/2;
    mid = mean(tmp);
    % mid = (max(tmp)+min(tmp))/2; % C.A.'s code : not the same
    % does not seem to change a lot the results
% end

% MARCHE BIEN ET RAPIDE - voir test.m
% Si la cellule est exclue (voir le code dans I_SHG.m)
if (isempty(B) && tmp(1) == 0)
    phi_opt = 1.05;
    amp_opt = 0;
    err_rel = 0;
    % Si la cellule est inclue
else
    
    t = x_phase'; % Vecteur phase
    
    % meth = 1; 0 for inv. of system, 1 for inv. of matrix then multiply
    if (isempty(B) && (use_invA || decomp_LU)) % % not input
        B = [sum(tmp); sum(tmp.*cos(t/180*pi)); -sum(tmp.*sin(t/180*pi))];
    end
    amp1 = 1;

    if (use_invA && ~decomp_LU)  % 1 for inv. of matrix, then multiply (FASTER
        % took 28 sec, example # 1  PPLN
        % # example 2 (cornea, 21.10.2015) 4.0 sec
        % (mu0.42pi, sigma0.47pi) ; (mu-0.58pi, sigma0.49pi)
        beta = inv_A*B;
    elseif decomp_LU
        % %         [L,U] = inv_A;
        % this is an easy, triangular solve
        Y = inv_A{1}\B; % L
        % this is another triangular solve
        beta = inv_A{2}\Y; % U
    else % 0 for inv. of system
        
        % took 34 sec, example # 1
        % # example 2 (cornea, 21.10.2015) 5.3 sec
        % (mu0.42pi, sigma0.47pi) ; (mu-0.58pi, sigma0.49pi)
        X = ones(length(t),3);
        X(:,2) = mid+amp*cos(t/180*pi)'; % creation de la base (Identité cos sin)
        X(:,3) = mid+amp*sin(t/180*pi)'; % creation de la base (Identité cos sin)
        beta = X\tmp; % solve the equation X*beta = tmp
        
        amp1 = amp;
    end
    
    switch algo_vib
        case 0 % only tilted algo
            mid = mean(tmp);
        case 1 % tilted with vib. algo    
            amp1  = mean(alpha_temp); % amp0*
            mid = 0;
    end
    amp_opt = amp1*sqrt(beta(2)^2+beta(3)^2); %/beta(1);            
    phi_opt = atan2(-beta(3),beta(2))/pi; % atan but in the interval -pi pi
    % putting a minus sign increases the error !!
    
    mid_model = mid - mean(amp_opt*cos(t/180*pi-phi_opt*pi));%trapz(test)/length(test); % mean(test) -
    model = mid_model+amp_opt*cos(t/180*pi-phi_opt*pi);
    
    ss_tot = sum((tmp-mean(tmp)).^2);
    % ss_reg = sum((model-mean(tmp)).^2);
    ss_res = sum((tmp-model).^2);
    err_rel = ss_res/ss_tot;
    % h = figure;
    % plot(t,tmp,'g')
    % hold on;
    % plot(t,model,'b')
    % close(h)
    
end

end
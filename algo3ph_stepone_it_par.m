function [ phase_opt, amp, err, Delta, alpha, beta ] = algo3ph_stepone_it_par( k, contr, end_j, algo_tilted, algo_vib, raw_shift, ...
    Delta, phase_opt, amp, err, deltax, deltay, kx, ky, dx, dy, nb_frame, alpha_row, beta_row, alpha_col, beta_col, ...
    alpha, beta, delta, V1, use_invA, decomp_LU, parallel_comput, EndOfLoop, dividerWaitbar, hh)
%algo3ph_stepone_it_par Summary of this function goes here
%   Detailed explanation goes here
%
% Maxime PINSARD

%% init of var, matrices

phase_opt_no_vib_no_tilt = phase_opt; amp_no_vib_no_tilt = amp; err_no_vib_no_tilt = err;

alpha_temp_no_vib = 1;
beta_temp_no_vib = 1;
coeff_A_no_vib = nb_frame;

alpha_temp_notilt = alpha_temp_no_vib;
beta_temp_notilt = beta_temp_no_vib;
coeff_A_no_notilt = coeff_A_no_vib;
inv_A_not_algotilted = zeros(9); % for nothing
% % alpha_temp = 0;

delta_current = delta(:,1) ; % for nothing

if ~algo_tilted % && ~algo_vib) % no vib., no tilt.
    
    if (algo_vib || use_invA)
        switch algo_vib
            case 1
                alpha_temp_notilt = squeeze(alpha_col(:, 1, k-1)); % it's k-1 here because it's calculated at previous iteration !
                beta_temp_notilt = squeeze(beta_col(:, 1, k-1));
                coeff_A_no_notilt = sum(alpha_temp_notilt.^2); %  in one example, = 2.3435e+09 whereas nb_frame = 24
            case 0
                alpha_temp_notilt = alpha_temp_no_vib;
                beta_temp_notilt = beta_temp_no_vib;
                coeff_A_no_notilt = coeff_A_no_vib;
        end
        delta_current = delta(:, k);
        % %         DD = delta(:, k); % it's k here because it's calculated at previous iteration, but with an initial value imposed !
        
        term_cos1 = sum(alpha_temp_notilt.*beta_temp_notilt.*cos(delta(:, k)/180*pi));
        term_sin1 = sum(alpha_temp_notilt.*beta_temp_notilt.*sin(delta(:, k)/180*pi));
        term_cos2 = sum(beta_temp_notilt.^2.*cos(delta(:, k)/180*pi).^2);
        term_sin2 = sum(beta_temp_notilt.^2.*sin(delta(:, k)/180*pi).^2);
        term_sincos = sum(beta_temp_notilt.^2.*cos(delta(:, k)/180*pi).*sin(delta(:, k)/180*pi));
        % %  t1 = tic;
        % %                 if cond(A) > fact_sing
        % %                     disp('Matrix A for phase is close to singular !')
        % %                 end
        % % toc(t1); % 0.000058/px : long !!
        
        A_not_algotilted = [coeff_A_no_notilt, term_cos1, term_sin1; ...
            term_cos1, term_cos2, term_sincos; ...
            term_sin1, term_sincos, term_sin2];
        
        if decomp_LU
            [L, U] = lu( A_not_algotilted); % [L, U]
            inv_A_not_algotilted = {L, U};
        elseif use_invA
            inv_A_not_algotilted = A_not_algotilted^(-1);
        else
            
        end
    end
end

size_big = size(contr,1)*size(contr,2);
beta_row_current = zeros(size_big, size(alpha_col, 1));
alpha_row_current = zeros(size_big, size(alpha_col, 1));
beta_col_current = zeros(size(alpha_col, 1), size_big);
alpha_col_current = zeros(size(alpha_col, 1), size_big);

switch algo_tilted
    case 1
        end_i = size(contr,1);
        if algo_vib
            beta_row_current = beta_row(:, :, k-1);
            alpha_row_current = alpha_row(:, :, k-1);
            beta_col_current = beta_col(:, :, k-1);
            alpha_col_current = alpha_col(:, :, k-1);
        end
        kx_current = kx(:, k-1);
        ky_current = ky(:, k-1);
        dx_current = dx(:, k-1);
        dy_current = dy(:, k-1);
        deltax_current = deltax(:, :, k-1);
        deltay_current = deltay(:, :, k-1);
        
    case 0
        alpha = zeros(1, size_big, size(alpha,3));
        beta = zeros(1, size_big, size(alpha,3));
        contr = reshape(contr,[1,size_big, size(contr,3)]);
        Delta = reshape(Delta,[1, size_big, size(Delta,3)]);
        kx_current = kx;  % for nothing
        ky_current = ky;  % for nothing
        dx_current = dx;  % for nothing
        dy_current = dy; % for nothing
        deltax_current = zeros(size(alpha,3), size_big);  % for nothing
        deltay_current = zeros(size_big, size(alpha,3)); % for nothing
        phase_opt = reshape(phase_opt,[1, size_big]);
        end_i = 1;
end

if ~parallel_comput
    parforArg = 0;
else % no parallel computing
    parforArg = Inf; % max. of workers available
end

%% loop(s) within the array
%par
parfor(j = 1:end_j, parforArg) % on x
    % for j = 1:end_j % on x
    for i = 1:end_i % on y
        
        % parfor : 7.6 sec for t0 (one phase mat)
        % without parfor 21sec 273x199
        
        %                     t1 = tic; % !!!
        % 0.00009 sec with assignements, same without
        
        switch algo_tilted
            case 0
                %                 ind_wait_bar = i;
                switch algo_vib
                    
                    case 0 % no vib., no tilt.
                        
                        [phase_opt1, amp1, err1] = fit_I_SHG_2(delta_current, V1(j, :), inv_A_not_algotilted, use_invA, decomp_LU, [], algo_vib, 0);
                        phase_opt_no_vib_no_tilt(j) = phase_opt1; amp_no_vib_no_tilt(j) = amp1; err_no_vib_no_tilt(j) = err1;
                        
                    case 1 % with vibrations
                        
                        tmp = transpose(V1(j, :));
                end
            case 1
                %                 ind_wait_bar = j;
                switch raw_shift % raw vect of phase shifts
                    case 1
                        
                        Delta(i, j, :) = (squeeze(deltax_current(:, j)) + transpose(squeeze(deltay_current(i, :)))); % column vect
                    case 0 % fit vect of phase shifts
                        Delta(i, j, :) = (kx_current*j + dx_current + ky_current*i + dy_current); % column vect
                end
                
                DD= squeeze(Delta(i, j, :));% only the last it. will be available, otherwise create a 4 dim matrix
                
                tmp = squeeze(contr(i, j, :)); % is a Mx1 vector
                
                switch algo_vib
                    
                    
                    case 1 % tilted with vib. algo
                        
                        beta_temp = (transpose(squeeze(beta_row_current(i, :))) + squeeze(beta_col_current(:, j)))/2; % column vect, avg of both
                        alpha_temp  = (transpose(squeeze(alpha_row_current(i, :))) + squeeze(alpha_col_current(:, j)))/2; % column vect, avg of both
                        
                        beta(i, j, :) = beta_temp;
                        alpha(i, j, :) = alpha_temp;
                        
                        coeff_A = sum(alpha_temp.^2);
                    case 0
                        alpha_temp  = alpha_temp_no_vib;
                        beta_temp = beta_temp_no_vib;
                        coeff_A = coeff_A_no_vib;
                        
                end
        end
        % %  t1 = tic;
        % %                 if cond(A) > fact_sing
        % %                     disp('Matrix A for phase is close to singular !')
        % %                 end
        % % toc(t1); % 0.000058/px : long !!
        if (algo_vib || algo_tilted) % not normal case
            
            switch algo_tilted
                case 1
                    term_cos1 = sum(alpha_temp.*beta_temp.*cos(DD/180*pi));
                    term_sin1 = sum(alpha_temp.*beta_temp.*sin(DD/180*pi));
                    term_cos2 = sum(beta_temp.^2.*cos(DD/180*pi).^2);
                    term_sin2 = sum(beta_temp.^2.*sin(DD/180*pi).^2);
                    term_sincos = sum(beta_temp.^2.*cos(DD/180*pi).*sin(DD/180*pi));
                    % %  t1 = tic;
                    
                    A = [coeff_A, term_cos1, term_sin1; ...
                        term_cos1, term_cos2, term_sincos; ...
                        term_sin1, term_sincos, term_sin2];
                case 0
                    alpha_temp = alpha_temp_notilt;
                    beta_temp = beta_temp_notilt;
                    DD = delta_current;
            end
            
            B = [sum(alpha_temp.*tmp); sum(beta_temp.*tmp.*cos(DD/180*pi)); sum(beta_temp.*tmp.*sin(DD/180*pi))];
            
            if (algo_tilted && decomp_LU)  % LU is safe and tested
                [L,U] = lu(A);
                inv_A = {L, U};
            elseif algo_tilted
                % normal inversion (not LU), or not algo_tilted
                % takes 0.000341 seconds per pixel for just A --> very long
                % % t1 = tic;
                inv_A = A^(-1);
            else
                inv_A = inv_A_not_algotilted;
            end
            
            % % toc(t1); % 0.000044 seconds /px : LU is faster
            
            [phase_opt1, amp1, err1] = fit_I_SHG_2(delta_current, tmp', inv_A, use_invA, decomp_LU, B, algo_vib, alpha_temp);
            phase_opt(i, j) = phase_opt1;
            amp(i,j) = amp1;
            err(i,j) = err1;
            
            % %         phase_opt(i, j) = atan2(-U_cap(3),U_cap(2))/pi; % atan but in the interval -pi pi
            % %
            % %         model = mid+amp(i,j)*cos(DD/180*pi + phase_opt(i, j)*pi);
            % %         ss_tot = sum((tmp-mean(tmp)).^2);
            % %         % ss_reg = sum((model-mean(tmp)).^2);
            % %         ss_res = sum((tmp-model).^2);
            % %         err(i,j) = ss_res/ss_tot;
        end
    end
    
    %         if ~parforArg
    %             if ((round(ind_wait_bar/dividerWaitbar)==ind_wait_bar/dividerWaitbar) && ~algo_tilted)
    %                 waitbar(ind_wait_bar/EndOfLoop,hh, sprintf('it. # %d : fit-I-SHG calculation : %.3g %%', k, ind_wait_bar/EndOfLoop*100));
    %             end
    %         end
    
    if ~parforArg
        ind_wait_bar = end_j-j;
        if ((round(ind_wait_bar/dividerWaitbar)==ind_wait_bar/dividerWaitbar))
            waitbar(ind_wait_bar/EndOfLoop,hh, sprintf('it. # %d : fit-I-SHG calculation : %.3g %%', k, ind_wait_bar/EndOfLoop*100));
        end
    end
    
end

if (~algo_vib && ~algo_tilted)
    phase_opt = phase_opt_no_vib_no_tilt;
    amp = amp_no_vib_no_tilt;
    err = err_no_vib_no_tilt;
end

end


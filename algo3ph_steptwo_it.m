function [deltax, kx, dx, delta, beta_col, alpha_col, diff_k, h, no_conv, eps_fact ] = algo3ph_steptwo_it( x_phase0, coeff_sum2, algo_tilted, k, raw_shift, phase00, deltay, ky, dy, contr, V1, decomp_LU ,...
    deltax, kx, dx, delta, algo_vib, beta_col, alpha_col, corr, epsilon, diff_k, coeff_Ap, verbose, crit2, corr_str, first_delta_to_zero)
% algo3ph_steptwo_it
%  step 2 of the algo 3 phases code
% calculate kx, dx, deltax if tilted, delta if not,
% beta_col, alpha_col if tilted + vib
% need ky, dy, deltay
%
% 2016.7.7, Maxime PINSARD


h = 0; % for output
eps_fact = 0; % for nothing

no_conv = zeros(1, length(x_phase0));

step_theo = x_phase0(2)- x_phase0(1); % in degree
% % step_theo = 0;

if (algo_tilted && k > 1)
    phase_opt1=zeros(size(contr,1), size(contr,2));
else % if it's the first iteration, deltay is not calculated so use the phase00
    % if not tilted algo, deltay is never calculated
    phase_opt1 = phase00;
end

delta1 = zeros(1, coeff_sum2);

for m = 1:length(x_phase0) % on all phase shifts
    
    %          t1 = tic;
    for jj = 1:coeff_sum2 % X, on columns of phase ( no loop if no tilt)
        
        switch algo_tilted
            case 1
                
                if k > 1
                    
                    switch raw_shift
                        case 1 % raw values of deltay considered
                            phase_opt1(:, jj) = phase00( :, jj) + squeeze(deltay( :, m, k-1))/180; % k-1 and not k (deltay was determined at repvious it.)
                        case 0 % fit values of deltay considered
                            phase_opt1(:, jj) = phase00( :, jj) + (ky(m, k-1).*(1:size(contr, 1))' + dy(m, k-1))/180 ;
                    end
                    
                    %                 else % k == 1, ky and dy do not exist yet
                    %                     phase_opt1(:, jj) = phase00( :, jj) + delta(m,1)/180 ;
                end
                % delta_y must be a row, dependent of m
                tmp2 = contr( :, jj, m);
                
            case 0
                tmp2 = V1(:, m);
        end
        
        term_cos1p = sum(cos(phase_opt1(:, jj)*pi)); % if no algo_tilt it's normal that A is calculated here and with j index because phase_opt1 is re-organized in one column
        term_sin1p = sum(sin(phase_opt1(:, jj)*pi));
        term_cos2p = sum(cos(phase_opt1(:, jj)*pi).^2);
        term_sin2p = sum(sin(phase_opt1(:, jj)*pi).^2);
        term_sincosp = sum(cos(phase_opt1(:, jj)*pi).*sin(phase_opt1(:, jj)*pi));
        
        Ap = [coeff_Ap, term_cos1p, term_sin1p,; ...
            term_cos1p, term_cos2p, term_sincosp; ...
            term_sin1p, term_sincosp, term_sin2p];
        % %             if cond(Ap) > fact_sing
        % %                 disp('Matrix Ap for deltax is close to singular !')
        % %             end
        
        Bp = [sum(tmp2); sum(tmp2.*cos(phase_opt1(:, jj)*pi)); sum(tmp2.*sin(phase_opt1(:, jj)*pi))];
        
        switch decomp_LU % LU is safe and tested
            case 1
                % % t1 = tic;
                [Lp,Up] = lu(Ap);
                Yp = Lp \ Bp; % this is an easy, triangular solve
                betap = Up \ Yp; % this is another triangular solve
                % % toc(t1); % 0.000015 seconds/px : LU is faster
                
            case 0 % normal inversion
                % takes 0.000341 seconds per pixel for just A --> very long
                % % t1 = tic;
                inv_Ap = Ap^(-1);
                %     val = 0; % ini
                betap = inv_Ap*Bp;
                % % toc(t1); % 0.000044 seconds /px : LU is faster
        end
        
        deltax(m, jj, k) = atan2(-betap(3),betap(2))/pi*180; % vector of phase shifts, in degree
        
        % %             if deltax(m, jj, k) < 0
        % %                 deltax(m, jj, k) = deltax(m, jj, k) + 360; % in interval [0, 2pi]
        % %             end
        
        switch algo_vib
            case 1 % with vib
                %                 if k > 1
                beta_col(m, jj, k) = (betap(2).^2 + betap(3).^2).^(1/2); % modulation fluctuations in columns
                alpha_col(m, jj, k) = betap(1); % background fluctuations in columns
                %                 end
                % beta is 1000 on example, alpha 1e3 --> normal ?
        end
        
        % "After calculating the x-directional phase shift ?xj for all
        % the columns in the jth interferogram, unwrapping it ... " in
        % article --> unwrapping is performed (with tilted) on
        % xy directions only
        
        switch corr % correction or not : IS DONE ON a direction (x or y)
            case 4 % corr
                
                if m > 1
                    
                    if (~algo_tilted && first_delta_to_zero)
                        %                     case 0 % because it's not done after the loop in this case
                        deltax(m, 1, k) = deltax(m, 1, k) - delta1;
                    end
                    
                    diff1 = deltax(m, jj, k) - deltax(m-1, jj, k);
                    
                    ct = 1; pp = 0;
                    while (abs(diff1 - sign(diff1)*step_theo + pp*180) > 90)
                        
                        ct=ct+1;
                        pp = abs(pp)*(-1)^ct;
                        if pp >= 0
                            pp=pp+1;
                        end
                    end
                    
                    if (m > 2 && sign(deltax(m-1, 1, k) - deltax(m-2, 1, k)) ~= sign(diff1))
                        [~, posmin] = min([abs(diff1  + pp*180), abs(diff1  - pp*180), abs(diff1  + (pp+1)*180), abs(diff1  + (pp-1)*180)]);
                        %                     pp = pp*(-1)^mod(posmin+1,2); % =1 if 1, -1 if 2
                        switch posmin
                            case 2
                                pp=-pp;
                            case 3
                                pp=pp+1;
                            case 4
                                pp=pp-1;
                        end
                    end
                    deltax(m, jj, k) = deltax(m, jj, k) + pp*180;
                    
                    %                     deltax_pp = deltax(m, jj, k);
                    %                     for pp=[-2, -1, 1, 2]
                    %                         if (abs(diff1 + pp*180) < abs(deltax_pp - deltax(m-1, jj, k)))
                    %                             deltax_pp = deltax(m, jj, k) + pp*180;
                    %                             % %                             break; % only one value is better ?
                    %                         end
                    %                     end
                    %                     deltax(m, jj, k) = deltax_pp; % - delta1(jj);
                else % m = 1
                    delta1(jj) = deltax(m, jj, k);
                    if (~algo_tilted && first_delta_to_zero)
                        deltax(1, 1, k) = 0;
                    end
                end
        end
        
    end % end of jj loop
    
    %         toc(t1); % 0.03 sec
    
    % % *********** UNWRAPPING + REG LIN ****************
    
    switch algo_tilted
        
        case 0 % normal algo
            
            delta(m, k+1) = deltax(m, 1, k);
            
            % little corr
            diff_k(m, k+1) = min( abs( -(delta(m, k+1) - delta(1, k+1)) - (delta(m, k) - delta(1, k))),...
            abs( (delta(m, k+1) - delta(1, k+1)) - (delta(m, k) - delta(1, k))));
            
%             diff_k(m, k+1) = sign(diff_k(m, k+1))*min([abs(diff_k(m, k+1)), abs(diff_k(m, k+1)-360), abs(diff_k(m, k+1)+360)]);
            
            fprintf('Remaining difference = %.5f (degree)\n', diff_k(m, end))
            eps_fact = epsilon*180/pi;
            
            criterion1 = diff_k(m, end) > eps_fact; % because epsilon is in rad
            
            if criterion1 % probing convergence
                
                no_conv(m) = 1; % no convergence
                
            else
                no_conv(m) = 0; % convergence for one phase shift
            end
    end
    
end % trying to correct phase also on m, need to do 2 loops

% this loop has to be separated from the other, otherwise the deltax values
% are considered after the subtraction of delta1 and not before
switch algo_tilted
    case 1 % tilted algo
        
        %             if corr == 2 % is done on X direction ??
        %                 deltax(m, :, k) = 180/pi*unwrap(deltax(m, :, k)/180*pi);
        %
        %             end
        for m = 1:length(x_phase0) % on all phase shifts
            
            if (m> 1 && first_delta_to_zero)
                deltax(m, :, k) = deltax(m, :, k)-delta1;
            end
            
            y_fit_col = squeeze(deltax(m, :, k));
            y_fit_col = y_fit_col(abs(y_fit_col-median(y_fit_col)) < max(90, step_theo)); % no aberration values for this fit
            
            xv = 1:length(y_fit_col);
            [p_fit, gof1] = polyfit(xv, y_fit_col, 1); % a*x+b
            kx(m, k) = p_fit(1); % col. tilt factor on mth interferogram
            dx(m, k) = p_fit(2); % col. translational value on mth interferogram
            
            
            err1 = gof1.normr; % error, must be as small as possible
            
            switch verbose
                case 1
                    fprintf(' err on fit columns = %.3g\n', err1);
            end
            
            if m == crit2
                if k > 1
                    try
                        close(h);
                    catch
                        disp('Reopening fig.');
                    end
                end
                h = figure(100); subplot(2,1,1); hold on;
                plot(1:size(contr, 2), kx(m, k)*(1:size(contr, 2)) + dx(m, k), 'r');
                plot(xv, y_fit_col, 'bx');
                xlabel('column #'); ylabel('delta_x (deg)'); title(sprintf('Corrected by %s', corr_str));
            end
        end
        
        if (k>2)
            diff1 = dx(m, k) - dx(m, k-1);
            
            ct = 1; pp = 0;
            while (abs(diff1 + pp*180) > 90)
                
                ct=ct+1;
                pp = abs(pp)*(-1)^ct;
                if pp >= 0
                    pp=pp+1;
                end
            end
            dx(m, k) = dx(m, k)+ pp*180;
        end
        
        if first_delta_to_zero
            dx(1, k) = 0;%dx(1, k) - mean(delta1);
            deltax(1, :, k) = deltax(1, :, k) - delta1;
        end
end


end


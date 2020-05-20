function [ deltay, ky, dy, beta_row, alpha_row, diff_k1, diff_k2, no_conv, eps_fact1, eps_fact2] = algo3ph_stepthree_it( x_phase0, k, raw_shift, phase00, deltay, ky, dy, contr, decomp_LU ,...
    deltax, kx, dx, algo_vib, beta_row, alpha_row, corr, verbose, crit2, corr_str, diff_k1, diff_k2, epsilon1, epsilon2, h, k_max, delta, first_delta_to_zero)
% algo3ph_stepthree_it
% calculate ky, dy, deltay, beta_row, alpha_row
% need kx, dx, deltax
%
% 2016.7.8 Maxime PINSARD


no_conv = zeros(1, length(x_phase0));

step_theo = x_phase0(2)- x_phase0(1); % in degree
% % step_theo = 0;

if k == 1
    phase_opt2 = phase00;
end

delta1 = zeros(1, size(phase00, 1));

for m = 1:length(x_phase0) % on all phase shifts
    if k > 1
        phase_opt2=zeros(size(contr,1), size(contr,2));
    end
    
    for ii=1:size(phase_opt2, 1) % Y, on rows
        
        tmp3 = contr(ii, :, m);
        
%         if k > 1 % on first iteration, the value of deltax is ignored to have deltay because deltay is ignored to have deltax
            switch raw_shift
                case 1 % raw values of deltay considered
                    phase_opt2(ii, :) = phase00(ii, :) + squeeze(deltax(m, :, k))/180; % must be k and not k-1 because the deltax was calculated at this step
                case 0 % fit values of deltay considered
                    phase_opt2(ii, :) = phase00(ii, :) + (kx(m, k).*(1:size(contr, 2)) + dx(m, k))/180 ; %
            end
            
%         else % k == 1, ky and dy do not exist yet
%              phase_opt2(ii, :) = phase00(ii, :) + delta(m,1)/180;
%          
%         end      
        
        term_cos1pp = sum(cos(phase_opt2(ii, :)*pi));
        term_sin1pp = sum(sin(phase_opt2(ii, :)*pi));
        term_cos2pp = sum(cos(phase_opt2(ii, :)*pi).^2);
        term_sin2pp = sum(sin(phase_opt2(ii, :)*pi).^2);
        term_sincospp = sum(cos(phase_opt2(ii, :)*pi).*sin(phase_opt2(ii, :)*pi));
        
        Mpp = [size(phase_opt2, 2),term_cos1pp, term_sin1pp,; ...
            term_cos1pp, term_cos2pp, term_sincospp; ...
            term_sin1pp, term_sincospp, term_sin2pp];
        
        % %                     if cond(Mpp) > fact_sing
        % %                         disp('Matrix Mpp for deltay is close to singular !')
        % %                     end
        
        Bpp = [sum(tmp3); sum(tmp3.*cos(phase_opt2(ii, :)*pi)); sum(tmp3.*sin(phase_opt2(ii, :)*pi))];
        
        switch decomp_LU % LU is safe and tested
            case 1
                % % t1 = tic;
                [Lpp,Upp] = lu(Mpp);
                Ypp = Lpp \ Bpp; % this is an easy, triangular solve
                betapp = Upp \ Ypp; % this is another triangular solve
                % % toc(t1); % 0.000015 seconds/px : LU is faster
                
            case 0 % normal inversion
                % takes 0.000341 seconds per pixel for just A --> very long
                % % t1 = tic;
                inv_Mpp = Mpp^(-1);
                betapp = inv_Mpp*Bpp;
                % % toc(t1); % 0.000044 seconds /px : LU is faster
        end
        
        deltay(ii, m, k) = atan2(-betapp(3),betapp(2))/pi*180; % vector of phase shifts, in degree
        
        % %                     if deltay(ii, m, k) < 0
        % %                         deltay(ii, m, k) = deltay(ii, m, k) + 360; % in interval [0, 2pi]
        % %                     end
        
        switch algo_vib
            case 1 % with vib
%                 if k > 1
                    beta_row(ii, m, k) = (betapp(2).^2 + betapp(3).^2).^(1/2); % modulation fluctuations in columns
                    alpha_row(ii, m,  k) = betapp(1); % background fluctuations in columns
%                 end
        end
        
        switch corr % correction or not : IS DONE ON a direction (x or y)
            case 4 % corr

                if m > 1
                    
                    diff1 = deltay(ii, m, k) - deltay(ii, m-1, k);
                    
                    ct = 1; pp = 0;    
                    while (abs(diff1 - sign(diff1)*step_theo + pp*180) > 90)
                        
                        ct=ct+1;
                        pp = abs(pp)*(-1)^ct;
                        if pp >= 0
                            pp=pp+1;
                        end
                    end
                    
                    if (m > 2 && sign(deltay(ii, m-1, k) - deltay(ii, m-2, k)) ~= sign(diff1))
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

                    deltay(ii, m, k) = deltay(ii, m, k) + pp*180;
%                     deltay_pp = deltay(ii, m, k);
%                     for pp=[-2, -1, 1, 2]
%                         if (abs(diff1 + pp*180) < abs(deltay_pp - deltay(ii, m-1, k)))
%                             deltay_pp = deltay(ii, m, k) + pp*180;
% % %                             break; % only one value is better ?
%                         end
%                     end
%                     deltay(ii, m, k) = deltay_pp;% - delta1;
                else % m = 1
                    delta1(ii) = deltay(ii, m, k);
                end
        end
        
        % "After calculating the x-directional phase shift ?xj for all
        % the columns in the jth interferogram, unwrapping it ... " in
        % article --> unwrapping is performed (with tilted) on
        % xy directions only
                
    end % end of ii loop
    
    %                 toc(t2); % 0.03 sec
    
end % end of mm loop

for m = 1:length(x_phase0) % on all phase shifts    
%     if corr == 2 % is done on direction
%         deltay(:, m, k) = 180/pi*unwrap(deltay(:, m, k)/180*pi);
%     end

    if (m> 1 && first_delta_to_zero)
        deltay(:, m, k) = deltay(:, m, k) - delta1';
    end
    
    y_fit_row = squeeze(deltay(:, m, k));
    y_fit_row = y_fit_row(abs(y_fit_row-median(y_fit_row)) < max(90, step_theo));
    
    yv = 1:length(y_fit_row);
    [p_fit2, gof2] = polyfit(yv', y_fit_row, 1); % a*x+b
    ky(m, k) = p_fit2(1); % row tilt factor on mth interferogram
    dy(m, k) = p_fit2(2); % row translational value on mth interferogram
    
    err2 = gof2.normr; % error, must be as small as possible
    switch verbose
        case 1
            fprintf(' err on fit rows = %.3g\n', err2);
    end
   
    if m == crit2
        try
            figure(h); 
        catch
            h = figure(100);
        end
        subplot(2,1,2); hold on;
        plot(1:size(contr, 1), ky(m, k)*(1:size(contr, 1)) + dy(m, k) , 'r');
        plot(yv', y_fit_row, 'bx');
        xlabel('row #'); ylabel('delta_y (deg)'); title(sprintf('Corrected by %s', corr_str));
    end
    
    if (k>2)
        diff1 = dy(m, k) - dy(m, k-1);
        
        ct = 1; pp = 0;
        while (abs(diff1 + pp*180) > 90)
            
            ct=ct+1;
            pp = abs(pp)*(-1)^ct;
            if pp >= 0
                pp=pp+1;
            end
        end
        dy(m, k) = dy(m, k)+ pp*180;
    end
    
    if first_delta_to_zero
        dy(1, k) = 0;%dy(1, k) - mean(delta1);
        deltay(:, 1, k) = deltay(:, 1, k) - delta1';
    end
        [diff_k1, diff_k2, no_conv, eps_fact1, eps_fact2 ] = verify_convergence_STP3_algo( k, m, dx, dy, kx, ky, diff_k1, diff_k2, epsilon1, epsilon2, no_conv, k_max, delta );
end

end


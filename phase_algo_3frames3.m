function [ phase_opt, amp, err, delta, diff_k, divers_vect ] = phase_algo_3frames3( contr, x_phase0, use_invA, ...
    k_max, epsilon, epsilon1, algo_tilted, algo_vib, corr, h_f, method, parallel_comput, decomp_LU, ...
    reit_prev_phase, prev_res, ramp_phshft )
% phase_algo_3frames : Calculate phase with the optimized algorithm in Wang
% et al., OL, 2004
%
% 2016.6.3 by Maxime PINSARD
%
% [ phase, amp, err, delta, diff_k ] = phase_algo_3frames( contr, x_phase0, use_invA, k_max, epsilon, algo_tilted, algo_vib, corr )
%
% pb alpha, beta : there are equal to ~ 10000 whereas they should be only
% multiplication factor of intensity (around 1.2)
% initial value of alpha, beta !!
% vib. algorithm is a combination of previous ones
% tried to calculate alpha, beta only at 2nd step, with different PPLN, to
% correct it (they are too big, like 1e.3 or 1e4)
% but did not improved
%
% 2017.6.28 : improvements on corr of deltax and deltay : use of +m*180 to
% corr, m in [-2,-1,1,2]
% at 1st iteration, use of delta_ini instead of nothing (and instead of
% deltay(k) that do not exist yet)
% try to put the 1st delta to 0, so the comparison in row/columns has the
% same reference (dangerous !)
%
% look at alpha, beta because they seem high ...

% x_phase0 = x_phase0/180*pi; % all is converted in rad

% supposing there will be at least 3 iterations

nb_frame = length(x_phase0);

if algo_vib
    fprintf(2, 'In vib. algo, alpha, beta are calculated only from 2nd it. !!\n');
end

verbose = 0; % 1 for ON

% epsilon2 = 2e-3; % rad, 1e-6 in article
% epsilon1 = 0.01; % rad, 1e-2; % 1e-5 in article
% corr = 1; % warning : change if not tilted algo

%fact_sing = 1e8; % threshold value to considered matrix singular (Inf if 100% singular)

if corr == 4
    corr_str = 'Outside only on phase-shifts corr.';
else
    fprintf(2, 'CAUTION : delta_x and y are NOT unwrapped !\n');
    corr_str = 'Nothing';
end

% Pinsard : with 1D matrices
% example 2 cornea : took 3.4 sec

raw_shift = 0; % 1 for raw

first_delta_to_zero = 1; % putting this to zero makes the algo to converge, but ther is a progressive shift of the delta which results in a very tilted phase map

% KEEP IT THERE for parfor
deltax = zeros( size(contr, 3), size(contr, 2)); % phase-shift of columns
deltay = zeros(  size(contr, 1), size(contr, 3)); % phase-shift of rows
dx = zeros(length(x_phase0)); % col. translational values
dy = zeros(length(x_phase0)); % row translational values
kx  = zeros(length(x_phase0)); % col. tilt factors
ky   = zeros(length(x_phase0)); % row tilt factors
beta_col = zeros( size(contr, 3), size(contr, 2)); % modulation fluctuations in columns
alpha_col = zeros( size(contr, 3), size(contr, 2)); % background fluctuations in columns

beta_row = zeros(size(contr, 1), size(contr, 3)); % modulation fluctuations in rows
alpha_row = zeros(size(contr, 1), size(contr, 3)); % background fluctuations in rows

beta_col(:,:,1) = ones(size(contr, 3), size(contr, 2));
alpha_col(:,:,1) = ones(size(contr, 3), size(contr, 2));
beta_row(:,:,1) = ones(size(contr, 1), size(contr, 3));
alpha_row(:,:,1) = ones(size(contr, 1), size(contr, 3));

alpha = zeros(size(contr,1), size(contr,2), size(contr,3)); % background fluctuations
beta = zeros(size(contr,1), size(contr,2), size(contr,3)); % modulation fluctuations
Delta  = zeros(size(contr,1), size(contr,2), size(contr,3));


if method == 1 % with 2D matrices
    
    phase_ini = zeros(size(contr,1), size(contr,2));
    EndOfLoop0=size(contr,1);
    end_j00 = size(contr,2);
    V1 = 0;
    amp=zeros(size(contr,1), size(contr,2));
    err = zeros(size(contr,1), size(contr,2));
    
else % with 1D matrices
    
    phase_ini = zeros(1, size(contr,1)*size(contr,2));
    V1 = reshape(contr, size(contr,1)*size(contr,2), size(contr,3));
    EndOfLoop0=size(V1,1);
    end_j00 = size(contr,1)*size(contr,2);
    amp=zeros(1, size(contr,2)*size(contr,1));
    err = zeros(1, size(contr,2)*size(contr,1));
    
    delta = x_phase0';
    
    switch raw_shift
        case 1 % 1 for raw
            fprintf(2, 'CAUTION : delta_x and y are input with their raw values !\n');
    end
    
    switch algo_tilted
        case 1
            phase_opt=zeros(size(contr,1), size(contr,2));
            amp=zeros(size(contr,1), size(contr,2));
            err = zeros(size(contr,1), size(contr,2));
            
            coeff_sum2 = size(contr, 2);
            
            coeff_Ap = size(contr, 1); % Y
            end_j = size(contr,2);
            EndOfLoop=size(contr,2);
            
            %         algo_vib = 1; % 1 for algo with vibrations
            diff_k1 = [];%zeros(length(x_phase0), 2); % error vector for tilted algo : on translational values
            diff_k2 = [];%zeros(length(x_phase0), 2); % error vector for tilted algo : on tilted factors
            diff_k = 0; % for nothing
            %         diff_k1(:, 1) = ones(length(x_phase0), 1);
            %         diff_k2(:, 1) = ones(length(x_phase0), 1);
            %         switch algo_vib
            %             case 1
            %                 %                 alpha_temp = ones(1, nb_frame);
            %                 %                 beta_temp = ones(1, nb_frame);
            %         end
            %
            
        case 0
            diff_k = zeros(length(x_phase0), 1); % error vector for normal algo
            amp=zeros(1, size(contr,2)*size(contr,1));
            err = zeros(1, size(contr,2)*size(contr,1));
            coeff_sum2 = 1;
            %         corr = 1;
            coeff_Ap = size(contr,2)*size(contr,1);
            end_j = size(contr,1)*size(contr,2);
            EndOfLoop=size(contr,1)*size(contr,2);
            %         algo_vib = 0; % stays to zero
    end
end


k = 1;
no_convergence = 1;

while no_convergence % there is a break later for the case of 1 iteration
    
    t11 = tic;
    fprintf('\n ---------\n Iteration # %d\n', k)
    
    if getappdata(0,'param_algo_changed')
        k_max = getappdata(0,'k_max');
        epsilon1 = getappdata(0 , 'epsilon_tilt_th');
        epsilon = getappdata(0 , 'epsilon_ph_th');
        setappdata(0 , 'param_algo_changed', 0);
    end
    if k > k_max % too many iterations
        
        switch k_max > 1
            case 1 % not 1st it.
                switch algo_tilted
                    case 0
                        diff_k_temp = diff_k;
                        %                         delta(:, k+1) = delta(:, iii+1);
                    case 1
                        diff_k_temp = diff_k1;
                        %                         [~, iii] = min(abs(mean(diff_k1(:,2:end), 1)));
                end
                [~, iii] = min(abs(mean(diff_k_temp(:,2:end), 1)));
                fprintf('\n Best iteration is # %d\n', iii) % real # of iteration is iii
                iii = iii+1; % because 1st was ignored as it is equal to zero
                k=iii-1;  % must be iii-1 because to calculate phase # iii, need deltax and y # iii-1
                
            case 0 % 1st it.
                delta(:,2) = x_phase0'; % !!
        end
        no_convergence = 0;
        fprintf('it. max reached (%d) \n', k_max);
        
        if iii == length(diff_k_temp)
            break
        end
        %                         k =1;
        %                         phase_opt = phase_ini; amp = amp_ini; err = err_ini;
        
        %     else % normal running
        %         ind = k-1; % because k has become k+1 since the last results
    end
    
    crit2 = randsample(size(contr, 3),1); % for the plot
    
    %% STEP ONE : classic step --> assuming background intensity + modulation will be constant frame to frame
    % CALCULATES phase from phase-shifts
    
    if ~algo_tilted
        phase_opt=zeros(1, size(contr,1)*size(contr,2));
    else
        phase_opt=zeros(size(contr,1), size(contr,2));
    end
    % t0 = tic;
    %     try
    
    % SPECIAL CASE OF STEP ONE IF IT'S THE FIRST ITERATION - IT'S THE ONLY STEP IF NO ITERATIVE ALGO
    if k == 1 % 1st it.
        
        if reit_prev_phase
            phase_ini = prev_res{1}; amp_ini = prev_res{2}; err_ini = prev_res{3};
            divers_vect = prev_res{4};
            if method ~= 1 % with 2D matrices
                phase_ini = reshape(phase_ini, 1, size(contr,1)*size(contr,2));
            end
            delta=[divers_vect{1}, divers_vect{1}]; deltax=divers_vect{2}; kx=divers_vect{3}; dx=divers_vect{4}; beta_col=divers_vect{5}; alpha_col=divers_vect{6};
            deltay=divers_vect{7}; ky=divers_vect{8}; dy=divers_vect{9}; beta_row=divers_vect{10}; alpha_row=divers_vect{11};
            
            re_calc=0;
            if algo_tilted
                if (~isempty(Delta(Delta~=0)) && ~isempty(deltax(deltax~=0)) && ~isempty(deltay(deltay~=0))...
                        && ~isempty(deltay(deltay~=0)) && ~isempty(kx(kx~=0)) && ~isempty(ky(ky~=0)) && ~isempty(dx(dx~=0)) && ~isempty(dy(dy~=0))...
                        && ~isempty(kx(kx~=0))) % non-zero
                    re_calc=0;
                else
                    re_calc=1;
                end
            end
            if algo_vib
                if (~isempty(alpha_row(alpha_row~=0)) && ~isempty(beta_row(beta_row~=0))...
                        && ~isempty(alpha_col(alpha_col~=0)) && ~isempty(beta_col(beta_col~=0)))
                    re_calc=0;
                else
                    re_calc=1;
                end
            end
            if ~re_calc
                %             k, contr, end_j, algo_tilted, algo_vib, raw_shift, ...
                %                     Delta, phase_opt, amp, err, deltax, deltay, kx, ky, dx, dy, nb_frame, alpha_row, beta_row, alpha_col, beta_col, ...
                %                     alpha, beta, delta
                k = k+1; continue; % not to check directly the convergence
            end
        else
            % it does parallel_comput case or not
            [phase_ini, amp_ini, err_ini ] = algo3ph_zero_it_par(k, contr, phase_ini, amp, err, x_phase0, use_invA, end_j00, ...
                method, V1, EndOfLoop0, decomp_LU, parallel_comput);
        end
        
        if algo_tilted
            phase_opt = reshape(phase_ini, size(contr,1),size(contr,2));
        else % can keep 1D matrices
            phase_opt = phase_ini;
        end
        %         phase_opt1 = phase_opt;
        
    else % k > 1
        
        switch parallel_comput
            
            case 1
                
                hh=waitbar(0,sprintf('it. # %d : fit-I-SHG calc. (parallel computing implies no view or progress) ...', k));
            case 0
                
                hh=waitbar(0,sprintf('it. # %d : fit-I-SHG calculation ...', k));
        end
        dividerWaitbar=10^(floor(log10(EndOfLoop))-1);
        
        [ phase_opt, amp, err, Delta, alpha, beta ] = algo3ph_stepone_it_par( k, contr, end_j, algo_tilted, algo_vib, raw_shift, ...
            Delta, phase_opt, amp, err, deltax, deltay, kx, ky, dx, dy, nb_frame, alpha_row, beta_row, alpha_col, beta_col, ...
            alpha, beta, delta, V1, use_invA, decomp_LU, parallel_comput, EndOfLoop, dividerWaitbar, hh);
        
    end
    
    %     catch
    %         r
    %     end
    % toc(t0); % 4.38 sec for 199*273*24, normal algo
    % 5.84 sec for k ==1 titlted algo
    % 13.4 sec for k > 1 titlted algo
    
    % phase is between -1 , 1
    
    if method < 2 % DO NOT PUT IT BEFORE
        amp = amp_ini;
        err = err_ini;
        break; % go outside 'while' loop on k if method no algorithm
    else
        try
            delete (hh);
        catch
            disp('Couldn`t delete bar');
        end
    end
    
    switch algo_tilted
        case 0
            phase_opt = phase_opt';
    end
    
    try
        axes(h_f);
    catch
        figure; h_f = axes;
    end
    imagesc(reshape(phase_opt, size(contr,1),size(contr,2))); colormap(hsv); colorbar; axis('image');
    
    %% STEP 2 : ON COLUMNS tilt phase shift in x-direction if algo_tilted : calculate deltax, kx, dx
    % otherwise just global phase-shift delta + convergence checking as in "Advanced iterative algorithm
    % for phase extraction of randomly phase-shifted interferograms"
    % IF VIB calculates beta_col, alpha_col
    
    [deltax, kx, dx, delta, beta_col, alpha_col, diff_k, h, no_conv, eps_fact ] = algo3ph_steptwo_it( x_phase0, coeff_sum2, algo_tilted, k, raw_shift, phase_opt, deltay, ky, dy, contr, V1, decomp_LU ,...
        deltax, kx, dx, delta, algo_vib, beta_col, alpha_col, corr, epsilon, diff_k, coeff_Ap, verbose, crit2, corr_str, first_delta_to_zero);
    
    %% STEP 3 : ON ROWS tilt phase shift in y-direction if algo_tilted  : calculate deltay, ky, dy
    % otherwise just convergence check as in "Advanced iterative algorithm
    % for phase extraction of randomly phase-shifted interferograms"
    % IF VIB calculates beta_row, alpha_row
    
    switch algo_tilted
        
        case 1
            [deltay, ky, dy, beta_row, alpha_row, diff_k1, diff_k2, no_conv, eps_fact1, eps_fact2] = algo3ph_stepthree_it( x_phase0, k, raw_shift, phase_opt, deltay, ky, dy, contr, decomp_LU ,...
                deltax, kx, dx, algo_vib, beta_row, alpha_row, corr, verbose, crit2, corr_str, diff_k1, diff_k2, epsilon1, epsilon, h, k_max, delta, first_delta_to_zero);
    end
    
    if sum(no_conv) ~= 0 % no convergence for at least one phase shift
        
        k = k+1;
    else
        % %         if k > 1 % must do at least 2 iterations to compare
        no_convergence = 0; % total convergence, stop
        switch algo_tilted
            case 1 % tilted algo
                fprintf('Threshold is %.3g (in degree) for transl. phase-shift\n', eps_fact1)
                fprintf('Threshold is %.3g (in degree) for tilt phase-shift\n', eps_fact2)
            case 0 % normal algo
                fprintf('Threshold is %.3g (in degree)\n', eps_fact)
        end
        % %         else % must do at least 2 iterations to compare
        % %             k = k+1;
        % %         end
    end
    
    t11 = toc(t11);
    fprintf('This iteration last %.3g sec. \n', t11);
end

if method ~= 1 % not with 2D matrices
    switch algo_tilted
        %         case 0 % normal algo
        case 1
            delta=squeeze(mean(mean(Delta,1),2));
            diff_k = diff_k1; % for output
            fprintf('\n%g\n', delta);
        case 0
            disp('delta   diff');
            dd=delta(:, k);
            disp([dd [0;diff(dd)]]);
            %             fprintf('\n%g\n', delta(:, k)); % do not put 'end' for cases when it does not converge well and the best value is used
    end
    
    phase_opt = reshape(phase_opt, size(contr,1),size(contr,2));
    amp = reshape(amp, size(contr,1),size(contr,2));
    err = reshape(err, size(contr,1),size(contr,2));
else
    delta=0; diff_k=0;
end

if ramp_phshft
    amp = amp/sinc((x_phase0(2) - x_phase0(1))/2);
end

divers_vect = {delta(:,end), deltax(:,:,end), kx(:,end), dx(:,end), beta_col(:,:,end), alpha_col(:,:,end), ...
    deltay(:,:,end), ky(:,end), dy(:,end), beta_row(:,:,end), alpha_row(:,:,end)};

end


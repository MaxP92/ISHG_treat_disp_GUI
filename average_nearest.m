function tmp = average_nearest( contr, coef_centre, coef_cross, coef_diag , mean_meth, sigma_med, real_sigma, skip_asking_coeff_median, answer_med_sz, circ_avg)
% average_nearest
%
% 2016.6.9 :  Maxime PINSARD
%
% To average on nearest neighbours
%
% tmp = average_nearest( contr, coef_centre, coef_cross, coef_diag, mean_meth, sigma_med )

% contr(isnan(
do_px_px_avg = 0; % init
if circ_avg
    addpath('C:\Users\pc\Documents\These\codes Matlab\Variety\CircStat2012a'); % fileexchange
    if ((mean_meth == 0 && sigma_med > 0) || (mean_meth > 1 && (coef_diag ~= coef_cross || coef_centre~= coef_cross || coef_centre~= coef_diag))) % weight
        fprintf(2, 'weighted circular avg is not implemented: will do not weighted ! \n');
    end    
    if circ_avg == 1 % in mulitple of pi [-1, 1]
        fact = pi;
    elseif circ_avg == 3 % already in deg
        fact = pi/180;
    elseif circ_avg == 2 % already in rad
        fact = 1;
    end
%         [M, N] = size(contr);
    contr = contr*fact;
    do_px_px_avg = 1;
        
else % no circular avg
    
    if mean_meth % FAST meth : mean
        tmp = imfilter(contr, [coef_diag, coef_cross, coef_diag; ...
            coef_cross, coef_centre, coef_cross; coef_diag, coef_cross, coef_diag])/(coef_centre + 4*coef_cross + 4*coef_diag);

    else % median, slower
        tmp = contr*0;

        for k = 1: size(contr, 3) % loop on third dimension

            %         sure_meth = 1; % malab method
            if sigma_med == 0 % % no weight
                %((coef_cross == 0 && coef_diag == 0) || sure_meth)
                % filter neighbors, with weights of 1 around
                cond = k <=1 && ~skip_asking_coeff_median;
                if cond
    % %                 prompt = {'Size of median filter along lines', 'Size of median filter along column'};
    % %                 dlg_title = 'Size of med. filter';
    % %                 num_lines = 1;
    % %                 % Valeurs par défaut
    % %                 def = {'3','3'};
    % %                 answer_med_sz = inputdlg(prompt,dlg_title,num_lines,def);
                    % Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
                    if ~isempty(answer_med_sz)
    %                     medx = str2double(cell2mat(answer(1)));
    %                     medy = str2double(cell2mat(answer(2)));
                    a=num2cell(str2num(answer_med_sz)); [medx, medy]=deal(a{:}); %#ok<ST2NM>
                    else
                        cond=0;
                    end
                end

                if ~cond % skip_asking_coeff_median
                    medx = 3; medy = 3;
                    if k <=1
                        fprintf(2, 'In average_nearest : median is taken on a 3x3 window with equal weight ! \n');
                    end
                end
                tmp(:,:,k) = medfilt2(contr(:,:,k), [medx, medy]);

            elseif  sigma_med > 0
                if k == 1
                    addpath('C:\Users\pc\Documents\These\codes Matlab\Codes_I-SHG\MP\jointWMF x64');
                    fprintf('Using WMF median filter (radius %d sigma %d)!\n', sigma_med, real_sigma);
                end
                tmp(:,:,k) = jointWMF(contr(:,:,k),contr(:,:,k),sigma_med,real_sigma,256,256,1,'exp');
                % jointWMF(I,F,r,sigma,nI,nF,iter,weightType,mask)
            end
        end

        if sigma_med == -1
            do_px_px_avg = 1;
        end
    end
end
if do_px_px_avg
    % SLOW
    hb = msgbox('Wait until the slow avg is being done !');
    tmp = contr*0;
    for posx=1:size(contr,1) % dimension Y
        for posy=1:size(contr,2) % dimension X
            %                 t1 = tic; %!!!
            % Moyennage sur une région 3X3 si le pixel n'est pas aux bords de l'image
            % MP : possiblement améliorable
            centre = squeeze(contr(posx,posy,:)); % if 3D, last dim will come to the 1st
%             if size(contr, 3) > 1; dim = 2; else; dim = 1; end
            dim = 2;
            if (posx~=1 && posx~=size(contr,1) && posy~=1 && posy~=size(contr,2)) % if not on the edges

                sides = [squeeze(contr(posx-1,posy,:)), squeeze(contr(posx+1,posy,:)) , squeeze(contr(posx,posy-1,:)), squeeze(contr(posx,posy+1,:))];
                corners = [squeeze(contr(posx-1,posy-1,:)), squeeze(contr(posx-1,posy+1,:)), squeeze(contr(posx+1,posy-1,:)), squeeze(contr(posx+1,posy+1,:))];
                if circ_avg
                   if  mean_meth == 0 % median
                       tmp(posx, posy,:) =circ_median([centre, sides, corners] , dim);
                   else % mean
                       tmp(posx, posy,:) =circ_mean([centre, sides, corners], [] ,dim);
                   end
                else
                    tmp(posx, posy,:) = median([repmat(sides, 1, coef_cross), repmat(corners, 1, coef_diag), repmat(centre, 1, coef_centre)],dim); % /(coef_centre + 4*coef_cross + 4*coef_diag)
                end
            else
                tmp(posx, posy,:) = centre; % tmp is a (180/diff_phase)x1 double vector
            end
        end
    end
    close(hb)
end
if circ_avg
    tmp = tmp/fact; %reshape(, M, N);
end
end


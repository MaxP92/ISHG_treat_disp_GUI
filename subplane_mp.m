function   [IM, x, y]=subplane_mp(Z, x, y, coplanar_points, h, was_unwrap, surf_fit_tilt_auto_bacth, ampflag)
%SUBtract selected PLANE, or remove slope (detilt) from data matrix
%
%Call:
%           IM=subplane(Z)
%Input:
%           Z = (double) data matrix
%Output:
%           IM = Z - plane     (the SUBtracted PLANE is defined by manually selected points)
%
%2017.3.28 Maxime PINSARD
%==============================

use_median = 0; % use a median to avg the pointed points did not give good results
nm ='';
switch coplanar_points
    
    %   Vassili Pastushenko	March	2005
    %
    % modified 2017.3.28 Maxime PINSARD
    
    case 1 % % coplanar_points
        
        if (sum(abs(x)) == 0 && sum(abs(y)) == 0)
            axes(h);
            imagesc(Z); axis image; colormap(hsv)
            % shg
            set(gca,'fontsize',15)
            %select at least three (x,y) points which define a plane in 3D
            title('Click at least three coplanar points')
            [x,y]=getpts(gca);
            
            if numel(x)<3
                title('More points please');
                [xx,yy]=getpts(gca);
                x=[x;xx];y=[y;yy];
            end
        end
        
        x=round(x);
        y=round(y);
        PONT=numel(x);
        ZIN=ones(PONT,1);
        M=[x y ZIN];
        for i=1:PONT
            switch use_median        
                case 1
                    if y(i)+1 <= size(Z, 1)
                        if y(i)-1 >= 1
                            ymax = y(i)+1;
                        else
                            if y(i) >= 1
                                ymax = y(i)+2;
                            else
                                ymax = 3;
                            end
                        end
                    else
                        if y(i) <= size(Z, 1)
                            ymax = y(i);
                        else
                            ymax = size(Z, 1);
                        end
                    end
                    
                    if x(i)+1 <= size(Z, 2)
                        if x(i)-1 >= 1
                            xmax = x(i)+1;
                        else
                            if x(i) >= 1
                                xmax = x(i)+2;
                            else
                                xmax = 3;
                            end
                        end
                    else
                        if x(i) <= size(Z, 2)
                            xmax = x(i);
                        else
                            xmax = size(Z, 2);
                        end
                    end
                    median([Z(ymax, xmax), Z(ymax - 2,xmax - 2), Z(ymax - 1,xmax ), Z(ymax ,xmax - 1), Z(ymax - 1,xmax - 1), Z(ymax - 1,xmax - 2), Z(ymax - 2,xmax - 1), Z(ymax ,xmax - 2), Z(ymax - 2,xmax )] );
                case 0
                    ZIN(i)=Z(y(i),x(i));
            end
        end
        
        V=M\ZIN;
        [VY,VX]=size(Z);
        [X,Y]=meshgrid(1:VX,1:VY);
        BAS=V(1)*X+V(2)*Y+V(3);
        IM=Z-BAS;
        
        axes(h);
        
    case 0 % tilt from a ROI that dictates the orientation

         addpath('C:\Users\pc\Documents\These\codes Matlab\Codes_I-SHG\MP\jointWMF x64');
        
        % averaging to help the fit of the surface
        % you need to specify your exp. width (generally 3 is enough)
        % also the center of one peak (e.g. -0.4 for -0.4pi), its width 
        % and the distance between the two peaks ( one peak will be canceled out)
        cstr_onedir = '0'; % dflt
        deg_fit = '1'; % dflt
        if ampflag; avgd = '0'; else; avgd = '1'; end
        if surf_fit_tilt_auto_bacth < 0 % short version
            surf_fit_tilt_auto_bacth = -surf_fit_tilt_auto_bacth;
            switch surf_fit_tilt_auto_bacth
                case 1;  cstr_onedir = '0'; 
                case 11; cstr_onedir = '1';
                case 12; cstr_onedir = '2';
                case 2;  deg_fit = '2';
                case 21; cstr_onedir = '1'; deg_fit = '2';
                case 22; cstr_onedir = '2'; deg_fit = '2';
            end
            surf_fit_tilt_auto_bacth = 1;
            short_vers = 1; ret_fit = 1; wrapfitmat_dflt = '0';
        else; short_vers = 0;ret_fit = 0; wrapfitmat_dflt = '0';
        end
        
        if (surf_fit_tilt_auto_bacth < -2 || surf_fit_tilt_auto_bacth > 2 ) % anything but 0,1 2 
            surf_fit_tilt_auto_bacth = 0;
            autosave = '1';
        else; autosave = '0'; 
        end
        def = {avgd, '3', '0.5', '0', '1', '0', '0', wrapfitmat_dflt , cstr_onedir, deg_fit, autosave, ''};
        if (surf_fit_tilt_auto_bacth >=1 &&  surf_fit_tilt_auto_bacth <= 2 ) % batch auto no ask
            def{end-2} = num2str(surf_fit_tilt_auto_bacth);
            def{end-1} = '1'; % dflt save
            if ~was_unwrap
                def(end+1) = {'0'};
            end
            answer = def;
        else
            prompt = {'Do you want to average with a gaussian filter ?' ; 'Sigma of gaussian filter (>0)' ; 'Center of the peak (in /pi, usually 0.5) ?' ; ...
            'width at bottom of this peak (in /pi, e.g. 1.0 for large, 0 for no !) ?'; 'Coeff to multiply the fit'; 'Distance between peaks (in /pi, usually 1.0, negative for exclude other peak, 0 for no !)'; ...
            'Crop to ROI?'; 'wrap fit mat in [-1, 1]'; 'Constrain fit surface to (X =1, Y =2) other no ?'; 'Parabola(2) or linear (1)'; 'autosave of fit fig.'; 'coeffs pred'};
            
            if ~was_unwrap
                prompt(end+1) = {'phase outside [-pi, pi] (unwrapped)'};
                def(end+1) = {'0'};
            end
            answer = inputdlg(prompt,'Params fit surface',1,def);
%             % Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
        end
        if (length(x) == 1 && x ~= 0); sx = x; % imposed
        else; sx = size(Z, 2);
        end
        if (length(y) == 1 && y ~= 0); sy = y;% imposed
        else; sy = size(Z, 1);
        end
        
        xx2 = 1:sx; % before crop
        yy2 = 1:sy;
        [X2, Y2] = meshgrid(xx2, yy2);
        if ~isempty(answer)
            avg_flag = str2double(cell2mat(answer(1)));
            sigma = str2double(cell2mat(answer(2)));
            center_peak = str2double(cell2mat(answer(3)));
            width_peak = str2double(cell2mat(answer(4)));
            coeff_fit_sup = str2double(cell2mat(answer(5)));
            dist_peak = str2double(cell2mat(answer(6)));
            crop_ROI = str2double(answer{7});
            wrap_fitmat = str2double(answer{8});
            constrain_surf_axis = str2double(answer{9});
            if (constrain_surf_axis < 1 || constrain_surf_axis >2)
                constrain_surf_axis  = 0;
            end
            degree_fit = str2double(answer{10});
            autosave_fit_fig = str2double(answer{11});
            coeffpred = answer{12};
            if (~isempty(coeffpred) || ~strcmp(coeffpred, '0'))
                coeffpred =str2num(coeffpred); %#ok<ST2NM>
            else
                coeffpred = [];
            end
            if ~was_unwrap
                was_unwrap= str2double(answer{end});
            end
        else
            avg_flag = 0; coeffpred = [];
        end
        if isempty(coeffpred)
            if crop_ROI
    %             axes(h);
                figure; h =axes;
                imagesc(Z); 
                if (size(Z, 1) > 10*size(Z, 2) || size(Z, 1) < 10*size(Z, 2))
                    axis auto;
                else
                    axis image;
                end
                title('Double-click to confirm your crop');
                colormap(hsv);  colorbar; set(gca,'fontsize',15)
                Z1=imcrop(h);
            else
                Z1 =Z;
            end

            if avg_flag
                Z_avg = jointWMF(Z1, Z1, sigma, sigma, 256, 256, 1, 'exp');
                if ~short_vers
                    figure(131); subplot(2,1,1); imagesc(Z1);  colormap(hsv); colorbar; axis image; title('Normal');%caxis( [-sat_value, sat_value]);
                    subplot(2,1,2); imagesc(Z_avg);  colormap(hsv); colorbar; axis image; title('gaussian avg'); %caxis( [-sat_value, sat_value]);
                end
            else
                Z_avg = Z1;
            end

            if (dist_peak~=0 && width_peak~=0) % otherwise no
                borne_sup = center_peak+width_peak/2;
                
                borne_inf = center_peak-width_peak/2;
                if dist_peak>=0 % do not exclude other peak
                    if borne_inf < -1
                        borne_inf = borne_inf + 2;
                    end
                    if borne_sup > 1
                        borne_sup = borne_sup - 2;
                    end
                end

        %         if borne_inf > borne_sup
        %             borne_inf_disp = -1; borne_sup_disp = 1;
        %         else
        %             borne_sup_disp = borne_sup;
        %             borne_inf_disp = borne_inf;
        %             if borne_sup + dist_peak > 1
        %                 borne_sup_disp = 1;
        %             end
        %             
        %         end
            end
            if was_unwrap
                im_calib =Z_avg; 
            else
                dat2=Z_avg; 
                if (dist_peak~=0 && width_peak~=0) % otherwise no
                    if dist_peak<0 % exclude other peak
                        correct=dat2(dat2 > borne_sup | dat2 < borne_inf);
                        dat2(dat2 <= borne_sup & dat2 >= borne_inf)= median(correct);
                    else
                        dat2(dat2 < borne_sup & dat2 > borne_inf) = dat2(dat2 < borne_sup & dat2 > borne_inf)+ dist_peak; 
                    end
                    dat2(dat2>1) = dat2(dat2>1) - 2 ; dat2(dat2<-1) = dat2(dat2<-1) + 2 ;
                end
                if ~short_vers
                    figure(66);imagesc(dat2); axis image;colorbar;   colormap(jet); set(gca, 'CLim', [min(min(dat2)), max(max(dat2))]); 
                    title('All data recast in one interval');
                end
                im_calib = dat2;
            end
            xx = 1:size(im_calib, 2);
            yy = 1:size(im_calib, 1);

            % [X, Y] = meshgrid(xx, yy);

            % Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2;
            
            if degree_fit <= 1 % linear

                % fit surface

                [XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im_calib);
                boundX = Inf; boundY = Inf;
                if constrain_surf_axis == 1 % const X
                    boundY = 0;
                elseif constrain_surf_axis == 2 % const Y
                    boundX = 0;
                end

                [fitobject, gof] = fit([XOut, YOut], ZOut, 'poly11',...
                        'Lower',[-Inf, -boundX, -boundY],...
                       'Upper',[Inf, boundX, boundY]) ;

                fprintf('R squared of fit double linear = %f \n', gof.adjrsquare);
                % % takes only predominant values
%                 nb = 50; %
%                 ZOut(ZOut < prctile(fliplr(ZOut),50-nb/2) | ZOut > prctile(fliplr(ZOut),50+nb/2)) = prctile(fliplr(ZOut),nb);
        %         figure;imagesc(ddd); colorbar;
        %         set(gca, 'CLim', [-1,1]); colormap(hsv);

                coefX =fitobject.p10; coefY = fitobject.p01;
                if constrain_surf_axis == 1 % const X
                    coefY = 0;
                    nm = 'X lin';
                elseif constrain_surf_axis == 2 % const Y
                    coefX = 0;
                    nm = 'Y lin';
                else
                    nm = '2D lin';
                end
                mat_fit = coefX*X2 + coefY*Y2 + fitobject.p00;
                figure(132); hh=axes; plot3(hh, XOut, YOut, ZOut, 'd');  hold on; plot(fitobject); colormap('hot');colorbar;title(sprintf('Every point in 3D and the surf fit (r2 %.2f) %s', gof.adjrsquare, nm)); % do not put hh for last
                 hold off;
            else % parabola
                nm = 'parabola';
                addpath('C:\Users\pc\Documents\These\codes Matlab\Variety');
    % %             xx = 1:size(im_calib, 2);
    % %             yy = 1:size(im_calib, 1);
                mat_fit = parabola_fit_function(im_calib, 0, xx2, yy2, 0); % 0 for not uint16 data, 0 for no pos val imposed
            end
        else % predef coeffs for building fit mat
            coefX = coeffpred(1);coefY=0; 
            if (length(coeffpred) > 1 && constrain_surf_axis == 0); coefY = coeffpred(2); end
            if constrain_surf_axis == 2;  coefY = coeffpred(1);coefX=0;end 
            mat_fit = coefX*X2 + coefY*Y2;
        end
        if ~ret_fit % return im corrected not fit
            IM = Z - coeff_fit_sup*mat_fit;
        else
            IM = mat_fit;
            
        end
        if wrap_fitmat
          mat_fit= wrapToPi(mat_fit*pi)/pi ;
        end
% %         if ~short_vers
            figure(55);imagesc(mat_fit); axis image; colorbar; title(sprintf('The fit matrix %s, in a map', nm));
            %         coeff_fit_sup = 2;
            if autosave_fit_fig
                savefig(55, 'fit_mat.fig','compact');
            end   
% %         end
end

%% plot

if ~ampflag
    IM(IM>1) = IM(IM>1) -2 ; IM(IM<-1) = IM(IM<-1) +2 ;
end
if ~short_vers
    figure(133);
    if size(Z, 2) > 16/9*size(Z, 1)
        subplot(2,1,1); 
    else
        subplot(1,2,1); 
    end
    imagesc(Z); colorbar;  axis image; title('sans corr');
    if size(Z, 2) > 4/3*size(Z, 1)
        subplot(2,1,2); 
    else
        subplot(1,2,2); 
    end
    imagesc(IM);  axis image; title('de-tilted');

    if ~ampflag
        colormap('hsv');
    else
        colormap('parula');
    end
    colorbar; 
    %         figure(103); imagesc(IM); colorbar; colormap(hsv); title('de-tilted');
    %                 set(gca, 'CLim', [-1,-0.2]); colormap(hsv);

    figure(134) ; subplot(1,2,1); 
    histogram(Z(:), linspace(-1,1,300));  title('sans corr');
    subplot(1,2,2); histogram(IM(:), linspace(-1,1,300));  title('de-tilted');
end


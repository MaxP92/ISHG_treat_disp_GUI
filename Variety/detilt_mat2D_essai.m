x=1:500; y=x; [X, Y] = meshgrid(x, y); 

%% perfect fit

RR1 = X*cos(33/180*pi) + Y*sin(33/180*pi);

        
%% noisy fit

p = 0.1; RR2 = X*cos(33/180*pi) + Y*sin(33/180*pi) + mean2(X)/10*cos((X*sin(33/180*pi) + Y*cos(33/180*pi))*2*pi);

%% plot

RR = RR2;

figure(1);imagesc(RR); colorbar;

[XOut, YOut, ZOut] = prepareSurfaceData(x, y, RR);
        
[fitobject, gof] = fit([XOut, YOut], ZOut, 'poly11') ;

fprintf('R squared of fit double linear = %f \n', gof.adjrsquare);

figure(2); plot3( XOut, YOut, ZOut, 'd');  hold on; plot(fitobject); % 

%% avg




% NN = 3;
% coef_diag = 3; coef_cross = 3; coef_centre = 3; matdiv = ones(NN, NN);
% 
% ddd = imfilter(dat2, matdiv)/size(matdiv, 1).^2;
%     figure(1);imagesc(ddd);         set(gca, 'CLim', [0.16,1]); colormap(hsv);
% colorbar;

        %
                addpath('C:\Users\pc\Documents\These\codes Matlab\weighted_median_filter\matlab_interface\x64');

        diff_circ = dat;
        
        prompt = {'Do you want to average with a guassian filter ?' ; 'Sigma of gaussian filter (>0)'};
        dlg_title = 'Size of med. filter';
        num_lines = 1;
        % Valeurs par défaut
        def = {'1','3'};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        % Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
        if ~isempty(answer)
            avg_flag = str2double(cell2mat(answer(1)));
            sigma = str2double(cell2mat(answer(2)));
        else
            avg_flag = 0;
        end
    
    if avg_flag
        diff_circ_avg = jointWMF(diff_circ, diff_circ, sigma, sigma,256,256,1,'exp');
        figure(121); subplot(2,1,2); imagesc(diff_circ_avg);  colormap(hsv); colorbar; axis image; %caxis( [-sat_value, sat_value]);
        subplot(2,1,1); imagesc(diff_circ);  colormap(hsv); colorbar; axis image; %caxis( [-sat_value, sat_value]);
    end
    
    dat2=diff_circ_avg; dat2(dat2<0.16 & dat2>-0.83) = dat2(dat2<0.16 & dat2>-0.83)+0.85; dat2(dat2>1) = dat2(dat2>1) -2 ; dat2(dat2<-1) = dat2(dat2<-1) +2 ;
        figure(1);imagesc(dat2); colorbar;   colormap(hsv); set(gca, 'CLim', [0.16,1]); 

    
    im_calib = dat2;
xx = 1:size(im_calib, 2);
        yy = 1:size(im_calib, 1);
        
        % [X, Y] = meshgrid(xx, yy);
        
        % Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2;
        
        % fit surface
        
        [XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im_calib);
        
        [fitobject, gof] = fit([XOut, YOut], ZOut, 'poly11') ;
        
        fprintf('R squared of fit double linear = %f \n', gof.adjrsquare);
        
        figure(1);imagesc(ddd); colorbar;
        set(gca, 'CLim', [-1,1]); colormap(hsv);
        figure(2); plot3( XOut, YOut, ZOut, 'd');  hold on; plot(fitobject); % 
        
        xx2 = 1:size(dat, 2);
        yy2 = 1:size(dat, 1);
        [X2, Y2] = meshgrid(xx2, yy2);
        
        mat_fit = fitobject.p10*X2 + fitobject.p01*Y2; % + fitobject.p00 
        figure(3);imagesc(mat_fit); colorbar;
        
        coeff = 2;
        
        IM = dat - coeff*mat_fit;
%         x=fitobject; y= gof;

%% plot
IM(IM>1) = IM(IM>1) -2 ; IM(IM<-1) = IM(IM<-1) +2 ;

        figure(101);  
        subplot(1,2,1); imagesc(dat); colorbar; colormap(hsv); title('sans corr'); 
        subplot(1,2,2); imagesc(IM); colorbar; colormap(hsv); title('de-tilted'); 
        
        figure(103); imagesc(IM); colorbar; colormap(hsv); title('de-tilted'); 
                set(gca, 'CLim', [-1,-0.2]); colormap(hsv);
        
figure(102) ; subplot(1,2,1); histogram(dat(:), linspace(-1,1,500));  title('sans corr'); 
        subplot(1,2,2); histogram(IM(:), linspace(-1,1,500));  title('de-tilted');
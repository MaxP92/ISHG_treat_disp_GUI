function [ x, y, complete, test, model, r2 ] = analyse_cos_ISHG( contr, img_3D, resx, resy, rect, crop, xTitle_dflt, yTitle_dflt, screensize, fact, left_offset_fig, top_offset_fig, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, short_version, x_phase, double_contrast, contr_mode, cmap_redgreen)
% [ x, y, complete, test, model, r2 ] = analyse_cos_ISHG( contr, img_3D, resx, resy, rect, start_phase, diff_phase, crop, screensize, fact, left_offset_fig, top_offset_fig, ...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, short_version, x_phase )
%
% 2013 CA Couture
%
% 2015-10 edited by Maxime PINSARD
%
%  To simply do the analyze in cosinus

x=0; y=0; complete=0; %test=0; model=0; r2=0; % for nothing
ask_meth_avg = 0;

if (~short_version && ask_meth_avg)
    
    val=menu('Measure by','Signal average','Abs of the signal average', 'Median');
else
    val = 3; % median
end

%     test=zeros(1, size(contr,3));
%
%     for i = 1:size(contr,3)
%         test(i) = mean(mean(contr(:,:,i)));
%         if val == 2
%             test(i) = abs(test(i));
%         end
%     end

% Affichage graphique (1ère image et région croppée si il y a et profil d'intensité vs la phase)
tit_cos_beh = 'Analyze of the cosinus behaviour';

% outside fig
if ~short_version
    
    hf = figure(20);
    set(hf,'outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
        screensize(3)*fact screensize(4)*fact]);
    
    h = axes;
%subplot(2,1,1); % plot of the fit, on the same figure
    
    A1 = img_3D(:, :, 1);
    A2 = img_3D(:, :, 2);
    complete = single(A2(:,:))-single(A1(:,:));
    x = linspace(0,size(A1,2)*resx,size(A1,2));
    y = linspace(0,size(A1,1)*resy,size(A1,1));
    cmap = cmap_redgreen;
    %imagesc(x,y,A1)
    %colormap(gray)
    
    draw_plots_ISHG( 0, 0, complete, x, y, h,...
        xTitle_dflt, yTitle_dflt, tit_cos_beh, cmap, 0, '', 0, 0, ...
        axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz );
    %  To plot the figures of the ISHG program : hist or image figure
    % see the help of this function for entry vars
    
    % Maximum = max(max(abs(contr(:,:,1))));
    % caxis([-1*contrast*Maximum contrast*Maximum])
    if crop==2
        rect(1)=rect(1)*resx;
        rect(3)=rect(3)*resx;
        rect(2)=rect(2)*resy;
        rect(4)=rect(4)*resy;
        rectangle('Position',rect,'LineWidth',5,'EdgeColor','y')
    end
end

if short_version
    hf = figure(19);
    set(hf, 'outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
        screensize(3)*fact screensize(4)*fact]);
    %     h=axes;  
    ax2=gca;
else
    
    crop_img = 2 - menu('Crop img ?', 'ROI rect', 'whole img', 'ROI 3Dimg', 'whole 3Dimg');

    if (abs(crop_img) == 1) % 2 or 0 is no crop
        title('Double-click to confirm' );
        [~, rect] =imcrop(hf);
        title(tit_cos_beh );
        xmin = rect(1); ymin = rect(2); width = rect(3); height = rect(4);
        contr = contr(round(ymin/resy): round(ymin/resy + height/resy - 1), ...
            round(xmin/resx): round(xmin/resx + width/resx - 1), :);
        rectangle(gca, 'Position',rect,'EdgeColor','y','LineWidth',3)
    end

    subplot(2,1,1,h)
    if (crop_img<0) 
        [x,y,z] = meshgrid(1:1000,1:200,1:12);
        figure(63);scatter3(z(:),x(:),y(:),15,contr(:));view(3,178);%h,
        axis image;dar=get(gca,'DataAspectRatio'); set(gca,'DataAspectRatio', [dar(1), dar(2), min(8,size(contr,1)/size(contr,3))]); 
        pbar=get(gca,'PlotBoxAspectRatio'); set(gca,'PlotBoxAspectRatio', [pbar(1),pbar(2),min(4,size(contr,1)/size(contr,3))]);
        colormap(cmap); colorbar; xlabel('# of ctr frame');
    end
    figure(hf);
    ax2=subplot(2,1,2); % create new subplot
end

test = zeros(size(contr, 3), 1);
for k = 1:size(contr, 3)
    contr1 = contr(:,:,k);
    % Valeur moyenne de la région
    test(k) = mean(contr1(~isnan(contr1)));
    if val == 2
        test(k) = abs(test);
    end
    if val == 3 % median
        test(k) =median(contr1(~isnan(contr1)));
    end
end

test = squeeze(test);
% OLD VERSION

plot(x_phase,test,'.',...
    'MarkerSize',15) % 'bx');
hold on
% Fit automatisé d'une fonction sinusoidale
% On mesure l'amplitude et le décalage (en y) du cos obtenu
amp = (max(test)-min(test))/2;
mid = mean(test) ; 
%mid = (max(test)+min(test))/2; % in C.A.'s code : depends on the
%length of test, but not the same !

t = x_phase';
X = ones(size(t,1),3);
X(:,2) = mid+amp*cos(t/180*pi);
X(:,3) = mid+amp*sin(t/180*pi);
beta = X\test;
new_amp = sqrt(beta(2)^2+beta(3)^2);
phi_opt = atan2(-beta(3),beta(2))/pi;

amp_opt = amp*new_amp;  % *new_amp

mid_model = mid - mean(amp_opt*cos(t/180*pi+phi_opt*pi));%trapz(test)/length(test); % mean(test) - 

model = mid_model + amp_opt*cos(t/180*pi+phi_opt*pi);


hh=plot(ax2,x_phase, model, 'LineWidth',2 ); %,'r')
ss_tot = sum((test-mean(test)).^2);
% ss_reg = sum((model-mean(tmp)).^2);
ss_res = sum((test-model).^2);
r2 = 1-ss_res/ss_tot;
if short_version
    dim = [.15+rand(1)*0.2, .4+rand(1)*0.2, .3, .3];
else
    dim = [.15+rand(1)*0.2, .1+rand(1)*0.2, .3, .3];
end
annotation(hf, 'textbox',dim, 'String', sprintf('r2 = %.3f', r2),'FitBoxToText','on', 'Color', hh.Color);

legend('Data', 'Model'); legend('off'); legend('show');
end


function [ ph_hist, second_hist, data_tot, ylabel_hist3, title_hist3, nbinsx ] = hist_tot_ISHG( phasee, img, int_x, int_y, phi_mat_default, Counts, Titre4, Titre5, Yaxis1,...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig )
% [ ph_hist, second_hist, data_tot, ylabel_hist3, title_hist3, nbinsx ] = hist_tot_ISHG( phasee, img, int_x, int_y, phi_mat_default, Counts, Titre4, Titre5, Yaxis1,...
% axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig )
%
% 2015-10 edited by Maxime PINSARD
%
%  To calculate and plot the histograms (2D and 3D)

%% Construction des vecteurs pour les histogrammes

[ ph_hist, second_hist ] = init_hist_ISHG( phasee, img );
% Calculate vectors for histograms

%% Histogramme 2D

offset_pi2 = 0; % no change at first

[~, ~, ~] = hist2D_ISHG( fact, left_offset_fig, ...
top_offset_fig, ph_hist, int_x, phi_mat_default, Counts, Titre4, ...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, screensize, offset_pi2 );

%     % Fit gaussien sur l'histo, encore en version bêta je crois, mais
%     % marchais OK dernière fois que je l'ai essayé. Voir à la toute fin,
%     % partie commentée pour un fit à 2 gaussiennes ...
%     [val pos] = max(hist_ydata);
%     [val2 pos2] = min(abs(hist_ydata(1:pos)-val/2));
%     [val3 pos3] = min(abs(hist_ydata(pos+1:end)-val/2));
%     FWHM = hist_xdata(pos+pos3)-hist_xdata(pos2);
%     a1 = val;
%     b1 = hist_xdata(pos);
%     c1 = FWHM/(2*sqrt(2*log(2)));
%
%     gauss = @(a, b, c, x) a*exp(-(x-b).^2/(2*c^2));
% %     a1=3842; b1=-0.1829; c1=0.06129;
% %     a2=3941; b2=0.02606; c2=0.1109;
% %     x=hist_xdata;
% %     gauss2 = @(a1,b1,c1,a2,b2,c2,x) a1*exp(-((x-b1)/c1).^2)+a2*exp(-((x-b2)/c2).^2);
%
%     g = fittype(gauss);
%     f1 = fit(hist_xdata',hist_ydata',g,'Startpoint',[a1 b1 c1]);
%
% %     c1=c1/sqrt(2);
% %     c2=c2/sqrt(2);
%
%     string = sprintf('%.2f',  round(f1.c*100)/100);
% %     string = sprintf('%.2f',  round(c1*100)/100);
% %     string2 = sprintf('%.2f',  round(c2*100)/100);
%
% %     hold on; plot(hist_xdata,gauss2,'Color','r','LineWidth',3);
% %     h=legend('Experimental histogram',['Gaussian fit - \boldmath$(\delta\varphi^2)^{1/2}= ' string '\pi$ ; \boldmath$(\delta\varphi^2)^{1/2}= ' string2 '\pi$']);
%     hold on; plot(hist_xdata,gauss(f1.a,f1.b,f1.c,hist_xdata),'Color','r','LineWidth',3);
%     h=legend('Experimental histogram',['Gaussian fit - \boldmath$(\delta\varphi^2)^{1/2}= ' string '\pi$']);
%     set(h,'Interpreter','Latex','Location','Best')


%% Histogramme 3D

[ data_tot, ylabel_hist3, title_hist3, nbinsx,~ ] = hist_3D_ISHG( ph_hist, second_hist, int_x, int_y, Titre5, Yaxis1, Counts, phi_mat_default,...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig, offset_pi2 );
% Histogram 3D function

end


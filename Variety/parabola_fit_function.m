function mat_fit2 = parabola_fit_function(im_calib, uint16_data, xx2, yy2, pos_val)
% mat_fit2 = parabola_fit_function(im_calib)
%
% 2017.2.23
%
% Maxime PINSARD
%
% Does a parabola 2D fit on a 2D array

if uint16_data
    im_calib = (im_calib - min(min(im_calib)))*(2^16-1)/(max(max(im_calib)) - min(min(im_calib)));
    fprintf(2, 'warning, I cast the data in uint16 ! \n')
end
% grid vectors

xx = 1:size(im_calib, 2);
yy = 1:size(im_calib, 1);

% % [X, Y] = meshgrid(xx, yy);
[X2, Y2] = meshgrid(xx2, yy2);
% Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2; 

% fit surface

[XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im_calib);

[fitobject, gof] = fit([XOut, YOut], ZOut, 'poly22') ;

fprintf('R squared of fit = %f \n', gof.adjrsquare);
nm = 'parabola';

if  pos_val >= 0 % <0 don;t plot
    figure(132); hh=axes; plot3(hh, XOut, YOut, ZOut, 'd'); hold on; plot(fitobject);  
    colormap('hot');colorbar;title(sprintf('Every point in 3D and the surf fit (r2 %.2f) %s', gof.adjrsquare,nm)); % do not put hh for last
    % do not put hh for last
    hold off;
end
mat_fit = fitobject.p00 + fitobject.p10*X2 + fitobject.p01*Y2 + fitobject.p20*X2.^2 + fitobject.p11.*X2.*Y2 + fitobject.p02*Y2.^2;

if  pos_val
    % mat_fit2(mat_fit2 < 0) = 0;
    % mat_fit2 = mat_fit2 + 1;
    mat_fit2 = mat_fit + abs(min(min(mat_fit)))+1;
else
    mat_fit2 = mat_fit;
end


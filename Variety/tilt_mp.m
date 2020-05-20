im_calib = aa(:, 300:400);

xx = 1:size(im_calib, 2);
yy = 1:size(im_calib, 1);

% [X, Y] = meshgrid(xx, yy);

% Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2; 

% fit surface

[XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im_calib);

[fitobject, gof] = fit([XOut, YOut], ZOut, 'poly11') ;

fprintf('R squared of fit = %f \n', gof.adjrsquare);

figure; plot3( XOut, YOut, ZOut, 'd'); hold on; plot(fitobject);

xx2 = 1:size(aa, 2);
yy2 = 1:size(aa, 1);
[X2, Y2] = meshgrid(xx2, yy2);

mat_fit = fitobject.p00 + fitobject.p10*X2 + fitobject.p01*Y2;
aa2 = aa - mat_fit;
figure;subplot(2,1,1); imagesc(aa); axis image; colormap('hsv')
subplot(2,1,2); imagesc(aa2); axis image; colormap('hsv')

figure; subplot(2,1,1); histogram(aa, 500); 
subplot(2,1,2); histogram(aa2, 500); 
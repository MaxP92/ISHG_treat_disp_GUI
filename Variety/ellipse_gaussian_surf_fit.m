

% im_calib = (im_calib - min(min(im_calib)))*(2^16-1)/(max(max(im_calib)) - min(min(im_calib)));
im1 = imread('C:\Users\pc\Desktop\FFT_of_CH01-F_PhiR2f8bMap_z00.0.tif');
im_calib = double(im1);

t = linspace(0,2*pi,50);

BW = imbinarize(im_calib, 110);
stats= regionprops(BW,'Centroid',...
    'MajorAxisLength','MinorAxisLength', 'Orientation');
stats2 = struct2table(stats);
centers = stats2.Centroid;
diameters = mean([stats2.MajorAxisLength stats2.MinorAxisLength],2);
[MM, i_M] = max(diameters);
[MM2, i_M2] = max(diameters(diameters~=MM));
ss = stats(i_M);
a = ss.MajorAxisLength/2;
b = ss.MinorAxisLength/2;
Xc = ss.Centroid(1);
Yc = ss.Centroid(2);
phi = deg2rad(-ss.Orientation);
x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);

linwdth = 2;
figure(1); imagesc(im_calib); axis image; colorbar; colormap gray; 
hold on
plot(x,y,'r','Linewidth',linwdth)
% viscircles(centers(i_M,:),radii(i_M));
hold off

figure(2); imagesc(BW); axis image;  colormap gray;axis off;
hold on
plot(x,y,'r','Linewidth',linwdth)
% viscircles(centers(i_M,:),radii(i_M));
hold off

disp(b/a)

%% gaussian fit

% grid vectors

xx = 1:size(im_calib, 2);
yy = 1:size(im_calib, 1);

[X, Y] = meshgrid(xx, yy);

% Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2; 

% fit surface

[XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im_calib);

[fitobject, gof] = fit([XOut, YOut], ZOut, 'poly22') ;

fprintf('R squared of fit = %f \n', gof.adjrsquare);

figure; plot3( XOut, YOut, ZOut, 'd'); hold on; plot(fitobject); 

mat_fit = fitobject.p00 + fitobject.p10*X + fitobject.p01*Y + fitobject.p20*X.^2 + fitobject.p11.*X.*Y + fitobject.p02*Y.^2;

mat_fit2 = mat_fit;
% mat_fit2(mat_fit2 < 0) = 0;
% mat_fit2 = mat_fit2 + 1;
mat_fit2 = mat_fit2 + abs(min(min(mat_fit2)))+1;
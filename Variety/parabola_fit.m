clearvars
clc
close all

%% create a parabola in 2D

xx = 1:100;
yy = xx;

[X, Y] = meshgrid(xx, yy);

Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2; 

%% fit surface

[XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, Z);

fitobject = fit([XOut, YOut], ZOut, 'poly22') ;

plot3( XOut, YOut, ZOut, 'd'); hold on; plot(fitobject);

mat_fit = fitobject.p00 + fitobject.p10*X + fitobject.p01*Y + fitobject.p20*X.^2 + fitobject.p11.*X.*Y + fitobject.p02*Y.^2;

results = Z./mat_fit;

%% try with real image

%pcomp = 'D:\';
pcomp = 'C:\Users\pc';
path1 = fullfile(pcomp, 'Documents\These\Biomedical\Menisques Saint-Hyacinthe\Correction mosaique');
cd(path1);

[fname, fld, fltr] = uigetfile('data.tif');

im_data = double(imread(fname));

im_data = (im_data - min(min(im_data)))*(2^16-1)/(max(max(im_data)) - min(min(im_data)));

% calib

[fname_calib, fld, fltr] = uigetfile('calib.tif');

if fltr ~= 1
    %fprintf(2, 
    msgbox('No calib. file: the intensity profile won`t be corrected')
end

warning('OFF', 'all'); % because of Tiff warning
info=imfinfo(fname_calib); % disp a warning 'Division by zero when processing YResolution'


if (info(1).Width ~= size(im_data, 2) || info(1).Height ~= size(im_data, 1))
    msgbox('Calib. file is not the same size as the image file')
else
    msgbox('You should verify that the physical size of calibration is the same as your data ...')
end

warning('ON', 'all'); % because of Tiff warning

im = double(imread(fname_calib));

im = (im - min(min(im)))*(2^16-1)/(max(max(im)) - min(min(im)));

% grid vectors

xx = 1:size(im, 1);
yy = 1:size(im, 2);

[X, Y] = meshgrid(xx, yy);

% Z = -(X - 50).^2 - (Y - 50).^2 + 100^2/2; 

% fit surface

[XOut, YOut, ZOut] = prepareSurfaceData(xx, yy, im);

[fitobject, gof] = fit([XOut, YOut], ZOut, 'poly22') ;

fprintf('R squared of fit = %f \n', gof.adjrsquare);

figure(101); plot3( XOut, YOut, ZOut, 'd'); hold on; plot(fitobject);

mat_fit = fitobject.p00 + fitobject.p10*X + fitobject.p01*Y + fitobject.p20*X.^2 + fitobject.p11.*X.*Y + fitobject.p02*Y.^2;

mat_fit2 = mat_fit;
% mat_fit2(mat_fit2 < 0) = 0;
% mat_fit2 = mat_fit2 + 1;
mat_fit2 = mat_fit2 + abs(min(min(mat_fit2)))+1;


im_corr = im_data./mat_fit2;

im_corr(im_corr > mean(mean(im_corr))*2^16/50) = mean(mean(im_corr));

im_corr = (im_corr - min(min(im_corr)))*(2^16-1)/(max(max(im_corr)) - min(min(im_corr)));

im_corr_16 = uint16(im_corr);

imwrite(im_corr_16, 'out.tif');
imwrite(uint16(im_data), 'in_compare.tif');
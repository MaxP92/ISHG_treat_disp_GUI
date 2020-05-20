function [img_phase_out2, img_calib_out2, img_out2] = realign_pattern(img_phase, img_calib, img0)
% 2016.10.16 : aligning the ref and phase patterns by correlation
% Stephane Bancelin & Maxime Pinsard

 %% find the center and re-align images of calibration 
corr_calib = fftshift(ifft2((fft2(img_calib)).*conj(fft2(img_calib(end:-1:1, end:-1:1)))));

[c,ci]=max(corr_calib);
[~,max_x_calib]=max(c);
max_y_calib=ci(max_x_calib);

img_calib_out = img_calib(max(1, max_y_calib - round(size(img_calib, 1)/2)):min(size(img_calib, 1), max_y_calib + round(size(img_calib, 1)/2)),...
    max(1, max_x_calib - round(size(img_calib, 2)/2)):min(size(img_calib, 2), max_x_calib + round(size(img_calib, 2)/2)));

%% find the center and re-align images of phase
corr = fftshift(ifft2((fft2(img_phase)).*conj(fft2(img_phase(end:-1:1, end:-1:1)))));

[c,ci]=max(corr);
[~,max_x]=max(c);
max_y=ci(max_x);

img_phase_out = img_phase(max(1, max_y - round(size(img_phase, 1)/2)):min(size(img_phase, 1), max_y + round(size(img_phase, 1)/2)),...
    max(1, max_x - round(size(img_phase, 2)/2)):min(size(img_phase, 2), max_x + round(size(img_phase, 2)/2)));


img_out =img0(max(1, max_y - round(size(img_phase, 1)/2)):min(size(img_phase, 1), max_y + round(size(img_phase, 1)/2)),...
    max(1, max_x - round(size(img_phase, 2)/2)):min(size(img_phase, 2), max_x + round(size(img_phase, 2)/2)));

%% resize the phase image accordingly to calib size


img_phase_out2 = img_phase_out(max(1,  round(size(img_phase_out,1)/2) - min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))+1):min(size(img_phase_out,1) ,  round(size(img_phase_out,1)/2) + min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))),...
    max(1,  round(size(img_phase_out,2)/2) - min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))+1):min(size(img_phase_out,2) ,  round(size(img_phase_out,2)/2) + min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))));

img_out2 = img_out(max(1,  round(size(img_phase_out,1)/2) - min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))+1):min(size(img_phase_out,1) ,  round(size(img_phase_out,1)/2) + min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))),...
    max(1,  round(size(img_phase_out,2)/2) - min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))+1):min(size(img_phase_out,2) ,  round(size(img_phase_out,2)/2) + min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))));
       
img_calib_out2 = img_calib_out(max(1,  round(size(img_calib_out,1)/2) - min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))+1):min(size(img_calib_out,1) ,  round(size(img_calib_out,1)/2) + min(round(size(img_phase_out,1)/2), round(size(img_calib_out,1)/2))),...
    max(1,  round(size(img_calib_out,2)/2) - min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))+1):min(size(img_calib_out,2) ,  round(size(img_calib_out,2)/2) + min(round(size(img_phase_out,2)/2), round(size(img_calib_out,2)/2))));

if size(img_phase_out2,1) ~= size(img_calib_out2,1)
    img_phase_out0=img_phase_out2;
    img_calib_out0=img_calib_out2;
    img_out0 = img_out2;
    img_phase_out2 = img_phase_out0(1:end-floor(size(img_phase_out0,1)/max(size(img_phase_out0,1), size(img_calib_out0,1))), :);
    img_calib_out2 = img_calib_out0(1:end-floor(size(img_calib_out0,1)/max(size(img_phase_out0,1), size(img_calib_out0,1))), :);
    img_out2 = img_out0(1:end-floor(size(img_phase_out0,1)/max(size(img_phase_out0,1), size(img_calib_out0,1))), :);
end

if size(img_phase_out2,2) ~= size(img_calib_out2,2)
    img_phase_out0=img_phase_out2;
    img_calib_out0=img_calib_out2;
    img_out0 = img_out2;
    img_phase_out2 = img_phase_out0(:, 1:end-floor(size(img_phase_out0,2)/max(size(img_phase_out0,2), size(img_calib_out0,2))));
    img_calib_out2 = img_calib_out0(:, 1:end-floor(size(img_calib_out0,2)/max(size(img_phase_out0,2), size(img_calib_out0,2))));
    img_out2 = img_out0(:, 1:end-floor(size(img_phase_out0,2)/max(size(img_phase_out0,2), size(img_calib_out0,2))));
end

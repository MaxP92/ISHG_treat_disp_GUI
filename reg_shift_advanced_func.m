function [shiftv, im] = reg_shift_advanced_func(im0, lim1, off_shift, shift00)
% 
% im0 image to treat (array of numbers), [] if load
% lim1 limit the ROI to X first lines (sizeY dflt)
% off_shift : final imposed offset (0 dflt)
% % !!! transpose the image to find shift on Y (here it's on X) !!!
% shift00 is [] dflt, unless to just shift


if isempty(im0)
    [FileName,PathName,FilterIndex]=uigetfile('.tif');
    if FilterIndex
        im0 = imread(fullfile(PathName,FileName));
    else; return
    end
end

if (~isempty(lim1) && lim1> 0)
    off= size(im0,1)-lim1;
else
    off=0;
end

A = im0(1:2:end-off, :);
B = im0(2:2:end-off, :);
if isempty(shift00)
    B1 = B;
    A1 = A;
    shiftv =0;
    shiftvector2(1)=1;
    ct=0;
    while (ct<20 && (abs(shiftvector2(1)) >= 1))
        ct=ct+1;
        shiftvector2 = findshift(A1, B1,'proj'); % DIPimage 2.9
        sx = shiftvector2(1);
        if sx < 0
            B1 = B1(:,abs(round())+1:end);
            A1 = A1(:, 1:end+round(sx));
        else
            B1 = B1(:,1:end-round(sx));
            A1 = A1(:, abs(round(sx))+1:end);
        end
        shiftv=shiftv+sx;
        disp(sx)
    end

    shiftv = shiftv + off_shift;
else
    shiftv = shift00;
end

if shiftv < 0
    B2 = B(:,abs(round(shiftv))+1:end);
    A2 = A(:, 1:end+round(shiftv));
else
    B2 = B(:,1:end-round(shiftv));
    A2 = A(:, abs(round(shiftv))+1:end);
end

im=reshape([A2(:) B2(:)]',2*size(A2,1), []);

if isempty(shift00)
    figure(11);
    subplot(2,1,1); imagesc(A2);
    subplot(2,1,2); imagesc(B2);

    figure(12);
    subplot(2,1,1);imagesc(im0);title('original');% subplot(2,1,1); 
    subplot(2,1,2); imagesc(im);title(sprintf('corrected %.1f', shiftv));
end

end
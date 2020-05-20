%%

imgLennaSpectrum = fftshift(fft2(phase_ppln));
mxLennaSize = size(phase_ppln);
	intRows = mxLennaSize(1);
	intCols = mxLennaSize(2);

% Part 2 - Lowpass Filter Design
    % a)
    LowpassFilter = fspecial('gaussian', [11 11],100);
    
    
    
    % b)  fft of lowpass filter
    LowpassFilterSpectrum = fftshift(fft2(LowpassFilter, intRows, intCols));
    
     LowpassFilterSpectrum = phase_ppln*0; size1=200;
     LowpassFilterSpectrum(round(size(LowpassFilterSpectrum, 1)/2-size1/2):round(size(LowpassFilterSpectrum, 1)/2+size1/2),...
 round(size(LowpassFilterSpectrum, 2)/2-size1/2):round(size(LowpassFilterSpectrum, 2)/2+size1/2))=1;
    
    % c)  plots of lowpass filter function in frequency domain
    figure('Name', 'Lowpass Filter', 'NumberTitle', 'off');
    colormap('gray');
    subplot(1,2,1);
    imagesc(log10(1+abs(LowpassFilterSpectrum)));
    title('Log Magnitude');
    subplot(1,2,2);
    imagesc(angle(LowpassFilterSpectrum));
    title('Phase');
    
    % d)  Filtering
    imgLennaFiltered = LowpassFilterSpectrum .* imgLennaSpectrum;
    figure('Name', 'Lowpass Filtered Image', 'NumberTitle', 'off');
    colormap('gray');
    subplot(2,2,1);
    imagesc(log10(1+abs(imgLennaFiltered)));
    title('Log Magnitude');
    subplot(2,2,2);
    imagesc(angle(imgLennaFiltered));
    title('Phase');

    % e) Inverse FFT
    subplot(2,2,3); 
    imgLennaLowpassFiltered = abs(ifft2(imgLennaFiltered));
    imgLennaLowpassFiltered = circshift(imgLennaLowpassFiltered,[-1.*floor(length(LowpassFilter)/2) -1.*floor(length(LowpassFilter)/2)]);
    imagesc(imgLennaLowpassFiltered);
    title('Inverse FFT of Lowpass Filtered Image'); colormap jet; colorbar;
    
    subplot(2,2,4); 
    
    imagesc(phase_ppln);
    title('Original'); colormap jet;colorbar;
    
    %-------------------------------------------------------------------------------
%% Part 3 - High Pass Filter Design
    % a) 
    HighpassFilter = fspecial('laplacian'); % Using default alpha value.
    
    % b) fft of highpass filter
    HighpassFilterSpectrum = fftshift(fft2(HighpassFilter, intRows, intCols));
    
     HighpassFilterSpectrum = phase_ppln.^0; size1=10;
     HighpassFilterSpectrum(round(size(HighpassFilterSpectrum, 1)/2-size1/2):round(size(HighpassFilterSpectrum, 1)/2+size1/2),...
 round(size(HighpassFilterSpectrum, 2)/2-size1/2):round(size(HighpassFilterSpectrum, 2)/2+size1/2))=0;
    
    % c) plots of highpass filter function in frequency domain
    figure('Name', 'Highpass Filter', 'NumberTitle', 'off');
    colormap(gray);
    subplot(1,2,1);
    imagesc(log10(1+abs(HighpassFilterSpectrum)));
    title('Log Magnitude');
    subplot(1,2,2);
    imagesc(angle(HighpassFilterSpectrum));
    title('Phase');
    
    % d) Filtering
    imgLennaFiltered = HighpassFilterSpectrum .* imgLennaSpectrum;
    figure('Name', 'Highpass Filtered Image', 'NumberTitle', 'off');
    colormap(gray);
    subplot(2,2,1);
    imagesc(log10(1+abs(imgLennaFiltered)));
    title('Log Magnitude');
    subplot(2,2,2);
    imagesc(angle(imgLennaFiltered));
    title('Phase');
    
    % e) Inverse FFT
    subplot(2,2,3);
    imgLennaHighpassFiltered = abs(ifft2(imgLennaFiltered));
    imgLennaHighpassFiltered = circshift(imgLennaHighpassFiltered,[-1.*floor(length(HighpassFilter)/2) -1.*floor(length(HighpassFilter)/2)]);
    imagesc(imgLennaHighpassFiltered);
    title('Inverse FFT of Highpass Filtered Image'); colormap jet; colorbar;
    
    subplot(2,2,4); 
    
    imagesc(phase_ppln);
    title('Original'); colormap jet;colorbar;

%% diff

coeff0 = pi/10;
aa=diff(phase_ppln); 
aa(aa>coeff0) = coeff0*1.05;aa(aa<-coeff0) = -coeff0*1.05;
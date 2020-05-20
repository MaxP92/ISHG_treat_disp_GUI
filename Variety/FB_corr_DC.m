% ESSAI CORRECTION DU F/B en presence de backscattering

% IDEE : On a B_reel = B + epsilon*F
% donc (B/F)_reel = (B/F) + epsilon
% Il faut alors filtrer la partie constante de B/F 

cc = corr_PMT/r_center_calib;

b_f = bwd_crop_avg_neighb./(fwd_crop_avg_neighb*cc);

Vfft = fft2(b_f); %tf

Vfft2 = fftshift(Vfft); 

figure;imagesc(abs(Vfft2));
Vfft2(floor(size(Vfft, 1)/2)+1, floor(size(Vfft, 2)/2)+1)=0;

figure;imagesc(abs(Vfft2));

Vreloaded=abs(ifft2(Vfft2));

figure; subplot(1,2,1); imagesc(b_f); colormap gray; 
subplot(1,2,2); imagesc(Vreloaded); colormap gray;

%% max of the hist of B/F map

figure; histogram(abs(Vreloaded));
[N,edges] = histcounts(abs(Vreloaded), 1000);

[~, II] = max(N);

ratio_FB_max_hist = 1/edges(II);

%% max of the hist of F/B map

vv=1./Vreloaded; vv(vv>50)=0;
figure; histogram(vv), round(max(max(vv )));
[N2,edges2] = histcounts( vv, 1000);

[~, II2] = max(N2);
ratio_FB_max_hist2 = edges2(II2);

%% max of the hist of F/B map without correction

figure; histogram((fwd_crop_avg_neighb*cc)./bwd_crop_avg_neighb, 1000);
[N_classic, edges_classic] = histcounts((fwd_crop_avg_neighb*cc)./bwd_crop_avg_neighb, 1000);

[~, II_classic] = max(N_classic);
ratio_FB_max_hist_classic = edges_classic(II_classic);

% F/B_corr = 2.1 against F/B_classic = 0.987
% with 11.01.2014, with 1" water column

% F/B_corr = 1.4 against F/B_classic = 0.72
% with 11.01.2014, with no water column

% --> difference  = la meme
function [phase_test, img, hfun2, AA] = map_spec_region_ISHG( x, y, phase_test, amp, xTitle_dflt, yTitle_dflt, Titre1, phi_mat_default, pb, ph,...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig1, hfun2, screensize, fact, left_offset_fig, top_offset_fig, cmap_brbgsat  )
% [phase_test_hist, second_hist, phase_test, hfun2] = map_spec_region_ISHG( x, y, phase, amp, img_shg, shg, xTitle_dflt, yTitle_dflt, Titre1, phi_mat_default, pb, ph,...
%  axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, undocked_fig1, hfun2, screensize, fact, left_offset_fig, top_offset_fig, cmap_brbgsat)
% 
% 2015 - 10 edited by Maxime PINSARD
% 
% To calculate and plot the histogram 2D and 3D for a specific range of
% intensity

    % nbinsy2 = ceil((ph-pb)/int_y);
   
%     if shg==1 % SHG mode chosen
%         
%         AA = find(img_shg(1:size(phase_test,1), 1:size(phase_test,2))>ph | img_shg(1:size(phase_test,1), 1:size(phase_test,2))<pb);
%         img_shg(AA) = max(max(img_shg))+1;
%         img =img_shg; 
% %         shg_test_hist = reshape(shg_test,1,size(shg_test,1)*size(shg_test,2));
% %         shg_test_hist(shg_test==max(max(img_shg))+1) = [];
% %         second_hist = shg_test_hist;
%     else % interferometric contrast chosen
        
        AA = find(amp(1:size(phase_test,1), 1:size(phase_test,2))>ph | amp(1:size(phase_test,1), 1:size(phase_test,2))<pb);
        amp(AA) = 0;
         img =amp; 
%         amp_test_hist = reshape(amp_test,1,size(amp_test,1)*size(amp_test,2));
%         amp_test_hist(amp_test==0) = [];
%         second_hist = amp_test_hist;
%     end
    
   
    phase_test(AA) = 1.05;
    
%     phase_test_hist = reshape(phase_test,1,size(phase_test,1)*size(phase_test,2));
%     phase_test_hist(phase_test_hist==1.05) = [];
    
%     
    %% Phase map uploaded
    
    Titre1 = sprintf('%s in [%.1f, %.1f]', Titre1, pb, ph);
    
    hfun2 = pl_phasemap_ISHG(undocked_fig1, screensize, fact, left_offset_fig, ...
top_offset_fig, hfun2, phase_test, x, y, ...
    xTitle_dflt, yTitle_dflt, Titre1, cmap_brbgsat, phi_mat_default, ...
axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz);
    
    set(hfun2,'CLim',[-1 1.05]);
   
    
end


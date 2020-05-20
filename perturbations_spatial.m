function im3D_out = perturbations_spatial(im3D, eta, modX, modY, freq_sinY, freq_sinX)
% % im3D_out = perturbations_spatial(im3D, eta, modX, modY, freq_sinY, freq_sinX)
% % 
% % purposely impose some perturbations (sin ...) modulation, additionnal
% (eta = 0) or multiplicative (eta = 1), or in-between.

% % eta = 1 : multiplicative 100%
% % norm_add_then_mult = 1; 
% % arr_norm3 = arr_norm2; 
% % NO min subtract x2
% % works
% % 
% % eta = 0 : additive 100%
% % norm_add_then_mult = 1; 
% % arr_norm3 = mean(im3D_out(:,: , 1:round(ss/2)), 3); 
% % min subtract x1 
% % works
% % or no norm

im3D_out = zeros(size(im3D));
[XX, YY] =meshgrid(1:size(im3D,2) , 1:size(im3D,1));

mod_sin = cell(1,2); mod_flag = [modX, modY]; freq=[freq_sinX, freq_sinY]; coord = {XX, YY};

for mm=1:2
    if mod_flag(mm)
        mod_sin{mm} = sin(coord{mm}*2*pi*freq(mm)/size(im3D,3-mm));
    else
        mod_sin{mm} = im3D_out(:,:,1)+1; % ones
    end
end

mask = im3D == 0;
im3D(mask) = nan;
for k =1:size(im3D,3)
    im3_k = im3D(:, :, k);
    if eta ~=0
        im3D_out(:, :, k) = eta*im3_k.*mod_sin{1}.*mod_sin{2};
    end
    if eta ~=1
        im3D_out(:, :, k) = im3D_out(:, :, k) + (1-eta)*(im3_k+ mean2(im3_k(~isnan(im3_k))).*mod_sin{1}.*mod_sin{2});
    end
end
im3D_out(mask) = 0;
% mm = min(min(min(im3D_out )));
% if mm < 0
%     im3D_out = im3D_out - mm;
% end

end

%% test 

% % im3D=ones(400,400);
% % [XX, YY] =meshgrid(1:size(im3D,2) , 1:size(im3D,1));
% % modX = 1; modY = 1; % % modulate along X and/or Y
% % eta = 1; % raio mult vs addition of perturb
% % freq_sinY = 10; freq_sinX = 5;
% % mod_sin = cell(1,2); mod_flag = [modX, modY]; freq=[freq_sinX, freq_sinY]; coord = {XX, YY};
% % for mm=1:2
% %     if mod_flag(mm)
% %         mod_sin{mm} = sin(coord{mm}*2*pi*freq(mm)/size(im3D,3-mm));
% %     else
% %         mod_sin{mm} = im3D_out(:,:,1)+1; % ones
% %     end
% % end
% % im3D_out = im3D.*mod_sin{1}.*mod_sin{2};
% % 
% % arr_norm = im3D_out*18; % %sqrt(sum(im3D.^2, 3));
% % % % !!! if you put sqrt and square sum, an additive perturbation will be
% % % visible !!!
% % arr_norm2 = arr_norm; arr_norm2(arr_norm2 == 0) = nan; % division by 0 will make nan
% % norm = mean2(arr_norm)./arr_norm2;
% % 
% % im3D_out2 = im3D_out.*norm;
% % figure;imagesc(im3D_out2);
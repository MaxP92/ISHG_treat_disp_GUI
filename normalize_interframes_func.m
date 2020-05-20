function im3D_out = normalize_interframes_func(im3D, near_avg, norm_add_then_mult, contrast_flag)

% % mm = min(min(min(im3D )));
% % if mm < 0
% %     im3D = im3D - mm;
% % %     arr_norm3 = arr_norm3 - mm;
% % end

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

ss = size(im3D, 3);
im3D_out = im3D;
arr_norm = mean(im3D, 3); % %sqrt(sum(im3D.^2, 3));
% % !!! if you put sqrt and square sum, an additive perturbation will be
% visible !!!
arr_norm2 = arr_norm; 
if ~contrast_flag
    arr_norm2(arr_norm2 == 0) = nan; % division by 0 will make nan
end

if norm_add_then_mult % two steps
    for k =1:ss
        im3D_out(:, :, k) = (im3D(:, :, k)-arr_norm2);
    end
    arr_norm3 =mean(im3D_out(:,: , 1:round(ss/2)), 3);
else
    arr_norm3 = arr_norm2; %mean(im3D_out(:,: , 1:round(ss/2)), 3);  %arr_norm2; %

end

for k =1:ss
    im3D_out(:, :, k) = im3D_out(:, :, k).*mean2(arr_norm(~isnan(arr_norm)))./arr_norm3;
end

% % im3D_out (isnan(im3D_out)) = 0; % division by 0 will make nan

near_avg = 0;
if ~near_avg % % otherwise some values might be very offsetted
    if ~contrast_flag
        im = im3D_out(im3D_out~= 0);
    else
        im = im3D_out(:);
    end
    nb = round(length(im)/10);
    [N,edges_hist] =histcounts(im, nb);
    coef = 0.99;
    x1=find(cumsum(N)>=length(im)*(1-coef)); if ~isempty(x1); x1 = x1(1); else; x1 = 1; end
    x2=find(cumsum(N)>=length(im)*coef); if ~isempty(x2); x2 = x2(1); else; x2 = length(N); end
    if x1~= x2
        im3D_out(im3D_out < edges_hist(x1)) = edges_hist(x1);  im3D_out(im3D_out > edges_hist(x2)) = edges_hist(x2);
    end
end

mm = min(min(min(im3D_out )));
if mm < 0
    im3D_out = im3D_out - mm;
%     arr_norm3 = arr_norm3 - mm;
end

end
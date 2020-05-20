function [new_img, num_images] = avg_whole_img_ISHG (avg_img, num_images, A1, select_only_few_frame_per_ps, fname1, ...
    axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, hh_plot)
%  [new_img, num_images] = avg_whole_img_ISHG (avg_img, num_images,  A1)
%
% 2017.3 , edited by Maxime PINSARD
%
% OLD version is % 2015.10.23 C.A COUTURE
%
%  Allow to do the average on the whole image
% Ex : 12 images, moyenné sur 3
% On suppose que les images 1 à avg_img sont identiques, puis que la phase change
% Je veux 1-4-7-10 donc 1:avg:12-(avg_img-1)
% Puis pour chaque i, ajouter 2 puis 3 et moyenner (2:avg_img)

%% avg
if avg_img > 1
    
    if select_only_few_frame_per_ps
        sel = menu('Save avg stack after ?', 'Y', 'N');
        array_frames = erase_img_avg_stck(A1, axes_font_size, xaxis_sz, yaxis_sz, title_sz, clrbr_tl_sz, hh_plot, avg_img, size(A1,3)/avg_img);
    else
        array_frames=[];
        sel= 0;
    end
   
    % Compteur
    new_cpt = 1;

    % MP : probably improvable in speed, with a method to read fast the
    % images at the beginning

    % MP : no method to average 3D matrices without loops, because 3D
    % multiplication is not possible
    new_img = zeros(size(A1,1), size(A1,2), size(A1,3)/avg_img);
    for i=1:avg_img:size(A1,3)-(avg_img-1)
        vect_range = i:i+avg_img-1;
        if length(array_frames(:)) ~= sum(~isnan(array_frames(:))) % some frame deleted
            vect_range2 = vect_range(~isnan(array_frames(:, (i-1)/avg_img + 1)));
        else
            vect_range2 = vect_range;
        end  
        
        new_img(:, :, new_cpt) = mean(A1(:, :, vect_range2), 3);
        new_cpt = new_cpt +1;
    end

    % OLD METHOD

    % %     for i = 1:avg_img:num_images-(avg_img-1)
    % %         new_img(:, :, i) = mean(A1(:, :, i:i+avg_img-1), 3);
    % % %         A1 = imread(fname, i);
    % %         tmp_tot = A1(:,:,i);
    % %         for j = 2:avg_img
    % % %             tmp = imread(fname, i+(j-1));
    % %             tmp = A1(:,:, i+(j-1));
    % %
    % %             tmp_tot = tmp_tot+tmp;
    % %             clear tmp
    % %         end
    % %         tmp_tot = tmp_tot./avg_img;
    % %         new_img(:,:,new_cpt) = tmp_tot;
    % %         new_cpt = new_cpt+1;
    % %         clear tmp_tot
    % %     end

    % Nouveau nombre d'images
    num_images = num_images/avg_img;
    
    if sel == 1
        outputFileName = sprintf('%s_%s_%d.tif', fname1(1:end-4), 'AVG', avg_img ); 

        for k = 1:num_images
            if k == 1
                w_mode = 'overwrite';
            else
                w_mode = 'append';
            end
            imwrite(uint16(new_img(:, :, k)), outputFileName, 'WriteMode', w_mode,  'Compression','none');
        end
    end
end


end
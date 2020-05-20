function [ A1, A2 ] = filtering_img_stef ( A1, A2 )
% filtering_img_stef  : [ A1, A2 ] = filtering_img_stef ( A1, A2 )
%  
% 2015.10.23 : Stephane BANCELIN, edited by Maxime PINSARD
% 
% To apply a line filter (if looking at lineic structures for instance)

            % directional sampling in 2-D
            Ndir=36; % number of directional samples for the 2-D oriented laplacian filter, must be even number
            seuil=0.8;
            % laplacian filter parameters
            th3 = 15;    % longueur d'un segment rectiligne
            im_temp(:,:,1)=A1;
            im_temp(:,:,2)=A2;
            
            for ii = 1:2
                
                im=im_temp(:,:,ii);
                % mask=single(ones(size(im))); % useless ??
                %     im = imfilter(im, 1/9*ones(3,3)).*imdilate(mask, strel('disk', th3));
                im = imfilter(im, 1/9*ones(2,2));
                %                     im = medfilt2(im, [2 2]);
                i_max_d = single(zeros(size(im)));
                i_dir = single(zeros(size(im)));
                i_d = single(zeros([size(im) Ndir]));
                d_sum = zeros(1,Ndir);
                
                for j = 1:Ndir ;% for loop ok
                    % Angle de l'élément structurant
                    th = (j-1)/Ndir*180;
                    i_d(:,:,j) = imopen(im, strel('line', th3, th));
                    % Orientation local
                    i_dir (i_d(:,:,j)>i_max_d) = th;
                    % Valeur de la fibrilles réhaussée
                    i_max_d = max(i_d(:,:,j), i_max_d);
                end
                
                sup_open = max(i_d, [], 3);
                mean_open = mean(i_d, 3);
                
                % orientation component
                for dd=1:Ndir % for loop compulsory because looking at 3rd dimension
                    temp(:,:) = max(0, i_d(:,:,dd)-mean_open); % for each element, take the max of (0, the element)
                    % i_max_d_enh(:,:,dd)=temp(:,:);
                    d_sum(dd) = sum(sum(temp(:,:)));
                end
                
                i_max_d = sup_open-seuil*mean_open;
                
                % i_max_d=medfilt2(i_max_d,[1 1]);
                
                % figure(20); imagesc(im); axis image; colorbar; colormap gray
                % figure(21); imagesc(i_max_d); axis image; colorbar;colormap gray
                
                % sauvergarde de l'image filtrée
                % outdir = [dir '\img_fil']; [~,~,~] = mkdir(outdir);
                % imwrite(uint8(i_max_d/max(max(i_max_d))*255), fullfile(outdir, ['2D_out_filt.' num2str(i) '.tif']));
                im_temp(:,:,ii)=i_max_d;
            end
            
            A1=single(im_temp(:,:,1));
            A2=single(im_temp(:,:,2));

                %         else
        %             A1 = squeeze(new_img(i-1,:,:));
        %             A2 = squeeze(new_img(i,:,:));
        %         end
        %
        %         B1 = medfilt2(A1, [2 2]);
        %         B1=imfilter(A1,1/9*ones(2,2));
        %         A1 = B1;
        
        %         B2 = medfilt2(A2, [2 2]);
        % %         B2=imfilter(A2,1/9*ones(2,2));
        %         A2 = B2;

end


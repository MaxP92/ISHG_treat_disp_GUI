function [ contr, img_shg ] = exclude_area_ISHG( contr, obj_test, crop, shg, img_3D, rect, img_shg, nbr_objects, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig  )
% [ contr, img_shg ] = exclude_area_ISHG( contr, obj_test, crop, shg, img_3D, rect, img_shg, nbr_objects, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig  )
%
% 2015.10 edited by Maxime PINSARD
%
% To exclude (selection = manually) zones from your ISHG images

% Les counts à 0 sont remplacés par 1 pour ne pas être exclus plus tard

% OLD VERSION
% % for k = 1:size(contr,3)
% %     for i = 1:size(contr,1)
% %         for j = 1:size(contr,2)
% %             if contr(i,j,k)==0
% %                 contr(i,j,k) = 1;
% %             end
% %         end
% %     end
% % end

contr(contr == 0) = 1;

% S'il y a des cellules à exclure de l'analyse

if obj_test==1
    A1 =img_3D(: , :, 1);
    
    if undocked_fig
        
        figure('outerposition',...
        [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
        screensize(3)*fact screensize(4)*fact]);
        % [left bottom width height]
        h = axes;
    else
        if h ~= 0
            axes(h);
        end
    end
    
    if crop==2
        imagesc(imcrop(A1,rect));
    else
        imagesc(A1);
    end
    axis(h, 'image');
    colormap(h, gray);
    xlabel(h, 'X (pixels)'); 
    ylabel(h, 'Y (pixels)');
    title(h, sprintf('Circle your %d areas to exclude with freeHand', nbr_objects))
    % On laisse l'utilisateur sélectionner les régions à exclure une après l'autre
    cells = cell(1,nbr_objects);
    
    for i = 1:nbr_objects
        hFH = imfreehand(); % selection of excluded zones
        
        binaryImage = hFH.createMask();
        cells{i} = binaryImage;
        if shg==1
            img_shg(cells{i}) = 0;
        end
        
        for k = 1:size(contr,3)
            
            tmp2 = contr(:,:,k);
            
            tmp2(cells{i}) = 0;
            contr(:,:,k) = tmp2;
            
        end
    end
    
    clear tmp2
     
    % OLD
% %     for i = 1:nbr_objects
% %         hFH = imfreehand(); % selection of excluded zones
% %         
% %         binaryImage = hFH.createMask();
% %         cells{i} = binaryImage;
% %     end
% %     
% %     % aire_obj = 0;
% %     % On met la valeur des images à zéro là où il y a des cellules ...
% %     for i = 1:size(contr,3)
% %         tmp = contr(:,:,i);
% %         for j = 1:nbr_objects
% %             tmp(cells{j}) = 0;
% %             if shg==1
% %                 img_shg(cells{j}) = 0;
% %             end
% %             % if i==1
% %             % aire_obj = aire_obj+size(find(cells{j}==1),1);
% %             % end
% %         end
% %         contr(:,:,i) = tmp;
% %     end
% %     
% %     clear tmp

end

end


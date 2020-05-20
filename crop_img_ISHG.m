function [ x, y, cropped_img, rect ] = crop_img_ISHG(img_3D, method_crop, resx, resy, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig )
% [ x, y, cropped_img, rect ] = crop_img_ISHG(img_3D, method_crop, resx, resy, undocked_fig, h, screensize, fact, left_offset_fig, top_offset_fig )
%
% 2015.10 edited by Maxime PINSARD
% edited 2016.6.8
% for frame by frame, see V13 or earlier
%
%   Crop or not the image according to the area defined by user

usr = 0;

while usr~=1 % Variable pour le contentement de l'usager
    
    if undocked_fig
        
        h1 = figure('outerposition',...
            [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
            screensize(3)*fact screensize(4)*fact]);
        % [left bottom width height]
        h = axes; % create axes
        
    end
    
    %method_crop = 2; % rectangle
    axes(h);
    % %             h = figure;
    imagesc(img_3D(:,:,1))
    axis(h, 'image')
    colormap(h, gray)
    
    switch method_crop
        
        case 1 % CROP BY FREE-HAND
            
            %                draw a line over the zone you want to select
            hFH = imfreehand();
            binaryImage = hFH.createMask();
            
            % Si on veut vraiment juste la région entourée -
            % Il y a probablement d'autres lignes à changer dans ce cas ...
            %tmp1(~binaryImage)=0;
            %tmp2(~binaryImage)=0;
            title(h, 'Select an ROI with free-hand')
            structBoundaries = bwboundaries(binaryImage);
            xy=structBoundaries{1}; % Get n by 2 array of x,y coordinates.
            xs = xy(:, 2); % Columns.
            ys = xy(:, 1); % Rows
            
            xmin = min(xs);
            ymin = min(ys);
            width = max(xs)-min(xs);
            height = max(ys)-min(ys);
            rect = [xmin ymin width height];
            
        case 2
            title(h, 'Select an ROI with the mouse (confirm = double-click)')
%             rect = getrect(h); % [xmin ymin width height];
            [~, rect] = imcrop(h); % [xmin ymin width height]; % 1st arg is the croped image
            
% %             disp('ah');
            
    end
    % [cropped_img1 rect] = imcrop(tmp1); % crop by
    % rectangle instead of free hand
    
    %cropped_img = imcrop(img_3D,rect);
    
    xmin = rect(1); ymin = rect(2); width = rect(3); height = rect(4);
   
    cropped_img = img_3D(max(1, round(ymin)): min(size(img_3D, 1), round(ymin + height - 1)), max(1, round(xmin)): min(size(img_3D, 2), round(xmin + width - 1)), :);
 
    if undocked_fig
        
        figure(h1);
        % [left bottom width height]
        %                 h = axes;
    end
    
    imagesc(cropped_img(:,:,1))
    colormap(gray)
    title('ROI selected')
    
    x = linspace(rect(1)*resx,(rect(1)+rect(3))*resx,size(cropped_img,2));
    y = linspace(rect(2)*resy,(rect(2)+rect(4))*resy,size(cropped_img,1));
    
    usr = menu('Région OK ?','Oui','Non');
    
end

if undocked_fig
    try
        close(h1);
    end
end

end


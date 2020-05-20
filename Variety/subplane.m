function   [IM, x, y]=subplane(Z)
%SUBtract selected PLANE, or remove slope (detilt) from data matrix
% 
%Call:
%           IM=subplane(Z)
%Input:
%           Z = (double) data matrix 
%Output:
%           IM = Z - plane     (the SUBtracted PLANE is defined by manually selected points)
%
%Vassili Pastushenko	March	2005
%==============================
%figure;
imagesc(Z); axis image; colormap(hsv)
% shg
set(gca,'fontsize',15)
%select at least three (x,y) points which define a plane in 3D
title('Click at least three coplanar points')
[x,y]=getpts(gca);

if numel(x)<3,
    title('More points please');
    [xx,yy]=getpts(gca);
    x=[x;xx];y=[y;yy];
end
    
x=round(x);
y=round(y);
PONT=numel(x);
ZIN=ones(PONT,1);
M=[x y ZIN];
for i=1:PONT
ZIN(i)=Z(y(i),x(i));
end

V=M\ZIN;
[VY,VX]=size(Z);
[X,Y]=meshgrid(1:VX,1:VY);
BAS=V(1)*X+V(2)*Y+V(3);
IM=Z-BAS;

IM(IM>1) = IM(IM>1) -2 ; IM(IM<-1) = IM(IM<-1) +2 ;

imagesc(IM); axis image; colormap(hsv)
title('Detilted data')


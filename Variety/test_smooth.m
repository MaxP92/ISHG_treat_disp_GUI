clear all
close all

phi = 0.3*pi;
a = 5;

t=0:15:525;
t=t';

x=a*cos(t/180*pi+phi);
% x=ones(size(t));
r = a/2*((rand(size(x))*2)-1);
%r=0;
% plot(t,x)

x=x+r;
% hold on; plot(t,x,'r')
plot(t,x)

amp = (max(x)-min(x))/2;
mid = (max(x)+min(x))/2;

X=ones(size(t,1),3);
X(:,2)=mid+amp*cos(t/180*pi)';
X(:,3)=mid+amp*sin(t/180*pi)';

% Résoud X*beta = x où beta donne 2 chiffres : 1 pour le cos, 1 pour le sin
beta=X\x;
amp_opt = sqrt(beta(2)^2+beta(3)^2); % Chiffre entre 0 et 1
% amp
new_amp = amp_opt*amp
phi_opt = atan2(-beta(3),beta(2))/pi

%err1=sqrt(sum((x-(mid+amp*cos(t/180*pi+phi_opt*pi))).^2)/size(t,1));
err2=sqrt(sum((x-(mid+new_amp*cos(t/180*pi+phi_opt*pi))).^2)/size(t,1));

% if err1<err2
%     err1
%     plot(t,x)
%     hold on; plot(t,mid+amp*cos(t/180*pi+phi_opt*pi),'r')
% elseif err1==err2
%     err1
%     err2
% else
%     err2
%     plot(t,x)
     hold on; plot(t,mid+new_amp*cos(t/180*pi+phi_opt*pi),'r')
% end

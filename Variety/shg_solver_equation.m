syms Ah(x,y,z) x y z kh wh c chi A0 k1 kh zR delta_k
zz = -100:0.01:150; %um
zR = 4; % um
L = 300; %um
chi0 = 1;
LL = 0:0.1:zR*100; %um
z= LL/2;
% zzR = 0:0.01:L*10; %um

q = @(z) z-1i*zR;

A1 = @(x, y, z) A0*(2/pi)^0.5*zR/(1i*q(z))*exp(1i*k1*(x^2+y^2)/(2*q(z)));

eqn = diff(Ah, x, 2)+ diff(Ah, y, 2) + 2*1i*kh*diff(Ah, z) == -wh*chi/c^2*Ah^3*exp(1i*(3*k1-kh)*z);

AhSol(x,y,z) = solve(eqn); % dsolve

f = @(z) 1i*chi*(1/q(z)-1/q(L/2));

A3 = @(x, y, z) f(z)*zR/(1i*q(z))*exp(1i*kh*(x^2+y^2)/(2*q(z)));

b= simplify(diff(A3(x, y, z), x, 2)+ diff(A3(x, y, z), y, 2) + 2*1i*kh*diff(A3(x, y, z), z));

a= simplify(-wh*chi/c^2*A3(x, y, z)^3*exp(1i*(3*k1-kh)*z));

a==b

syms chi(z) I$2 chi0
assume(L > 0); assume(chi0 > 0);
chi(z) = piecewise(z<-L/2, 0, z>L/2, 0, chi0);
% assume(z>=L/2)
int(chi(z)/(z-1i*zR)^2,-Inf, z, 'IgnoreSpecialCases', true);

int_THG = int(chi(z)/(z-1i*zR)^3,-Inf, z, 'IgnoreSpecialCases', true);

int_SHG = int(chi(z)/(z-1i*zR)^1,-Inf, z, 'IgnoreSpecialCases', true);

pretty(int_SHG)


%% function of Z, L >> zR

plot(zz, abs(-chi0*(log(- L - zR*2i) - log(2*zz - zR*2i)))); % module

ASHG1 = (-chi0*(log(abs((- L - zR*2i)./ (2*zz - zR*2i))))).^2;
ATHG1 = abs(-chi0*(1./(- L - zR*2i)-1./ (2*zz - zR*2i))).^(2);
As = ASHG1/max(ASHG1); At = ATHG1/max(ATHG1);

figure;
plot(zz, As); % module
hold on; plot(zz, At); legend('SHG', 'THG')


%% function of Z, L << zR
L = zR/100; %um
figure;
% plot(zz, -real(-chi0*(log(- L - zR*2i) - log(2*zz - zR*2i))));
plot(zz, abs(-chi0*(log(- L - zR*2i) - log(2*zz - zR*2i)))); hold on; 
for ii=[1/100, 1/10, 1/2, 1, 2, 3,  5, 8, 10, 100]
    L = zR*ii; %um
%     plot(zz, abs(-chi0*(log(- L - zR*2i) - log(2*zz - zR*2i)))); hold on;
    plot(zz, abs(-chi0*(log(abs((- L - zR*2i )./(2*zz - zR*2i))))));hold on;
end

plot(LL, abs(-chi0*(log(- LL - zR*2i) - log(2*z - zR*2i))));

plot(LL, abs(-chi0*(log(- LL ) - log(2*z - zR*2i))));

plot(LL, -chi0*(log(abs((- LL )./(2*z - zR*2i)))));

figure;
shgLIM = (chi0*abs((log((- LL + zR*2i)./(LL + zR*2i))))).^2; GSHG = shgLIM/max(shgLIM);
thgLIM =(chi0*(LL)./((LL/2).^2 + zR.^2)).^2; GTHG = thgLIM/max(thgLIM);
plot(LL/zR, GSHG); hold on;
plot(LL/zR, GTHG);
legend('SHG', 'THG')

%% phase mismatch
syms z
assume(delta_k>0); assume(z>0);
% assume(L > 0); assume(chi0 > 0); assume(zR > 0);
int_SHG_mismatch = int(chi0*exp(1i*delta_k*z)/(z-1i*zR),z, -Inf, Inf, 'IgnoreSpecialCases', true);

int_THG_mismatch = int(chi(z)*exp(1i*delta_k*z)/(z-1i*zR).^2, z,-Inf, Inf, 'IgnoreSpecialCases', true);

zR = 4;% um
L = 500;% um
delta_n = -0.2;
lambda = 0.8; % um
delta_k = 4*pi/lambda*delta_n; % um-1
deltaZ = 0.1;

for ii=1:round((2*500+L)/deltaZ)
    z0(ii) = -500+ii*deltaZ;
    z = z0(ii):deltaZ:z0(ii)+L;
    A(ii) = sum(abs(exp(1i*delta_k*z)./(1+1i*z/zR))*deltaZ).^2;
    ATHG(ii) = sum(abs(exp(1i*delta_k*z)./(1+1i*z/zR).^2)*deltaZ).^2;
end
figure;plot(z0,A); hold on; plot(z0,ATHG);

z = zz; %-50:deltaZ:-50+L;
 ASHG = (abs(exp(1i*delta_k*z)./(1+1i*z/zR))*deltaZ).^2;
  
 ATHG = (abs(exp(1i*delta_k*z)./(1+1i*z/zR).^2)*deltaZ).^2;
 As2=ASHG/max(ASHG); At2 = ATHG/max(ATHG);
 figure; 
 plot(z,   As2); hold on; plot(z,  At2 ); 
 legend('SHG', 'THG', 'SHG phase-mismatched', 'THG phase-mismatched')


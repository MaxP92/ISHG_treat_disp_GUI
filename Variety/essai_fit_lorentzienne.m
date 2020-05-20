% essai fit lorentzienne

%% obtain plot from hist

clear all; clc; close all;
x = randn(10000,1);
nbins = 40;
% histogram(x,nbins); hold on;
[N,edges] = histcounts(x,nbins );
xvector = edges(1:end-1) + mean(-edges(1:end-1)+edges(2:end))/2;

xvector_2 = linspace(min(xvector)-1, max(xvector)+1, 200000);

% lorentzian

ff = 'a1/(1 + ((x-b1)/c1)^2)';
fit_f = fittype(ff,'independent',{'x'},...
    'coefficients',{'a1','b1','c1'}); % works
[curve2, ~] = fit (xvector', N', fit_f, 'Lower' ,[0; 0; 0]);
yfit_lor = (curve2.a1)./(1 + ((xvector_2-curve2.b1)./(curve2.c1)).^2) ;

% gaussian
ff2 = 'a1*exp(-((x-b1)/(c1))^2)';

fit_f2 = fittype(ff2,'independent',{'x'},...
    'coefficients',{'a1','b1', 'c1'});
% % opt_fit1 = fitoptions(fit_f2);
% % % opt_fit1.Method= 'NearestInterpolant';
% % % opt_fit1.StartPoint = [758; 0; 1.42]; % critical !!
% % opt_fit1.Lower = [0; 0; 0]; % critical !!
[curve22, ~] = fit (xvector', N', fit_f2, 'Lower' ,[0; 0; 0]);
[curve222, ~] = fit (xvector', N', 'gauss1');

yfit_gauss = curve22.a1.*exp(-((xvector_2-curve22.b1)/(curve22.c1)).^2);

yfit222 = curve222.a1.*exp(-((xvector_2-curve222.b1)/(curve222.c1)).^2);

% product of lorentzian/gaussian

% manual

y_prod = yfit_lor.*yfit_gauss/max(yfit_lor.*yfit_gauss)*max(max(yfit_gauss), max(yfit_lor));

% fit prod
ff_prod = 'a1/(1 + ((x-b1)/c1)^2)*a2*exp(-((x-b2)/(c2)).^2)';
fit_prod = fittype(ff_prod,'independent',{'x'},...
    'coefficients',{'a1','b1', 'c1','a2','b2', 'c2'});
[curve_prod, ~] = fit (xvector', N', fit_prod, 'Lower' ,[0; 0; 0; 0; 0; 0]);
y_prod_fit = curve_prod.a1./(1 + ((xvector_2-curve_prod.b1)./curve_prod.c1).^2).*curve_prod.a2.*exp(-((xvector_2-curve_prod.b2)./(curve_prod.c2)).^2);

%% ------ Voigt -------

voigt_met = 3; % 3 works

switch voigt_met
    
    case 1
        % Voigt NOT WORKING
        % syms x t c0 b0 a1 b1 c1
        fun_v = @(t,c0,b0,a1,b1,c1,x) exp(-c0.*(x-t-b0).^2).*a1./(1 + ((t-b1)./c1).^2);
        ff_v = @(c0,b0,a1,b1,c1,x) integral(@(t)fun_v(t,c0,b0,a1,b1,c1,x),-Inf, Inf);
        fit_f = fittype(ff_v, 'coefficients',{'c0','b0','a1','b1','c1'}, 'independent',{'x'}); %
        [curve_voigt, ~] = fit (xvector', N', fit_f, 'Lower' ,[0; 0; 0; 0; 0]);
        yfit_lor = (curve_voigt.a1)./(1 + ((xvector_2-curve_voigt.b1)./(curve_voigt.c1)).^2) ;
    case 2
        % Voigt FFT NOT WORKING
        syms x y a0 c0 b0 a1 b1 c1
        yfit_voigt3 = ifourier((fourier(a0/(1 + ((x-b0)/c0)^2), x, y).*fourier(a1*exp(-((x-b1)/(c1))^2), x, y)), y, x);
        yfit_voigt3_func = @(a0,c0,b0,a1,b1,c1, x) yfit_voigt3(a0,c0,b0,a1,b1,c1, x);
        fit_f = fittype(yfit_voigt3_func, 'coefficients',{'a0', 'c0','b0','a1','b1','c1'}, 'independent',{'x'}); %
        [curve_voigt, ~] = fit (xvector', N', fit_f, 'Lower' ,[0; 0; 0; 0; 0]);
        yfit_lor = (curve_voigt.a1)./(1 + ((xvector_2-curve_voigt.b1)./(curve_voigt.c1)).^2) ;
    case 3
        % Voigt sum
        eta = 0.834; % best value
        yfit_voigt = (1-eta)*yfit_gauss + eta*yfit_lor;
        
        
        ff_voigt = sprintf('a1.*((1-%g)/(1 + ((x-b1)/c1)^2) + %g*exp(-((x-b2)/(c2))^2))', eta, eta);
        fit_f_voigt = fittype(ff_voigt,'independent',{'x'},...
            'coefficients',{'a1','b1','c1','b2', 'c2'});
        [curve_prod_voigt, ~] = fit (xvector', N', fit_f_voigt, ...
            'Lower' ,[0; 0; 0; 0; 0], ...
            'StartPoint', []);
        vect_fit = curve_prod_voigt.a1.*((1-eta)./(1 + ((xvector_2-curve_prod_voigt.b1)./curve_prod_voigt.c1).^2) + eta.*exp(-((xvector_2-curve_prod_voigt.b2)./(curve_prod_voigt.c2)).^2));
        %             meth_fit = 'Voigt sum';
        
end

%% Von Mises

fit_bad = 1;
goodness_value = 0.9;
max_it = 5;
ff = 'a1*exp(c1*cos(x-b1))/(2*3.14*besseli(0,c1))'; % 
coeff_fit = {'a1','b1','c1'};
low_vect = [ 0;-0.1; 0]; % indeed taken
up_vect = [10000; 0.1; 8*log(2)/2^2]; %
start_vect = [];
%'Lower' ,[0; 0; 0; 0; 0], ...
fit_f = fittype(ff,'independent',{'x'},...
    'coefficients',coeff_fit);
while (fit_bad >= 1)
    [fit_struct, gof] = fit (xvector', N', fit_f, 'Lower' ,low_vect, 'Upper' ,up_vect, 'StartPoint', start_vect);
    if (gof.adjrsquare >= goodness_value || fit_bad >= max_it)  % good perf fit
        fit_bad0 = fit_bad;
        fit_bad = 0;
    else
        fit_bad = fit_bad + 1;
    end
    
end
meth_fit = 'Von Mises';
vect_fit = fit_struct.a1.*exp(fit_struct.c1.*cos(xvector_2-fit_struct.b1))./(2*3.14*besseli(0,fit_struct.c1));
% c1 = 2*sigma^2, FWHM = 2*sqrt(2*log(2))*sigma for
% GAUSS
% 1/a1 = sigma^2
% so FWHM = 2*sqrt(2*log(2)/a1)
fit_struct
b(1) = fit_struct.b1;
c(1) = 2*sqrt(2*log(2)/fit_struct.a1); % real FWHM
fprintf('gof = %.3g \n', gof.adjrsquare);

% plot
histogram(x,nbins); hold on;
plot(xvector, N, 'x'); % edges

hold on; plot(xvector_2, vect_fit, 'k', 'LineWidth', 2) % Von Mises
%%
hold on; plot(xvector_2, yfit_lor, 'k', 'LineWidth', 2) % lor

hold on;plot(xvector_2, yfit_gauss, 'g-') % gauss manu

hold on;plot(xvector_2, yfit222, 'b') % gauss auto
hold on;plot(xvector_2,y_prod  , 'm') % prod direct

hold on;plot(xvector_2,y_prod_fit  , 'r') % prod fit

hold on;plot(xvector_2,yfit_voigt  , 'y') % Voigt sum dir

hold on;plot(xvector_2, vect_fit  , 'x') % Voigt sum fit

legend('Hist rand', 'Edges', 'Lorentzian', 'Gauss manu', 'Gauss auto', 'Prod direct', 'Prod fit', 'Voigt sum dir', 'Voigt sum fit')

%% fft

yfit_voigt1 = fftshift(ifft(fft(yfit_gauss).*fft(yfit_lor)));
yfit_voigt2 = conv(yfit_gauss, yfit_lor);
% yfit_voigt3 = ifourier(fourier(yfit_gauss).*fourier(yfit_lor));

figure; hold on; plot(yfit_gauss/max(yfit_gauss));
plot(yfit_lor/max(yfit_lor)); plot(yfit_voigt1/max(yfit_voigt1));
plot(yfit_voigt2/max(yfit_voigt2));
% plot(yfit_voigt3/max(yfit_voigt3))
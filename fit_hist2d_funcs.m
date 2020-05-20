function func_hdl = fit_hist2d_funcs
func_hdl.fit_hist2d_ishg = @fit_hist2d_ishg;
func_hdl.norm_peaks_hist = @norm_peaks_hist;
% func_hdl.load_stack_plot_ISHG_func = @load_stack_plot_ISHG_func;
end

function [ stringb1, stringc1, stringb2, stringc2, cmap1, vect_fit_struct, fit_warning ] = fit_hist2d_ishg( hist2_xdata, hist2_ydata, h_ax, Legendehisto, offset, fit_chosen, M,fact_norm, baseline_hist)
%fit_hist2d_ishg :  fit (gaussian) for the hist in 2D
%
% edited 2016.2.11 by Maxime PINSARD
%   [ stringb1, stringc1, stringb2, stringc2, cmap1 ] = fit_hist2d_ishg( hist2_xdata, h_ax, Legendehisto, offset, fit_chosen, M )

%% ****** PARAMETERS *******
goodness_value = 0.88; % value min for R square of fit
% % M = 0.5; % for Pearson VII : M = 1 = Lorentz, M = Inf = Gauss, M < 1 = super-Lorentz
max_it = 5; % nb of try for doing fit expecting a good Rsquared, after 'max_it' give up
div_fact = 6; % the width of the window used to remove the first peak, in order to select well the 2nd one
% see 'div_fact' in the following for better understanding
diff_peak_max = 1.5; % frac of pi, maximum ratio between the 2 phase peaks accepted
% change this value if bad results with peaks two different in height (decrease it)

if baseline_hist>0
    hist2_ydata = hist2_ydata-baseline_hist;
end

%% ****** ELSE *******

warning('off','curvefit:fit:noStartPoint')
xfine = -1:0.001:1;
xmin = zeros(1,2); xmax = zeros(1,2);
width = zeros(1,2); a = zeros(1,2); b = zeros(1,2); c = zeros(1,2);
fit_bad = 1; % counter

hh_l = legend(h_ax);
str_l = get(hh_l, 'String');
nb = numel(str_l);
if (isempty(str_l) || strcmp(str_l{1}, 'contr'))  % first plot
    str_l = {Legendehisto};
    clr{1} = 'r'; clr{2} = 'g';
else
    if nb>=14
        menu( 'Too many plot !', 'OK');
        stringb1=0; stringc1=0; stringb2=0; stringc2=0;
        return;
    end
    clr{1} = [1-0.3-0.1*(nb-1)/2, 0.2*(nb-3)/2, 0.2*(nb-3)/2];
    clr{2} = [0.2*(nb-3)/2, 1-0.3-0.1*(nb-1)/2, 0.2*(nb-3)/2];
end

Nmax = 6; % % number of fit MANU

%% fit with 1 curve (auto)

if fit_chosen == 0 % 1 fit gaussian
    cmap1 = zeros(length(xfine),3);
    
    [~,II] = max(hist2_ydata);
    selection = max(II-length(hist2_xdata)/4, 1): min(II+length(hist2_xdata)/4, length(hist2_xdata));
    
    fitx = hist2_xdata(selection)';
    fity = hist2_ydata(selection)'; % indeed taken
    a = max(fity);
    b = sum(fity.*fitx)/sum(fity);
    c = sqrt(sum(fity.*(fitx-b).^2)/sum(fity));
    
    %                 if ii == 1
    %                     stringb1 = sprintf('%.2f', round(b(ii)*100)/100 - offset);
    %                     stringc1 = sprintf('%.2f', round(c(ii)*100)/100 - offset);
    %                 else
    %                     stringb2 = sprintf('%.2f', round(b(ii)*100)/100 - offset);
    %                     stringc2 = sprintf('%.2f', round(c(ii)*100)/100 - offset);
    %                 end
    
    vect_fit = a*exp(-(xfine-b).^2/(2*c^2));
    vect_fit_struct = vect_fit;
    meth_fit = 'Gaussian';
    
    if baseline_hist>0; vect_fit = vect_fit+baseline_hist; end
    
    c = c*2*sqrt(2*log(2)) ; % real FWHM
    c(2) =c;   b(2) =b;
    axes(h_ax); hold on;
    plot(h_ax, xfine, vect_fit,'Color', clr{1},'LineWidth',2.5);
    plot(h_ax, xfine, vect_fit,'Color', clr{2},'LineWidth',2.5);
    
elseif fit_chosen <= Nmax % % fit MANU !!
    cmap1 = zeros(length(xfine),3);
    for ii=1:2
        fit_bad = 1; % counter
        switch ii
            
            case 1
                
                edge = 'left';
                
            case 2
                
                edge = 'right';
                %                 vect_fit1 = vect_fit;
        end
        
        menu(['Choose region for ' edge ' fit in hist 2D'], '     OK     ');
        if ii == 1
            axes(h_ax); hold on; %#ok<LAXES>
        end
        
        rect = getrect(h_ax); % [xmin ymin width height]
        xmin(ii) = rect(1);width(ii)=rect(3);
        xmax(ii) = ((xmin(ii)+width(ii))+1)*(length(hist2_xdata)-1)/2 + 1; % convert to pos in array
        xmin(ii) = ((xmin(ii))+1)*(length(hist2_xdata)-1)/2 + 1; % convert to pos in array
        
        % ******** Fit manuel à 2 gaussiennes pour l'histo ***********
        
        % Il est nécessaire d'ajuster l'intervalle "selection" pour obtenir un bon fit
        %     fact1 = 3;
        %     offset_sel1 = 5 + fact1; % in number of elements (units)
        %     size_fit1 = round(length(hist2_xdata)/2) - 25 - fact1; % in number of elements
        %     selection = 1 + offset_sel1:size_fit1+offset_sel1; % !!!
        selection = max(1, round(xmin(ii))):min(length(hist2_xdata), round(xmax(ii)));
        
        fitx = hist2_xdata(selection)';
        fity = hist2_ydata(selection)'; % indeed taken
        
        switch fit_chosen
            
            case 1 % gaussian
                
                a(ii) = max(fity);
                b(ii) = sum(fity.*fitx)/sum(fity);
                c(ii) = sqrt(sum(fity.*(fitx-b(ii)).^2)/sum(fity));
                fprintf('In fit2d l55 : a=%f, b=%f, c=%f \n', a(ii),b(ii),c(ii));
                %                 if ii == 1
                %                     stringb1 = sprintf('%.2f', round(b(ii)*100)/100 - offset);
                %                     stringc1 = sprintf('%.2f', round(c(ii)*100)/100 - offset);
                %                 else
                %                     stringb2 = sprintf('%.2f', round(b(ii)*100)/100 - offset);
                %                     stringc2 = sprintf('%.2f', round(c(ii)*100)/100 - offset);
                %                 end
                
                vect_fit = a(ii)*exp(-(xfine-b(ii)).^2/(2*c(ii)^2));
                meth_fit = 'Gaussian';
                fit_bad0 = 1;
                gof.adjrsquare = 1;
                c(ii) = c(ii)*2*sqrt(2*log(2)) ; % real FWHM
                intb1(ii) = 0;
                intc1(ii) = 0;
            case 2 % lorentzian
                
                %             figure; plot(fitx, fity);
                ff = 'a1/(1 + ((x-b1)/c1)^2)'; % FWHM = 2*c1
                
                coeff_fit = {'a1','b1', 'c1'};
                low_vect = [0; rect(1) + width(ii)/2-0.01 ; (width(ii)-0.01)/4]; % indeed taken
                up_vect = [max(fity)*1.2; rect(1) + width(ii); (width(ii)+0.01)/4];
                start_vect = [];
                
                meth_fit = 'Lorentzian';
                
            case 3 % prod lorentzian/gaussian
                
                ff = 'a1/(1 + ((x-b1)/c1)^2)*a2*exp(-((x-b2)/(c2)).^2)';
                coeff_fit = {'a1','b1', 'c1','a2','b2', 'c2'};
                low_vect = [0; rect(1) + width(ii)/2-0.01 ; (width(ii)-0.01)/4; 0; rect(1) + width(ii)/2-0.01 ; (width(ii)-0.01)/2];
                up_vect = [max(fity)*1.2; rect(1) + width(ii); (width(ii)+0.01)/4; max(fity)*1.2; rect(1) + width(ii); (width(ii)+0.01)/2];
                start_vect = [];
                
                meth_fit = 'Gauss*Lorentz';
                
            case 4 % voigt
                % TO KNOW : Voigt = conv of Gaussian/Lorentzian
                % But did not succeed to do conv into 'fit' operator (use of conv, fft, or int ...)
                % Voigt can be approximated by the sum of gauss and lorentz : best param is 0.834
                
                eta = 0.834; % best value
                ff = sprintf('a1.*((1-%g)/(1 + ((x-b1)/c1)^2) + %g*exp(-((x-b2)/(c2))^2))', eta, eta);
                coeff_fit = {'a1','b1','c1','b2', 'c2'};
                low_vect =[0; rect(1) + width(ii)/2-0.01; (width(ii)-0.01)/4; rect(1) + width(ii)/2-0.01; (width(ii)-0.01)/2];
                up_vect = [max(fity)*1.2; rect(1) + width(ii); (width(ii)+0.01)/4; rect(1) + width(ii); (width(ii)+0.01)/2];
                start_vect = [];
                %'Lower' ,[0; 0; 0; 0; 0], ...
                
                meth_fit = 'pseudo-Voigt';
                
            case 5 % Pearson VII
                
                ff = sprintf('a1/(1 + ((x-b1)/c1)^2)^%g', M); % FWHM = 2*c1
                coeff_fit = {'a1','b1', 'c1'};
                low_vect = [0; rect(1) + width(ii)/2-0.01 ; (width(ii)-0.01)/4]; % indeed taken
                up_vect = [max(fity)*1.2; rect(1) + width(ii); (width(ii)+0.01)/4];
                start_vect = [];
                %'Lower' ,[0; 0; 0; 0; 0], ...
                
                meth_fit = sprintf('Pearson VII (M = %.2g)', M);
                
            case 6 % Von Mises
                
                ff = 'a1*exp(c1*cos(x-b1))/(2*3.14*besseli(0,c1))'; %
                coeff_fit = {'a1','b1','c1'};
                low_vect = [0; rect(1) + width(ii)/2-0.01; 0]; % indeed taken
                up_vect = [10000; rect(1) + width(ii); 8*log(2)/0.5^2];
                start_vect = [];
                %'Lower' ,[0; 0; 0; 0; 0], ...
                
                meth_fit = 'Von Mises';
        end
        
        if fit_chosen ~= 1 % gaussian is done another way
            intervals = confint(fit_struct); % % a1 b1 c1 a2 b2 c2
            
            fit_f = fittype(ff,'independent',{'x'},...
                'coefficients',coeff_fit);
            while (fit_bad >= 1)
                [fit_struct, gof] = fit (fitx, fity, fit_f, 'Lower' ,low_vect, 'Upper' ,up_vect, 'StartPoint', start_vect);
                if (gof.adjrsquare >= goodness_value || fit_bad >= max_it)  % good perf fit
                    fit_bad0 = fit_bad;
                    fit_bad = 0;
                else
                    fit_bad = fit_bad + 1;
                end
                
            end
            %             fit_struct.b1
            switch fit_chosen
                case 2 % lorentzian
                    
                    vect_fit = (fit_struct.a1)./(1 + ((xfine-fit_struct.b1)./(fit_struct.c1)).^2);
                    b(ii) = fit_struct.b1;
                    c(ii) = 2*fit_struct.c1; % real FWHM
                    intb1(ii) = -1*(intervals(1,2)-intervals(2,2))*0.5;
                    intc1(ii) = -2*(intervals(1,3)-intervals(2,3))*0.5;
                case 3 % prod lorentzian/gaussian
                    
                    vect_fit = fit_struct.a1./(1 + ((xfine-fit_struct.b1)./fit_struct.c1).^2).*fit_struct.a2.*exp(-((xfine-fit_struct.b2)./(fit_struct.c2)).^2);
                    b(ii) = sign(fit_struct.b1)*sqrt(abs(fit_struct.b1 * fit_struct.b2));
                    c(ii) = sqrt(2*fit_struct.c1*fit_struct.c2); % real FWHM
                    intb1(ii) = -1*(intervals(1,2)-intervals(2,2))*0.5;
                    intc1(ii) = -coef*(intervals(1,3)-intervals(2,3))*0.5;
                case 4 % voigt
                    
                    vect_fit = fit_struct.a1*((1-eta)./(1 + ((xfine-fit_struct.b1)./fit_struct.c1).^2) + eta.*exp(-((xfine-fit_struct.b2)./(fit_struct.c2)).^2));
                    
                    b(ii) = (fit_struct.b2);
                    c(ii) = (1-eta)*2*fit_struct.c1 + (eta)*sqrt(log(2)*fit_struct.c2); % real FWHM
                    intb1(ii) = -1*(intervals(1,2)-intervals(2,5))*0.5;
                    intc1(ii) = -(1-eta)*2*(intervals(1,3)-intervals(2,3)) - (eta)*sqrt(log(2)*(intervals(1,6)-intervals(2,6)));
                case 5 % Pearson VII
                    
                    vect_fit = (fit_struct.a1)./(1 + ((xfine-fit_struct.b1)./(fit_struct.c1)).^2).^M;
                    b(ii) = fit_struct.b1;
                    coef=sqrt(4*abs(2^(1/M)-1));
                    c(ii) = fit_struct.c1*coef; % real FWHM
                    intb1(ii) = -1*(intervals(1,2)-intervals(2,2))*0.5;
                    intc1(ii) = -coef*(intervals(1,3)-intervals(2,3))*0.5;
                case 6 % Von Mises
                    
                    vect_fit = fit_struct.a1.*exp(fit_struct.c1.*cos(xfine-fit_struct.b1))./(2*3.14*besseli(0,fit_struct.c1));
                    % c1 = 2*sigma^2, FWHM = 2*sqrt(2*log(2))*sigma for
                    % GAUSS
                    % 1/a1 = sigma^2
                    % so FWHM = 2*sqrt(2*log(2)/a1)
                    b(ii) = fit_struct.b1;
                    c(ii) = 2*sqrt(2*log(2)/fit_struct.a1); % real FWHM
                    intb1(ii) = -1*(intervals(1,2)-intervals(2,2))*0.5;
                    intc1(ii) = -2*sqrt(2*log(2)/(intervals(1,1)-intervals(2,1))*0.5);
            end
        end
        
        if baseline_hist>0; vect_fit = vect_fit+baseline_hist; end
        
        plot(h_ax, xfine, vect_fit,'Color', clr{ii},'LineWidth',2.5);
        %     legend(Legendehisto,['Fit Gaussien - \mu = ' stringb1 '\pi ; \sigma = ' stringc1 '\pi']);
        
        %     set(h,'Interpreter','Latex','Location','Best') %taille texte légende ~16 à 20
        % ******** END fit manuel à 2 gaussiennes pour l'histo ***********
        fprintf('Number of fit attempt = %d \n R-square adj. = %.3g (->1 for good fit)\n\n', fit_bad0, gof.adjrsquare);
        if fit_bad0 == max_it
            fprintf('Tried %d times, but no good fit\n', max_it);
        end
        
        vect_fit2 = zeros(1, length(vect_fit));
        vect_fit2(1:end-round(offset/2*length(vect_fit))) = vect_fit(round(offset/2*length(vect_fit))+1:end);
        vect_fit2(end-round(offset/2*length(vect_fit))+1:end) = vect_fit(1:round(offset/2*length(vect_fit)));
        
        vect_fit_struct{ii} = vect_fit2;
        
        % %         whos cmap1
        % %         v=vect_fit_struct{1};
        % %         whos v
        cmap1(:,ii) = (vect_fit_struct{ii}'-min(vect_fit_struct{ii}))/(max(vect_fit_struct{ii})-min(vect_fit_struct{ii}));
        % %    cmap1(:,1) = (vect_fit_struct{1}'-min(vect_fit_struct{1}))/(max(vect_fit_struct{1})-min(vect_fit_struct{1}));
        
    end
    
    %% fit with 2 curves
    
else % fit 2 curves
    
    fitx = hist2_xdata';
    fity = hist2_ydata';
    % indeed use, not always set by func
    norm_imposed = 0; erase = 0;norm_vect=0;
    if fact_norm<=0; fact_norm = 1;
    else
        ch=get(h_ax, 'children');
        for k = length(ch):-1:1
            if isa(ch(k), 'matlab.graphics.chart.primitive.Histogram')
                break
            elseif k==1 % not found
                erase = 1;
            end
        end
        if erase
            cla(h_ax);
        else
            if (strcmp(ch(k).FaceColor,'auto') || sum(ch(k).FaceColor ~= [0 0 0.5000])) % dflt histogram function
                disp('no clear hist');
            else;cla(h_ax);  erase = 1;
            end
        end
        norm_imposed = 1;
    end
    %     if erase
    [norm_vect, max_increase_with_array_index, fity, max1, max2]=norm_peaks_hist(hist2_ydata, div_fact, diff_peak_max,fact_norm); % fact_norm
    
    if norm_imposed % %h_ax.XLim(1) == 0 % was cla, otherwise -1, so normalize imposed !
        if max(fity) > 2; fity=fity/max(fity);end
        %             if max(fity) < 2; fitydisp = fity*max1;
        %             else;fitydisp = fity;
        %             end
        if erase % was erased
            histogram(h_ax, 'BinEdges',[hist2_xdata, 2*hist2_xdata(end)-hist2_xdata(end-1)],'BinCounts',abs(fity)*max1);hold on; clear fitydisp;
        end
    end
    %     end
    if fit_chosen == Nmax+1 % fit 2 gaussians
        
        low_vect = [0; 0 ; 0; 0; -1 ; 0];
        up_vect = [max(fity)*1.2; 1 ; 1; max(fity)*1.2; 0 ; 1];
        start_vect = [];
        
        while (fit_bad >= 1 && fit_bad <= max_it)
            [fit_struct, gof] = fit (fitx, fity, 'gauss2', 'Lower' ,low_vect, ...
                'Upper' ,up_vect, 'StartPoint', start_vect);
            if (gof.adjrsquare >= goodness_value || fit_bad >= max_it)  % good perf fit
                fit_bad0 = fit_bad;
                fit_bad = 0;
            else
                fit_bad = fit_bad + 1;
            end
        end
        vect_fit = fit_struct.a1.*exp(-((xfine-fit_struct.b1)./fit_struct.c1).^2) + fit_struct.a2.*exp(-((xfine-fit_struct.b2)./fit_struct.c2).^2);
        
        meth_fit = 'Gauss';
        % c1^2 = 2*sigma^2, FWHM = 2*sqrt(2*log(2))*sigma
        % FWHM = 2*sqrt(log(2))*c1
        coef = 2*sqrt(log(2));
        b(1) = fit_struct.b2; b(2) = fit_struct.b1; % inverted
        c(1) = coef*fit_struct.c2;
        c(2) = coef*fit_struct.c1;
        intervals = confint(fit_struct);
        intb1 = -(intervals(1,2)-intervals(2,2))*0.5;
        intb2 = -(intervals(1,5)-intervals(2,5))*0.5;
        intc1 = -coef*(intervals(1,3)-intervals(2,3))*0.5;
        intc2 = -coef*(intervals(1,6)-intervals(2,6))*0.5;
        
    else %% 2 curves not gaussian
        
        switch fit_chosen
            
            case Nmax+2 % 2 lorentzian
                
                ff = 'a1/(1 + ((x-b1)/c1)^2) + a2/(1 + ((x-b2)/c2)^2)'; % FWHM = 2*c1
                coeff_fit = {'a1','b1','c1','a2','b2','c2'};
                low_vect = [0; 0 ; 0; 0; -1 ; 0];
                up_vect = [max(fity)*1.2; 1 ; 1; max(fity)*1.2; 0 ; 1];
                if max(fity) > 1 % more difficult
                    start_vect = [max(fity(round(length(fity)/2)+1:end))*0.8; 0.5 ; 0.5; max(fity(1:round(length(fity)/2)))*0.8; -0.5 ; 0.5];
                else; start_vect = [];
                end
                
            case Nmax+3 % 2 prod lorentzian/gaussian
                
                ff = 'a1/(1 + ((x-b1)/c1)^2)*a2*exp(-((x-b2)/(c2)).^2)*a11/(1 + ((x-b11)/c11)^2)*a22*exp(-((x-b22)/(c22)).^2)';
                coeff_fit = {'a1','b1', 'c1','a2','b2', 'c2','a11','b11', 'c11','a22','b22', 'c22'};
                low_vect = [0; 0 ; 0; 0; 0 ; 0; ...
                    0; -1 ; 0; 0; -1 ; 0];
                up_vect = [max(fity)*1.2; 1 ; 1; max(fity)*1.2; 1 ; 1; ...
                    max(fity)*1.2; 0 ; 1; max(fity)*1.2; 0 ; 1];
                start_vect = [];
                
            case Nmax+4 % 2 voigt
                % TO KNOW : Voigt = conv of Gaussian/Lorentzian
                % But did not succeed to do conv into 'fit' operator (use of conv, fft, or int ...)
                % Voigt can be approximated by the sum of gauss and lorentz : best param is 0.834
                
                eta = 0.834; % best value
                ff = sprintf('a1.*((1-%g)/(1 + ((x-b1)/c1)^2) + %g*exp(-((x-b2)/(c2))^2))+ a11.*((1-%g)/(1 + ((x-b11)/c11)^2) + %g*exp(-((x-b22)/(c22))^2))', eta, eta, eta, eta);
                coeff_fit = {'a1','b1','c1','b2', 'c2','a11','b11','c11','b22', 'c22'};
                low_vect = [0; 0 ; 0; 0 ; 0; ...
                    0; -1 ; 0; -1 ; 0];
                up_vect = [max(fity)*1.2; 1 ; 1; 1 ; 1; ...
                    max(fity)*1.2; 0 ; 1; 0 ; 1];
                start_vect = [];
                %'Lower' ,[0; 0; 0; 0; 0], ...
                
            case Nmax+5 % 2 Pearson VII
                % lorentzienne puissance
                
                ff = sprintf('a1/(1 + ((x-b1)/c1)^2)^%g + a2/(1 + ((x-b2)/c2)^2)^%g', M, M); % FWHM = 2*c1
                coeff_fit = {'a1','b1','c1','a2','b2','c2'};
                low_vect = [0; 0 ; 0; 0; -1 ; 0];
                up_vect = [max(fity)*1.2; 1 ; 1; max(fity)*1.2; 0 ; 1];
                start_vect = [];
                
            case Nmax+6 % 2 Von Mises
                
                ff = 'a1*exp(c1*cos(x-b1))/(2*3.14*besseli(0,c1)) + a2*exp(c2*cos(x-b2))/(2*3.14*besseli(0,c2))'; %
                coeff_fit = {'a1','b1','c1', 'a2','b2','c2'};
                low_vect = [0; 0; 0; 0;-1;0 ]; % indeed taken
                up_vect = [max(fity)*1.2; 1; 10*8*log(2)/0.5^2; max(fity)*1.2; 0; 10*8*log(2)/0.5^2];
                start_vect = [];
                %'Lower' ,[0; 0; 0; 0; 0], ...
                
                meth_fit = 'Von Mises';
                
        end
        
        fit_f = fittype(ff,'independent',{'x'},...
            'coefficients',coeff_fit);
        
        while (fit_bad >= 1 && fit_bad <= max_it)
            [fit_struct, gof] = fit (fitx, fity, fit_f, ...
                'Lower' ,low_vect, ...
                'Upper' ,up_vect, 'StartPoint', start_vect);
            if (gof.adjrsquare >= goodness_value || fit_bad >= max_it)  % good perf fit
                fit_bad0 = fit_bad;
                fit_bad = 0;
            else
                fit_bad = fit_bad + 1;
                if length(coeff_fit)==6
                    start_vect = [max(fity(round(length(fity)/2)+1:end))*0.8; fit_struct.b1 ; fit_struct.c1; max(fity(1:round(length(fity)/2)))*0.8; fit_struct.b2 ; fit_struct.c2];
                end
            end
        end
        
        intervals = confint(fit_struct); % % a1 b1 c1 a2 b2 c2
        
        switch fit_chosen
            
            case Nmax+2 % 2 lorentzian
                
                vect_fit = (fit_struct.a1)./(1 + ((xfine-fit_struct.b1)./(fit_struct.c1)).^2) + ...
                    (fit_struct.a2)./(1 + ((xfine-fit_struct.b2)./(fit_struct.c2)).^2);
                
                meth_fit = 'Lorentzian';
                b(1) = fit_struct.b2; b(2) = fit_struct.b1; % inversed !
                c(1) = 2*fit_struct.c2; c(2) = 2*fit_struct.c1; % inversed !
                intb1 = -1*(intervals(1,2)-intervals(2,2))*0.5;
                intb2 = -1*(intervals(1,5)-intervals(2,5))*0.5;
                intc1 = -(intervals(1,3)-intervals(2,3));
                intc2 = -(intervals(1,6)-intervals(2,6));
                
            case Nmax+3 % 2 prod lorentzian/gaussian
                
                vect_fit = fit_struct.a1./(1 + ((xfine-fit_struct.b1)./fit_struct.c1).^2).*fit_struct.a2.*exp(-((xfine-fit_struct.b2)./(fit_struct.c2)).^2).*fit_struct.a11./(1 + ((xfine-fit_struct.b11)./fit_struct.c11).^2).*fit_struct.a22.*exp(-((xfine-fit_struct.b22)./(fit_struct.c22)).^2);
                meth_fit = 'Gauss*Lorentz';
                
                % \\\\\\\\ CODE for finding the exact FWHM : applicable at
                % other functions ////////
                [max_vectfit1, I_maxg] = max(vect_fit);
                vectfit_g = [vect_fit(1:max(round(I_maxg-length(xfine)/8), 1)), vect_fit(min(round(I_maxg+length(xfine)/8), length(vect_fit)):end)];
                [max_vectfit2, I_maxd] = max(vectfit_g);
                vectfit_d = [vect_fit(1:max(round(I_maxd-length(xfine)/8), 1)), vect_fit(min(round(I_maxd+length(xfine)/8), length(vect_fit)):end)];
                
                %                     if I_maxm > I_maxM
                %                         I_maxm = I_maxm +  length(xfine)/4;
                %                     end
                
                [~,I_maxM1] =  min(abs(vectfit_d - max_vectfit1/2));
                [~,I_maxm1] =  min(abs([vectfit_d(1:max(round(I_maxM1-length(xfine)/12), 1)), vectfit_d(min(round(I_maxM1+length(xfine)/12), length(vect_fit)):end)] - max_vectfit1/2));
                if I_maxm1 > I_maxM1
                    I_maxm1 = I_maxm1 +  length(xfine)/12;
                end
                [~,I_maxM2] =  min(abs(vectfit_g - max_vectfit2/2));
                [~,I_maxm2] =  min(abs([vectfit_g(1:max(round(I_maxM2-length(xfine)/12), 1)), vectfit_g(min(round(I_maxM2+length(xfine)/12), length(vect_fit)):end)] - max_vectfit2/2));
                if I_maxm2 > I_maxM2
                    I_maxm2 = I_maxm2 +  length(xfine)/12;
                end
                % \\\\\\\\ END of CODE for finding the exact FWHM : applicable at
                % other functions /////////////
                
                b(2) = fit_struct.b1; %sign(fit_struct.b1)*sqrt(abs(fit_struct.b1 * fit_struct.b2));
                c(2) = abs(I_maxM2-I_maxm2)/length(xfine)*2;%fit_struct.c1; %sqrt(fit_struct.c1*sqrt(fit_struct.c2*log(2))); % real FWHM
                b(1) = fit_struct.b11; %sign(fit_struct.b11)*sqrt(abs(fit_struct.b11 * fit_struct.b22));
                c(1) = abs(I_maxM1-I_maxm1)/length(xfine)*2; %sqrt(fit_struct.c11*sqrt(fit_struct.c22*log(2))); % real FWHM
                
                intb1 = -(intervals(1,2)-intervals(2,2))*0.5;
                intb2 = -(intervals(1,5)-intervals(2,5))*0.5;
                intc1 = -1*(intervals(1,3)-intervals(2,3))*0.5;
                intc2 = -1*(intervals(1,6)-intervals(2,6))*0.5;
                
            case Nmax+4 % 2 voigt
                
                vect_fit = fit_struct.a1.*((1-eta)./(1 + ((xfine-fit_struct.b1)/fit_struct.c1).^2) + eta*exp(-((xfine-fit_struct.b2)./(fit_struct.c2)).^2))+ ...
                    fit_struct.a11.*((1-eta)./(1 + ((xfine-fit_struct.b11)./fit_struct.c11).^2) + eta.*exp(-((xfine-fit_struct.b22)./(fit_struct.c22)).^2));
                meth_fit = 'pseudo-Voigt';
                b(2) = (fit_struct.b2); % + fit_struct.b2);
                c(2) = (1-eta)*2*fit_struct.c1 + (eta)*sqrt(log(2)*fit_struct.c2); % + 2*fit_struct.c2); % real FWHM
                b(1) = fit_struct.b22; % + fit_struct.b11);
                c(1) = (1-eta)*2*fit_struct.c11 + (eta)*sqrt(log(2)*fit_struct.c22);%2*fit_struct.c11; % + 2*fit_struct.c22); % real FWHM
                
                intb1 = -(intervals(1,9)-intervals(2,9))*0.5;
                intb2 = -(intervals(1,4)-intervals(2,4))*0.5;
                intc1 = -(1-eta)*2*(intervals(1,3)-intervals(2,3))*0.5-...
                    (eta)*sqrt(log(2)*(intervals(1,5)-intervals(2,5))*0.5);
                intc2 =-(1-eta)*2*(intervals(1,7)-intervals(2,7))*0.5-...
                    (eta)*sqrt(log(2)*(intervals(1,9)-intervals(2,9))*0.5);
                
            case Nmax+5 % 2 Pearson VII
                
                vect_fit = (fit_struct.a1)./(1 + ((xfine-fit_struct.b1)./(fit_struct.c1)).^2).^M + ...
                    (fit_struct.a2)./(1 + ((xfine-fit_struct.b2)./(fit_struct.c2)).^2).^M;
                
                meth_fit = sprintf('Pearson VII (M = %.2g)', M);
                b(1) = fit_struct.b2; b(2) = fit_struct.b1; % inverted !
                coef = sqrt(4*abs(2^(1/M)-1));
                c(1) = fit_struct.c2*coef;
                c(2) = fit_struct.c1*coef; % inverted !
                
                intb1 = -(intervals(1,2)-intervals(2,2))*0.5;
                intb2 = -(intervals(1,5)-intervals(2,5))*0.5;
                intc1 = -coef*(intervals(1,3)-intervals(2,3))*0.5;
                intc2 = -coef*(intervals(1,6)-intervals(2,6))*0.5;
                
            case Nmax+6 % Von Mises
                
                vect_fit = fit_struct.a1.*exp(fit_struct.c1.*cos(xfine-fit_struct.b1))./(2*3.14*besseli(0,fit_struct.c1)) + ...
                    fit_struct.a2.*exp(fit_struct.c2.*cos(xfine-fit_struct.b2))./(2*3.14*besseli(0,fit_struct.c2));
                % c1 = 2*sigma^2, FWHM = 2*sqrt(2*log(2))*sigma for
                % GAUSS
                % 1/a1 = sigma^2
                % so FWHM = 2*sqrt(2*log(2)/a1)
                
                b(1) = fit_struct.b1;
                c(1) = 2*sqrt(2*log(2)/fit_struct.a1); % real FWHM
                b(2) = fit_struct.b2;
                c(2) = 2*sqrt(2*log(2)/fit_struct.a2); % real FWHM
                intb1 = -(intervals(1,2)-intervals(2,2))*0.5;
                intb2 = -(intervals(1,5)-intervals(2,5))*0.5;
                intc1 = -2*sqrt(2*log(2)/((intervals(1,1)-intervals(2,1))*0.5));
                intc2 = -2*sqrt(2*log(2)/((intervals(1,4)-intervals(2,4))*0.5));
        end
    end
    
    [~, I_maxM] = max(vect_fit);
    [~, I_maxm] = max([vect_fit(1:max(round(I_maxM-length(xfine)/div_fact), 1)), vect_fit(min(round(I_maxM+length(xfine)/div_fact), length(vect_fit)):end)]);
    if I_maxm > I_maxM
        I_maxm = I_maxm +  length(xfine)/(div_fact/2);
    end
    
    % %     % !!!!!!
    % %     figure; plot([vect_fit(1:max(round(I_maxM-length(xfine)/div_fact), 1)), vect_fit(min(round(I_maxM+length(xfine)/div_fact), length(vect_fit)):end)])
    % %
    % %     if I_maxM <= round(length(vect_fit)/2)
    % %
    % %         [~, I_maxm] = max(vect_fit(round(length(vect_fit)/2)+1:end));
    % %     else
    % %         [~, I_maxm] = max(vect_fit(1:round(length(vect_fit)/2)));
    % %     end
    mid = round((I_maxm + I_maxM)/2);
    
    vect_fit_struct{1} = vect_fit(1:mid);
    vect_fit_struct{2} =  vect_fit(mid+1:end);
    
    % to correct the plot when norm has been applied
    if (norm_vect || norm_imposed)
        if norm_imposed; max_first =max1;max_nd=max1;
        else% normal
            if max_increase_with_array_index
                max_first = max1;
                max_nd = max2;
            else % 2nd max is before first one
                max_first = max2;
                max_nd = max1;
            end
        end
        vect_fit_struct{1} = vect_fit_struct{1}*max_first;
        vect_fit_struct{2} = vect_fit_struct{2}*max_nd;
    end
    if baseline_hist>0
        vect_fit_struct{1} = vect_fit_struct{1}+baseline_hist;
        vect_fit_struct{2} = vect_fit_struct{2}+baseline_hist;
    end
    
    axes(h_ax); hold on;
    plot(h_ax, xfine(1:mid), vect_fit_struct{1},'Color', clr{1},'LineWidth',2.5);
    plot(h_ax, xfine(mid+1:end), vect_fit_struct{2},'Color', clr{2},'LineWidth',2.5);
    
    fprintf('Number of fit attempt = %d \n R-square adj. = %.3g (->1 for good fit)\n\n', fit_bad0, gof.adjrsquare);
    if fit_bad0 == max_it
        fprintf('Tried %d times, but no good fit\n', max_it);
    end
    %     legend(Legendehisto,['Fit ' meth_fit ' 1 - \mu = ' stringb1 '\pi, \sigma = ' stringc1 '\pi' ' ; ' 'Fit ' meth_fit ' 2 - \mu = ' stringb2 '\pi, \sigma = ' stringc2 '\pi'], 'Location', 'NorthEast');
    
    vect_fit2 = zeros(1, length(vect_fit));
    vect_fit2(1:end-round(offset/2*length(vect_fit))) = vect_fit(round(offset/2*length(vect_fit))+1:end);
    vect_fit2(round(offset/2*length(vect_fit))+1:end) = vect_fit(1:end-round(offset/2*length(vect_fit)));
    
    cmap1 = zeros(length(vect_fit2),3); cmap1(:,1) = (vect_fit2'-min(vect_fit2(1:mid)))/(max(vect_fit2(1:mid))-min(vect_fit2(1:mid))); % defining red values
    cmap1(:,2) = (vect_fit2'-min(vect_fit2(mid+1:end)))/(max(vect_fit2(mid+1:end))-min(vect_fit2(mid+1:end)));
    cmap1(max(1, mid+1-round(offset/2*length(vect_fit))):end - round(offset/2*length(vect_fit)),1) = 0;
    cmap1(1:max(1, mid-round(offset/2*length(vect_fit))), 2) = 0;
    cmap1(end - round(offset/2*length(vect_fit)):end, 2) = 0;
    
end
%% plot

stringb1 = sprintf('%.4f', round(b(1)*10000)/10000 - offset);
stringc1 = sprintf('%.4f', round(c(1)*10000)/10000 - offset);
stringb2 = sprintf('%.4f', round(b(2)*10000)/10000 - offset);
stringc2 = sprintf('%.4f', round(c(2)*10000)/10000 - offset);
if length(intb1) > 1; intb2=intb1(end);intb1=intb1(1);end
if length(intc1) > 1; intc2=intc1(end);intc1=intc1(1);end
disp1= @(intb1, intc1) ['Fit ' meth_fit ' 1 : \mu = ' stringb1 intb1 '\pi, FWHM = ' stringc1 intc1 '\pi'];
disp2=@(intb2, intc2) ['Fit ' meth_fit ' 2 : \mu = ' stringb2 intb2 '\pi, FWHM = ' stringc2 intc2 '\pi'];
disp(disp1(['±', num2str(round(intb1,4))], ['±', num2str(round(intc1,4))]));
disp(disp2(['±', num2str(round(intb2,4))], ['±', num2str(round(intc2,4))]));

hh_l = legend(h_ax, [str_l, disp1('',''), disp2('','')]);

set(hh_l,  'Location', 'NorthEast');
set(hh_l,  'FontSize', round(16 - nb/2));
warning('on','curvefit:fit:noStartPoint')

cmap1(end,:) = [1 1 1]; % last value is for sat

if (2*str2double(stringc1) > (1-abs(str2double(stringb1))) || 2*str2double(stringc2) > (1-abs(str2double(stringb2))))
    fit_warning = 'verify fit !';
else
    fit_warning = 'ok';
end
if isa(vect_fit_struct, 'cell')
    if fit_chosen <= Nmax % % fit MANU
        vect_fit_struct{end+1} = vect_fit_struct{1}+vect_fit_struct{2};
    else % auto
        vect_fit_struct{end+1} = [vect_fit_struct{1}, vect_fit_struct{2}];
    end
end
% x = -32:1:31; r = exp(-(x+16).^2/(2*4^2));g = exp(-(x-16).^2/(2*4^2));
% cmap1 = zeros(64,3);cmap1(:,1) = r;cmap1(:,2) = g;
%
% figure; imagesc(magic(25)); colormap(cmap1); colorbar;
end

function [norm_vect, max_increase_with_array_index, fity, max1, max2]=norm_peaks_hist( hist2_ydata, div_fact, diff_peak_max, fact_norm)
% code to normalize the peaks if necessary
fity = hist2_ydata';
norm_vect = 0; max_increase_with_array_index = 0; % dflt
[max1, I_maxM] = max(hist2_ydata);
vect_trou = [hist2_ydata(1:max(round(I_maxM-length(hist2_ydata)/div_fact), 1)), hist2_ydata(min(round(I_maxM+length(hist2_ydata)/div_fact), length(hist2_ydata)):end)];
[max2, I_maxm] = max(vect_trou);
if I_maxm > I_maxM
    %             I_maxm = I_maxm +  length(hist2_ydata)/(div_fact/2);
    max_increase_with_array_index = 1;
end
% %         figure; plot(vect_trou)
max1 = fact_norm*max1;
% ALWAYS : max2< max1
if max(max1, max2) > diff_peak_max*min(max1, max2)
    fity=fity';
    fity(:, :) = [fity(1:max(round(I_maxM-length(hist2_ydata)/div_fact), 1))/max2, ...
        fity(max(round(I_maxM-length(hist2_ydata)/div_fact), 1)+1:min(round(I_maxM+length(hist2_ydata)/div_fact), length(hist2_ydata))-1)/max1, ...
        fity(min(round(I_maxM+length(hist2_ydata)/div_fact), length(hist2_ydata)):end)/max2];
    fity=fity';
    norm_vect = 1;
end
% END OF code to normalize the peaks if necessary
end
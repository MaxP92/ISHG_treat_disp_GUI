function phcorr = osc_corr_movemean(ph0)
% osc_corr_movemean
% % Max Pinsard 2019.10.01

prompt = {'method:movemean, sgolayfilt' , 'Center of the peak to exclude(in /pi, usually 0.5) ?' ; ...
            'width at bottom of this peak (in /pi, e.g. 1.0 for large, 0 for NO', 'size of kernel for avg'};
def = {'movemean', '0.5', '1', '10'};
answer = inputdlg(prompt,'Params corr. oscillations',1, def);
meth =answer{1} ; center=str2double(answer{2}); width=str2double(answer{3}); kernelsz=str2double(answer{4});
if width~=0
    ph=ph0; 
    borne_sup = center+width/2;            
    borne_inf = center-width/2;
%     correct=ph(ph > borne_sup | ph < borne_inf);
    ph(ph <= borne_sup & ph >= borne_inf)= ph(ph <= borne_sup & ph >= borne_inf)-1 ; %median(correct);
    ph(ph<-1)=ph(ph<-1)+2; ph(ph>1)=ph(ph>1)-2;
     %ph(ph>= center-abs(width)/2 | ph<= center+abs(width)/2)=ph(ph>= center-abs(width)/2 | ph<= center+abs(width)/2)
else; ph =ph0;
end
switch meth
    case 'movemean'
        mm=movmean(ph,kernelsz);
        figure(66); subplot(2,1,1); imagesc(ph);  colormap(jet); colorbar; axis image; title('same range');%caxis( [-sat_value, sat_value]);
subplot(2,1,2); imagesc(mm);axis image;colormap(jet); colorbar; title('movemean-ed');

        phcorr=ph0-mm; % > R2016a
    case ''
        ss2=sgolayfilt(double(ph4),2,kernelsz+(1-mod(kernelsz,2)),[],2);phcorr=sgolayfilt(ss2,2,kernelsz+(1-mod(kernelsz,2)),[],1);
end

figure(135); subplot(2,1,1); imagesc(ph0);  colormap(hsv); colorbar; axis image; title('00');%caxis( [-sat_value, sat_value]);
subplot(2,1,2); imagesc(phcorr);  colormap(hsv); colorbar; axis image; title('corr');
end
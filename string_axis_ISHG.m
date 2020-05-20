function [ Titre1, Titre2, Titre3, Titre4, Titre5, Titre6, Counts, Yaxis1, Yaxis2, Legendehisto ] = string_axis_ISHG( langue )
% [ Titre1, Titre2, Titre3, Titre4, Titre5, Titre6, Counts, Yaxis1, Yaxis2, Legendehisto ] = string_axis_ISHG( langue )
% 
% 2015-10 edited by Maxime PINSARD
% 
%  The strings for xlabel, title, lengends etc. of axis

% Le texte qui servira � la construction des graphiques est ennregistr� ici

if langue==2 % fr
    Titre1 = 'Phase relative';
    Titre2 = 'Contraste interf�rom�trique';
    Titre3 = 'Erreur relative';
    Titre4 = 'Distribution de la phase';
    Titre5 = sprintf('Distribution en fonction de la phase \n et de l''intensit� GSH'); 
    Titre6 = sprintf('Distribution en fonction de la phase \n et du contraste interf�rom�trique');
    Counts = 'Comptes (nombre de pixels)';
    Yaxis1 = 'Intensit� GSH (a.u.)';
    Yaxis2 = 'Contraste interf�rom�trique (a.u.)';
    Legendehisto = 'Histogramme des donnees experimentales';
else % English
    Titre1 = 'Relative phase';
    Titre2 = 'Interferometric contrast';
    Titre3 = 'Relative error';
    Titre4 = 'Phase distribution';
    Titre5 = 'Phase and SHG intensity distribution';
    Titre6 = 'Phase and interferometric contrast distribution';
    Counts = 'Counts (number of pixels)';
    Yaxis1 = 'SHG intensity (a.u.)';
    Yaxis2 = 'Interferometric contrast (a.u.)';
    Legendehisto = 'Experimental data histogram';
end

end


% refraction of air
%
% 2016.4.28 Maxime PINSARD

clearvars; close all; clc

path00 = 'C:\Users\Maxime\Documents\';
addpath(fullfile(path00, 'These\codes Matlab\Variety'));

prompt = {'Methode (0 for simple Ciddor, 1 for index table, 2 for ciddor)',...
'Max lambda (um)','Min lambda (um)', 'Temperature (deg C)', 'Pressure (Pa)', 'CO2 (ppm)', 'Humidity (%)'};
dlg_title = 'Paramètres';
num_lines = 1;
% Valeurs par défaut
def = {'2','0.8','0.4','20', '1e5', '450', '20'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
% Les réponses en caractères sont converties en chiffres qui sont enregistrés dans des variables.
meth = str2double(cell2mat(answer(1)));
% 0 for simple, 1 for table, 2 for ciddor

range_l = [str2double(cell2mat(answer(2))) str2double(cell2mat(answer(3)))]; % in um
% range_l = [1.064 0.532]; % in um

P = str2double(cell2mat(answer(5))); % Pa
T = 273.15 + str2double(cell2mat(answer(4))); %K
co2 = str2double(cell2mat(answer(6))); % Pppm
H = str2double(cell2mat(answer(7)));% in %

diff_marche_corr_20X_mm = 2.5; %   mm, totale
diff_marche_corr_40X_mm = 2.0;%  mm, totale

Ps = 101325; %Pa
Ts = 288.15; %K

ct = 1;

tabl = [0.2300,1.0003080029552;
    0.2446,1.0003030361805;
    0.2592,1.0002991101754;
    0.2738,1.000295937727;
    0.2884,1.0002933286388;
    0.3030,1.000291151532;
    0.3176,1.000289312581;
    0.3322,1.000287742941;
    0.3468,1.0002863909483;
    0.3614,1.0002852170881;
    0.3760,1.0002841906396;
    0.3906,1.0002832873764;
    0.4052,1.0002824879531;
    0.4198,1.000281776749;
    0.4344,1.0002811410247;
    0.4490,1.0002805702966;
    0.4636,1.0002800558675;
    0.4782,1.0002795904676;
    0.4928,1.000279167979;
    0.5074,1.0002787832208;
    0.5220,1.0002784317799;
    0.5366,1.0002781098769;
    0.5512,1.0002778142591;
    0.5658,1.0002775421136;
    0.5804,1.0002772909975;
    0.5950,1.0002770587802;
    0.6096,1.0002768435966;
    0.6242,1.0002766438078;
    0.6388,1.0002764579689;
    0.6534,1.0002762848018;
    0.6680,1.0002761231724;
    0.6826,1.0002759720715;
    0.6972,1.0002758305987;
    0.7118,1.0002756979483;
    0.7264,1.0002755733975;
    0.7410,1.0002754562966;
    0.7556,1.0002753460601;
    0.7702,1.0002752421591;
    0.7848,1.0002751441149;
    0.7994,1.0002750514935;
    0.8140,1.0002749639005;
    0.8286,1.0002748809769;
    0.8432,1.0002748023951;
    0.8578,1.0002747278561;
    0.8724,1.0002746570862;
    0.8870,1.0002745898346;
    0.9016,1.000274525871;
    0.9162,1.0002744649838;
    0.9308,1.0002744069782;
    0.9454,1.0002743516746;
    0.9600,1.0002742989071;
    0.9746,1.0002742485227;
    0.9892,1.0002742003794;
    1.004,1.0002741537295;
    1.018,1.0002741114819;
    1.033,1.0002740681302;
    1.048,1.0002740266464;
    1.062,1.0002739895199;
    1.077,1.0002739513533;
    1.091,1.0002739171561;
    1.106,1.0002738819616;
    1.121,1.0002738481835;
    1.135,1.0002738178689;
    1.150,1.0002737866205;
    1.164,1.0002737585479;
    1.179,1.0002737295823;
    1.194,1.0002737017104;
    1.208,1.0002736766354;
    1.223,1.0002736507265;
    1.237,1.0002736273968;
    1.252,1.0002736032705;
    1.267,1.0002735800022;
    1.281,1.0002735590237;
    1.296,1.0002735373021;
    1.310,1.0002735177028;
    1.325,1.0002734973936;
    1.340,1.000273477767;
    1.354,1.000273460038;
    1.369,1.0002734416467;
    1.383,1.0002734250218;
    1.398,1.0002734077641;
    1.413,1.0002733910563;
    1.427,1.000273375938;
    1.442,1.0002733602287;
    1.456,1.000273346005;
    1.471,1.0002733312162;
    1.486,1.0002733168752;
    1.500,1.0002733038787;
    1.515,1.0002732903537;
    1.529,1.0002732780896;
    1.544,1.0002732653197;
    1.559,1.0002732529184;
    1.573,1.000273241664;
    1.588,1.000273229936;
    1.602,1.0002732192871;
    1.617,1.0002732081845;
    1.632,1.0002731973879;
    1.646,1.0002731875774;
    1.661,1.0002731773414;
    1.675,1.0002731680358;
    1.690,1.0002731583221];

for lambda = range_l
    
    switch meth
        case 0 % simple formula
            ns = 1 + 0.05792105/(238.0185-lambda^-2) + 0.00167917/(57.362-lambda^-2);
            % too little precision
            lambdareal = lambda;
             n = 1 + (ns - 1) * (P / Ps) * (Ts / T);
        case 1 % table from refractive index site
            tmp = abs(tabl(:,1)-lambda);
            [~, idx] = min(tmp); %index of closest value
            ns = tabl(idx,2);
            lambdareal = tabl(idx,1);
            n = 1 + (ns - 1) * (P / Ps) * (Ts / T);
            
        case 2 % ciddor equation, take humidity into account
            lambdareal = lambda;
            n = ciddor(lambdareal, co2, T-273.15, P, H);
            % OUTPUTS:
            % * n = refractive index (dimensionless)s
            % INPUTS:
            % * w = vacuum wavelength of laser (micron) valid from 0.3 to 1.7
            % * c = CO2 Level (ppm mole) valid from 0 to 2000, typically 450
            % * t = Air Temperature (C) valid from -40 to 100
            % * p = Air Pressure (Pa) valid from 1E4 to 14E4 typically 101325
            % * h = Relative Humidity (%) valid from 0 to 100
            % * All inputs may be either scalars or column vectors
            % * All inputs must have the same dimensions
    end
    
    struct(ct).lambda = lambdareal;
    struct(ct).n = n;
    
    ct = ct + 1;
end

delta_phase_lineique = 2*pi*1e4*(struct(2).n - struct(1).n)/min(struct(1).lambda, struct(2).lambda); % rad/cm

dist_deph_2pi_cm = 2*pi/delta_phase_lineique;

prop_corr_20X_m = 1e-3*diff_marche_corr_20X_mm/(struct(2).n - struct(1).n); %   mm, totale
prop_corr_40X_m = 1e-3*diff_marche_corr_40X_mm/(struct(2).n - struct(1).n);%  mm, totale

fprintf('Deph. lineique entre %g nm et sa SHG a %g nm = %g rad/cm\n Pour dephaser de 2pi, propagation de %g cm\n Corriger dispersion 20X, propager %.4g m\n Corriger dispersion 40X, propager %.4g m\n',...
max(struct(2).lambda, struct(1).lambda), min(struct(2).lambda, struct(1).lambda), delta_phase_lineique, dist_deph_2pi_cm, prop_corr_20X_m, prop_corr_40X_m);



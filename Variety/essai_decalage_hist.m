% essai decalage hist

clc; clear all; close all

cd('C:\Users\Maxime\Documents\These\codes Matlab\');

load('C:\Users\Maxime\Documents\These\codes Matlab\exemple_hist_phase.mat');

int_x = 0.025;
int_y = 100;

hhist1 = axes;

hist2_xdata = -1+int_x/2:int_x:1-int_x/2;
 Counts = 'Counts';

phi_mat_default = '\boldmath$\phi_{mat}$';

 draw_plots_ISHG( 2, 0, phase_test_hist, hist2_xdata, 0, hhist1, ...
        phi_mat_default, Counts, '', 0, 0, '', 0, 0, ...
        16, 14, 14, 16, 14 );
    
    % offset
    offset = 0.5; % offset of pi/2
    
    phase_test_hist = phase_test_hist + offset;
    
   phase_test_hist(phase_test_hist>1) = phase_test_hist(phase_test_hist>1)-2;
   
   figure; hhist2= axes;
    draw_plots_ISHG( 2, 0, phase_test_hist, hist2_xdata, 0, hhist2, ...
        phi_mat_default, Counts, '', 0, 0, '', 0, 0, ...
        16, 14, 14, 16, 14 );
      mbversion = version; % query the MatLab version
    if (str2double(mbversion(end-3)) > 0 && str2double(mbversion(end-2)) >= 4) % Matlab releases posterior to R2014a
        set(hhist2,'XTickLabel',{'\pi /2' '-\pi' '-\pi /2' '0' '\pi /2'}, 'TickLabelInterpreter', 'tex');
    else % because of Matalb f***ing update !
        set(hhist2,'XTickLabel','p /2 | -p | -p /2 | 0 | p /2', 'fontname', 'symbol');
    end
    
    
% Save in PNG the plots of I = f(U) for different times
%
% it reads text files TAB delimited
%
% 207.10.27 Maxime PINSARD

clearvars
close all
clc

path00 = 'C:\Users\pc\Desktop\cybelle';

cd(path00);

%% opening the file, path

[fname, folder_name, FILTERINDEX] = uigetfile('*.txt','Select your .txt files !','file.txt','MultiSelect','on');

if (~FILTERINDEX)
    error('File not chosen ! Program ends.');
end

% path_save = fullfile(folder_name, fname(1:end-4));
date_now = char(datetime('now','TimeZone','local','Format','d_MMM_y_HH.mm.ss'));
path_save = fullfile(folder_name, date_now);
mkdir(folder_name, date_now); % creates a new folder for results

%% parameters

prompt = {'Input sampling (ms)', 'Output saampling (ms)'};
dlg_title = 'Parameters';
num_lines = 1;
% Valeurs par défaut
def = {'1.5', '18'};
answer = inputdlg(prompt,dlg_title,num_lines,def);

if ~isempty(answer)
    sampling_input = str2double(answer{1});
    sampling_output = str2double(answer{2});
else % empty, default
    sampling_input = str2double(def{1});
    sampling_output = str2double(def{2});
end

if mod(sampling_output, sampling_input)
    warning('The output sampling is not a multiple of input sampling, I will approximate !');
end

%% reading the files

% time_corr_s{num_file} = cell(1, length(fname));
% current_A{num_file} = cell(1, length(fname));
% potential_V{num_file} = cell(1, length(fname));
% time_dwn_ms{num_file} = cell(1, length(fname));
% current_dwn_mA{num_file} = cell(1, length(fname));
% potential_dwn_mV{num_file} = cell(1, length(fname));
num_lines = 667;
time_corr_s = zeros(num_lines, length(fname));
current_A = zeros(num_lines, length(fname));
potential_V = zeros(num_lines, length(fname));

num_lines_dwn = round(num_lines/sampling_output*sampling_input);
time_dwn_ms = zeros(num_lines_dwn, length(fname));
current_dwn_mA = zeros(num_lines_dwn, length(fname));
potential_dwn_mV = zeros(num_lines_dwn, length(fname));


for num_file = 1:length(fname)
    
    content = dlmread(fullfile(folder_name, fname{num_file}), '\t', 1, 0);
    % it reads from the 2nd row (1st contains text)
    
    time_corr_s(:,  num_file) = content(:, 1);
    current_A(:, num_file) = content(:, 2);
    potential_V(:, num_file) = content(:, 3);
    
    %% downsampling
    
    %     time_dwn_ms{num_file} = time_corr_s{num_file}(1:round(sampling_output/sampling_input):end)*1000;
    %     current_dwn_mA{num_file} = current_A{num_file}(1:round(sampling_output/sampling_input):end)*1000;
    %     potential_dwn_mV{num_file} = potential_V{num_file}(1:round(sampling_output/sampling_input):end)*1000;
    time_dwn_ms(:, num_file) = time_corr_s(1:round(sampling_output/sampling_input):end,  num_file)*1000;
    current_dwn_mA(:, num_file) = current_A(1:round(sampling_output/sampling_input):end,  num_file)*1000;
    potential_dwn_mV(:, num_file) = potential_V(1:round(sampling_output/sampling_input):end,  num_file)*1000;
    
end
%% plot

screensize = get( 0, 'Screensize' ); % to get screen size
fact = 4/5;
left_offset_fig = 40;
top_offset_fig = 100;

h_1 = figure(1);
set(1, 'outerposition',...
    [min(screensize(3)*(1-fact), left_offset_fig) min(screensize(4)*(1-fact), top_offset_fig) ...
    screensize(3)*fact screensize(4)*fact]);
h_axes = axes;

for ii = 1:length(time_dwn_ms(:, 1))
    plot(h_axes, potential_dwn_mV(ii, :), current_dwn_mA(ii, :), 'x','MarkerSize', 12);
    xlabel('Potential (mV)'); ylabel('Current (mA)');
    title(sprintf('time = %.2f ms', time_dwn_ms(ii, 1)));
    set(h_axes, 'FontSize', 14);
    
    name_full = fullfile(path_save, sprintf('curr_vs_pot_%.2f_ms.png', time_dwn_ms(ii)));
    
    saveas(h_1, name_full);
    
end

disp('Saving is finsihed ! ;)');









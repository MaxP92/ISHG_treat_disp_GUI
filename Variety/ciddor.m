function n = ciddor(w, c, t, p, h)
% ciddor: Calculates refractive index using ciddor equation
% Dr Jody Muelaner 2015 jody@muelaner.com
% Department of Mechanical Engineering, The University of Bath
% 
% [n] = ciddor(w, c, t, p, h)
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
% 
% Sub-functions:
% * humidity_as_mole_fraction function
% * Saturation_Vapor_Pressure
% * Saturation_Vapor_Pressure_ice
%
% Validaton:
% Validated against the NIST Engineering Metrology Toolbox
% to 9 decimal places i.e n=1.000296813
% emtoolbox.nist.gov/

%% Validate inputs (this section can be removed if not required)

%Find size of w
sizew = size(w); rowsw = sizew(1,1); colsw = sizew(1,2);
%Find size of c
sizec = size(c); rowsc = sizec(1,1); colsc = sizec(1,2);
%Find size of t
sizet = size(t); rowst = sizet(1,1); colst = sizet(1,2);
%Find size of p
sizep = size(p); rowsp = sizep(1,1); colsp = sizep(1,2);
%Find size of h
sizeh = size(h); rowsh = sizeh(1,1); colsh = sizeh(1,2);

%Check all inputs are scalar or column vectors
if colsw == 1 && colsc == 1 && colst == 1 && colsp == 1 && colsh == 1
    %Do nothing
else
     errordlg('Error in ciddor function to calculate refractive index: All inputs must be either scalars or column vectors')
    RETURN
end   

%Check all inputs are same size
if rowsc==rowsw && rowst==rowsw && rowsp==rowsw && rowsh==rowsw
    %Do nothing
else
    errordlg('Error in ciddor function to calculate refractive index: All inputs must be the same size')
    return
end

%Check wavelength is within valid range
lessthan = w <= 1.7;
lessthan = sum(lessthan);
greaterthan = w >= 0.3;
greaterthan = sum(greaterthan);

if greaterthan == rowsw && lessthan == rowsw
    %Do nothing
else
    errordlg('Refractive index not accurate: Ciddor function requires that wavelength (w) must be between 0.3 and 1.7 (microns)')
end

%Check CO2 level is within valid range
lessthan = c <= 2000;
lessthan = sum(lessthan);
greaterthan = c >= 0;
greaterthan = sum(greaterthan);

if greaterthan == rowsc && lessthan == rowsc
    %Do nothing
else
    errordlg('Refractive index not accurate: Ciddor function requires that CO2 level (c) must be between 0 and 2000 (ppm mole)')
end

%Check temperature is within valid range
lessthan = t <= 100;
lessthan = sum(lessthan);
greaterthan = t >= -40;
greaterthan = sum(greaterthan);

if greaterthan == rowst && lessthan == rowst
    %Do nothing
else
    errordlg('Refractive index not accurate: Ciddor function requires that temperature (t) must be between -40 and 100 (deg C)')
end

%Check pressure is within valid range
lessthan = p <= 140000;
lessthan = sum(lessthan);
greaterthan = p >= 10000;
greaterthan = sum(greaterthan);

if greaterthan == rowsp && lessthan == rowsp
    %Do nothing
else
    errordlg('Refractive index not accurate: Ciddor function requires that Pressure (p) must be between 1E4 to 14E4 (pa)')
end


%Check Relative Humidity is within valid range
lessthan = h <= 100;
lessthan = sum(lessthan);
greaterthan = h >= 0;
greaterthan = sum(greaterthan);

if greaterthan == rowsh && lessthan == rowsh
    %Do nothing
else
    errordlg('Refractive index not accurate: Ciddor function requires that Relative Humidity (h) must be between 0 and 100 (%)')
end

%% Assign Values to Constants
w0 = 295.235;
w1 = 2.6422;
w2 = -0.03238;
w3 = 0.004028;
k0 = 238.0185;
k1 = 5792105;
k2 = 57.362;
k3 = 167917;
a0 = 1.58123 * (10 ^ (-6));
a1 = -2.9331 * (10 ^ (-8));
a2 = 1.1043 * (10 ^ (-10));
b0 = 5.707 * (10 ^ (-6));
b1 = -2.051 * (10 ^ (-8));
c0 = 1.9898 * (10 ^ (-4));
c1 = -2.376 * (10 ^ (-6));
d = 1.83 * (10 ^ (-11));
e = -0.765 * (10 ^ (-8));
pr1 = 101325;
tr1 = 288.15;
za = 0.9995922115;
r = 8.314472;
mv = 0.018015;

%% Perform Intermediate Calculations
s = 1 ./ (w .^ 2);
ras = (10 ^ (-8)) .* ((k1 ./ (k0 - s)) + (k3 ./ (k2 - s)));
rvs = 1.022 * (10 ^ (-8)) * (w0 + (w1 .* s) + (w2 .* (s .^ 2)) + (w3 .* (s .^ 3)));
ma = 0.0289635 + (1.2011 * (10 ^ (-8)) .* (c - 400));
raxs = ras .* (1 + (5.34 * (10 ^ (-7)) .* (c - 450)));
tk = t + 273.15;
% Calls sub-function to calculate relative humidity as mole fraction
x = humidity_as_mole_fraction(t, p, h);

zm1 = (p ./ tk);
zm2 = a0 + (a1 .* t) + (a2 .* (t .^ 2));
zm3 = (b0 + (b1 .* t)) .* x;
zm4 = (c0 + (c1 .* t)) .* (x .^ 2);
zm5 = ((p ./ tk) .^ 2) .* (d + (e .* (x .^ 2)));

zm = 1 - (zm1 .* (zm2 + zm3 + zm4)) + zm5;

% Calculate Zw
pref2 = 1333;
tref2 = 293.15;
zw1 = (pref2 / tref2);
zw2 = a0 + (a1 * (tref2 - 273.15)) + (a2 * ((tref2 - 273.15) ^ 2));
zw3 = (b0 + (b1 * (tref2 - 273.15))) * 1;
zw4 = (c0 + (c1 * (tref2 - 273.15))) * (1 ^ 2);
zw5 = ((pref2 / tref2) ^ 2) * (d + (e * (1 ^ 2)));

zw = 1 - (zw1 * (zw2 + zw3 + zw4)) + zw5;

rhoaxs = (pr1 .* ma) ./ (za * r * tr1);
rhov = (x .* p .* ma) ./ (zm * r .* tk) .* (1 - 1 .* (1 - mv ./ ma));
rhoa = ((1 - x) .* p .* ma) ./ (zm * r .* tk);
rhovs = (pref2 .* ma ./ (zw * r * tref2)) .* (1 - 1 .* (1 - mv ./ ma));


%% Final Result (the refractive index)
n = 1 + ((rhoa ./ rhoaxs) .* raxs) + ((rhov ./ rhovs) .* rvs);

end

%% Sub Funcition used to cacluate the humidity as a mole fraction
function [hmol] = humidity_as_mole_fraction(t, p, h)
% Calculates humidity_as_mole_fraction
% Dr Jody Muelaner 2015, University of Bath, jody@muelaner.com
% [hmol] = humidity_as_mole_fraction(t, p, h)
% OUTPUTS:
% * hmol = humidity as mole fraction (% mole)
% INPUTS:
% * t = Air Temperature (C)
% * p = Air Pressure (Pa)
% * h = Relative Humidity (%)
%
% Requires sub-functions:
% * Saturation_Vapor_Pressure
% * Saturation_Vapor_Pressure_ice

%% Assign Values to Constants
alpha = 1.00062;
beta = 3.14 * (10 ^ (-8));
gama = 5.6 * (10 ^ (-7));

%% Intermediate Calculations
f = alpha + (beta .* p) + (gama .* (t .^ 2));

%% Determine whether over water or ice and call relevant fuction to calculate the saturation vapor preassure
% Note: This section not vectorized
psv = zeros(max(size(t)),1); %preallocate psv before looping
for i=1:size(t)
if t > 0
    psv(i) = Saturation_Vapor_Pressure(t(i));
else
    psv(i) = Saturation_Vapor_Pressure_ice(t(i));
end %End if else statement
end %End loop through each temperature value

%% Final Result (the humidity as mole fraction)
hmol = ((h ./ 100) .* f .* psv) ./ p;
end

%% Sub Funcition used to cacluate the saturation vapour pressure over water
function [svp] = Saturation_Vapor_Pressure(t)
% Calculates the saturation vapour pressure over water using IAPWS formula
% Dr Jody Muelaner 2015, University of Bath, jody@muelaner.com
% [svp] = Saturation_Vapor_Pressure_ice(t)
% OUTPUTS:
% * svp = saturation vapour pressure (Pa)
% INPUTS:
% * t = Air Temperature (C)

%% Assign Values to Constants
k1 = 1.16705214528 * (10 ^ 3);
k2 = -7.24213167032 * (10 ^ 5);
k3 = -1.70738469401 * 10;
k4 = 1.20208247025 * (10 ^ 4);
k5 = -3.23255503223 * (10 ^ 6);
k6 = 1.49151086135 * 10;
k7 = -4.82326573616 * (10 ^ 3);
k8 = 4.05113405421 * (10 ^ 5);
k9 = -2.38555575678 * (10 ^ -1);
k10 = 6.50175348448 * (10 ^ 2);

%% Intermediate Calculations
tk = t + 273.15;
omega = tk + (k9 / (tk - k10));
a = (omega ^ 2) + (k1 * omega) + k2;
b = (k3 * (omega ^ 2)) + (k4 * omega) + k5;
c = (k6 * (omega ^ 2)) + (k7 * omega) + k8;
x = -b + sqrt((b ^ 2) - (4 * a * c));

%% Final Result (the saturation vapour pressure)
svp = (10 ^ 6) * (((2 * c) / x) ^ 4);
end

%% Sub Funcition used to cacluate the saturation vapour pressure over ice
function [svp] = Saturation_Vapor_Pressure_ice(t)
% Calculates the saturation vapour pressure over ice using IAPWS formula
% Dr Jody Muelaner 2015, University of Bath, jody@muelaner.com
% [svp] = Saturation_Vapor_Pressure_ice(t)
% OUTPUTS:
% * svp = saturation vapour pressure (?)
% INPUTS:
% * t = Air Temperature (C)

%% Assign Values to Constants
a1 = -13.928169;
a2 = 34.7078238;

%%Intermediate Calculations
tk = t + 273.15;
thi = tk / 273.16;
y = (a1 * (1 - (thi ^ (-1.5)))) + (a2 * (1 - (thi ^ (-1.25))));

%% Final Result (the saturation vapour pressure)
svp = 611.657 * ((exp(1)) ^ y);
end
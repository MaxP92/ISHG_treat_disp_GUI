% angle plaque de verre optgeo

y1 = 398;
x1 = 359;
x4 = 335;
x2 = 499;
x3 = 476;
y4 = 264;
y3 = 236;

r = (y1-y4)/(x1-x4); % 5.58

y2 = r*(x2-x3) + y3;

%% autre r
r = 0.5;
y1 = r*(x1-x4) + y4;
% r = (y1-y4)/(x1-x4); 
y2 = r*(x2-x3) + y3;

r2 = (y2-y1)/(x2-x1);

y3 = y4 + r2*(x3-x4);

% (y2-y3)/(x2-x3)

%% vitesse acquisition phase
N_phi = 5;
N_px = 500*500;
t_int = 1e-6; % s
f = 2;

t_motor = 0:1e-9:1e-3;

classic = N_phi*(N_px*t_int*f + t_motor);
px_px = N_px*(N_phi*(t_int + t_motor) +f*t_int);

figure; h = axes;
plot(t_motor, classic/60, 'b'); hold on;
plot(t_motor, px_px/60, 'r'); 
xlabel('t motor (s)')
ylabel('t TOT (min)')
legend('classic', 'px by px');
xscale('log');
set(h, 'XScale', 'log');
set(h, 'XLim', [0 3e-6]);
set(h, 'Fontsize', 16);

%% 

figure; x=0:1e-2:2.6;% en mm/s
plot(x, x./(4) + (800e-6)./x); % min a 0.06 mm/s
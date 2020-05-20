clearvars; clc;

lambda = sym('lambda');

% for calcite

% % n_e = sqrt(1+1.0856*lambda.^2/(lambda.^2 - 0.07897^2)+...
% % 0.0988*lambda.^2/(lambda.^2 - 0.142^2)+...
% % 0.317*lambda.^2/(lambda.^2 - 11.468^2)); % lambda in um
% % 
% % n_o = sqrt(1+0.8559*lambda.^2/(lambda.^2 - 0.0588^2)+...
% % 0.8391*lambda.^2/(lambda.^2 - 0.141^2)+...
% % 0.0009*lambda.^2/(lambda.^2 - 0.197^2)+...
% % 0.6845*lambda.^2/(lambda.^2 - 7.005^2)); % lambda in um

n_o = sqrt(1+ 0.73358749+ 0.96464345*lambda.^2/(lambda.^2-1.94325203e-2)+ 1.82831454*lambda.^2/(lambda.^2-120));
n_e = sqrt(1+ 0.35859695+ 0.82427830*lambda.^2/(lambda.^2-1.06689543e-2)+ 0.14429128*lambda.^2/(lambda.^2-120));

% for YVO4

n_o = sqrt(3.779 + 0.0705/(lambda.^2 - 0.04573^2) - 0.0097*lambda.^2); % lambda in um
n_e = sqrt(4.607 + 0.108/(lambda.^2 - 0.0525^2) - 0.0143*lambda.^2); % lambda in um


% for BK7
n_o = (1+1.04*(lambda/1000)^2/((lambda/1000)^2-0.006)+0.232*(lambda/1000)^2/((lambda/1000)^2-0.02)+1.01*(lambda/1000)^2/((lambda/1000)^2-103.56))^0.5;

dn_dlambda = diff(n_o);
ddn_dlambda = diff(dn_dlambda);

eval(subs(dn_dlambda, 810))*1e3;
eval(subs(ddn_dlambda, 810))*1e6;
eval(subs(n_o, 810))
eval(subs(n_o, 405))
% for N-SF11 SCHOTT

nn = sqrt(1 + 1.7375969*lambda.^2/(lambda.^2 - 0.01319)+...
0.3137*lambda.^2/(lambda.^2 - 0.0623)+...
1.899*lambda.^2/(lambda.^2 - 155.236));

% for SF11 SCHOTT

nn8 = sqrt(1 + 1.73845969*lambda.^2/(lambda.^2 - 0.013609)+...
0.31117*lambda.^2/(lambda.^2 - 0.06153)+...
1.1749*lambda.^2/(lambda.^2 - 121.923));

% for F2 SCHOTT

nn2 = sqrt(1 + 1.3453*lambda.^2/(lambda.^2 - 0.00997)+...
0.209*lambda.^2/(lambda.^2 - 0.0470)+...
0.937*lambda.^2/(lambda.^2 - 111.886));

% for SF10 SCHOTT

nn3 = sqrt(1 + 1.616*lambda.^2/(lambda.^2 - 0.0127)+...
0.259*lambda.^2/(lambda.^2 - 0.0582)+...
1.077*lambda.^2/(lambda.^2 - 116.601));

% for NSF6 SCHOTT

nn45 = sqrt(1 + 1.779*lambda.^2/(lambda.^2 - 0.0134)+...
0.339*lambda.^2/(lambda.^2 - 0.062)+...
2.087*lambda.^2/(lambda.^2 - 174.02));

% for NSF66 SCHOTT

nn4 = sqrt(1 + 2.024*lambda.^2/(lambda.^2 - 0.0147)+...
0.470*lambda.^2/(lambda.^2 - 0.0693)+...
2.599*lambda.^2/(lambda.^2 - 161.82));

% for BAF10 SCHOTT

nn44 = sqrt(1 + 1.585*lambda.^2/(lambda.^2 - 0.0093)+...
0.144*lambda.^2/(lambda.^2 - 0.0424)+...
1.085*lambda.^2/(lambda.^2 - 105.61));
n_e = nn4;

% for NSF5 SCHOTT

nn4 = sqrt(1 + 1.525*lambda.^2/(lambda.^2 - 0.01127)+...
0.187*lambda.^2/(lambda.^2 - 0.0588)+...
1.427*lambda.^2/(lambda.^2 - 129.14));
n_e = nn4;

% for LakL21 

nn5 = sqrt(1 + 1.227*lambda.^2/(lambda.^2 - 0.00602)+...
0.4207*lambda.^2/(lambda.^2 - 0.0197)+...
1.013*lambda.^2/(lambda.^2 - 88.44));

% for LakL22 

nn55 = sqrt(1 + 1.142*lambda.^2/(lambda.^2 - 0.00586)+...
0.535*lambda.^2/(lambda.^2 - 0.0199)+...
1.041*lambda.^2/(lambda.^2 - 100.83));

% for N-LaF35

nn9 = sqrt(1 + 1.5169*lambda.^2/(lambda.^2 - 0.00751)+...
0.4558*lambda.^2/(lambda.^2 - 0.0260)+...
1.0747*lambda.^2/(lambda.^2 - 80.59));

d_nn = diff(nn);
d_nn2 = diff(d_nn);

% for TiH6

nn98 = sqrt(1 + 1.7723*lambda.^2/(lambda.^2 - 0.0131)+...
0.346*lambda.^2/(lambda.^2 - 0.0614)+...
2.408*lambda.^2/(lambda.^2 - 200.75));

% for BaH11

nn99 = sqrt(2.718 - 0.00866*lambda.^2 + 0.021*lambda.^-2 +...
0.00015*lambda.^(-4)+ 0.0000327*lambda.^(-6) + 6.603e-8*lambda.^(-8));

% for SK11

nn991 = sqrt(1 + 1.1796*lambda.^2/(lambda.^2 - 0.0068)+...
0.2298*lambda.^2/(lambda.^2 - 0.02197)+...
0.936*lambda.^2/(lambda.^2 - 101.51));
n_e = nn4;

d_nn = diff(nn);
d_nn2 = diff(d_nn);

%% calc
d_n_e = diff(n_e);
d_n_o = diff(nn4);
dd2_n_e = diff(d_n_e);
dd2_n_o = diff(d_n_o);

d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));
d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
dd2_n_o_2w = eval(subs(dd2_n_o,lambda,0.405));
dd2_n_o_w = eval(subs(dd2_n_o,lambda,0.81));
dd2_n_e_2w = eval(subs(dd2_n_e,lambda,0.405));
dd2_n_e_w = eval(subs(dd2_n_e,lambda,0.81));

n_e_w = eval(subs(n_e,lambda,0.81));
n_e_2w = eval(subs(n_e,lambda,0.405));
n_o_2w = eval(subs(nn4,lambda,0.405));
n_o_w = eval(subs(nn4,lambda,0.81));
disp([n_o_w, n_o_2w, d_n_o_w,  d_n_o_2w, dd2_n_o_w, dd2_n_o_2w])
disp([n_e_w, n_e_2w, d_n_e_w,  d_n_e_2w, dd2_n_e_w, dd2_n_e_2w])

% n0 sqrt(1+ 0,73358749+ 0,96464345*$C9^2/($C9^2-1,94325203e-2)+ 1,82831454*$C9^2/($C9^2-120))
% ne sqrt(1+ 0,35859695+ 0,82427830*$C9^2/($C9^2-1,06689543e-2)+ 0,14429128*$C9^2/($C9^2-120))
% ((2058499170265015*$C9)/(562949953421312*($C9^2 - 120)) + (8688735763930779*$C9)/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)) - (2058499170265015*$C9^3)/(562949953421312*($C9^2 - 120)^2) - (8688735763930779*$C9^3)/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)^2))/(2*((2058499170265015*$C9^2)/(1125899906842624*$C9^2 - 135107988821114880) + (8688735763930779*$C9^2)/(9007199254740992*$C9^2 - 2800521317822387/16) + 3903691986989077/2251799813685248)^(1/2))
% ((1856109722364793*$C9)/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)) + (5198641238726495*$C9)/(18014398509481984*($C9^2 - 120)) - (1856109722364793*$C9^3)/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)^2) - (5198641238726495*$C9^3)/(18014398509481984*($C9^2 - 120)^2))/(2*((5198641238726495*$C9^2)/(36028797018963968*$C9^2 - 4323455642275676160) + (1856109722364793*$C9^2)/(2251799813685248*$C9^2 - 6150233422068845/256) + 1529644179441673/1125899906842624)^(1/2))
% (2058499170265015/(562949953421312*($C9^2 - 120)) - (10292495851325075*$C9^2)/(562949953421312*($C9^2 - 120)^2) + (2058499170265015*$C9^4)/(140737488355328*($C9^2 - 120)^3) + 8688735763930779/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)) - (43443678819653895*$C9^2)/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)^2) + (8688735763930779*$C9^4)/(1125899906842624*($C9^2 - 2800521317822387/144115188075855872)^3))/(2*((2058499170265015*$C9^2)/(1125899906842624*$C9^2 - 135107988821114880) + (8688735763930779*$C9^2)/(9007199254740992*$C9^2 - 2800521317822387/16) + 3903691986989077/2251799813685248)^(1/2)) - ((2058499170265015*$C9)/(562949953421312*($C9^2 - 120)) + (8688735763930779*$C9)/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)) - (2058499170265015*$C9^3)/(562949953421312*($C9^2 - 120)^2) - (8688735763930779*$C9^3)/(4503599627370496*($C9^2 - 2800521317822387/144115188075855872)^2))^2/(4*((2058499170265015*$C9^2)/(1125899906842624*$C9^2 - 135107988821114880) + (8688735763930779*$C9^2)/(9007199254740992*$C9^2 - 2800521317822387/16) + 3903691986989077/2251799813685248)^(3/2))
% (1856109722364793/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)) - (9280548611823965*$C9^2)/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)^2) + (1856109722364793*$C9^4)/(281474976710656*($C9^2 - 6150233422068845/576460752303423488)^3) + 5198641238726495/(18014398509481984*($C9^2 - 120)) - (25993206193632475*$C9^2)/(18014398509481984*($C9^2 - 120)^2) + (5198641238726495*$C9^4)/(4503599627370496*($C9^2 - 120)^3))/(2*((5198641238726495*$C9^2)/(36028797018963968*$C9^2 - 4323455642275676160) + (1856109722364793*$C9^2)/(2251799813685248*$C9^2 - 6150233422068845/256) + 1529644179441673/1125899906842624)^(1/2)) - ((1856109722364793*$C9)/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)) + (5198641238726495*$C9)/(18014398509481984*($C9^2 - 120)) - (1856109722364793*$C9^3)/(1125899906842624*($C9^2 - 6150233422068845/576460752303423488)^2) - (5198641238726495*$C9^3)/(18014398509481984*($C9^2 - 120)^2))^2/(4*((5198641238726495*$C9^2)/(36028797018963968*$C9^2 - 4323455642275676160) + (1856109722364793*$C9^2)/(2251799813685248*$C9^2 - 6150233422068845/256) + 1529644179441673/1125899906842624)^(3/2))
% 

%% for ZnSe

n = sqrt(4+1.9*lambda^2/(lambda^2-0.113));

dn_dlambda = diff(n);

dn_shg = eval(subs(dn_dlambda,lambda,0.405));

dn_fond = eval(subs(dn_dlambda,lambda,0.810));

% for LiNbO3

nn10 = sqrt(1 + 2.6734*lambda.^2/(lambda.^2 - 0.01764)+...
1.229*lambda.^2/(lambda.^2 - 0.05914)+...
12.61*lambda.^2/(lambda.^2 - 474.6));

nn10e = sqrt(1 + 2.9804*lambda.^2/(lambda.^2 - 0.02047)+...
0.5981*lambda.^2/(lambda.^2 - 0.0666)+...
8.954*lambda.^2/(lambda.^2 - 416.08));

d_n_e = diff(nn10e);
d_n_o = diff(nn10);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));

dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));

n_e_2w = eval(subs(nn10e,lambda,0.405));
n_e_w = eval(subs(nn10e,lambda,0.81));
n_o_w = eval(subs(nn10,lambda,0.81));
n_o_2w = eval(subs(nn10,lambda,0.405));

% for ADP

nn10 = sqrt( 2.30+ 15.1*lambda.^2/(lambda.^2 - 400)+...
0.011/(lambda.^2 - 0.0133));

nn10e = sqrt( 2.16+ 5.92*lambda.^2/(lambda.^2 - 400)+...
0.0096/(lambda.^2 - 0.013));

d_n_e = diff(nn10e);
d_n_o = diff(nn10);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));

dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));

n_e_2w = eval(subs(nn10e,lambda,0.405));
n_e_w = eval(subs(nn10e,lambda,0.81));
n_o_w = eval(subs(nn10,lambda,0.81));
n_o_2w = eval(subs(nn10,lambda,0.405));

% for LTA

nn10 = sqrt( 4.468 - 0.0179*lambda.^2 +...
0.10*lambda.^(-2)+...
2.85e-3*lambda.^(-4)+...
7.87e-4*lambda.^(-6));

nn10e = sqrt( 4.49 - 0.023*lambda.^2 +...
0.0966*lambda.^(-2)+...
1.05e-3*lambda.^(-4)+...
6.3e-4*lambda.^(-6));

d_n_e = diff(nn10e);
d_n_o = diff(nn10);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));

dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));

n_e_2w = eval(subs(nn10e,lambda,0.405));
n_e_w = eval(subs(nn10e,lambda,0.81));
n_o_w = eval(subs(nn10,lambda,0.81));
n_o_2w = eval(subs(nn10,lambda,0.405));

% for KTP

nn10a = sqrt( 3.29100 + 0.04140/(lambda.^2 - 0.03978) + 9.35522/(lambda.^2 - 31.45571));
nn10b = sqrt( 3.45018 + 0.04340/(lambda.^2 - 0.04597) + 16.98825/(lambda.^2 - 39.43799));
nn10g = sqrt( 4.59423 + 0.06206/(lambda.^2 - 0.04763) + 110.80672/(lambda.^2 - 86.12171));


d_n_e = diff(nn10b);
d_n_o = diff(nn10a);
d_n_g = diff(nn10g);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);
dd_n_g = diff(d_n_g);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));
d_n_g_w = eval(subs(d_n_g,lambda,0.81));
d_n_g_2w = eval(subs(d_n_g,lambda,0.405));

dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));
dd_n_g_w = eval(subs(dd_n_g,lambda,0.81));
dd_n_g_2w = eval(subs(dd_n_g,lambda,0.405));

n_e_2w = eval(subs(nn10b,lambda,0.405));
n_e_w = eval(subs(nn10b,lambda,0.81));
n_o_w = eval(subs(nn10a,lambda,0.81));
n_o_2w = eval(subs(nn10a,lambda,0.405));
n_g_w = eval(subs(nn10g,lambda,0.81));
n_g_2w = eval(subs(nn10g,lambda,0.405));

% for RTP

nn10a = sqrt( 1.68 + 1.43*lambda.^2/(lambda.^2 - 0.0325) -0.0119*lambda.^2);
nn10b = sqrt( 2.04 + 1.09*lambda.^2/(lambda.^2 - 0.044) -0.009*lambda.^2);
nn10g = sqrt( 2.29 + 1.13*lambda.^2/(lambda.^2 - 0.056) -0.0188*lambda.^2);

d_n_e = diff(nn10b);
d_n_o = diff(nn10a);
d_n_g = diff(nn10g);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);
dd_n_g = diff(d_n_g);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));
d_n_g_w = eval(subs(d_n_g,lambda,0.81));
d_n_g_2w = eval(subs(d_n_g,lambda,0.405));

dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));
dd_n_g_w = eval(subs(dd_n_g,lambda,0.81));
dd_n_g_2w = eval(subs(dd_n_g,lambda,0.405));

n_e_2w = eval(subs(nn10b,lambda,0.405));
n_e_w = eval(subs(nn10b,lambda,0.81));
n_o_w = eval(subs(nn10a,lambda,0.81));
n_o_2w = eval(subs(nn10a,lambda,0.405));
n_g_w = eval(subs(nn10g,lambda,0.81));
n_g_2w = eval(subs(nn10g,lambda,0.405));

% for KDP

nno = (2.258+0.0101/(lambda^2-0.0142)+1.762*lambda^2/(lambda^2-57.898))^0.5;
nne = (2.1295+0.0097/(lambda^2-0.014)+0.758*lambda^2/(lambda^2-127.0535))^0.5;


% for KD*P(DKDP)

nn10a = sqrt( 2.24 + 0.0097/(lambda.^2 - 0.0156) + 2.247*lambda.^2/(lambda.^2 - 127));
nn10b = sqrt( 2.126 + 0.0086/(lambda.^2 - 0.0120) + 0.784*lambda.^2/(lambda.^2 - 123.4));

d_n_e = diff(nn10b);
d_n_o = diff(nn10a);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));


dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));


n_e_2w = eval(subs(nn10b,lambda,0.405));
n_e_w = eval(subs(nn10b,lambda,0.81));
n_o_w = eval(subs(nn10a,lambda,0.81));
n_o_2w = eval(subs(nn10a,lambda,0.405));

% for MgF2
syms lambda
nn10 = sqrt(1 + 0.488*lambda.^2/(lambda.^2 - 0.04.^2)+...
0.398*lambda.^2/(lambda.^2 - 0.095.^2)+...
2.31*lambda.^2/(lambda.^2 - 23.8.^2));
nn10e = sqrt(1 + 0.413*lambda.^2/(lambda.^2 - 0.036.^2)+...
0.505*lambda.^2/(lambda.^2 - 0.0907.^2)+...
2.49*lambda.^2/(lambda.^2 - 23.77.^2));

d_n_e = diff(nn10e);
d_n_o = diff(nn10);

dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_e_w = eval(subs(d_n_e,lambda,0.81));
d_n_e_2w = eval(subs(d_n_e,lambda,0.405));
d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));


dd_n_e_w = eval(subs(dd_n_e,lambda,0.81));
dd_n_e_2w = eval(subs(dd_n_e,lambda,0.405));
dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));


n_e_2w = eval(subs(nn10e,lambda,0.405));
n_e_w = eval(subs(nn10e,lambda,0.81));
n_o_w = eval(subs(nn10,lambda,0.81));
n_o_2w = eval(subs(nn10,lambda,0.405));


% for quartz
syms lambda
nn10 = sqrt(1 + 0.696*lambda.^2/(lambda.^2 - 0.068.^2)+...
0.4081*lambda.^2/(lambda.^2 - 0.116.^2)+...
0.897*lambda.^2/(lambda.^2 - 9.892.^2));

% d_n_e = diff(nn10);
d_n_o = diff(nn10);

% dd_n_e = diff(d_n_e);
dd_n_o = diff(d_n_o);

d_n_o_w = eval(subs(d_n_o,lambda,0.81));
d_n_o_2w = eval(subs(d_n_o,lambda,0.405));

dd_n_o_w = eval(subs(dd_n_o,lambda,0.81));
dd_n_o_2w = eval(subs(dd_n_o,lambda,0.405));

% % n_e_2w = eval(subs(nn10,lambda,0.405));
% % n_e_w = eval(subs(nn10,lambda,0.81));
n_o_w = eval(subs(nn10,lambda,0.81));
n_o_2w = eval(subs(nn10,lambda,0.405));

% % for LaK22

% % HLAK67

n2=sqrt(2.7410828-0.016285022*lambda.^2 + 0.015866919*lambda.^(-2) + 0.0010819231*lambda.^(-4) - 7.7131038e-5*lambda.^(-6) + 4.033333e-6*lambda.^(-8));
d_n_o = diff(n2);
dd_n_o = diff(d_n_o);
n_o_w = eval(subs(n2,lambda,0.81));
dn_o_w = eval(subs(d_n_o,lambda,0.81));
ddn_o_w = eval(subs(dd_n_o,lambda,0.81));

% (2,7410828-0,016285022*($I$2/1000)^2 + 0,015866919*($I$2/1000)^(-2) + 0,0010819231*($I$2/1000)^(-4) - 7,7131038e-5*($I$2/1000)^(-6) + 4,033333e-6*($I$2/1000)^(-8))^0,5
% -((4693838016698901*(I$2/1000))/144115188075855872 + 2286664015869371/(72057594037927936*(I$2/1000)^3) + 4989489633283617/(1152921504606846976*(I$2/1000)^5) - 17073798217506795/(36893488147419103232*(I$2/1000)^7) + 4761719143363019/(147573952589676412928*(I$2/1000)^9))/(2*(2286664015869371/(144115188075855872*(I$2/1000)^2) - (4693838016698901*(I$2/1000)^2)/288230376151711744 + 4989489633283617/(4611686018427387904*(I$2/1000)^4) - 5691266072502265/(73786976294838206464*(I$2/1000)^6) + 4761719143363019/(1180591620717411303424*(I$2/1000)^8) + 3086184869167919/1125899906842624)^(1/2))
% (6859992047608113/(72057594037927936*(I$2/1000)^4) + 24947448166418085/(1152921504606846976*(I$2/1000)^6) - 119516587522547565/(36893488147419103232*(I$2/1000)^8) + 42855472290267171/(147573952589676412928*(I$2/1000)^10) - 4693838016698901/144115188075855872)/(2*(2286664015869371/(144115188075855872*(I$2/1000)^2) - (4693838016698901*(I$2/1000)^2)/288230376151711744 + 4989489633283617/(4611686018427387904*(I$2/1000)^4) - 5691266072502265/(73786976294838206464*(I$2/1000)^6) + 4761719143363019/(1180591620717411303424*(I$2/1000)^8) + 3086184869167919/1125899906842624)^(1/2)) - ((4693838016698901*(I$2/1000))/144115188075855872 + 2286664015869371/(72057594037927936*(I$2/1000)^3) + 4989489633283617/(1152921504606846976*(I$2/1000)^5) - 17073798217506795/(36893488147419103232*(I$2/1000)^7) + 4761719143363019/(147573952589676412928*(I$2/1000)^9))^2/(4*(2286664015869371/(144115188075855872*(I$2/1000)^2) - (4693838016698901*(I$2/1000)^2)/288230376151711744 + 4989489633283617/(4611686018427387904*(I$2/1000)^4) - 5691266072502265/(73786976294838206464*(I$2/1000)^6) + 4761719143363019/(1180591620717411303424*(I$2/1000)^8) + 3086184869167919/1125899906842624)^(3/2))
% 

%% EOM for phase modulation

L = 160; % mm
d = 2; % mm, d is the electrode separation
lambda_810 = 0.81; 
lambda_405 = 0.405; 

% LTA
r = 20; % tensor element, = 26.8 pm/V for KD*P, 20pm/V for LTA ; 31.5 pm/V for LiNbO3
n_0_810 = 2.15; 
n_0_405 = 2.277;

% LiNbO3
r = -4.7; % tensor element, = 26.8 pm/V for KD*P, 20pm/V for LTA ; 31.5 pm/V for LiNbO3 (d36 not d33); d_36 = 0.37 pm/V
n_0_810 = 2.25; 
n_0_405 = 2.43;

V = linspace(0, 500, 100);

% 
phase_shift_810 = pi*r*1e-6*V*L/d*(n_0_810^3/lambda_810)/pi; % in frac of pi
phase_shift_405 = pi*r*1e-6*V*L/d*(n_0_810^3/lambda_405)/pi; % in frac of pi

phase_shift = -phase_shift_810 + phase_shift_405; % in frac of pi


figure; plot(V, phase_shift, '+g'); hold on; plot(V, phase_shift_810, 'r'); hold on; plot(V, phase_shift_405, 'b');

%% grating compressor

m = 1;
m_shg = 2;
lambda = 0.8E-3; % mm
lambda_shg = 0.4E-3; % mm
lines_g = 1200; % lines/mm
Lg = 50; % mm, separation between gratings
theta_inc = 20; % deg
c = 3e5; % mm/s

phi2_w = -m^2*lambda^3*Lg/(2*pi*(c*lines_g)^2)*(1-(-m*lambda/lines_g-sin(theta_inc*pi/180))^2)^(-3/2)*10^30; % fs2

phi2_2w = -m_shg^2*lambda_shg^3*Lg/(2*pi*(c*lines_g)^2)*(1-(-m_shg*lambda_shg/lines_g-sin(theta_inc*pi/180))^2)^(-3/2)*10^30; % fs2

diff_shg_fund_phi2 = phi2_2w-phi2_w; % want it to be < -5000 fs2
ratio_shg_fund_phi2 = phi2_2w/phi2_w;% want it to abs > 1

%%

syms theta e

n = 1.51; 
n2 = 1.53;
lambda = 0.810e-3 ; % mm
% e = 1.5 ; % mm

assume(e >0)

phi = @(theta, e) 180/pi*4*pi*e/lambda*(sqrt(n2^2-sin(theta).^2)-sqrt(n^2-sin(theta).^2) - (n2 - n));

dy = @(theta, e) sin(theta).*e.*(1-1./sqrt(n2^2-sin(theta).^2));

eqn = phi(theta, e) == 3*180;

[solx, params, conds] =solve(eqn, theta,'ReturnConditions', true) ;

ee = [0.5, 1, 1.5, 2, 3, 5, 10];
theta_vect = linspace(-25, 25, 500);
 
 str_leg= cell(1, length(ee));
 figure;
 for i =1:length(ee)
    plot (theta_vect, phi(theta_vect*pi/180, ee(i))); hold on;
    str_leg{i} =  num2str(ee(i));
 end
 legend(str_leg);
 
 figure;
plot (theta_vect, dy(theta_vect*pi/180, 1.5))

for i =1:length(ee)
    eqn = phi(theta, ee(i)) == 3*180;

    solx =eval(solve(eqn, theta)) ;
    sol(i) = abs(solx(1))*180/pi;
end

figure;
plot (ee, dy(sol*pi/180, ee))

%% delta_Tg function of lambda BK7

syms lambda

c=3e8;
L=10e-3; % m

% for BK7
n_o = @(lambda) (1+1.04*(lambda./1000).^2./((lambda/1000).^2-0.006)+0.232*(lambda/1000).^2/((lambda/1000).^2-0.02)+1.01*(lambda/1000).^2./((lambda/1000).^2-103.56)).^0.5;

dn_dlambda = @(lambda) diff(n_o(lambda));

% vg = @(lambda) c./(n_o(lambda) - lambda.*dn_dlambda(lambda));
vg =  c./(n_o - lambda.*dn_dlambda);

Delta_Tg = @(lambda) L*(1./vg(lambda)-1./vg(lambda./2));

lambda_vect = linspace(700, 1500, 200); % nm
L_vect = [0.1, 1, 2, 5, 10, 20, 50, 100];

figure;
for L=L_vect %1:100
    Delta_Tg_vect = L*1e-3*(1./eval(subs(vg, lambda_vect))-1./eval(subs(vg, lambda_vect/2)));

    plot(lambda_vect, Delta_Tg_vect*1e12); hold on;
end

%% n THETA

syms c s no(lambda) ne(lambda)

nn = (c/no(lambda)^2+s/ne(lambda)^2)^-0.5;

diff (nn, lambda)


%((2*COS($S$4/180*3,14)^2*G4)/O4^3 + (2*SIN($S$4/180*3,14)^2*H4)/P4^3)/(2*(SIN($S$4/180*3,14)^2/P4^2 + COS($S$4/180*3,14)^2/O4^2)^(3/2))



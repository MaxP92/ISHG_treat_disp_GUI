function [diff_k1, diff_k2, no_conv, eps_fact1, eps_fact2 ] = verify_convergence_STP3_algo( k, m, dx, dy, kx, ky, diff_k1, diff_k2, ...
    epsilon1, epsilon2, no_conv, k_max, delta )
%  verify_convergence_STP3_algo
%
% 2016.7.13 Maxime PINSARD

if k > 1
    diff_k1(m, k+1) = abs( ((dx(m, k) + dy(m, k)) - (dx(1, k) + dy(1, k))) - ((dx(m, k-1) + dy(m, k-1)) - (dx(1, k-1) + dy(1, k-1))) );
    diff_k2(m, k+1) = abs( ((kx(m, k) - kx(1, k)) - (kx(m, k-1) - kx(1, k-1)))) + abs(((ky(m, k) - ky(1, k)) - (ky(m, k-1) - ky(1, k-1))) );
else % k=1
    diff_k1(m, k+1) = abs( dx(m, k) + dy(m, k) - (dx(1, k) + dy(1, k)) )/2- delta(m);
    diff_k2(m, k+1) = abs( kx(m, k) - kx(1, k)) + abs(ky(m, k) - ky(1, k));
end

% diff_k1(m, k+1) = sign(diff_k1(m, k+1))*min([abs(diff_k1(m, k+1)), abs(diff_k1(m, k+1)-360), abs(diff_k1(m, k+1)+360)]);

diff_k1(m, k+1) = sign(diff_k1(m, k+1))*wrapTo360(abs(diff_k1(m, k+1))); % put the value between 0 and 360 deg

if k == k_max
    
    diff_k1(m, end) = diff_k1(m, k+1); % ind = iii-1 but there is a column of zero so ind+2
    diff_k2(m, end) = diff_k2(m, k+1);
end

% diff_k1_end=sign(diff_k1(m, end))*min(abs(diff_k1(m, end)), 180-abs(diff_k1(m, end)), abs(diff_k1(m, end))-180);
diff_k1_end = diff_k1(m, end); if diff_k1_end > 180; diff_k1_end=diff_k1_end-360;end
if diff_k1_end < -180; diff_k1_end=diff_k1_end+360;end
fprintf('Remaining difference tilt gradient = %.5g (degree)\n', diff_k2(m, end));
fprintf('Remaining difference transl. = %.5g (degree)\n', diff_k1_end);
eps_fact1 = epsilon1*180/pi;
eps_fact2 = epsilon2*180/pi;
criterion1 = abs(diff_k1_end) > eps_fact1 || abs(diff_k2(m, end)) > eps_fact2; % because epsilon is in rad

if criterion1 % probing convergence
    
    no_conv(m) = 1; % no convergence
    
else
    no_conv(m) = 0; % convergence for one phase shift
end

end


function [err_polar, err_equator, slope, pass] = verify_dipole(epsilon_0, eps_s, domain_L, verbose)
%VERIFY_DIPOLE  Analytical cross-check for a z-directed electric dipole.
%
%  Three independent checks are performed:
%
%  1. Polar-axis value  (θ = 0, z-axis):
%       E_z(0,0,r) = 2p / (4πε₀ r³)
%
%  2. Equatorial-plane value  (θ = π/2, x-axis):
%       E_z(r,0,0) = -p / (4πε₀ r³),  E_x = E_y = 0
%
%  3. 1/r³ decay fit in the near field via log-log linear regression.
%     Ideal slope = -3.0; tolerance ±0.05 is considered passing.
%
% -----------------------------------------------------------------------
% Inputs
%   epsilon_0  - Vacuum permittivity  [F/m]
%   eps_s      - Softening factor     [m]  (must match simulation value)
%   domain_L   - Domain half-size     [m]  (sets the sampling range)
%   verbose    - true to print table; false for silent (default: true)
%
% Outputs
%   err_polar   - Max relative error on polar axis     (dimensionless)
%   err_equator - Max relative error on equatorial plane
%   slope       - Fitted log-log slope  (should be ≈ -3.0)
%   pass        - Logical: true when all three checks pass
% -----------------------------------------------------------------------

if nargin < 4, verbose = true; end

p_mag = 1.0;                         % normalised dipole moment
p_vec = [0, 0, p_mag];
r_src = [0, 0, 0];
coeff = p_mag / (4*pi*epsilon_0);

%% --- Sampling points: 5% to 90% of domain to stay in near-field ---
N   = 300;
r_s = linspace(0.05*domain_L, 0.90*domain_L, N)';   % column vector

%% --- Check 1: Polar axis (z-axis) ---
Xp = zeros(N,1);  Yp = zeros(N,1);  Zp = r_s;
[~, ~, Ez_num, ~] = electric_dipole_field(Xp, Yp, Zp, p_vec, r_src, eps_s, epsilon_0);
Ez_ana = 2*coeff ./ r_s.^3;
err_polar = max(abs(Ez_num - Ez_ana) ./ abs(Ez_ana));

%% --- Check 2: Equatorial plane (x-axis) ---
Xe = r_s;  Ye = zeros(N,1);  Ze = zeros(N,1);
[Ex_num, Ey_num, Ez_num_eq, ~] = electric_dipole_field(Xe, Ye, Ze, p_vec, r_src, eps_s, epsilon_0);
Ez_ana_eq = -coeff ./ r_s.^3;
Ex_ana_eq = zeros(N,1);

err_equator_Ez = max(abs(Ez_num_eq - Ez_ana_eq) ./ abs(Ez_ana_eq));
err_equator_Ex = max(abs(Ex_num));        % should be zero
err_equator    = max(err_equator_Ez, err_equator_Ex);

%% --- Check 3: log-log slope via linear regression ---
% Sample E_mag along polar axis
[~, ~, ~, E_num] = electric_dipole_field(Xp, Yp, Zp, p_vec, r_src, eps_s, epsilon_0);
log_r = log10(r_s);
log_E = log10(E_num + eps);
% Least-squares fit: log(E) = slope*log(r) + intercept
A_mat = [log_r, ones(N,1)];
coeffs = A_mat \ log_E;
slope  = coeffs(1);

%% --- Pass/fail thresholds ---
tol_field = 1e-3;      % 0.1% field error
tol_slope = 0.05;      % |slope + 3| < 0.05
pass_polar   = err_polar   < tol_field;
pass_equator = err_equator < tol_field;
pass_slope   = abs(slope + 3) < tol_slope;
pass = pass_polar && pass_equator && pass_slope;

%% --- Optional verbose output ---
if verbose
    sep = repmat('─', 1, 52);
    fprintf('\n%s\n', sep);
    fprintf('  偶极子场数值校验报告 (Verification Report)\n');
    fprintf('%s\n', sep);
    fprintf('  极轴 (polar) 相对误差 :  %.3e  %s\n', ...
            err_polar,   status_str(pass_polar));
    fprintf('  赤道 (equator) 相对误差:  %.3e  %s\n', ...
            err_equator, status_str(pass_equator));
    fprintf('  log-log 衰减斜率      :  %.4f  (理论: -3.000)  %s\n', ...
            slope,       status_str(pass_slope));
    fprintf('%s\n', sep);
    if pass
        fprintf('  结论: 全部校验通过 ✓\n');
    else
        fprintf('  结论: 存在不通过项，请检查 eps_s 或网格密度 ✗\n');
    end
    fprintf('%s\n\n', sep);
end
end

% -----------------------------------------------------------------------
function s = status_str(flag)
    if flag, s = '[PASS]'; else, s = '[FAIL]'; end
end

function [Hx, Hy, Hz, H_mag, Bx, By, Bz] = magnetic_dipole_field(X, Y, Z, m_vec, r_src, eps_s, mu_0, mode, k)
%MAGNETIC_DIPOLE_FIELD  Vectorized H and B fields of a magnetic dipole.
%
%  Static:
%    H(r) = (1/4π)   · [3r̂(m·r̂) - m] / r³
%    B(r) = (μ₀/4π)  · [3r̂(m·r̂) - m] / r³
%
%  Time-harmonic (dual of Hertzian electric dipole, e^{jωt}):
%    H(r) = (e^{-jkr}/4π) · { (k²/r)·m⊥  +  [(1+jkr)/r³]·m∥ }
%
%  Structural identity with electric_dipole_field; prefactor differs:
%    Electric E:  1/(4πε₀)  with p
%    Magnetic H:  1/(4π)    with m
%
% -----------------------------------------------------------------------
% Inputs
%   X,Y,Z      - Observation-point arrays [m] (same shape; from meshgrid)
%   m_vec      - [mx, my, mz]  magnetic dipole moment  [A·m²]
%   r_src      - [x0, y0, z0]  source location         [m]
%   eps_s      - Softening length                       [m]
%   mu_0       - Vacuum permeability                    [H/m]
%   mode       - 'static' (default) | 'timeharmonic'
%   k          - Wavenumber [rad/m]  (ignored for 'static')
%
% Outputs
%   Hx,Hy,Hz   - Magnetic field intensity  [A/m]
%   H_mag      - |H|  [A/m]
%   Bx,By,Bz   - Magnetic flux density    [T]  (B = μ₀H)
% -----------------------------------------------------------------------

if nargin < 8, mode = 'static'; end
if nargin < 9, k    = 0;        end

coeff_H = 1.0 / (4.0 * pi);      % H-field prefactor  (no μ₀)
mx = m_vec(1);  my = m_vec(2);  mz = m_vec(3);

%% --- Displacement and stabilised distance ---
dx = X - r_src(1);
dy = Y - r_src(2);
dz = Z - r_src(3);

r2 = dx.^2 + dy.^2 + dz.^2 + eps_s^2;
r  = sqrt(r2);
r3 = r2 .* r;
r5 = r2 .* r2 .* r;

m_dot_d = mx.*dx + my.*dy + mz.*dz;

%% --- H-field calculation (identical structure to electric dipole E) ---
switch lower(mode)

    case 'static'
        Hx = coeff_H .* (3.*dx.*m_dot_d - mx.*r2) ./ r5;
        Hy = coeff_H .* (3.*dy.*m_dot_d - my.*r2) ./ r5;
        Hz = coeff_H .* (3.*dz.*m_dot_d - mz.*r2) ./ r5;

    case 'timeharmonic'
        phase = exp(-1j .* k .* r);
        % e^{jωt} convention; reduces to static at k→0
        A = (1 + 1j.*k.*r) ./ r3;            % reactive/near-field:  (1+jkr)/r³
        B = k^2 ./ r;                        % radiative/far-field:  k²/r

        mt_x = mx - m_dot_d.*dx./r2;
        mt_y = my - m_dot_d.*dy./r2;
        mt_z = mz - m_dot_d.*dz./r2;

        ml_x = 3.*m_dot_d.*dx./r2 - mx;
        ml_y = 3.*m_dot_d.*dy./r2 - my;
        ml_z = 3.*m_dot_d.*dz./r2 - mz;

        Hx = coeff_H .* phase .* (B.*mt_x + A.*ml_x);
        Hy = coeff_H .* phase .* (B.*mt_y + A.*ml_y);
        Hz = coeff_H .* phase .* (B.*mt_z + A.*ml_z);

    otherwise
        error('magnetic_dipole_field: mode must be ''static'' or ''timeharmonic''.');
end

H_mag = sqrt(abs(Hx).^2 + abs(Hy).^2 + abs(Hz).^2);

%% --- B-field via constitutive relation B = μ₀H ---
Bx = mu_0 .* Hx;
By = mu_0 .* Hy;
Bz = mu_0 .* Hz;
end

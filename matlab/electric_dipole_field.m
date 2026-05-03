function [Ex, Ey, Ez, E_mag] = electric_dipole_field(X, Y, Z, p_vec, r_src, eps_s, epsilon_0, mode, k)
%ELECTRIC_DIPOLE_FIELD  Vectorized electric field of a point electric dipole.
%
%  Static near-field:
%    E(r) = (1/4πε₀) · [3r̂(p·r̂) - p] / r³
%
%  Time-harmonic (Hertzian dipole), e^{jωt} convention:
%    E(r) = (e^{-jkr}/4πε₀) · { (k²/r)·p⊥  +  [(1+jkr)/r³]·p∥ }
%    where  p⊥ = p - (p·r̂)r̂   (transverse, radiating)
%           p∥ = 3(p·r̂)r̂ - p  (longitudinal, reactive)
%
%  Singularity at r→0 is regularised: r² → r² + ε_s²
%  All computation is in Cartesian coordinates (no spherical mapping).
%
% -----------------------------------------------------------------------
% Inputs
%   X,Y,Z      - Observation-point coordinate arrays [m] (arbitrary shape,
%                must be the same size; produced by meshgrid / ndgrid)
%   p_vec      - [px, py, pz]  dipole moment vector  [C·m]
%   r_src      - [x0, y0, z0]  source location       [m]
%   eps_s      - Softening length for singularity     [m]
%                (Recommended: ~1e-4 × domain_half_size)
%   epsilon_0  - Vacuum permittivity                  [F/m]
%   mode       - 'static' (default) | 'timeharmonic'
%   k          - Wavenumber [rad/m]  (ignored for 'static')
%
% Outputs
%   Ex,Ey,Ez   - Field components  (real for static; complex for harmonic)
%   E_mag      - Field amplitude   sqrt(|Ex|²+|Ey|²+|Ez|²)
% -----------------------------------------------------------------------

if nargin < 8, mode = 'static'; end
if nargin < 9, k    = 0;        end

coeff = 1.0 / (4.0 * pi * epsilon_0);
px = p_vec(1);  py = p_vec(2);  pz = p_vec(3);

%% --- Displacement vectors and stabilised distance ---
dx = X - r_src(1);
dy = Y - r_src(2);
dz = Z - r_src(3);

r2 = dx.^2 + dy.^2 + dz.^2 + eps_s^2;   % regularised |r|²
r  = sqrt(r2);                             % regularised |r|
r3 = r2 .* r;                             % r^3
r5 = r2 .* r2 .* r;                       % r^5

% p · d  (unnormalised inner product; r factors handled by denominator)
p_dot_d = px.*dx + py.*dy + pz.*dz;

%% --- Field calculation ---
switch lower(mode)

    % ------------------------------------------------------------------
    case 'static'
        % E_i = (coeff/r⁵) · [3·d_i·(p·d) - p_i·r²]
        Ex = coeff .* (3.*dx.*p_dot_d - px.*r2) ./ r5;
        Ey = coeff .* (3.*dy.*p_dot_d - py.*r2) ./ r5;
        Ez = coeff .* (3.*dz.*p_dot_d - pz.*r2) ./ r5;

    % ------------------------------------------------------------------
    case 'timeharmonic'
        phase = exp(-1j .* k .* r);          % retardation phase

        % Weight factors  (e^{jωt} convention; reduces to static at k→0)
        A = (1 + 1j.*k.*r) ./ r3;            % reactive/near-field:  (1+jkr)/r³
        B = k^2 ./ r;                        % radiative/far-field:  k²/r

        % Transverse component of p relative to r̂:  p - (p·d/r²)·d
        pt_x = px - p_dot_d.*dx./r2;
        pt_y = py - p_dot_d.*dy./r2;
        pt_z = pz - p_dot_d.*dz./r2;

        % Longitudinal (aligned with r̂):  3(p·d/r²)·d - p
        pl_x = 3.*p_dot_d.*dx./r2 - px;
        pl_y = 3.*p_dot_d.*dy./r2 - py;
        pl_z = 3.*p_dot_d.*dz./r2 - pz;

        Ex = coeff .* phase .* (B.*pt_x + A.*pl_x);
        Ey = coeff .* phase .* (B.*pt_y + A.*pl_y);
        Ez = coeff .* phase .* (B.*pt_z + A.*pl_z);

    % ------------------------------------------------------------------
    otherwise
        error('electric_dipole_field: mode must be ''static'' or ''timeharmonic''.');
end

%% --- Amplitude ---
E_mag = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
end

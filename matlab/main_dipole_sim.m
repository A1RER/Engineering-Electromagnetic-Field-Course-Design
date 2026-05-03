%% ========================================================================
%  电（磁）偶极子空间场强数值仿真计算
%  Electric / Magnetic Dipole Spatial Field Strength Numerical Simulation
%
%  运行环境: MATLAB R2019b 及以上
%  依赖函数: electric_dipole_field.m  magnetic_dipole_field.m  verify_dipole.m
%
%  ▶ 直接按 F5 / Run 运行，可视化结果将自动弹出
%  ▶ 修改 Section 1 中的参数控制仿真内容
% =========================================================================
clear; clc; close all;

% =========================================================================
%  Module 1 ─ 数理逻辑层：物理常数 & 仿真参数
%  Mathematical Logic Layer
% =========================================================================

%% 1.1  基本物理常数
epsilon_0 = 8.854187817e-12;    % 真空介电常数  [F/m]
mu_0      = 4*pi*1e-7;          % 真空磁导率    [H/m]
c_light   = 1/sqrt(epsilon_0*mu_0);  % 光速       [m/s]

%% 1.2  偶极子类型选择
%   'electric' → 计算 E 场，输入 p [C·m]
%   'magnetic' → 计算 H/B 场，输入 m [A·m²]
DIPOLE_TYPE = 'electric';

%% 1.3  偶极子矩矢量（支持任意方向）
p_vec = [0, 0, 1e-30];   % z 方向电偶极子矩  [C·m]  (量级≈分子偶极矩)
m_vec = [0, 0, 1e-20];   % z 方向磁偶极子矩  [A·m²]

%% 1.4  源点位置
r_src = [0, 0, 0];       % 偶极子置于坐标原点

%% 1.5  时谐场参数（仅 mode='timeharmonic' 时有效）
SIM_MODE = 'static';     % 'static' | 'timeharmonic'
freq     = 1e9;          % 工作频率  [Hz]
omega    = 2*pi*freq;
k_wave   = omega / c_light;   % 波数  [rad/m]
lambda   = c_light / freq;    % 波长  [m]

% =========================================================================
%  Module 2 ─ 空间离散化引擎：计算网格
%  Computational Mesh Engine
% =========================================================================

%% 2.1  计算域参数
L  = 5e-9;      % 计算域半径  [m]  (5 nm，适用于静态近场展示)
N2 = 120;       % 2D 网格每边点数（精度与速度权衡）
N3 = 40;        % 3D 网格每边点数（体渲染，点数不宜过大）

eps_s = L * 1e-4;   % 奇异点软化因子 ε_s [m]（约 0.5 pm，远小于最近采样距离）

%% 2.2  二维切面网格（xz 平面，y = 0）
lin2 = linspace(-L, L, N2);
[X2, Z2] = meshgrid(lin2, lin2);
Y2 = zeros(size(X2));

fprintf('[网格] 2D: %d×%d = %d 点  |  3D: %d³ = %d 点\n', ...
        N2, N2, N2^2, N3, N3^3);
fprintf('[参数] 软化因子 ε_s = %.2e m  |  计算域 ±%.2e m\n', eps_s, L);

% =========================================================================
%  Module 3 ─ 计算内核：矢量化场强计算
%  Computation Kernel  (Vectorized — no nested for-loops over grid points)
% =========================================================================
tic;

if strcmp(DIPOLE_TYPE, 'electric')
    [Fx, Fy, Fz, F_mag] = electric_dipole_field( ...
        X2, Y2, Z2, p_vec, r_src, eps_s, epsilon_0, SIM_MODE, k_wave);
    field_label  = '|E| (V/m)';
    dipole_label = '电偶极子 (Electric Dipole)';

else   % 'magnetic'
    [Fx, Fy, Fz, F_mag, ~, ~, ~] = magnetic_dipole_field( ...
        X2, Y2, Z2, m_vec, r_src, eps_s, mu_0, SIM_MODE, k_wave);
    field_label  = '|H| (A/m)';
    dipole_label = '磁偶极子 (Magnetic Dipole)';
end

t_calc = toc;
fprintf('[计算] %s %s 场，耗时 %.3f s\n', dipole_label, SIM_MODE, t_calc);

% =========================================================================
%  Module 4 ─ 数据后处理与可视化
%  Post-Processing & Visualization
% =========================================================================

%% 4.1  场强模值热力图 + 矢量方向 + 衰减曲线
fig1 = figure('Name', [dipole_label ' ─ 场强分布'], ...
              'Position', [60, 60, 1400, 480], 'Color', 'w');

% ── 子图 A：log 幅值热力图 ──────────────────────────────────────────────
ax1 = subplot(1,3,1);
F_log = log10(F_mag + eps);
contourf(X2*1e9, Z2*1e9, F_log, 48, 'LineStyle', 'none');
cb = colorbar;  cb.Label.String = ['log_{10}(' field_label ')'];
colormap(ax1, 'turbo');
hold on;
plot(0, 0, 'w+', 'MarkerSize', 10, 'LineWidth', 2);  % 源点标记
title('场强幅值 (xz 平面)', 'FontSize', 12);
xlabel('x  (nm)');  ylabel('z  (nm)');
axis equal tight;  grid off;

% ── 子图 B：归一化矢量箭头图 (Quiver) ──────────────────────────────────
ax2 = subplot(1,3,2);
step = max(1, floor(N2/20));        % 降采样密度，避免箭头过密
idx  = 1:step:N2;
Qx = X2(idx,idx)*1e9;  Qz = Z2(idx,idx)*1e9;
QFx = real(Fx(idx,idx));  QFz = real(Fz(idx,idx));
QFm = F_mag(idx,idx) + eps;
quiver(Qx, Qz, QFx./QFm, QFz./QFm, 0.45, ...
       'Color', [0.15 0.45 0.85], 'LineWidth', 0.8);
hold on;
contour(X2*1e9, Z2*1e9, F_log, 12, 'k', 'LineWidth', 0.4);
plot(0, 0, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
title('场矢量方向 + 等幅线', 'FontSize', 12);
xlabel('x  (nm)');  ylabel('z  (nm)');
axis equal tight;  grid on;  set(gca, 'GridAlpha', 0.3);

% ── 子图 C：极轴衰减特性（log-log）──────────────────────────────────────
ax3 = subplot(1,3,3);
r_axis = linspace(L*0.05, L*0.95, 400);
x0 = zeros(size(r_axis));

if strcmp(DIPOLE_TYPE, 'electric')
    [~,~,~, E_axis] = electric_dipole_field( ...
        x0, x0, r_axis, p_vec, r_src, eps_s, epsilon_0);
    % 解析参考：E_z = 2p/(4πε₀r³)
    coeff_ref = 2 * norm(p_vec) / (4*pi*epsilon_0);
    E_ref = coeff_ref ./ r_axis.^3;
else
    [~,~,~, E_axis] = magnetic_dipole_field( ...
        x0, x0, r_axis, m_vec, r_src, eps_s, mu_0);
    coeff_ref = 2 * norm(m_vec) / (4*pi);
    E_ref = coeff_ref ./ r_axis.^3;
end

loglog(r_axis*1e9, E_axis, 'b-',  'LineWidth', 2,   'DisplayName', '数值计算');
hold on;
loglog(r_axis*1e9, E_ref,  'r--', 'LineWidth', 1.5, 'DisplayName', '1/r³ 解析');
legend('Location', 'southwest', 'FontSize', 9);
title('极轴方向场强衰减', 'FontSize', 12);
xlabel('r  (nm)');  ylabel(field_label);
grid on;  set(gca, 'FontSize', 10);

sgtitle([dipole_label '  |  计算模式: ' upper(SIM_MODE)], 'FontSize', 14);

%% 4.2  流线图 (Streamline) — 更直观地展示偶极子场拓扑
fig2 = figure('Name', [dipole_label ' ─ 场线拓扑'], ...
              'Position', [80, 80, 700, 600], 'Color', 'w');

% 流线用归一化场（streamslice 自动归一化）
streamslice(X2*1e9, Z2*1e9, real(Fx), real(Fz), 2.5);
hold on;
contour(X2*1e9, Z2*1e9, F_log, 10, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', '偶极子');
colormap('parula');
title([dipole_label '  —  场线图 (xz 平面)'], 'FontSize', 13);
xlabel('x  (nm)');  ylabel('z  (nm)');
axis equal tight;  grid on;
legend('Location', 'northeast');

%% 4.3  三维体渲染（可选，点数较少时运行）
if N3 <= 50
    lin3 = linspace(-L, L, N3);
    [X3,Y3,Z3] = meshgrid(lin3, lin3, lin3);

    if strcmp(DIPOLE_TYPE, 'electric')
        [~,~,~, Fm3] = electric_dipole_field( ...
            X3, Y3, Z3, p_vec, r_src, eps_s, epsilon_0);
    else
        [~,~,~, Fm3] = magnetic_dipole_field( ...
            X3, Y3, Z3, m_vec, r_src, eps_s, mu_0);
    end

    fig3 = figure('Name', [dipole_label ' ─ 3D 切片渲染'], ...
                  'Position', [100,100,700,600], 'Color', 'w');
    logFm3 = log10(Fm3 + eps);
    slice(X3*1e9, Y3*1e9, Z3*1e9, logFm3, 0, 0, 0);   % 三正交切面
    shading interp;  colorbar;  colormap('turbo');
    alpha(0.7);
    title([dipole_label '  —  三维正交切面场强'], 'FontSize', 13);
    xlabel('x (nm)');  ylabel('y (nm)');  zlabel('z (nm)');
    axis equal;  grid on;  view(35, 25);
end

% =========================================================================
%  Module 5 ─ 校验与收敛性分析
%  Verification & Convergence Analysis
% =========================================================================

if strcmp(DIPOLE_TYPE, 'electric')
    [err_p, err_e, slope, ok] = verify_dipole(epsilon_0, eps_s, L, true);

    % 绘制衰减拟合图
    fig4 = figure('Name', '收敛性分析', 'Position', [120,120,600,420], 'Color','w');
    r_chk = linspace(L*0.05, L*0.90, 200);
    z0 = zeros(size(r_chk));
    [~,~,~, Ev] = electric_dipole_field(z0,z0,r_chk, p_vec, r_src, eps_s, epsilon_0);
    loglog(r_chk*1e9, Ev, 'b-', 'LineWidth', 2, 'DisplayName', '数值结果');
    hold on;
    % 拟合参考线
    r0 = r_chk(50);  E0 = Ev(50);
    E_fit = E0 .* (r0./r_chk).^abs(slope);
    loglog(r_chk*1e9, E_fit, 'r--', 'LineWidth', 1.5, ...
           'DisplayName', sprintf('拟合斜率 = %.3f', slope));
    legend;
    title(sprintf('场强衰减 log-log 拟合  (斜率 = %.3f，理论 = -3.000)', slope));
    xlabel('r (nm)');  ylabel('|E| (V/m)');
    grid on;
end

% =========================================================================
%  数据导出 ─ Origin 兼容 CSV 格式
%  Data Export
% =========================================================================

out_file = fullfile(fileparts(mfilename('fullpath')), '..', ...
                    'dipole_field_data.csv');
fid = fopen(out_file, 'w');
fprintf(fid, 'x_nm,y_nm,z_nm,Fx,Fy,Fz,F_mag\n');
fclose(fid);

out_mat = [X2(:)*1e9, Y2(:)*1e9, Z2(:)*1e9, ...
           real(Fx(:)), real(Fy(:)), real(Fz(:)), F_mag(:)];
writematrix(out_mat, out_file, 'WriteMode', 'append');

fprintf('[导出] 数据已保存至: %s\n', out_file);

% =========================================================================
%  图片导出
% =========================================================================
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'output');
if ~isfolder(out_dir), mkdir(out_dir); end

tag = sprintf('%s_%s', lower(DIPOLE_TYPE(1:3)), lower(SIM_MODE(1:3)));
exportgraphics(fig1, fullfile(out_dir, [tag '_distribution.png']), 'Resolution', 150);
exportgraphics(fig2, fullfile(out_dir, [tag '_streamline.png']),   'Resolution', 150);
if exist('fig3', 'var')
    exportgraphics(fig3, fullfile(out_dir, [tag '_3d_slice.png']), 'Resolution', 150);
end
if exist('fig4', 'var')
    exportgraphics(fig4, fullfile(out_dir, [tag '_convergence.png']), 'Resolution', 150);
end
fprintf('[导出] 图片已保存至: %s\n', out_dir);
fprintf('[完成] 仿真结束 ─ 共 %d 个观测点\n', numel(F_mag));

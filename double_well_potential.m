%% Double-well potential: V(x) and gradient-descent trajectories
%
%   Potential extracted from MSS_figs.mlx (double Gaussian, mu = 0.5):
%
%       V(x) = -d_L * exp(-(x - x_L)^2 / (2*sig^2))
%             - d_R * exp(-(x - x_R)^2 / (2*sig^2))
%             + c * x^2
%
%   ODE (gradient descent on V):
%       dx/dt = -dV/dx
%
%   Equilibria: zeros of dV/dx, classified by sign of d²V/dx²
%     • Stable   (local min, d²V/dx² > 0) – solid red   filled circle
%     • Unstable (local max, d²V/dx² < 0) – solid black filled circle
%
%   Figures produced:
%     Left  panel – V(x) with equilibria marked on the curve
%     Right panel – x(t) trajectories from 12 initial conditions,
%                   coloured by basin of attraction
%
%   Saves: plots/double_well_potential.pdf

%% ── Parameters ──────────────────────────────────────────────────────────

x_L  = -1.0;           % left  Gaussian centre
x_R  =  1.0;           % right Gaussian centre
sig  =  0.55;          % Gaussian width
s2   =  sig^2;         % sig^2 = 0.3025 — appears in all derivative formulas
d_L  =  0.5;           % left  well depth  (mu = 0.5 => d_L = 1 - mu)
d_R  =  0.5;           % right well depth  (mu = 0.5 => d_R = mu)
c    =  0.18;          % quadratic confining coefficient

tspan = [0, 12];       % integration window

%% ── Potential and its derivatives (anonymous functions) ─────────────────

% Gaussian kernel factors
eL = @(x) exp(-(x - x_L).^2 ./ (2*s2));
eR = @(x) exp(-(x - x_R).^2 ./ (2*s2));

% Potential  V(x)
V   = @(x)  -d_L.*eL(x) - d_R.*eR(x) + c.*x.^2;

% First derivative  dV/dx  (chain rule on each Gaussian + quadratic term)
dV  = @(x)  (d_L/s2).*(x - x_L).*eL(x) ...
           + (d_R/s2).*(x - x_R).*eR(x) ...
           + 2*c.*x;

% Second derivative  d²V/dx²  (product rule: d/dx[(x-xi)*exp(...)])
d2V = @(x)  (d_L/s2).*(1 - (x - x_L).^2./s2).*eL(x) ...
           + (d_R/s2).*(1 - (x - x_R).^2./s2).*eR(x) ...
           + 2*c;

% ODE right-hand side: gradient descent
odefun = @(t, x) -dV(x);

%% ── ODE solver options ──────────────────────────────────────────────────

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ── Equilibrium analysis ────────────────────────────────────────────────

solOpts = optimoptions('fsolve', 'Display', 'off', ...
                       'TolFun', 1e-12, 'TolX', 1e-12);

% 15 seeds across the domain to catch all roots of dV/dx = 0
eq_seeds = linspace(-2.2, 2.2, 15);

eq_pts  = [];    % x* values
eq_stab = [];    % 1 = stable (local min), 0 = unstable (local max)

for k = 1 : numel(eq_seeds)
    [xstar, fval, flag] = fsolve(dV, eq_seeds(k), solOpts);
    if flag > 0 && abs(fval) < 1e-8
        % Reject duplicates (within tolerance 1e-4)
        if isempty(eq_pts) || min(abs(eq_pts - xstar)) > 1e-4
            eq_pts(end+1)  = xstar;                       %#ok<AGROW>
            eq_stab(end+1) = double(d2V(xstar) > 0);     %#ok<AGROW>
        end
    end
end

% Sort left to right
[eq_pts, ord] = sort(eq_pts);
eq_stab       = eq_stab(ord);

idx_s = eq_stab == 1;   % stable   mask
idx_u = eq_stab == 0;   % unstable mask

fprintf('\n── Equilibria ──────────────────────────────────────────\n');
for k = 1 : numel(eq_pts)
    if eq_stab(k); tag = 'stable (local min)'; else; tag = 'unstable (local max)'; end
    fprintf('  x* = %+.6f   V(x*) = %+.6f   [%s]\n', ...
            eq_pts(k), V(eq_pts(k)), tag);
end
fprintf('────────────────────────────────────────────────────────\n\n');

%% ── Trajectories ────────────────────────────────────────────────────────

% 12 ICs spanning [-2.2, 2.2] — 6 per basin (linspace with even N avoids x0=0)
IC_vals = linspace(-2.2, 2.2, 12);

% Basin colours (matching MSS_figs.mlx palette)
col_right = [0.5059, 0.6078, 0.4039];   % green  — right well (x0 >= 0)
col_left  = [0.7137, 0.6353, 0.4118];   % brown  — left  well (x0 <  0)

traj_t   = cell(numel(IC_vals), 1);
traj_x   = cell(numel(IC_vals), 1);
traj_col = zeros(numel(IC_vals), 3);

for k = 1 : numel(IC_vals)
    [tt, xx]   = ode45(odefun, tspan, IC_vals(k), odeOpts);
    traj_t{k}  = tt;
    traj_x{k}  = xx;
    if IC_vals(k) >= 0
        traj_col(k,:) = col_right;
    else
        traj_col(k,:) = col_left;
    end
end

%% ── Shared publication style ─────────────────────────────────────────────

FS      = 22;                    % base font size (pt)
LW      = 2.0;                   % trajectory line width
FN      = 'Helvetica';           % axis font
blue    = [0.13 0.47 0.71];      % potential curve colour
red_col = [0.85 0.07 0.07];      % stable equilibrium marker
blk_col = [0    0    0   ];      % unstable equilibrium marker
MS      = 10;                    % marker size (pt)

%% ── Figure ──────────────────────────────────────────────────────────────

fig      = figure('Units', 'centimeters', 'Position', [3 3 26 10], 'Color', 'w');
x_range  = linspace(-2.2, 2.2, 500);

% ─── Left panel: V(x) ────────────────────────────────────────────────────
ax1 = subplot(1, 2, 1);

plot(x_range, V(x_range), '-', 'Color', blue, 'LineWidth', LW);
hold on;

% Stable equilibria — solid red filled circles on the curve
plot(eq_pts(idx_s), V(eq_pts(idx_s)), 'o', ...
     'Color',           red_col, ...
     'MarkerFaceColor', red_col, ...
     'MarkerSize', MS, 'LineWidth', 1.5);

% Unstable equilibria — solid black filled circles on the curve
plot(eq_pts(idx_u), V(eq_pts(idx_u)), 'o', ...
     'Color',           blk_col, ...
     'MarkerFaceColor', blk_col, ...
     'MarkerSize', MS, 'LineWidth', 1.5);

xlabel('$x$',    'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('$V(x)$', 'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN, ...
       'Rotation', 0, 'HorizontalAlignment', 'right');

xlim([-2.2, 2.2]);

set(ax1, 'FontSize', FS, 'FontName', FN, 'LineWidth', 0.8, ...
         'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off');

% ─── Right panel: x(t) trajectories ──────────────────────────────────────
ax2 = subplot(1, 2, 2);
hold on;

for k = 1 : numel(IC_vals)
    plot(traj_t{k}, traj_x{k}, '-', ...
         'Color', traj_col(k,:), 'LineWidth', LW);
end

% Dashed horizontal lines at stable equilibrium positions
x_stab = eq_pts(idx_s);
for k = 1 : numel(x_stab)
    yline(x_stab(k), '--', 'Color', red_col, 'LineWidth', 1.0);
end

xlabel('Time, $t$', 'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('$x(t)$',    'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN, ...
       'Rotation', 0, 'HorizontalAlignment', 'right');

xlim([0,   12  ]);
ylim([-2.4, 2.4]);

set(ax2, 'FontSize', FS, 'FontName', FN, 'LineWidth', 0.8, ...
         'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off');

%% ── Save figure ──────────────────────────────────────────────────────────

if ~exist('plots', 'dir')
    mkdir('plots');
end

exportgraphics(fig, 'plots/double_well_potential.pdf', 'ContentType', 'vector');

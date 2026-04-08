%% 2D bistable ecosystem potential: surface and gradient-descent trajectories
%
%   Phenomenological double-well potential (cusp catastrophe normal form):
%
%       V(x, y) = (x² - 1)² + y²
%
%   This is the normal form of the fold/cusp catastrophe widely used in
%   ecology to model regime shifts and bistability (Zeeman 1976, May 1977,
%   Scheffer et al. 2001).  In the parametric form
%
%       V(x, y; b) = x⁴/4 - b·x²/2 + y²
%
%   the bifurcation parameter b controls the number of stable states:
%     b > 0  (here b = 2)  →  two stable states at x* = ±√b  (bistable)
%     b = 0                →  fold bifurcation (single degenerate state)
%     b < 0                →  one stable state at x* = 0 (monostable)
%
%   x = normalised ecosystem state (e.g. vegetation cover, fish biomass)
%   y = environmental fluctuation (e.g. soil moisture, turbidity)
%
%   Equilibria:
%     Stable (local minima):  (±1, 0),  V = 0
%     Unstable (saddle):       (0,  0),  V = 1
%
%   ODE system (gradient descent on V):
%       dx/dt = -∂V/∂x = 4x(1 - x²)
%       dy/dt = -∂V/∂y = -2y
%
%   Figure:
%     3D surface  V(x, y)  with colour-coded gradient-descent trajectories:
%     green  = trajectories converging to the right stable state  (x* = +1)
%     brown  = trajectories converging to the left  stable state  (x* = -1)
%
%   References:
%     Zeeman (1976) Catastrophe theory. Sci. Am. 234, 65-83.
%     May (1977) Thresholds and breakpoints. Nature 269, 471-477.
%     Scheffer et al. (2001) Catastrophic shifts. Nature 413, 591-596.
%
%   Saves: plots/double_well_potential.pdf

%% ── Parameters ──────────────────────────────────────────────────────────

b = 2;          % bifurcation parameter  (bistable for b > 0)

x_vec = linspace(-2, 2, 300);   % state-variable grid
y_vec = linspace(-2, 2, 300);   % environmental-variable grid

tspan = [0, 5];                  % integration window

%% ── Potential and partial derivatives (anonymous functions) ──────────────

% Potential  V(x, y; b) = x⁴/4 - b·x²/2 + y²  =  (x² - b/2·(x²) + y²)
% Written in terms of b for easy parameter sweeps:
V    = @(x, y)   (x.^2 - 1).^2  +  y.^2;          % b = 2 hardwired

% ∂V/∂x = 4x(x² - 1)
dVdx = @(x, y)   4.*x.*(x.^2 - 1);

% ∂V/∂y = 2y
dVdy = @(x, y)   2.*y;

% 2D ODE right-hand side:  [dx/dt; dy/dt] = -∇V
odefun = @(t, u) [ -dVdx(u(1), u(2)); ...
                   -dVdy(u(1), u(2)) ];

%% ── ODE solver options ───────────────────────────────────────────────────

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ── Equilibria ───────────────────────────────────────────────────────────
% Analytical: ∇V = 0 → 4x(x²-1) = 0 AND 2y = 0
%   Stable   (local minima, ∂²V/∂x² > 0): (±1, 0), V = 0
%   Unstable (saddle,       ∂²V/∂x² < 0): ( 0, 0), V = 1

eq_stable   = [ 1, 0;  -1, 0 ];    % [x, y] rows
eq_unstable = [ 0, 0           ];

fprintf('\n── Equilibria ──────────────────────────────────────────────\n');
for k = 1 : size(eq_stable, 1)
    fprintf('  Stable   (x*, y*) = (%+.1f, %+.1f)   V = %.4f\n', ...
            eq_stable(k,1), eq_stable(k,2), V(eq_stable(k,1), eq_stable(k,2)));
end
for k = 1 : size(eq_unstable, 1)
    fprintf('  Unstable (x*, y*) = (%+.1f, %+.1f)   V = %.4f\n', ...
            eq_unstable(k,1), eq_unstable(k,2), V(eq_unstable(k,1), eq_unstable(k,2)));
end
fprintf('────────────────────────────────────────────────────────────\n\n');

%% ── Trajectories ─────────────────────────────────────────────────────────

% 6 × 4 grid of ICs — linspace with even counts avoids x0 = 0 exactly
x0_vals = linspace(-1.8, 1.8, 6);   % 6 values: -1.8 -1.08 -0.36 +0.36 +1.08 +1.8
y0_vals = linspace(-1.5, 1.5, 4);   % 4 values: -1.5 -0.5 +0.5 +1.5

% Basin colours (matching MSS_figs.mlx palette)
col_right = [0.5059, 0.6078, 0.4039];   % green  — right basin (x0 > 0)
col_left  = [0.7137, 0.6353, 0.4118];   % brown  — left  basin (x0 < 0)

traj_all = {};

for ix = 1 : numel(x0_vals)
    for iy = 1 : numel(y0_vals)
        u0 = [x0_vals(ix); y0_vals(iy)];
        [~, uu] = ode45(odefun, tspan, u0, odeOpts);

        tr.x   = uu(:, 1);
        tr.y   = uu(:, 2);
        tr.V   = V(uu(:,1), uu(:,2));
        if x0_vals(ix) >= 0
            tr.col = col_right;
        else
            tr.col = col_left;
        end
        traj_all{end+1} = tr;    %#ok<AGROW>
    end
end

%% ── Shared publication style ─────────────────────────────────────────────

FS  = 22;             % base font size (pt)
LW  = 2.0;            % trajectory line width
FN  = 'Helvetica';    % axis font

red_col = [0.85 0.07 0.07];    % stable equilibrium marker
blk_col = [0    0    0   ];    % unstable equilibrium marker
MS  = 10;                      % marker size (pt)

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'centimeters', 'Position', [3 3 16 13], 'Color', 'w');
ax  = axes('Parent', fig);

% ─ 3D potential surface ───────────────────────────────────────────────────
[X_grid, Y_grid] = meshgrid(x_vec, y_vec);
V_surf = V(X_grid, Y_grid);

surf(X_grid, Y_grid, V_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.88);
shading interp
colormap(ax, flipud(parula));
hold on;

% ─ Gradient-descent trajectories ─────────────────────────────────────────
lift = 0.06;   % vertical offset so curves sit visibly above the surface

for k = 1 : numel(traj_all)
    tr = traj_all{k};
    plot3(tr.x, tr.y, tr.V + lift, '-', ...
          'Color', tr.col, 'LineWidth', LW);
    % Mark initial condition
    plot3(tr.x(1), tr.y(1), tr.V(1) + lift, 'o', ...
          'Color', tr.col, 'MarkerFaceColor', tr.col, ...
          'MarkerSize', 5, 'LineWidth', 1.0);
end

% ─ Equilibria on the surface ─────────────────────────────────────────────
% Stable — solid red filled circles at V = 0
plot3(eq_stable(:,1), eq_stable(:,2), V(eq_stable(:,1), eq_stable(:,2)) + lift, 'o', ...
      'Color', red_col, 'MarkerFaceColor', red_col, ...
      'MarkerSize', MS, 'LineWidth', 1.5);

% Unstable — solid black filled circles at V = 1
plot3(eq_unstable(:,1), eq_unstable(:,2), V(eq_unstable(:,1), eq_unstable(:,2)) + lift, 'o', ...
      'Color', blk_col, 'MarkerFaceColor', blk_col, ...
      'MarkerSize', MS, 'LineWidth', 1.5);

% ─ Axes labels and formatting ─────────────────────────────────────────────
xlabel('Ecosystem state, $x$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('Environmental variable, $y$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
zlabel('Potential, $V$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);

view([-40, 25]);
xlim([-2,  2]);
ylim([-2,  2]);

ax.Box       = 'on';
ax.FontSize  = FS;
ax.FontName  = FN;
ax.LineWidth = 0.8;
ax.GridAlpha = 0.2;

%% ── Save figure ──────────────────────────────────────────────────────────

if ~exist('plots', 'dir')
    mkdir('plots');
end

exportgraphics(fig, 'plots/double_well_potential.pdf', 'ContentType', 'vector');

%% 2D bistable ecosystem potential: surface and gradient-descent trajectories
%
%   Phenomenological double-well potential (cusp catastrophe normal form):
%
%       V(x, y; b) = (x² - b)² / (4b)  +  c · y²
%
%   This is the normal form of the fold/cusp catastrophe widely used in
%   ecology to model regime shifts and bistability (Zeeman 1976, May 1977,
%   Scheffer et al. 2001).  The bifurcation parameter b controls the
%   number and position of stable states:
%     b > 0  →  two stable states at x* = ±√b  (bistable)
%     b = 0  →  fold bifurcation (one degenerate equilibrium)
%     b < 0  →  one stable state at x* = 0 (monostable)
%
%   Here b = 4, c = 0.5:
%       V(x, y) = (x² - 4)² / 16  +  y² / 2
%
%   x = normalised ecosystem state (e.g. vegetation cover, fish biomass)
%   y = environmental fluctuation  (e.g. soil moisture, turbidity)
%
%   Equilibria:
%     Stable (local minima):  (±2, 0),  V = 0
%     Unstable (saddle):       ( 0, 0),  V = 1
%
%   ODE system (gradient descent on V):
%       dx/dt = -∂V/∂x =  x - x³/4   [≡ x(1 - x²/4)]
%       dy/dt = -∂V/∂y = -y
%
%   References:
%     Zeeman (1976) Catastrophe theory. Sci. Am. 234, 65-83.
%     May (1977) Thresholds and breakpoints. Nature 269, 471-477.
%     Scheffer et al. (2001) Catastrophic shifts. Nature 413, 591-596.
%
%   Saves: plots/double_well_potential.pdf

%% ── Parameters ──────────────────────────────────────────────────────────

b = 4;      % bifurcation parameter — stable states at x* = ±√b = ±2
c = 0.5;    % y-curvature coefficient

x_vec = linspace(-2.6, 2.6, 300);   % state-variable grid
y_vec = linspace(-1.6, 1.6, 300);   % environmental-variable grid

tspan = [0, 6];                      % integration window

%% ── Potential and partial derivatives (anonymous functions) ──────────────

% V(x,y) = (x²-b)²/(4b) + c·y²
V    = @(x, y)   (x.^2 - b).^2 ./ (4*b)  +  c .* y.^2;

% ∂V/∂x = x(x²-b)/b  [= x³/b - x]
dVdx = @(x, y)   x .* (x.^2 - b) ./ b;

% ∂V/∂y = 2c·y
dVdy = @(x, y)   2*c .* y;

% 2D ODE:  [dx/dt; dy/dt] = -∇V
odefun = @(t, u) [ -dVdx(u(1), u(2)); ...
                   -dVdy(u(1), u(2)) ];

%% ── ODE solver options ───────────────────────────────────────────────────

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ── Equilibria ───────────────────────────────────────────────────────────
% Analytical: ∇V = 0  →  x(x²-b)/b = 0  AND  2c·y = 0
%   Stable   (local minima, ∂²V/∂x²>0): (±√b, 0) = (±2, 0),  V = 0
%   Unstable (saddle,       ∂²V/∂x²<0): (  0,  0),            V = 1

eq_stable   = [ sqrt(b), 0 ; -sqrt(b), 0 ];   % [x, y] rows
eq_unstable = [ 0,       0 ];

fprintf('\n── Equilibria ──────────────────────────────────────────────\n');
for k = 1 : size(eq_stable, 1)
    fprintf('  Stable   (x*, y*) = (%+.4f, %+.4f)   V = %.4f\n', ...
            eq_stable(k,1), eq_stable(k,2), ...
            V(eq_stable(k,1), eq_stable(k,2)));
end
for k = 1 : size(eq_unstable, 1)
    fprintf('  Unstable (x*, y*) = (%+.4f, %+.4f)   V = %.4f\n', ...
            eq_unstable(k,1), eq_unstable(k,2), ...
            V(eq_unstable(k,1), eq_unstable(k,2)));
end
fprintf('────────────────────────────────────────────────────────────\n\n');

%% ── Trajectories ─────────────────────────────────────────────────────────

% 6 × 4 grid of ICs — even counts keep x0 = 0 out of the grid
x0_vals = linspace(-2.4, 2.4, 6);    % six ICs spanning both basins
y0_vals = linspace(-1.2, 1.2, 4);    % four ICs in y

% Basin colours (matching MSS_figs.mlx palette)
col_right = [0.5059, 0.6078, 0.4039];   % green  — right basin (x0 > 0)
col_left  = [0.7137, 0.6353, 0.4118];   % brown  — left  basin (x0 < 0)

traj_all = {};

for ix = 1 : numel(x0_vals)
    for iy = 1 : numel(y0_vals)
        u0       = [x0_vals(ix); y0_vals(iy)];
        [~, uu]  = ode45(odefun, tspan, u0, odeOpts);

        tr.x   = uu(:, 1);
        tr.y   = uu(:, 2);
        tr.V   = V(uu(:,1), uu(:,2));
        tr.col = col_right * (x0_vals(ix) >= 0) + col_left * (x0_vals(ix) < 0);
        traj_all{end+1} = tr;    %#ok<AGROW>
    end
end

%% ── Shared publication style ─────────────────────────────────────────────

FS  = 22;             % base font size (pt)
LW  = 2.0;            % trajectory line width
FN  = 'Helvetica';    % axis font

red_col = [0.85 0.07 0.07];   % stable equilibrium marker
blk_col = [0    0    0   ];   % unstable equilibrium marker
MS  = 10;                     % marker size (pt)

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'centimeters', 'Position', [3 3 18 13], 'Color', 'w');
ax  = axes('Parent', fig);

% ─ 3D potential surface ───────────────────────────────────────────────────
[X_grid, Y_grid] = meshgrid(x_vec, y_vec);
V_surf = V(X_grid, Y_grid);

surf(X_grid, Y_grid, V_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.88);
shading interp
colormap(ax, flipud(parula));
hold on;

% ─ Trajectories with directional arrows ──────────────────────────────────
lift      = 0.04;    % vertical offset so curves sit above the surface
arrow_len = 0.28;    % arrow length in data units (uniform across all curves)

for k = 1 : numel(traj_all)
    tr = traj_all{k};
    n  = length(tr.x);

    % Trajectory line
    plot3(tr.x, tr.y, tr.V + lift, '-', ...
          'Color', tr.col, 'LineWidth', LW);

    % IC marker
    plot3(tr.x(1), tr.y(1), tr.V(1) + lift, 'o', ...
          'Color', tr.col, 'MarkerFaceColor', tr.col, ...
          'MarkerSize', 5, 'LineWidth', 1.0);

    % Directional arrow — placed at ~35 % along the trajectory
    if n < 6; continue; end
    idx = max(3, round(n * 0.35));
    i0  = max(1,   idx - 2);
    i1  = min(n,   idx + 2);

    dx_a = tr.x(i1) - tr.x(i0);
    dy_a = tr.y(i1) - tr.y(i0);
    dV_a = tr.V(i1) - tr.V(i0);
    mag  = sqrt(dx_a^2 + dy_a^2 + dV_a^2);
    if mag < 1e-10; continue; end

    u = dx_a / mag * arrow_len;
    v = dy_a / mag * arrow_len;
    w = dV_a / mag * arrow_len;

    quiver3(tr.x(idx), tr.y(idx), tr.V(idx) + lift, u, v, w, 0, ...
            'Color', tr.col, 'LineWidth', LW, 'MaxHeadSize', 2.5);
end

% ─ Equilibria on the surface ─────────────────────────────────────────────
plot3(eq_stable(:,1),   eq_stable(:,2),   ...
      V(eq_stable(:,1),   eq_stable(:,2))   + lift, 'o', ...
      'Color', red_col, 'MarkerFaceColor', red_col, ...
      'MarkerSize', MS, 'LineWidth', 1.5);

plot3(eq_unstable(:,1), eq_unstable(:,2), ...
      V(eq_unstable(:,1), eq_unstable(:,2)) + lift, 'o', ...
      'Color', blk_col, 'MarkerFaceColor', blk_col, ...
      'MarkerSize', MS, 'LineWidth', 1.5);

% ─ Axes labels ───────────────────────────────────────────────────────────
xlabel('Ecosystem state, $x$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('Environmental variable, $y$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
zlabel('Potential, $V$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);

% ─ View and axis scaling ──────────────────────────────────────────────────
view([-40, 25]);
xlim([-2.6,  2.6]);
ylim([-1.6,  1.6]);
zlim([ 0,    2.0]);     % clip high edges; saddle at V=1, wells at V=0

% Flatten the z-axis visually so the double-well shape reads clearly
pbaspect([2.6, 1.6, 0.9]);

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

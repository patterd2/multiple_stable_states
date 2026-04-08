%% Double-well potential: 2D surface and gradient-descent trajectories
%
%   Potential extracted from MSS_figs.mlx (double Gaussian):
%
%       V(x,μ) = -(1-μ) · exp(-(x-x_L)²/(2σ²))
%               -    μ  · exp(-(x-x_R)²/(2σ²))
%               +  c · x²
%
%   where x  = ecological state   (x ∈ [-2, 2])
%         μ  = bifurcation param  (μ ∈ [ 0, 1])
%
%   ODE system (gradient descent on V):
%       dx/dt = -∂V/∂x
%       dμ/dt = -∂V/∂μ
%
%   Integration terminates when μ reaches 0 or 1 (ODE event function).
%
%   Figure:
%     3D surface  V(x,μ)  with colour-coded gradient-descent trajectories
%     overlaid — green for trajectories from x₀ > 0 (right basin),
%     brown for x₀ < 0 (left basin).
%
%   Saves: plots/double_well_potential.pdf

%% ── Parameters ──────────────────────────────────────────────────────────

x_L = -1.0;          % left  Gaussian centre
x_R =  1.0;          % right Gaussian centre
sig =  0.55;         % Gaussian width
s2  =  sig^2;        % s² = 0.3025 — appears in all derivative formulas
c   =  0.18;         % quadratic confining coefficient

x_vec  = linspace(-2, 2, 300);    % ecological-state grid
mu_vec = linspace( 0, 1, 300);    % bifurcation-parameter grid

tspan = [0, 8];      % integration window (event function caps μ ∈ [0,1])

%% ── Potential and partial derivatives (anonymous functions) ──────────────

% Gaussian kernel factors (functions of x only — μ enters via depth coefficients)
eL = @(x, ~)   exp(-(x - x_L).^2 ./ (2*s2));
eR = @(x, ~)   exp(-(x - x_R).^2 ./ (2*s2));

% Depth coefficients
dL = @(mu)     1 - mu;      % left  well depth
dR = @(mu)     mu;          % right well depth

% Potential  V(x, μ)
V  = @(x, mu)  -dL(mu).*eL(x,mu) - dR(mu).*eR(x,mu) + c.*x.^2;

% ∂V/∂x  (chain rule on each Gaussian + quadratic)
dVdx = @(x, mu)  (dL(mu)./s2).*(x - x_L).*eL(x,mu) ...
                + (dR(mu)./s2).*(x - x_R).*eR(x,mu) ...
                + 2*c.*x;

% ∂V/∂μ  = d/dμ[-(1-μ)·eL - μ·eR] = eL - eR
dVdmu = @(x, mu)  eL(x,mu) - eR(x,mu);

% 2D ODE: gradient descent  [dx/dt; dμ/dt] = -∇V
odefun = @(t, u) [ -dVdx(u(1), u(2)); ...
                   -dVdmu(u(1), u(2)) ];

%% ── ODE solver options ───────────────────────────────────────────────────

% Event function terminates integration when μ reaches 0 or 1
odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10, 'Events', @mu_bounds);

%% ── Surface mesh ─────────────────────────────────────────────────────────

% meshgrid convention matches MSS_figs.mlx Fig 1:
%   [X_grid, MU_grid] = meshgrid(x_vec, mu_vec)
%   → surf(MU_grid, X_grid, V_surf): μ on x-axis, x on y-axis
[X_grid, MU_grid] = meshgrid(x_vec, mu_vec);
V_surf = V(X_grid, MU_grid);

%% ── Trajectories ─────────────────────────────────────────────────────────

% 4 × 4 grid of initial conditions (x₀, μ₀) in the domain interior
x0_vals  = linspace(-1.5, 1.5, 4);
mu0_vals = linspace(0.15, 0.85, 4);

% Basin colours (matching MSS_figs.mlx palette)
col_right = [0.5059, 0.6078, 0.4039];   % green  — right basin (x₀ ≥ 0)
col_left  = [0.7137, 0.6353, 0.4118];   % brown  — left  basin (x₀ <  0)

traj_all = {};   % cell array of trajectory structs

for ix = 1 : numel(x0_vals)
    for im = 1 : numel(mu0_vals)
        u0 = [x0_vals(ix); mu0_vals(im)];
        [tt, uu] = ode45(odefun, tspan, u0, odeOpts);

        tr.t   = tt;
        tr.x   = uu(:, 1);
        tr.mu  = uu(:, 2);
        tr.V   = V(uu(:,1), uu(:,2));
        if x0_vals(ix) >= 0
            tr.col = col_right;
        else
            tr.col = col_left;
        end
        traj_all{end+1} = tr;             %#ok<AGROW>
    end
end

%% ── Shared publication style ─────────────────────────────────────────────

FS  = 22;              % base font size (pt)
LW  = 2.0;             % trajectory line width
FN  = 'Helvetica';     % axis font

%% ── Figure ───────────────────────────────────────────────────────────────

fig = figure('Units', 'centimeters', 'Position', [3 3 16 13], 'Color', 'w');
ax  = axes('Parent', fig);

% ─ 3D potential surface (style from MSS_figs.mlx Figure 1) ───────────────
surf(MU_grid, X_grid, V_surf, 'EdgeColor', 'none', 'FaceAlpha', 0.9);
shading interp
colormap(ax, flipud(parula));
hold on;

% ─ Gradient-descent trajectories overlaid on the surface ─────────────────
lift = 0.05;    % vertical offset so curves sit visibly above the surface

for k = 1 : numel(traj_all)
    tr = traj_all{k};

    % Trajectory curve
    plot3(tr.mu, tr.x, tr.V + lift, '-', ...
          'Color', tr.col, 'LineWidth', LW);

    % Mark initial condition with a small filled circle
    plot3(tr.mu(1), tr.x(1), tr.V(1) + lift, 'o', ...
          'Color',           tr.col, ...
          'MarkerFaceColor', tr.col, ...
          'MarkerSize', 5, 'LineWidth', 1.0);
end

% ─ Axes labels (matching MSS_figs.mlx Figure 1) ──────────────────────────
xlabel('Bifurcation parameter $\mu$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('Ecological state $x$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
zlabel('Potential $V$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);

% ─ View and axis formatting ───────────────────────────────────────────────
view([-40, 24]);
xlim([0,  1]);
ylim([-2, 2]);

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

%% ── Local functions ──────────────────────────────────────────────────────

function [val, isterm, dir] = mu_bounds(~, u)
%MU_BOUNDS  ODE event: terminate when μ reaches 0 or 1.
    val    = [u(2);   1 - u(2)];   % zero when μ = 0 or μ = 1
    isterm = [1;      1       ];   % stop integration at both events
    dir    = [-1;     1       ];   % μ approaching 0 from above, 1 from below
end

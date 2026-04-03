%% Two-species competition ODE system – publication-quality plots
%
%   dN1/dt = N1*(1 - N1 - c1*N2) + gamma1
%   dN2/dt = N2*(1 - c2*N1 - N2) + gamma2
%
%  Figures produced:
%    (1) Time courses for a single initial condition
%    (2) Phase plane: normalised vector field, trajectories, and fixed points
%        • Stable fixed points   – solid red  filled circle
%        • Unstable fixed points – solid black filled circle

%% ── Parameters ──────────────────────────────────────────────────────────
gamma1 = 0.1;
gamma2 = 0.1;
c1     = 2;     % competition coefficient: effect of N2 on N1
c2     = 2;     % competition coefficient: effect of N1 on N2

tspan = [0 50];

%% ── ODE right-hand side ─────────────────────────────────────────────────
odefun = @(t, N) [ ...
    N(1).*(1 - N(1) - c1.*N(2)) + gamma1; ...
    N(2).*(1 - c2.*N(1) - N(2)) + gamma2 ];

odeOpts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

%% ── Shared publication style ─────────────────────────────────────────────
FS   = 22;            % base font size (pt)
LW   = 2.0;           % trajectory line width
FN   = 'Helvetica';   % axis font  ('Times New Roman' for serif)
blue = [0.13 0.47 0.71];

%% ── Fixed-point analysis ─────────────────────────────────────────────────
F = @(N) [ N(1).*(1 - N(1) - c1.*N(2)) + gamma1; ...
           N(2).*(1 - c2.*N(1) - N(2)) + gamma2 ];

J = @(N) [ 1 - 2*N(1) - c1*N(2),   -c1*N(1); ...
           -c2*N(2),                 1 - c2*N(1) - 2*N(2) ];

nSeed    = 12;
sv       = linspace(0, 1.2, nSeed);
[S1, S2] = meshgrid(sv, sv);
seeds    = [S1(:), S2(:)];

fp_pts  = [];
fp_stab = [];

solOpts = optimoptions('fsolve', 'Display', 'off', ...
                       'TolFun', 1e-12, 'TolX', 1e-12);

for k = 1 : size(seeds, 1)
    [xstar, fval, flag] = fsolve(F, seeds(k,:).', solOpts);
    if flag > 0 && norm(fval) < 1e-8
        if isempty(fp_pts) || ...
                min(sqrt(sum((fp_pts - xstar.').^2, 2))) > 1e-4
            fp_pts(end+1, :) = xstar.';                     %#ok<AGROW>
            ev = eig(J(xstar));
            if all(real(ev) < 0)
                fp_stab(end+1) = 1;
            else
                fp_stab(end+1) = 0;
            end
        end
    end
end

idx_s = fp_stab == 1;
idx_u = fp_stab == 0;

fprintf('\n── Fixed points ────────────────────────────────\n');
for k = 1 : size(fp_pts, 1)
    if fp_stab(k); tag = 'stable'; else; tag = 'unstable'; end
    fprintf('  (N1*, N2*) = (%.5f, %.5f)  [%s]\n', ...
            fp_pts(k,1), fp_pts(k,2), tag);
end
fprintf('────────────────────────────────────────────────\n\n');

%% ── Figure 1: Time courses ───────────────────────────────────────────────
IC_tc = [0.2, 0.1];   % initial condition for time-course plot
[t, Nsol] = ode45(odefun, tspan, IC_tc.', odeOpts);

fig1 = figure('Units', 'centimeters', 'Position', [3 3 14 9], 'Color', 'w');

plot(t, Nsol(:,1), '-',  'Color', [0.13 0.47 0.71], 'LineWidth', 1.5);
hold on;
plot(t, Nsol(:,2), '--', 'Color', [0.84 0.15 0.16], 'LineWidth', 1.5);

xlabel('Time, $t$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('Population density, $N_i(t)$', ...
       'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
title('Time courses — initial condition 1', ...
      'FontSize', FS, 'FontWeight', 'bold', 'FontName', FN);
legend({'$N_1$', '$N_2$'}, ...
       'Interpreter', 'latex', 'FontSize', FS - 2, ...
       'Location', 'best', 'Box', 'off');

set(gca, 'FontSize', FS, 'FontName', FN, 'LineWidth', 0.8, ...
         'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off');
grid on;

%% ── Figure 2: Phase plane ────────────────────────────────────────────────
N1_lim = [0, 1.2];
N2_lim = [0, 1.2];
nGrid  = 22;

[N1g, N2g] = meshgrid( ...
    linspace(N1_lim(1), N1_lim(2), nGrid), ...
    linspace(N2_lim(1), N2_lim(2), nGrid));

dN1 = N1g.*(1 - N1g - c1.*N2g) + gamma1;
dN2 = N2g.*(1 - c2.*N1g - N2g) + gamma2;
mag = sqrt(dN1.^2 + dN2.^2);
mag(mag == 0) = 1;

fig2 = figure('Units', 'centimeters', 'Position', [3 3 14 13], 'Color', 'w');

% Normalised vector field
quiver(N1g, N2g, dN1./mag, dN2./mag, 0.5, ...
       'Color', [0.70 0.70 0.70], 'LineWidth', 0.6);
hold on;

% Grid of initial conditions covering the phase space
nSide = 7;
sv_ic = linspace(0.05, 1.15, nSide);
[IC1, IC2] = meshgrid(sv_ic, sv_ic);
ICs  = [IC1(:), IC2(:)];
nIC  = size(ICs, 1);

% Solve and plot all trajectories
for k = 1 : nIC
    N0 = ICs(k,:).';
    [~, Ntraj] = ode45(odefun, tspan, N0, odeOpts);
    plot(Ntraj(:,1), Ntraj(:,2), '-', 'Color', blue, 'LineWidth', LW);
end

% Stable fixed points – solid red filled circle
if any(idx_s)
    plot(fp_pts(idx_s, 1), fp_pts(idx_s, 2), 'o', ...
         'Color',           [0.85 0.07 0.07], ...
         'MarkerFaceColor', [0.85 0.07 0.07], ...
         'MarkerSize', 10, 'LineWidth', 1.5);
end

% Unstable fixed points – solid black filled circle
if any(idx_u)
    plot(fp_pts(idx_u, 1), fp_pts(idx_u, 2), 'o', ...
         'Color',           [0 0 0], ...
         'MarkerFaceColor', [0 0 0], ...
         'MarkerSize', 10, 'LineWidth', 1.5);
end

% Axes formatting – no title, no legend, no box
xlabel('$N_1$', 'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN);
ylabel('$N_2$', 'Interpreter', 'latex', 'FontSize', FS, 'FontName', FN,...
      'Rotation', 0);
axis([N1_lim(1) N1_lim(2) N2_lim(1) N2_lim(2)]);

set(gca, 'FontSize', FS, 'FontName', FN, 'LineWidth', 0.8, ...
         'TickDir', 'out', 'TickLength', [0.015 0.015], 'Box', 'off');


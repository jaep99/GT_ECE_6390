function d

mu = 1.327124400e11; % Sun's GM
dt = 21600; % time step [s]

% Initial values
AU = 149597870.7; % 1 au [km]
a_au = 17.83; % semi-major axis [au]
e_use = 0.967; % eccentricity
a_km = a_au * AU;

r = a_km * (1 - e_use); % perihelion distance [km]
theta = 0; % start angle
V_r0 = 0.0;
V_theta0 = -sqrt( mu*(1+e_use) / (a_km*(1 - e_use)) );

% Reset delta
delta_r = V_r0 * dt; % [km]
delta_theta = (V_theta0 / r) * dt; % [rad]

% Simulation Time & History
v2 = V_r0^2 + V_theta0^2; % speed^2 [km^2/s^2]
eps = 0.5*v2 - mu/r; % energy
h = r * V_theta0; % angular momentum
e = sqrt(1 + (2*eps*h^2)/(mu^2));

if eps < 0 && e < 1
    a = -mu/(2*eps);
    T = 2*pi*sqrt(a^3/mu);
    T_yr = T / (365.25*86400);
    fprintf('Elliptic orbit: a=%.3f km, e=%.5f, T=%.2f yr\n', a, e, T_yr);
else
    a = NaN; T = NaN;
    fprintf('Open orbit: e=%.5f (no period)\n', e);
end

orbits_to_draw = 1;        
if ~isnan(T)
    sim_time = orbits_to_draw * T;
else
    sim_time = 60000;
end
n_steps = ceil(sim_time / dt);

r_history = zeros(n_steps+1, 1);
theta_history = zeros(n_steps+1, 1);
r_history(1) = r;
theta_history(1) = theta;

for n = 1:n_steps
    % Position update
    r_next = r + delta_r;
    theta_next = theta + delta_theta;

    % Delta update
    r_mid = r + 0.5*delta_r;
    delta_r_next = delta_r + (r_mid*(delta_theta^2) - (mu/(r^2))*dt^2);
    delta_theta_next = delta_theta - (2*delta_r*delta_theta / r_mid);

    % Status update
    r = r_next;
    theta = theta_next;
    delta_r = delta_r_next;
    delta_theta = delta_theta_next;

    % Save
    r_history(n+1) = r;
    theta_history(n+1) = theta;
end

x = r_history.*cos(theta_history);
y = r_history.*sin(theta_history);

figure('Color', 'w');
plot(x, y, 'LineWidth', 1.4);
hold on
plot(0, 0, 'ko', 'MarkerFaceColor','k', 'MarkerSize', 7);
text(0, 0, ' Sun', 'VerticalAlignment','bottom', 'FontSize', 9);
axis equal; grid on;
xlabel('x [km]'); ylabel('y [km]');

ax = gca;

if ~isnan(T)
    infoLines = { ...
        sprintf('$e=%.5f$', e), ...
        sprintf('$T=%.2f\\,\\mathrm{yr}$', T_yr), ...
        sprintf('$a=%.2f\\,\\mathrm{au}$', a/AU) };
else
    infoLines = { ...
        sprintf('$e=%.5f$', e), ...
        'Open orbit (no period)'};
end

text(ax, 0.02, 0.98, infoLines, ...
    'Units','normalized', ...
    'Interpreter','latex', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'BackgroundColor','w', ...
    'EdgeColor','k', ...
    'Margin',6, ...
    'FontName','Times', 'FontSize',10);

end

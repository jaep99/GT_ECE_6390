function result = minimum_azimuth_lookangles()
% Choose Δv_t so that the 1-year mean azimuth drift (slope) ≈ 0.
% Choose Δλ0 to minimize the 1-year (solar-pressure) harmonic amplitude.
% Run a high-fidelity simulation and plot 1-year look angles/metrics.
%
% Output:
% Console: Δv_t, Δλ0, slope, annual amplitude, detrended/annual-removed Max/RMS
% Figure: (1) Raw azimuth, (2) Linear detrended, (3) After removing annual term
% result struct: dv_t_mps, dlambda0_deg, metrics

%% Constants
mu = 398600.4418; % μ=GM [km^3/s^2]
r_earth = 6378.137; % Earth radius [km]
w_earth = 2*pi/86164.0905; % Earth sidereal rotation [rad/s]
sidereal_year = 365.25636*86400; % Sidereal year [s]
sidereal_days = 365.25636; % Sidereal year [day]
w_yr = 2*pi/sidereal_year; % Sun-direction rotation [rad/s]
alpha_r = 0.5; % Effective reflectivity (given)
rho_A = 0.5; % Area density [kg/m^2] (fixed)
r0 = 42164.17; % GEO radius [km]

% Station (Van Leer, Atlanta) — latitude, longitude [rad]
station_lat = deg2rad(33.7758);
station_lon = deg2rad(-84.39738);

%% Find Δv_t such that mean drift slope ≈ 0
dt_A = 120; % [s] 2 min (tighter slope estimation)

slope_fn = @(dv_t) slope_deg_per_day( ...
    dv_t, 0.0, dt_A, mu, r_earth, w_earth, sidereal_year, w_yr, ...
    alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days);

% Grid to bracket a sign change for fzero
cand = linspace(-1, 1, 81); % Δv_t in [-1, +1] m/s
svals = arrayfun(slope_fn, cand);
ix = find(sign(svals(1:end-1)).*sign(svals(2:end)) <= 0, 1, 'first');

if ~isempty(ix)
    bracket = [cand(ix), cand(ix+1)];
else
    cand2 = linspace(-5, 5, 81);
    svals2 = arrayfun(slope_fn, cand2);
    ix2 = find(sign(svals2(1:end-1)).*sign(svals2(2:end)) <= 0, 1, 'first');
    if ~isempty(ix2)
        bracket = [cand2(ix2), cand2(ix2+1)];
    else
        bracket = [-5, 5];
    end
end

dv_t_star = fzero(slope_fn, bracket); % [m/s]

%% Minimize 1-year amplitude by setting Δλ0
dt_B = dt_A;
amp_fn = @(d_lambda0) annual_amp_deg( ...
    dv_t_star, d_lambda0, dt_B, mu, r_earth, w_earth, sidereal_year, w_yr, ...
    alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days);

% Search range ±15 deg
d_lambda0_star = fminbnd(amp_fn, deg2rad(-15), deg2rad(15));

%% Simulation & plots
dt_F = 60; % [s]

[t_hist, az_deg_raw, az_deg_detr_lin, resid_noannual_deg, metrics] = ...
    sim_and_metrics( ...
        dv_t_star, d_lambda0_star, dt_F, mu, r_earth, w_earth, ...
        sidereal_year, w_yr, alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days);

fprintf('\nMin-azimuth design (rho_A = %.3f kg/m^2)\n', rho_A);
fprintf('Δv_t (slope≈0) = %+7.3f m/s\n', dv_t_star);
fprintf('Δλ0 (min annual amp) = %+7.3f deg\n', rad2deg(d_lambda0_star));
fprintf('Mean drift slope = %10.6f deg/day\n', metrics.slope_deg_per_day);
fprintf('Annual amp (Az) = %7.3f deg (peak-to-peak ≈ %7.3f deg)\n', ...
        metrics.annual_amp_deg, 2*metrics.annual_amp_deg);
fprintf('Max |Az detrended| = %7.3f deg; RMS = %7.3f deg\n', ...
        metrics.max_abs_detrended_deg, metrics.rms_detrended_deg);
fprintf('After removing annual: Max = %7.3f deg; RMS = %7.3f deg\n', ...
        metrics.max_abs_after_annual_deg, metrics.rms_after_annual_deg);

% Plots
days = t_hist/86400;
figure('Name','Az over sidereal year','Color','w');
subplot(3,1,1);
plot(days, az_deg_raw, 'LineWidth', 1); grid on;
xlabel('Days'); ylabel('Az [deg]');
title('Raw azimuth');

subplot(3,1,2);
plot(days, az_deg_detr_lin, 'LineWidth', 1); grid on;
xlabel('Days'); ylabel('Az - linear [deg]');
title(sprintf('Linear detrended (slope = %.4f deg/day)', metrics.slope_deg_per_day));

subplot(3,1,3);
plot(days, resid_noannual_deg, 'LineWidth', 1); grid on;
xlabel('Days'); ylabel('Residual [deg]');
title(sprintf('After removing annual term (amp = %.3f deg)', metrics.annual_amp_deg));

result = struct('dv_t_mps', dv_t_star, ...
                'dlambda0_deg', rad2deg(d_lambda0_star), ...
                'metrics', metrics);
end

%%
function slope = slope_deg_per_day(dv_t_mps, d_lambda0, dt, ...
        mu, r_earth, w_earth, T, w_yr, alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days)
% Estimate mean azimuth drift slope (deg/day) via 4-term regression:
% [1, t(day), sin(2πt/1y), cos(2πt/1y)] to separate the annual component.

[t_hist, az_deg] = forward_az( ...
    dv_t_mps, d_lambda0, dt, mu, r_earth, w_earth, T, w_yr, ...
    alpha_r, rho_A, r0, station_lat, station_lon);

tt_day = t_hist/86400;
y_deg = unwrap(deg2rad(az_deg))*180/pi; % unwrap in rad, then convert to deg
X = [ones(size(tt_day)), tt_day, ...
          sin(2*pi*tt_day/ sidereal_days), cos(2*pi*tt_day/ sidereal_days)];
b = X\y_deg;
slope = b(2); % deg/day
end

function A_deg = annual_amp_deg(dv_t_mps, d_lambda0, dt, ...
        mu, r_earth, w_earth, T, w_yr, alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days)
% Return annual amplitude (deg) from the sin/cos coefficients of the 4-term fit.

[t_hist, az_deg] = forward_az( ...
    dv_t_mps, d_lambda0, dt, mu, r_earth, w_earth, T, w_yr, ...
    alpha_r, rho_A, r0, station_lat, station_lon);

tt_day = t_hist/86400;
y_deg = unwrap(deg2rad(az_deg))*180/pi;
X = [ones(size(tt_day)), tt_day, ...
          sin(2*pi*tt_day/ sidereal_days), cos(2*pi*tt_day/ sidereal_days)];
b = X\y_deg;
A_deg = hypot(b(3), b(4));
end

function [t_hist, az_deg, az_deg_detr_lin, resid_noannual_deg, m] = ...
    sim_and_metrics(dv_t_mps, d_lambda0, dt, ...
        mu, r_earth, w_earth, T, w_yr, alpha_r, rho_A, r0, station_lat, station_lon, sidereal_days)
% Final simulation and metric computation.

[t_hist, az_deg] = forward_az( ...
    dv_t_mps, d_lambda0, dt, mu, r_earth, w_earth, T, w_yr, ...
    alpha_r, rho_A, r0, station_lat, station_lon);

tt_day = t_hist/86400;
y_deg  = unwrap(deg2rad(az_deg))*180/pi;

% Remove linear trend
P = polyfit(tt_day, y_deg, 1);
trend = polyval(P, tt_day);
az_deg_detr_lin = y_deg - trend;

% Remove annual term (sin/cos at 1/year)
X = [sin(2*pi*tt_day/ sidereal_days), cos(2*pi*tt_day/ sidereal_days)];
ab = X\az_deg_detr_lin;
annual = X*ab;
resid = az_deg_detr_lin - annual;

m = struct();
m.slope_deg_per_day = P(1);
m.annual_amp_deg = hypot(ab(1), ab(2));
m.max_abs_detrended_deg = max(abs(az_deg_detr_lin));
m.rms_detrended_deg = sqrt(mean(az_deg_detr_lin.^2));
m.max_abs_after_annual_deg = max(abs(resid));
m.rms_after_annual_deg = sqrt(mean(resid.^2));

resid_noannual_deg = resid;
end

function [t_hist, az_deg] = forward_az( ...
        dv_t_mps, d_lambda0, dt, mu, r_earth, w_earth, T, w_yr, ...
        alpha_r, rho_A, r0, station_lat, station_lon)
% Earth-centered initial (ECI)
% Earth-centered - Earth-fixed (ECEF)
% Planar ECI -> ECEF -> East North Up -> azimuth time-series over T seconds.
% East North Up (ENU)

% Propagate in Cartesian with RK4 (planar, solar pressure rotating with sun)
[s_hist, t_hist] = propagate_rk4( ...
    init_state(dv_t_mps, d_lambda0, r0, mu, station_lon), ...
    dt, T, mu, alpha_r, rho_A, w_yr);

% ECI -> ECEF
x = s_hist(:,1);  y = s_hist(:,2);
theta_e = w_earth * t_hist;
c = cos(theta_e); s = sin(theta_e);
xe = c.*x + s.*y;
ye = -s.*x + c.*y;
ze = zeros(size(xe)); % planar z = 0

% Station ECEF
[xs, ys, zs] = station_ecef(r_earth, station_lat, station_lon);

% Line-of-sight in ECEF
dx = xe - xs; dy = ye - ys; dz = ze - zs;

% ECEF -> East North Up
[E, N, U] = ecef_to_enu(dx, dy, dz, station_lat, station_lon); %#ok<ASGLU>

% Look angles
azimuth = atan2(E, N);            % [rad], 0..2π
azimuth = mod(azimuth, 2*pi);

az_deg = rad2deg(azimuth);
end

function s0 = init_state(dv_t_mps, d_lambda0, r0, mu, station_lon)
% Initial planar ECI state from Δv_t and Δλ0.

v_circ = sqrt(mu/r0); % [km/s]
v_t = v_circ + dv_t_mps/1000; % [km/s]
theta0 = station_lon + d_lambda0; % start longitude offset

x0 = r0*cos(theta0);
y0 = r0*sin(theta0);
vx0 = -v_t*sin(theta0);
vy0 =  v_t*cos(theta0);

s0 = [x0; y0; vx0; vy0];
end

function [S, t_hist] = propagate_rk4(s0, dt, T, mu, alpha_r, rho_A, w_yr)
% 4th-order Runge–Kutta propagation (planar 2D + rotating solar-pressure accel).

n = floor(T/dt);
t_hist = (0:n)'*dt;

S = zeros(n+1, 4);
S(1,:) = s0.';
a_solar = 9.08e-9 * alpha_r / rho_A; % [km/s^2] (from 9.08e-6 m/s^2)

for k = 1:n
    tk = t_hist(k);
    s = S(k,:).';

    k1 = f(tk,           s);
    k2 = f(tk + 0.5*dt,  s + 0.5*dt*k1);
    k3 = f(tk + 0.5*dt,  s + 0.5*dt*k2);
    k4 = f(tk + dt,      s + dt*k3);

    S(k+1,:) = (s + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)).';
end

    function ds = f(ti, s)
        x=s(1); y=s(2); vx=s(3); vy=s(4);
        r = hypot(x,y);
        aG = -mu/r^3 * [x; y]; % gravity
        alpha = w_yr * ti; % sun direction angle
        aS = a_solar * [cos(alpha); sin(alpha)]; % solar-pressure accel
        ds = [vx; vy; aG(1)+aS(1); aG(2)+aS(2)];
    end
end

function [xs, ys, zs] = station_ecef(r_earth, station_lat, station_lon)
% Spherical Earth station to ECEF.
xs = r_earth*cos(station_lat)*cos(station_lon);
ys = r_earth*cos(station_lat)*sin(station_lon);
zs = r_earth*sin(station_lat);
end

function [E, N, U] = ecef_to_enu(dx, dy, dz, station_lat, station_lon)
% ECEF delta -> local East North Up.
sl = sin(station_lon); 
cl = cos(station_lon);
sp = sin(station_lat); 
cp = cos(station_lat);
E = -sl.*dx +  cl.*dy;
N = -sp.*cl.*dx - sp.*sl.*dy + cp.*dz;
U = cp.*cl.*dx + cp.*sl.*dy + sp.*dz;
end

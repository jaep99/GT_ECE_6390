function minimum_area_density()
% Finds the minimum area density rho_A that yields a target azimuth drift
% (1°, 5°, 15°) over one sidereal year as seen from Atlanta (Van Leer).

%% Constants
mu = 398600.4418; % μ=GM [km^3/s^2]
r_earth = 6378.137; % [km]
w_earth = 2*pi/86164.0905; % Earth sidereal rotation [rad/s]
sidereal_year = 365.25636*86400; % Sidereal year [s]
w_yr = 2*pi/sidereal_year; % Sun-direction rotation [rad/s]
alpha_r = 0.5; % Effective reflectivity (Given in problem)

% Station (Van Leer, Atlanta) latitute, longitude [rad]
station_lat = deg2rad(33.7758);
station_lon = deg2rad(-84.39738);

% GEO initial conditions (circular, equatorial)
r0 = (mu / w_earth^2)^(1/3); % [km]
vtheta0 = r0 * w_earth;
theta0 = 0.0;
lon0 = station_lon;
vr0 = 0.0;

% Integrator
% Test for 1, 0.5
dt = 0.5; % [s]

%% Solve for rho_A at three drift targets
targets = [1, 5, 15]; % degrees
rho_out = zeros(size(targets));
for k = 1:numel(targets)
    rho_out(k) = rho_for_target_deg(targets(k));
end

fprintf('\nMinimum area density over 1 sidereal year\n');
for k = 1:numel(targets)
    fprintf('Target drift %2d° : rho_A = %.3f kg/m^2\n', targets(k), rho_out(k));
end

%%
    function rho = rho_for_target_deg(target_deg)
        % Find rho_A such that |ΔAz| ≈ target_deg after one sidereal year.
        % Drift decreases as rho_A increases.

        % Initial bracket [rho_lo, rho_hi] where drift_lo > target > drift_hi
        
        rho_lo = 0.5; 
        rho_hi = 200; % broad bracket [kg/m^2]

        drift_lo = drift_for_rho(rho_lo);
        drift_hi = drift_for_rho(rho_hi);

        it = 0;
        while ~(drift_lo>target_deg && drift_hi<target_deg) && it < 50
            if drift_lo <= target_deg, rho_lo = rho_lo/2; drift_lo = drift_for_rho(rho_lo); end
            if drift_hi >= target_deg, rho_hi = rho_hi*2; drift_hi = drift_for_rho(rho_hi); end
            it = it + 1;
        end
        % Bisection
        for i=1:50
            mid = 0.5*(rho_lo+rho_hi);
            d = drift_for_rho(mid);
            if d > target_deg
                rho_lo = mid;
            else
                rho_hi = mid;
            end
        end
        rho = 0.5*(rho_lo+rho_hi);
    end

    function drift_deg = drift_for_rho(rho_A)
        [r_hist,theta_hist,tt] = propagate(r0,theta0,vr0,vtheta0, ...
            mu,alpha_r,rho_A,w_yr,dt,sidereal_year);
        lsh = wrapToPi( (theta_hist - w_earth*tt) + lon0 ); % ECI->ECEF subpoint longitude
        azimuth_0 = az_from_station(station_lat,station_lon,0,lsh(1),r_earth);
        azimuth_f = az_from_station(station_lat,station_lon,0,lsh(end),r_earth);
        drift_deg = abs(rad2deg(angdiff(azimuth_0,azimuth_f)));
    end
end

%% Propagator with solar-pressure
function [r_hist,theta_hist,t_hist] = propagate(r0,theta0,vr0,vtheta0,...
    mu,alpha_r,rho_A,w_yr,dt,simTime)

n = ceil(simTime/dt);
r_hist = zeros(n+1,1); 
t_hist = (0:n)'*dt;
r_hist(1) = r0;
theta_hist = zeros(n+1,1); 
theta_hist(1) = theta0;

% Delta variables store Δ over one step
delta_r = vr0*dt;
delta_theta = (vtheta0/r0)*dt;

% Solar acceleration magnitude [km/s^2]
a_solar = (9.08e-6*alpha_r/rho_A)/1000;  % m/s^2 to km/s^2

r = r0; theta = theta0;
for k = 1:n
    t = (k-1)*dt;

    % Position update
    r_next = r + delta_r;
    theta_next = theta + delta_theta;

    % Mid-step radius
    r_mid = r + 0.5*delta_r;
    theta_mid = theta + 0.5*delta_theta;

    % Sun direction and components
    alpha = w_yr*(t + 0.5*dt); % Sun angle (equatorial)
    psi = alpha - theta_mid; % Relative angle
    a_r = a_solar*cos(psi);
    a_theta = a_solar*sin(psi);

    % Modified updates
    dr_next  = delta_r + ( r_mid*(delta_theta^2) - (mu/(r_mid^2))*dt^2 + a_r*dt^2 );
    dth_next = delta_theta - (2*delta_r*delta_theta)/r_mid + (a_theta/r_mid)*dt^2;

    % Roll
    r = r_next;
    theta = theta_next;
    delta_r = dr_next;
    delta_theta = dth_next;

    r_hist(k+1) = r;
    theta_hist(k+1) = theta;
end
end

%% Azimuth from station (station_lat, station_lon) to sub-satellite point (target_lat, target_lon)
function azimuth = az_from_station(station_lat, station_lon, target_lat, target_lon, earth_radius)

    % Spherical law of cosines for central angle between station and target
    cos_gamma = sin(target_lat)*sin(station_lat) + ...
                cos(target_lat)*cos(station_lat)*cos(target_lon - station_lon);
    central_angle = acos( min(1, max(-1, cos_gamma)) );

    % Magnitude of the azimuth offset (an acute angle in [0, π/2])
    sin_azimuth_abs = abs( sin(abs(station_lon - target_lon)) * cos(target_lat) / ...
                           max(1e-12, sin(central_angle)) );
    azimuth_offset = asin( min(1, max(0, sin_azimuth_abs)) ); % in [0, π/2]

    % Quadrant resolution (for GEO with target_lat = 0 and station_lat > 0, target is to the south)
    if target_lon > station_lon % South & East
        azimuth = pi - azimuth_offset;
    else % South & West
        azimuth = pi + azimuth_offset;
    end

    azimuth = mod(azimuth, 2*pi);
end

%%
function x = wrapToPi(a)
    x = mod(a+pi, 2*pi) - pi;
end

function d = angdiff(a,b)
% minimal signed difference b-a in [-pi,pi]
    d = atan2(sin(b-a), cos(b-a));
end

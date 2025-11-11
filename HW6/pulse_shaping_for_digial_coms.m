clear; clc; close all;

%% --------------------- Global parameters ----------------------------
fc       = 2.45e9;            % Carrier frequency at ISM center (Hz)
f_lo     = 2.40e9;            % ISM lower edge (Hz)
f_hi     = 2.50e9;            % ISM upper edge (Hz)
mask_dB  = -50;               % Required OOB attenuation (dB)
spanSym  = 4;                 % Time support in symbols (must be ≤ 4)
Fs       = 10e9;              % Simulation sampling rate (Hz) -> Nyquist > 2*2.5GHz
dt       = 1/Fs;              % Time step
alphas   = 0.15:0.05:0.90;    % Roll-off search grid (avoid alpha=0 to prevent singularities)
Rb_step  = 1e6;               % 1 Mbps search step
padSym   = 64;                % Zero-padding window length in symbols for FFT resolution
rng(1);

% Helper: One-sided ISM bandwidth (Hz)
B_ISM = (f_hi - f_lo)/2;      % = 50 MHz

%% --------------------- Search best (alpha, Rb) -----------------------
best.Rb     = 0;
best.alpha  = NaN;
best.margin = -Inf;  % (dB) = min attenuation outside band, relative to peak in-band

fprintf('Searching for maximum Rb with span=%d symbols and mask %d dB ...\n', ...
        spanSym, mask_dB);

for alpha = alphas
    % Theoretical RRC bandwidth (one-sided) is B = (1+alpha)*Rb/2.
    % To live INSIDE 2.40–2.50 GHz around fc, require B <= B_ISM.
    Rb_th = floor((2*B_ISM)/(1+alpha)); % upper bound from band limits (Hz)

    % Start at theoretical bound; step down until mask is met.
    for Rb = Rb_th:-Rb_step:1e6
        Tb = 1/Rb;

        % Build 4-Tb RRC pulse at baseband (unit-energy)
        sps = max(8, round(Fs*Tb));           % samples/symbol (ensure at least 8)
        [p, t_p] = rrc_pulse(alpha, spanSym, sps, Tb);
        p = p / sqrt(sum(p.^2));              % unit-energy normalization

        % Place pulse in a longer zero-padded window for fine FFT
        % Total window length (seconds):
        T_win = (spanSym + padSym)*Tb;
        N_win = 2^nextpow2( max(numel(p), ceil(T_win*Fs)) );
        t = (-N_win/2:N_win/2-1).' * dt;      % symmetric time axis
        x = zeros(N_win,1);
        % Center the 4Tb pulse in the window
        idx0 = floor(N_win/2) - floor(numel(p)/2) + (1:numel(p));
        x(idx0) = p;

        % AM modulation to fc: p_M(t) = p(t)*cos(2π f_c t)
        x_mod = x .* cos(2*pi*fc*t);

        % FFT (scaled and shifted per handout)
        [f, Xmod] = spec_fft(x_mod, dt);

        % Measure in-band peak and OOB maximum (positive frequencies only)
        pos = f >= 0;
        inb = pos & (f >= f_lo) & (f <= f_hi);
        oob = pos & ~inb;  % everything nonnegative frequency outside ISM

        peak_inband = max(abs(Xmod(inb)) + eps);
        peak_oob    = max(abs(Xmod(oob)) + eps);
        att_dB      = 20*log10(peak_oob/peak_inband);  % must be ≤ mask_dB

        % If mask is met, record best
        if att_dB <= mask_dB
            if Rb > best.Rb || (Rb == best.Rb && att_dB < best.margin)
                best.Rb     = Rb;
                best.alpha  = alpha;
                best.margin = att_dB;
                best.f      = f;
                best.Xmod   = Xmod;
                best.t      = t;
                best.x      = x;
                best.x_mod  = x_mod;
                best.p      = p;
                best.t_p    = t_p;
                best.Tb     = Tb;
                best.sps    = sps;
            end
            break; % for this alpha we cannot do better than current Rb (descending)
        end
    end

    fprintf(' alpha=%.2f -> feasible Rb up to %.1f Mbps (att margin %.1f dB)\n', ...
            alpha, best.Rb/1e6, best.margin);
end

if best.Rb <= 0
    error('No (alpha, Rb) combination met the -50 dB mask with span=%d.', spanSym);
end

%% --------------------- Report results --------------------------------
fprintf('\Best design \n');
fprintf('  Roll-off alpha  : %.2f\n', best.alpha);
fprintf('  Bit rate Rb     : %.3f Mb/s\n', best.Rb/1e6);
fprintf('  Symbol period Tb: %.3f ns\n', best.Tb*1e9);
B_theory = (1+best.alpha)*best.Rb/2;  % theoretical one-sided baseband bandwidth
fprintf('  Theoretical baseband B = (1+alpha)*Rb/2 = %.1f MHz\n', B_theory/1e6);
fprintf('  OOB attenuation : %.1f dB (≤ %d dB required) --> PASS\n', best.margin, mask_dB);

%% --------------------- Plots -----------------------------------------
% Baseband spectrum |P(f)| needs Xbb, f_bb first
x_bb = zeros(size(best.x));
idx0 = floor(numel(best.x)/2) - floor(numel(best.p)/2) + (1:numel(best.p));
x_bb(idx0) = best.p;
[f_bb, Xbb] = spec_fft(x_bb, dt);

% In-band peak for mask line
peak_inband = max(abs(best.Xmod(best.f>=f_lo & best.f<=f_hi)));

% One window with 3 vertically stacked subplots
fig = figure('Name','Result','Position',[100 100 960 900]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% (1) Time-domain baseband pulse p(t)
nexttile;
plot(best.t_p*1e9, best.p, 'LineWidth', 1.4); grid on;
xlabel('Time (ns)'); ylabel('Amplitude');
title(sprintf('4T_b RRC pulse (\\alpha=%.2f, T_b=%.2f ns, R_b=%.2f Mb/s)', ...
      best.alpha, best.Tb*1e9, best.Rb/1e6));

% (2) Baseband spectrum |P(f)| around DC
nexttile;
semilogy(f_bb/1e6, abs(Xbb)+eps, 'LineWidth', 1.2); grid on;
xlim([0 200]); % show up to 200 MHz
xlabel('Frequency (MHz)'); ylabel('|P(f)|');
title('|P(f)| of 4T_b pulse (semilog-y)');

% (3) Passband spectrum |P_M(f)| with ISM band and -50 dB mask
nexttile;
semilogy(best.f/1e9, abs(best.Xmod)+eps, 'LineWidth', 1.2); hold on; grid on;
yl = ylim;
% Shade ISM band [2.40, 2.50] GHz
patch([f_lo f_hi f_hi f_lo]/1e9, [yl(1) yl(1) yl(2) yl(2)], [0.95 0.95 0.95], ...
      'EdgeColor','none','FaceAlpha',0.35);
% -50 dB line referenced to in-band peak
plot([min(best.f) max(best.f)]/1e9, peak_inband*10^(mask_dB/20)*[1 1], 'k--', 'LineWidth', 1.0);

xlabel('Frequency (GHz)'); ylabel('|P_M(f)|');
title(sprintf('|P_M(f)| with ISM band & %d dB OOB mask', mask_dB));
legend({'|P_M(f)|','ISM band','-50 dB level'}, 'Location','southoutside');

% Title
sgtitle(sprintf('R_b=%.2f Mb/s, \\alpha=%.2f, span=%dT_b', ...
        best.Rb/1e6, best.alpha, spanSym));


function [p, t] = rrc_pulse(alpha, spanSym, sps, Tb)

    T = Tb;
    % Symmetric time axis covering 'spanSym' symbols
    t = (-spanSym/2 : 1/sps : spanSym/2).' * T;
    t = t(:);
    p = zeros(size(t));
    tn = t/T; % normalized time

    % Indices for special cases
    iz0  = abs(t) < 1e-15;
    iz1  = (abs(abs(tn) - 1/(4*alpha)) < 1e-8);

    % General formula (from standard RRC definition)
    num  = sin(pi*(1-alpha)*tn) + 4*alpha*tn.*cos(pi*(1+alpha)*tn);
    den  = pi*tn.*(1 - (4*alpha*tn).^2);
    p(~(iz0|iz1)) = (1/T) * (num(~(iz0|iz1))./den(~(iz0|iz1)));

    % t = 0 limit
    p(iz0) = (1/T) * (1 + alpha*(4/pi - 1));

    % t = ±T/(4α) limit
    if any(iz1)
        % Widely used closed-form limit at t = ±T/(4α)
        p(iz1) = (1/T) * (alpha/sqrt(2)) * ...
                 ( (1 + 2/pi)*sin(pi/(4*alpha)) + (1 - 2/pi)*cos(pi/(4*alpha)) );
    end

    % Normalize to unit energy
    p = p / sqrt(sum(p.^2));
end

function [f, X] = spec_fft(x, dt)
    N = numel(x);
    X = fftshift(fft(x)) * dt; % scale by dt for CTFT units
    df = 1/(N*dt);
    f = (-N/2:N/2-1).' * df;
end

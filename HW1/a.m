function minimum_area_density

% --- Step 1: 초기 설정 및 정지궤도(GEO) 조건 적용 ---
mu = 398600.4418; % 지구 중력 변수 [km^3/s^2]
dt = 600; % 시간 간격 [s], 1년은 길기 때문에 간격을 늘려줍니다 (예: 10분)

% 정지궤도(GEO) 초기 조건 (이전 답변에서 계산한 값)
r0 = 42164.14;       % 초기 거리 [km]
theta0 = -84.39738 * (pi/180); % 초기 각도 [rad], 애틀랜타 경도
vr0 = 0;             % 초기 레이디얼 속도 [km/s]
omega_geo = 7.2921e-5; % 지구 자전 각속도 [rad/s]
v_theta0 = r0 * omega_geo; % 초기 접선 속도 [km/s]

% 시뮬레이션 시간 설정 (1년)
one_year_in_seconds = 365.25 * 24 * 3600;
sim_time = one_year_in_seconds;
n_steps = ceil(sim_time / dt);

% 변수 초기화
r = r0;
theta = theta0;
delta_r = vr0 * dt;
delta_theta = (v_theta0 / r) * dt;

% 결과 저장을 위한 history 배열 초기화
r_history = zeros(n_steps+1, 1);
theta_history = zeros(n_steps+1, 1);
r_history(1) = r;
theta_history(1) = theta;


% --- Step 2: 태양압 모델 변수 추가 ---
rho_A = 20; % 위성의 면적 밀도 [kg/m^2] (과제의 첫 질문을 풀기 위해 이 값을 바꿔가며 테스트)
alpha_r = 0.5; % 유효 반사율 (가정)
solar_const = 9.08e-9; % 태양압 상수 [kN/m^2] 또는 [km*kg/m^2/s^2]

% 가속도 계산 상수 (단위: km/s^2)
% a_solar = (solar_const * alpha_r) / rho_A; 이 계산은 루프 안에서 수행합니다.
% 1 kN/kg = 1 km/s^2 이므로 단위 변환이 간단합니다.

% 1년에 대한 각속도 (태양의 공전 효과)
omega_yr = 2*pi / one_year_in_seconds; % [rad/s]

% --- Step 3: 메인 루프 수정 ---
for n = 1:n_steps
    % 현재 시간 계산
    t = n * dt;

    % --- 태양압 계산 시작 ---
    a_solar = (solar_const * alpha_r) / rho_A; % 태양압 가속도 크기 [km/s^2]
    
    alpha_sun = omega_yr * t; % 현재 시간에서의 태양의 각도 [rad] (초기 각도 0으로 가정)
    psi = alpha_sun - theta;  % 위성과 태양 사이의 상대 각도 [rad]
    
    a_r = a_solar * cos(psi);     % 레이디얼(radial) 성분 가속도
    a_theta = a_solar * sin(psi); % 접선(tangential) 성분 가속도
    % --- 태양압 계산 종료 ---

    % Position update (기존과 동일)
    r_next = r + delta_r;
    theta_next = theta + delta_theta;
    
    % Delta update (★여기가 핵심 수정 부분★)
    r_mid = r + 0.5*delta_r;
    
    % 기존 업데이트 식에 태양압 가속도 항(a_r) 추가
    delta_r_next = delta_r + (r_mid*(delta_theta^2) - (mu/(r^2))*dt^2 + a_r*dt^2);
    
    % 기존 업데이트 식에 태양압 가속도 항(a_theta) 추가
    delta_theta_next = delta_theta - (2*delta_r*delta_theta / r_mid) + (a_theta/r_mid)*dt^2;

    % Status update (기존과 동일)
    r = r_next;
    theta = theta_next;
    delta_r = delta_r_next;
    delta_theta = delta_theta_next;
    
    % Save (기존과 동일)
    r_history(n+1) = r;
    theta_history(n+1) = theta;
end

% --- Step 4: 결과 분석 - 관측각 계산 ---
% 지상국(ES) 위치 - Van Leer 좌표 (라디안으로 변환)
L_e = 33.7758 * (pi/180); % 지상국 위도
l_e = -84.39738 * (pi/180); % 지상국 경도

% 위성 직하점(SSP) 위치 계산
% 정지궤도이므로 위성 위도 L_s는 항상 0 입니다.
L_s = 0; 
% 위성 경도 l_s는 시뮬레이션의 theta와 같습니다.
l_s_history = theta_history;

% 방위각(Azimuth) 저장을 위한 배열
azimuth_history = zeros(n_steps+1, 1);

for i = 1:(n_steps+1)
    l_s = l_s_history(i);
    
    % 중심각 gamma 계산
    cos_gamma = sin(L_s)*sin(L_e) + cos(L_s)*cos(L_e)*cos(l_s - l_e);
    gamma = acos(cos_gamma);

    % 초기 방위각 alpha 계산
    sin_alpha_num = sin(abs(l_e - l_s)) * cos(L_s);
    sin_alpha = sin_alpha_num / sin(gamma);
    % MATLAB의 asin은 -pi/2 ~ pi/2 범위의 값을 반환하므로 주의
    alpha = asin(sin_alpha);

    % 사분면 판단 및 최종 방위각 Az 계산
    % 위성은 항상 적도(L_s=0)에 있으므로 지상국(북반구)보다 남쪽에 있습니다.
    % 위성 경도(l_s)가 지상국 경도(l_e)보다 크면 동쪽, 작으면 서쪽입니다.
    if l_s > l_e  % 남동쪽 (South & East)
        Az = pi - alpha;
    else % 남서쪽 (South & West)
        Az = pi + alpha;
    end
    azimuth_history(i) = Az * (180/pi); % 결과를 도(degree)로 저장
end

% 1년간의 총 방위각 변화량 계산
initial_azimuth = azimuth_history(1);
final_azimuth = azimuth_history(end);
azimuth_drift = final_azimuth - initial_azimuth;

fprintf('Area Density (rho_A): %.2f kg/m^2\n', rho_A);
fprintf('Initial Azimuth: %.4f degrees\n', initial_azimuth);
fprintf('Final Azimuth: %.4f degrees\n', final_azimuth);
fprintf('Total Azimuth Drift over 1 year: %.4f degrees\n', azimuth_drift);

% 이제 그래프는 궤도 대신 시간-방위각 그래프를 그릴 수 있습니다.
figure;
time_axis = (0:n_steps) * dt / (24 * 3600); % 시간을 일(day) 단위로
plot(time_axis, azimuth_history);
grid on;
xlabel('Time [days]');
ylabel('Azimuth [degrees]');
title(['Azimuth Drift for \rho_A = ', num2str(rho_A), ' kg/m^2']);

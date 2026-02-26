clear; clc; close all;

%% ============================================================
% ACTIVE CAMBER FLUTTER SUPPRESSION
% Clean, Verified, HD-Ready Implementation
% =============================================================

%% ===================== PARAMETERS ===========================

p.m     = 35;       % kg/m
p.I     = 3.2;      % kg*m^2
p.kh    = 8e4;      % N/m
p.ka    = 2e4;      % N*m/rad
p.ch    = 200;      
p.ca    = 100;

p.rho   = 1.225;
p.S     = 0.6;
p.cbar  = 0.5;

p.CL_a  = 5.5;
p.CL_d  = 0.8;
p.CM_a  = -1.2;
p.CM_d  = -0.5;

%% ================= EIGENVALUE SWEEP =========================

V_range = linspace(1,80,400);
eigvals = zeros(4,length(V_range));

for i = 1:length(V_range)
    [A,~] = aeroelastic_model(V_range(i),p);
    eigvals(:,i) = eig(A);
end

figure
plot(V_range, real(eigvals),'LineWidth',1.5)
xlabel('Airspeed (m/s)')
ylabel('Real Part of Eigenvalues')
title('Open-Loop Eigenvalue Migration')
grid on

max_real = max(real(eigvals));
[~,idx] = min(abs(max_real));
V_flutter = V_range(idx);

fprintf('\nEstimated Flutter Speed ≈ %.2f m/s\n\n',V_flutter);

%% ================= OPEN LOOP RESPONSE =======================

V = 1.1 * V_flutter;
[A,B] = aeroelastic_model(V,p);

sys = ss(A,B,eye(4),0);

t = 0:0.01:10;

figure
initial(sys,[0.01 0 0.01 0],t)
title('Open-Loop Flutter Response')
grid on

%% ================= CONTROLLABILITY CHECK ====================

Co = ctrb(A,B);
fprintf('Controllability Rank = %d\n',rank(Co));

%% ================= LQR DESIGN ===============================

Q = diag([200 1 400 1]);
R = 1;

K = lqr(A,B,Q,R);

Acl = A - B*K;
sys_cl = ss(Acl,B,eye(4),0);

figure
initial(sys_cl,[0.01 0 0.01 0],t)
title('Closed-Loop Stabilised Response')
grid on

%% ================= CONTROLLED FLUTTER SWEEP =================

eigvals_cl = zeros(4,length(V_range));

for i = 1:length(V_range)
    [A_temp,B_temp] = aeroelastic_model(V_range(i),p);
    K_temp = lqr(A_temp,B_temp,Q,R);
    eigvals_cl(:,i) = eig(A_temp - B_temp*K_temp);
end

figure
plot(V_range, real(eigvals),'--','LineWidth',1.2)
hold on
plot(V_range, real(eigvals_cl),'LineWidth',1.5)
xlabel('Airspeed (m/s)')
ylabel('Real Part of Eigenvalues')
legend('Open Loop','Closed Loop')
title('Flutter Boundary Extension')
grid on

%% ================= OBSERVABILITY CHECK ======================

C = [1 0 0 0;     % measure plunge
     0 0 1 0];    % measure torsion

Ob = obsv(A,C);
fprintf('Observability Rank = %d\n',rank(Ob));

%% ================= OBSERVER DESIGN ==========================

observer_poles = [-8 -9 -10 -11];
L = place(A',C',observer_poles)';

A_aug = [A -B*K;
         L*C A-B*K-L*C];

B_aug = [B;
         B];

sys_aug = ss(A_aug,B_aug,eye(8),0);

figure
initial(sys_aug,[0.01 0 0.01 0 zeros(1,4)],t)
title('Output-Feedback Implementation')
grid on

%% ================= SATURATION TEST ==========================

x0 = [0.01 0 0.01 0];
dt = 0.001;
t_sat = 0:dt:10;

x = x0';
delta_limit = 0.05;

x_store = zeros(4,length(t_sat));

for k = 1:length(t_sat)
    delta = -K*x;
    delta = max(min(delta,delta_limit),-delta_limit);
    dx = A*x + B*delta;
    x = x + dx*dt;
    x_store(:,k) = x;
end

figure
plot(t_sat,x_store(3,:),'LineWidth',1.5)
xlabel('Time (s)')
ylabel('\alpha (rad)')
title('Closed-Loop with Actuator Saturation')
grid on

%% ============================================================
% =================== MODEL FUNCTION ==========================
% =============================================================

function [A,B] = aeroelastic_model(V,p)

q = 0.5 * p.rho * V^2;

L_alpha = q * p.S * p.CL_a;
L_delta = q * p.S * p.CL_d;

M_alpha = q * p.S * p.cbar * p.CM_a;
M_delta = q * p.S * p.cbar * p.CM_d;

A = [ 0 1 0 0;
     -p.kh/p.m -p.ch/p.m -L_alpha/p.m 0;
      0 0 0 1;
      0 0 -(p.ka + M_alpha)/p.I -p.ca/p.I ];

B = [0;
     -L_delta/p.m;
      0;
      M_delta/p.I];
end
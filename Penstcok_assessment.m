%% Intermediate Penstock Problem
%% Material Properties for Steel A and Steel B
close all; clear all ;
% To Define material properties for Steel A and Steel B

% Steel A
steelA.yield_strength            = 700e06; % Yield strength of Steel A in Pa
steelA.ultimate_tensile_strength = 800e06; % Ultimate tensile strength of Steel A in Pa
steelA.fracture_toughness        = 100e06; % Fracture Toughness of Steel A in Pa*m^(1/2)
steelA.elongation                = 18    ; % Elongation in percentage
steelA.price                     = 965   ; % Price of the rolled plate in Pounds per 1000 Kg
steelA.density                   = 9750  ; % Density of Steel A in Kg/m^3

% Steel B
steelB.yield_strength            = 350e06; % Pa
steelB.ultimate_tensile_strength = 500e06; % Pa
steelB.fracture_toughness        = 130e06; % Pa*m^(1/2)
steelB.elongation                = 20    ; % in percentage
steelB.price                     = 525   ; % Pounds per 1000 Kg
steelB.density                   = 9750  ; % Kg/m^3
%% Input Parameters
% Define input parameters
crack_length          = 3e-03            ; % a, Crack depth in m
crack                 = 15e-03           ; % 2c, crack limit in m
crack_size            = crack / 2        ; % a, Crack size in m
pipe_diameter         = 2.6              ; % Pipe diameter in m
pipe_radius           = pipe_diameter/2  ; % Pipe radius in m
static_head           = 529              ; % Static pressure of water in m
water_hammer          = 720              ; % Transient pressure of water in m
rho_water             = 1000             ; % Density of water in Kg/m^3
g                     = 9.81             ; % Gravitational Constant in N
residual_stress       = 70e06            ; % Residual Stress in Pa
%% Overload assessment
% To calculate hoop stress due to static head and water hammer
static_pressure = static_head * rho_water* g;
water_hammer_pressure = water_hammer * rho_water * g;

% To determine the Limiting hoop stress
Max_pressure = max(static_pressure, water_hammer_pressure);

% To calculate the allowable hoop stress for steel A and B (ie. 60 % of the yield strength of Steel A and B)
allowable_hoop_stress_A = 0.6 * steelA.yield_strength;
allowable_hoop_stress_B = 0.6 * steelB.yield_strength;

% To calculate the Wall thickness for steel A and Steel B
wall_thickness_A = (Max_pressure*(pipe_diameter/2))/allowable_hoop_stress_A;
wall_thickness_B = (Max_pressure*(pipe_diameter/2))/allowable_hoop_stress_B;

KIC_A = steelA.fracture_toughness ;
KIC_B = steelB.fracture_toughness ;

fprintf('For Steel A:\n')
fprintf('Wall thickness is %.4f mm\n', wall_thickness_A * 1e03)
fprintf('Hoop stress is %.2f MPa\n', allowable_hoop_stress_A/1e6)
fprintf('For Steel B:\n')
fprintf('Wall thickness is %.4f mm\n', wall_thickness_B * 1e03)
fprintf('Hoop stress is %.2f MPa\n', allowable_hoop_stress_B/1e6)

%% Fracture Assessment
% To calculate m for Steel  A and Steel B
m_A = sqrt(1 + (0.526 *(2*crack_size)^2 / (wall_thickness_A * pipe_diameter)));
m_B = sqrt(1 + (0.526 *(2*crack_size)^2 / (wall_thickness_B * pipe_diameter)));

% To calculate Lr for Steel A and Steel B
Lr_A = ((( wall_thickness_A / crack_length)-( 1 / m_A ))/(( wall_thickness_A / crack_length ) - 1 )) * ( allowable_hoop_stress_A / steelA.yield_strength);
Lr_B = ((( wall_thickness_B / crack_length)-( 1 / m_B ))/(( wall_thickness_B / crack_length ) - 1 )) * ( allowable_hoop_stress_B / steelB.yield_strength);

fprintf('limit load ratio of Steel A at reference length is %.2f \n', Lr_A)
fprintf('limit load ratio of Steel B at reference length is %.2f \n', Lr_B)

% To determine the Maximum Lr for Steel A and Steel B
Lr_max_A = ( steelA.yield_strength + steelA.ultimate_tensile_strength ) / ( 2 * steelA.yield_strength ) ;
Lr_max_B = ( steelB.yield_strength + steelB.ultimate_tensile_strength ) / ( 2 * steelB.yield_strength ) ;

fprintf('Maximum limit load ratio of Steel A is %.4f \n', Lr_max_A)
fprintf('Maximum limit load ratio of Steel B is %.4f \n', Lr_max_B)

% To calculate the Geometric Calibration factor for Steel A and Steel B
Y = ( 1.13 - (0.09 * (crack_length / crack_size ))) / sqrt( 1 + (1.464 * ( crack_length / crack_size )^1.65) ) ;
Y_best= ( 1.13 - (0.09 * (0 / crack_size ))) / sqrt( 1 + (1.464 * ( 0 / crack_size )^1.65) ) ;

fprintf('Geometric Calibration Factor for crack is %.2f \n', Y);
fprintf('Geometric Calibration Factor for best case is %.2f \n', Y_best);

% To Calculate local stress for Steel A and Steel B
Local_stress_A = allowable_hoop_stress_A + residual_stress ;
Local_stress_B = allowable_hoop_stress_B + residual_stress ;

% To calculate local stress for steel A and Steel B with no residual stress
Local_stress_A_NR = allowable_hoop_stress_A ;
Local_stress_B_NR = allowable_hoop_stress_B ;

% To calculate KI for Steel A and Steel B with no residual stress
KI_A_NR =  Y * Local_stress_A_NR * sqrt( pi * crack_length) ;
KI_B_NR =  Y * Local_stress_B_NR * sqrt( pi * crack_length) ;

% To calculate KI for Steel A and Steel B
KI_A_R = Y * Local_stress_A * sqrt( pi * crack_length) ;
KI_B_R = Y * Local_stress_B * sqrt( pi * crack_length) ;

% To calculate  primary KI for Steel A and Steel B
KI_A_P = Y * residual_stress * sqrt( pi * crack_length) ;
KI_B_P = Y * residual_stress * sqrt( pi * crack_length) ;

% To caculate Kr for Steel A and Steel B with no residual
Kr_A_NR = KI_A_NR / KIC_A ;
Kr_B_NR = KI_B_NR / KIC_B ;

fprintf('Stress intensity factor for Steel A without residual stress is %.4f \n', Kr_A_NR)
fprintf('Stress intensity factor for Steel B without residual stress is %.4f \n', Kr_B_NR)

% To caculate Kr for Steel A and Steel B with residual stress
Kr_A_R = KI_A_R / KIC_A ;
Kr_B_R = KI_B_R / KIC_B ;

fprintf('Stress intensity factor for Steel A with residual stress is %.4f \n', Kr_A_R)
fprintf('Stress intensity factor for Steel B with residual stress is %.4f \n', Kr_B_R)

% To caculate Kr for Steel A and Steel B with no residual
Kr_A_P = KI_A_P / KIC_A ;
Kr_B_P = KI_B_P / KIC_B ;

%% Fracture Assessment Diagram
% For Steel A
Lr = linspace(0,3,100);
Kr = (1 + (0.5 .* (Lr).^2)).^(-1/2) .* (0.3 + (0.7 .* exp(-0.6 .* (Lr).^6)));

plot(Lr, Kr,'m');

hold on

% Define the coordinates and slope
Kr_A_axis = Kr_A_R - Kr_A_NR;

x_A_R_axis = 0;
y_A_R_axis = Kr_A_axis;

x_A_R = Lr_A;
y_A_R = Kr_A_R;

slope_A = (y_A_R - y_A_R_axis) / (x_A_R - x_A_R_axis);

x_A_NR = Lr_A;
y_A_NR = Kr_A_NR;

slope_A_NR = (y_A_NR ) / x_A_NR;

x_ext_A_R = linspace(0, max(Lr));
y_ext_A_R = (slope_A * x_ext_A_R) + Kr_A_axis;
[xint_A_R, yint_A_R] = polyxpoly(x_ext_A_R, y_ext_A_R, Lr, Kr);

x_ext_A_NR = linspace(0, max(Lr));
y_ext_A_NR = (slope_A_NR * x_ext_A_NR) ;
[xint_A_NR, yint_A_NR] = polyxpoly(x_ext_A_NR, y_ext_A_NR, Lr, Kr);

Lr_cutoff = Lr_max_A;
[~, Lr_index] = min(abs(Lr-Lr_cutoff));

x1 = Lr_cutoff;
y1 = 0;

[x_intersect_max, y_intersect_max] = polyxpoly([x1, x1], [y1, max(Kr)], Lr(Lr_index:end), Kr(Lr_index:end));

% plot

plot(x_A_R, y_A_R, 'r.','MarkerSize', 15, 'LineWidth', 2);
plot(x_A_NR, y_A_NR, 'b.','MarkerSize', 15, 'LineWidth', 2);


plot(xint_A_R, yint_A_R, 'kx','linewidth',1, 'MarkerSize', 10);
plot([0, xint_A_R], [Kr_A_axis, yint_A_R],'b--');

plot(xint_A_NR, yint_A_NR, 'k+','linewidth',1, 'MarkerSize', 10);
plot([0, xint_A_NR], [0, yint_A_NR],'r--');

plot([x1, x_intersect_max], [y1, y_intersect_max],'linestyle', '-.','LineWidth', 1, 'Color', 'green');
plot(x_intersect_max, y_intersect_max, 'k*','MarkerSize', 10, 'LineWidth', 1);

title("Steel A")
xlabel("Lr_A","FontSize",15);
ylabel("Kr_A","FontSize",15);
legend('L_r - K_r curve ', '(Lr_A,res(Kr_A))','(Lr_A,nores(Kr_A))',' Residual line intersect ',' Residual local stress',' Zero residual line Intersect',' Zero residual local stress','Cutoff (Lr_A) ','Max(Lr_A)', 'FontSize', 10)

grid on
hold off

% For Steel B
Lr = linspace(0,3,100);
Kr = (1 + (0.5 .* (Lr).^2)).^(-1/2) .* (0.3 + (0.7 .* exp(-0.6 .* (Lr).^6)));

plot(Lr, Kr,'m');
hold on;
Kr_B_axis = Kr_B_R - Kr_B_NR;

x_B_R_axis = 0;
y_B_R_axis = Kr_B_axis;

x_B_R = Lr_B;
y_B_R = Kr_B_R;

slope_B_R = (y_B_R - y_B_R_axis) / x_B_R;

x_B_NR = Lr_B;
y_B_NR = Kr_B_NR;

slope_B_NR = (y_B_NR ) / x_B_NR;

Lr_cutoff = Lr_max_B;
[~, Lr_index] = min(abs(Lr-Lr_cutoff));

x1 = Lr_cutoff;
y1 = 0;

[x_intersect_max, y_intersect_max] = polyxpoly([x1, x1], [y1, max(Kr)], Lr(Lr_index:end), Kr(Lr_index:end));

x_ext_B_R = linspace(0, max(Lr));
y_ext_B_R = (slope_B_R * x_ext_B_R) + Kr_B_axis;
[xint_B_R, yint_B_R] = polyxpoly(x_ext_B_R, y_ext_B_R, Lr, Kr);

x_ext_B_NR = linspace(0, max(Lr));
y_ext_B_NR = (slope_B_NR * x_ext_B_NR) ;
[xint_B_NR, yint_B_NR] = polyxpoly(x_ext_B_NR, y_ext_B_NR, Lr, Kr);

% plot
plot(x_B_R, y_B_R, 'r.','MarkerSize', 15, 'LineWidth', 2);
plot(x_B_NR, y_B_NR, 'b.','MarkerSize', 15, 'LineWidth', 2);

plot(xint_B_R, yint_B_R, 'kx','linewidth',1, 'MarkerSize', 10);
plot([0, xint_B_R], [Kr_B_axis, yint_B_R],'b--');

plot(xint_B_NR, yint_B_NR, 'k+','linewidth',1, 'MarkerSize', 10);
plot([0, xint_B_NR], [0, yint_B_NR],'r--');

plot([x1, x_intersect_max], [y1, y_intersect_max],'linestyle', '-.','LineWidth', 1, 'Color', 'green');
plot(x_intersect_max, y_intersect_max, 'k*','MarkerSize', 10, 'LineWidth', 1);

title("Steel B");
xlabel("Lr_B","FontSize",15);
ylabel("Kr_B","FontSize",15);
legend(' L_r - K_r curve ', '(Lr_B,res(Kr_B))','(Lr_B,nores(Kr_B))',' Residual line intersect ',' Residual local stress',' Zero residual line Intersect',' Zero residual local stress','Cutoff (Lr_B) ','Max(Lr_B)', 'FontSize', 10)
grid on;
hold off;

%% Factor of Safety
o_a_A = sqrt((Lr_A^2)+(Kr_A_NR^2));
o_b_A = sqrt((xint_A_NR^2)+(yint_A_NR^2));
F_A = (o_b_A)/(o_a_A);

o_a_B = sqrt((Lr_B^2)+(Kr_B_NR^2));
o_b_B = sqrt((xint_B_NR^2)+(yint_B_NR^2));
F_B = (o_b_B)/(o_a_B);
fprintf(' factor of safety for Steel A is %.2f  \n', F_A);
fprintf(' factor of safety for Steel B is %.2f  \n', F_B);

%% Lifetime Assessment & Inspection Interval 
% To find Critical 
m = 3; % Paris law exponent
C = 1e-11; % Paris law constant
N = 7500;

% Determine the inspection interval and critical crack size
a_crit_A = fsolve(@(ac) KIC_A - Y * allowable_hoop_stress_A * sqrt(pi*ac), crack); % Critical crack length of steel A (m)
a_crit_B = fsolve(@(ac) KIC_B - Y * allowable_hoop_stress_B * sqrt(pi*ac), crack); % Critical crack length of steel B (m)

fprintf('Critical length for Steel A is %.4f \n', a_crit_A/1e-03)
fprintf('Critical length for Steel B is %.4f \n', a_crit_B/1e-03)


% Static pressure for wall thickness of steel A and steel B
hoopstatic_wt_A = (static_pressure * pipe_diameter /2) / wall_thickness_A;
hoopstatic_wt_B = (static_pressure * pipe_diameter /2) / wall_thickness_B;

delta_A = hoopstatic_wt_A/(1e06);
delta_B = hoopstatic_wt_B/(1e06);

% number of cycle until failuer for steel A
Nf_A = (((a_crit_A)^(-0.5)) - ((crack_length)^(-0.5))) / (C * (-0.5) * ((Y * delta_A * (sqrt(pi)))^(3)));
% number of cycle until failuer for steel B
Nf_B = (((a_crit_B)^(-0.5)) - ((crack_length)^(-0.5))) / (C * (-0.5) * ((Y * delta_B * (sqrt(pi)))^(3)));

Tf_a = Nf_A / N;
Tf_b = Nf_B / N;

fprintf(' Total fatigue life for Steel A is %.f  \n', Tf_a);
fprintf(' Total fatigue life for Steel B is %.f  \n', Tf_b);

a_i = crack_length; % initial crack size
sigma_range_A = delta_A; % range of stress intensity factor
cycles_A = 0:7500; % number of cycles for 50 years

% Calculate the crack growth rate for each cycle
da_a = zeros(size(cycles_A));
a_a = a_i;
for i = 2:length(cycles_A)
    if a_a <= a_crit_A
        delta_K_a = Y * sigma_range_A^m * sqrt(pi^1.5 *(a_a + da_a(i-1)))/(1*10^6);
        da_a(i) = C *((delta_K_a)^m) *cycles_A(i);
        a_a = a_a + cumsum(da_a(i));
    else
        break
    end
end

% Calculate the new values of a for each cycle
cycles_A = 0:length(da_a)-1;
a_a = a_i + cumsum(da_a);

plot(cycles_A, a_a,'r--')
xlabel('Number of cycles')
ylabel('Crack size (mm)')

a_i = crack_length; % initial crack size
sigma_range_B = delta_B; % range of stress intensity factor
cycles_B = 0:7500; % number of cycles for 50 years
% Calculate the crack growth rate for each cycle
da_b = zeros(size(cycles_B));
a_b = a_i;
for i = 2:length(cycles_B)
    if a_b <= a_crit_B
        delta_K_b = (Y * sigma_range_B^m *sqrt(pi^1.5*(a_b + da_b(i-1))))/(1*10^6);
        da_b(i) = C * (delta_K_b)^m *cycles_B(i);
        a_b = a_b + cumsum(da_b(i));
    else
        break;
    end
end
cycles_B = 0:length(da_b)-1;
a_b = a_i + cumsum(da_b);

% Plot the crack growth curve

plot(cycles_B, a_b,'r--')
xlabel('Number of cycles')
ylabel('Crack size (mm)')
%% Revised Wall Thickness
%Steel A
wall_thickness_A_new = wall_thickness_A * 1.26;
Hoop_stress_new = (static_pressure*(pipe_diameter/2))/wall_thickness_A_new;
m_A_new = (1+0.526*((2*crack_length)^2/wall_thickness_A_new*pipe_diameter))^0.5;
Lr_A_new = ((wall_thickness_A_new/crack_length)-(1/m_A_new))/((wall_thickness_A_new/crack_length)-1)*(Hoop_stress_new/steelA.yield_strength);
Lr_max_A_new = (steelA.yield_strength +steelA.ultimate_tensile_strength)/(2*steelA.yield_strength);
Kr_A_NR_new = (Y*Local_stress_A_NR*sqrt(pi*crack_length))/(KIC_A);
Kr_A_R_new = ((Y*Local_stress_A*sqrt(pi*crack_length))/KIC_A);

Lr = linspace(0, 3, 100);
Kr = (1 + 0.5*Lr.^2).^(-0.5) .* (0.3 + 0.7*exp(-0.6*Lr.^6));
plot(Lr, Kr);
hold on;

Kr_A_axis_new = Kr_A_R_new - Kr_A_NR_new;

x_A_R_axis_new = 0;
y_A_R_axis_new = Kr_A_axis_new;

x_A_R_new = Lr_A_new;
y_A_R_new = Kr_A_R_new;

slope_A_R_new = (y_A_R_new - y_A_R_axis_new) / x_A_R_new;

x_A_NR_new = Lr_A_new;
y_A_NR_new = Kr_A_NR_new;

slope_A_NR_new = (y_A_NR_new) / x_A_NR_new;

Lr_cutoff = Lr_max_A_new;
[~, Lr_index] = min(abs(Lr-Lr_cutoff));

x1 = Lr_cutoff;
y1 = 0;

[x_intersect_max, y_intersect_max] = polyxpoly([x1, x1], [y1, max(Kr)], Lr(Lr_index:end), Kr(Lr_index:end));
x_ext_A_R_new = linspace(0, max(Lr));
y_ext_A_R_new = (slope_A_R_new * x_ext_A_R_new) + Kr_A_axis_new;
[xint_A_R_new, yint_A_R_new] = polyxpoly(x_ext_A_R_new, y_ext_A_R_new, Lr, Kr);

x_ext_A_NR_new = linspace(0, max(Lr));
y_ext_A_NR_new = (slope_A_NR_new * x_ext_A_NR_new) ;
[xint_A_NR_new, yint_A_NR_new] = polyxpoly(x_ext_A_NR_new, y_ext_A_NR_new, Lr, Kr);

% plot

plot(x_A_R_new, y_A_R_new, 'r.','MarkerSize', 15, 'LineWidth', 2);
plot(x_A_NR_new, y_A_NR_new, 'b.','MarkerSize', 15, 'LineWidth', 2);

plot(xint_A_R_new, yint_A_R_new, 'kx','linewidth',1, 'MarkerSize', 10);
plot([0, xint_A_R_new], [Kr_A_axis_new, yint_A_R_new],'b--');

plot(xint_A_NR_new, yint_A_NR_new, 'k+','linewidth',1, 'MarkerSize', 10);
plot([0, xint_A_NR_new], [0, yint_A_NR_new],'r--');

plot([x1, x_intersect_max], [y1, y_intersect_max],'linestyle', '-.','LineWidth', 1, 'Color', 'green');
plot(x_intersect_max, y_intersect_max, 'k*','MarkerSize', 10, 'LineWidth', 1);

title("Revised Steel A");
xlabel("Lr_A(new)","FontSize",15);
ylabel("Kr_A(new)","FontSize",15);
legend(' L_r - K_r curve ', '(Lr_A_n,res(Kr_A_n))','((Lr_A_n,nores(Kr_A_n))',' Residual line intersect ',' Residual local stress',' Zero residual line Intersect',' Zero residual local stress','Cutoff (Lr_A_n) ','Max(Lr_A_n)', 'FontSize', 10)
grid on
hold off
o_a_A_new = sqrt((Lr_A_new^2)+(Kr_A_NR_new^2));
o_b_A_new = sqrt((xint_A_NR_new^2)+(yint_A_NR_new^2));
F_A_new = (o_b_A_new)/(o_a_A_new);

%critical crack length
% steel A:
a_crit_A_new = fsolve(@(ac) KIC_A - Y * Hoop_stress_new* sqrt(pi*ac), crack);

%Lifetime assessment

delta_a_new = Hoop_stress_new/(1*10^6);

Nf_A_new = (2*(crack_length^-0.5 - a_crit_A_new^-0.5))/( C*(Y^m * delta_a_new^m * pi^1.5 ));

Tf_a_new = Nf_A_new/N;

a_i = crack_length; % initial crack size
sigma_range_a = delta_a_new; % range of stress intensity factor
cycles_a_new = 0:7500; % number of cycles

% Calculate the crack growth rate for each cycle
da_a_new = zeros(size(cycles_a_new));
a_a_new = a_i;
for i = 2:length(cycles_a_new)
    if a_a_new <= a_crit_A_new
        delta_K_a_new = Y*sigma_range_a^m*sqrt(pi^1.5*(a_a_new + da_a_new(i-1)))/(1*10^6);
        da_a_new(i) = C*(delta_K_a_new)^m *cycles_a_new(i);
        a_a_new = a_a_new + da_a_new(i);
    else
        break
    end
end

% Calculate the new values of a for each cycle
cycles_a_new = 0:length(da_a_new)-1;
a_a_new = a_i + cumsum(da_a_new);

%% Plot the crack growth curve
plot(cycles_a_new, a_a_new,'r--')
grid on
title('Revised Steel A')
xlabel('Number of cycles')
ylabel('Crack size (mm)')
fprintf('Revised Steel A:\n')
fprintf('Revised Wall thickness is %.4f mm\n', wall_thickness_A_new * 1e03)
fprintf('Revised Hoop stress is %.2f MPa\n', Hoop_stress_new/1e6)
fprintf('Revised limit load ratio of Steel A at reference length is %.2f \n', Lr_A_new)
fprintf('Revised Maximum limit load ratio of Steel A is %.4f \n', Lr_max_A_new)
fprintf('Revised Stress intensity factor for Steel A without residual stress is %.4f \n', Kr_A_NR_new)
fprintf('Revised Stress intensity factor for Steel A with residual stress is %.4f \n', Kr_A_R_new)
fprintf('Revised critical length for Steel A is %.4f \n', a_crit_A_new/1e-03)
fprintf('Revised factor of safety for Steel A is %.2f  \n', F_A_new)
fprintf('Revised Total fatigue life for Steel A is %.4f  \n', Tf_a_new)


% Steel B
wall_thickness_B_new = wall_thickness_B * 1.26;
Hoop_stress_B_new = (static_pressure*(pipe_diameter/2))/wall_thickness_B_new;
m_B_new = (1+0.526*((2*crack_length)^2/wall_thickness_B_new*pipe_diameter))^0.5;
Lr_B_new = ((wall_thickness_B_new/crack_length)-(1/m_B_new))/((wall_thickness_B_new/crack_length)-1)*(Hoop_stress_B_new/steelB.yield_strength);
Lr_max_B_new = (steelB.yield_strength +steelB.ultimate_tensile_strength)/(2*steelB.yield_strength);
Kr_B_NR_new = (Y*Local_stress_B_NR*sqrt(pi*crack_length))/(KIC_B);
Kr_B_R_new = ((Y*Local_stress_B*sqrt(pi*crack_length))/KIC_B);

Lr = linspace(0, 3, 100);
Kr = (1 + 0.5*Lr.^2).^(-0.5) .* (0.3 + 0.7*exp(-0.6*Lr.^6));
plot(Lr, Kr);
hold on;

Kr_B_axis_new = Kr_B_R_new - Kr_B_NR_new;

x_B_R_axis_new = 0;
y_B_R_axis_new = Kr_B_axis_new;

x_B_R_new = Lr_B_new;
y_B_R_new = Kr_B_R_new;

slope_B_R_new = (y_B_R_new - y_B_R_axis_new) / x_B_R_new;

x_B_NR_new = Lr_B_new;
y_B_NR_new = Kr_B_NR_new;

slope_B_NR_new = (y_B_NR_new) / x_B_NR_new;

Lr_cutoff = Lr_max_B_new;
[~, Lr_index] = min(abs(Lr-Lr_cutoff));

x1 = Lr_cutoff;
y1 = 0;

[x_intersect_max, y_intersect_max] = polyxpoly([x1, x1], [y1, max(Kr)], Lr(Lr_index:end), Kr(Lr_index:end));
x_ext_B_R_new = linspace(0, max(Lr));
y_ext_B_R_new = (slope_B_R_new * x_ext_B_R_new) + Kr_B_axis_new;
[xint_B_R_new, yint_B_R_new] = polyxpoly(x_ext_B_R_new, y_ext_B_R_new, Lr, Kr);

x_ext_B_NR_new = linspace(0, max(Lr));
y_ext_B_NR_new = (slope_B_NR_new * x_ext_B_NR_new) ;
[xint_B_NR_new, yint_B_NR_new] = polyxpoly(x_ext_B_NR_new, y_ext_B_NR_new, Lr, Kr);

% plot
plot(x_B_R_new, y_B_R_new, 'r.','MarkerSize', 15, 'LineWidth', 2);
plot(x_B_NR_new, y_B_NR_new, 'b.','MarkerSize', 15, 'LineWidth', 2);

plot(xint_B_R_new, yint_B_R_new, 'kx','linewidth',1, 'MarkerSize', 10);
plot([0, xint_B_R_new], [Kr_B_axis_new, yint_B_R_new],'b--');

plot(xint_B_NR_new, yint_B_NR_new, 'k+','linewidth',1, 'MarkerSize', 10);
plot([0, xint_B_NR_new], [0, yint_B_NR_new],'r--');

plot([x1, x_intersect_max], [y1, y_intersect_max],'linestyle', '-.','LineWidth', 1, 'Color', 'green');
plot(x_intersect_max, y_intersect_max, 'k*','MarkerSize', 10, 'LineWidth', 1);

title("Revised Steel B");
xlabel("Lr_B(new)","FontSize",15);
ylabel("Kr_B(new)","FontSize",15);
legend(' L_r - K_r curve ', '(Lr_B_n,res(Kr_B_n))','((Lr_B_n,nores(Kr_B_n))',' Residual line intersect ',' Residual local stress',' Zero residual line Intersect',' Zero residual local stress','Cutoff (Lr_B_n) ','Max(Lr_B_n)', 'FontSize', 10)
grid on;
hold off;
o_a_B_new = sqrt((Lr_B_new^2)+(Kr_B_NR_new^2));
o_b_B_new = sqrt((xint_B_NR_new^2)+(yint_B_NR_new^2));
F_B_new = (o_b_B_new)/(o_a_B_new);

%%critical crack length
% steel A:
a_crit_B_new = fsolve(@(ac) KIC_B - Y * Hoop_stress_B_new* sqrt(pi*ac), crack);

%Lifetime assessment

delta_b_new = Hoop_stress_B_new/(1*10^6);

Nf_B_new = (2*(crack_length^-0.5 - a_crit_B_new^-0.5))/( C*(Y^m * delta_b_new^m * pi^1.5 ));

Tf_b_new = Nf_B_new/N;

a_i = crack_length; % initial crack size
sigma_range_b = delta_b_new; % range of stress intensity factor
cycles_b_new = 0:7500; % number of cycles

% Calculate the crack growth rate for each cycle
da_b_new = zeros(size(cycles_b_new));
a_b_new = a_i;
for i = 2:length(cycles_b_new)
    if a_b_new <= a_crit_B_new
        delta_K_b_new = Y*sigma_range_b^m*sqrt(pi^1.5*(a_b_new + da_b_new(i-1)))/(1*10^6);
        da_b_new(i) = C*(delta_K_b_new)^m *cycles_b_new(i);
        a_b_new = a_b_new + da_b_new(i);
    else
        break
    end
end

% Calculate the new values of a for each cycle
cycles_b_new = 0:length(da_b_new)-1;
a_b_new = a_i + cumsum(da_b_new);

%% Plot the crack growth curve
plot(cycles_b_new, a_b_new,'r--')
grid on
title('Revised Steel B')
xlabel('Number of cycles')
ylabel('Crack size (mm)')
fprintf('Revised Steel B:\n')
fprintf('Revised Wall thickness is %.4f mm\n', wall_thickness_B_new * 1e03)
fprintf('Revised Hoop stress is %.2f MPa\n', Hoop_stress_B_new/1e6)
fprintf('Revised limit load ratio of Steel B at reference length is %.2f \n', Lr_B_new)
fprintf('Revised Maximum limit load ratio of Steel B is %.4f \n', Lr_max_B_new)
fprintf('Revised Stress intensity factor for Steel B without residual stress is %.4f \n', Kr_B_NR_new)
fprintf('Revised Stress intensity factor for Steel B with residual stress is %.4f \n', Kr_B_R_new)
fprintf('Revised critical length for Steel B is %.4f \n', a_crit_B_new/1e-03)
fprintf('Revised factor of safety for Steel B is %.2f  \n', F_B_new)
fprintf('Revised Total fatigue life for Steel B is %.4f  \n', Tf_b_new)

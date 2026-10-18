function saturacion = checkpoint13_saturacion_actuador( ...
    estados, control, region12, tau_limites)
%CHECKPOINT 13 - SATURACION DEL ACTUADOR
%Repite el barrido no lineal del Checkpoint 12 limitando
%el par realmente aplicado a la planta.
%Compara la region ideal con la region obtenida para
%diferentes limites de par total.

fprintf('\n');
fprintf('=============================================\n');
fprintf('CHECKPOINT 13 - SATURACION DEL ACTUADOR\n');
fprintf('=============================================\n\n');

if ~control.candidato_encontrado
    error('No existe un candidato estable del Checkpoint 10.');
end

% Ganancias
Kp_theta = control.ganancias(1);
Kd_theta = control.ganancias(2);
Kp_v = control.ganancias(3);
Ki_v = control.ganancias(4);
K = [Ki_v Kp_v Kp_theta Kd_theta];

fprintf('Ganancias utilizadas:\n');
fprintf('Kp_theta = %.6f\n',Kp_theta);
fprintf('Kd_theta = %.6f\n',Kd_theta);
fprintf('Kp_v     = %.6f\n',Kp_v);
fprintf('Ki_v     = %.6f\n\n',Ki_v);

% Modelo no lineal ya parametrizado en Checkpoint 12
f_num = matlabFunction( region12.f_parametros, ...
    'Vars',{estados.x, estados.tau});

% Datos y criterios del Checkpoint 12
theta0_deg_vec = region12.theta0_deg;
theta_tol = region12.theta_tol;
omega_tol = region12.omega_tol;
v_tol = region12.v_tol;
dr_tol = region12.dr_tol;
ventana_final = region12.ventana_final;
t_final = region12.referencia.t(end);

n_angulos = length(theta0_deg_vec);
n_limites = length(tau_limites);

recupera = false(n_angulos,n_limites);
max_tau_cmd = NaN(n_angulos,n_limites);
max_tau_real = NaN(n_angulos,n_limites);
max_dr = NaN(n_angulos,n_limites);
max_v = NaN(n_angulos,n_limites);
max_omega = NaN(n_angulos,n_limites);
tiempo_saturado = NaN(n_angulos,n_limites);
porcentaje_saturado = NaN(n_angulos,n_limites);

limite_bilateral = zeros(n_limites,1);
casos_recuperados = zeros(n_limites,1);

opciones = odeset( ...
    'RelTol',1e-7, ...
    'AbsTol',1e-9, ...
    'Events',@evento_caida);

fprintf('Region ideal del Checkpoint 12:\n');
fprintf('[-%.2f , +%.2f] grados\n\n', ...
    region12.limite_bilateral, ...
    region12.limite_bilateral);

% Barrido para cada limite de par
for j = 1:n_limites

    tau_lim = tau_limites(j);

    fprintf('Probando limite total: %.2f N*m\n',tau_lim);

    for i = 1:n_angulos

        theta0 = theta0_deg_vec(i)*pi/180;
        x0 = [0; 0; theta0; 0];

        try
            [t,x] = ode45( ...
                @(t,x) dinamica_saturada(t,x,f_num,K,tau_lim), ...
                [0 t_final],x0,opciones);

            dr = x(:,1);
            v = x(:,2);
            theta = x(:,3);
            omega_p = x(:,4);

            tau_cmd = -(x*K');
            tau_real = min(max(tau_cmd,-tau_lim),tau_lim);

            sat = abs(tau_cmd) > tau_lim;

            max_tau_cmd(i,j) = max(abs(tau_cmd));
            max_tau_real(i,j) = max(abs(tau_real));
            max_dr(i,j) = max(abs(dr));
            max_v(i,j) = max(abs(v));
            max_omega(i,j) = max(abs(omega_p));

            if length(t) > 1
                tiempo_saturado(i,j) = trapz(t,double(sat));
                porcentaje_saturado(i,j) = ...
                    100*tiempo_saturado(i,j)/t(end);
            else
                tiempo_saturado(i,j) = 0;
                porcentaje_saturado(i,j) = 0;
            end

            recupera(i,j) = evaluar_recuperacion( ...
                t,x,t_final,ventana_final, ...
                theta_tol,omega_tol,v_tol,dr_tol);

        catch
            recupera(i,j) = false;
        end
    end

    casos_recuperados(j) = sum(recupera(:,j));

    limite_bilateral(j) = obtener_limite_bilateral( ...
        theta0_deg_vec,recupera(:,j));

    fprintf('Region recuperada: [-%.2f , +%.2f] grados\n\n', ...
        limite_bilateral(j),limite_bilateral(j));
end

% Resumen
tau_motor_lim = tau_limites(:)/2;

tabla_resumen = table( ...
    tau_limites(:), ...
    tau_motor_lim, ...
    casos_recuperados, ...
    limite_bilateral, ...
    'VariableNames',{ ...
    'TauLimTotal_Nm', ...
    'TauLimMotor_Nm', ...
    'CasosRecuperados', ...
    'LimiteBilateral_deg'});

fprintf('---------------------------------------------\n');
fprintf('RESUMEN\n');
fprintf('---------------------------------------------\n\n');

disp(tabla_resumen);

% Grafica limite de par vs region recuperable
figure;
plot(tau_limites,limite_bilateral,'o-','LineWidth',1.3);
grid on;
xlabel('\tau_{lim,total} (N*m)');
ylabel('|\theta_0| maximo recuperable (grados)');
title('Region de recuperacion vs. limite de par');

% Comparacion ideal
figure;
plot(tau_limites,limite_bilateral,'o-','LineWidth',1.3);
hold on;
plot(tau_limites, ...
    region12.limite_bilateral*ones(size(tau_limites)), ...
    '--','LineWidth',1.2);
grid on;
xlabel('\tau_{lim,total} (N*m)');
ylabel('|\theta_0| maximo recuperable (grados)');
title('Recuperacion con saturacion vs. caso ideal');
legend('Con saturacion','Limite ideal');
hold off;

% Guardar resultados
saturacion.K = K;
saturacion.tau_limites = tau_limites;
saturacion.tau_motor_lim = tau_motor_lim;
saturacion.theta0_deg = theta0_deg_vec;
saturacion.recupera = recupera;
saturacion.max_tau_cmd = max_tau_cmd;
saturacion.max_tau_real = max_tau_real;
saturacion.max_dr = max_dr;
saturacion.max_v = max_v;
saturacion.max_omega = max_omega;
saturacion.tiempo_saturado = tiempo_saturado;
saturacion.porcentaje_saturado = porcentaje_saturado;
saturacion.casos_recuperados = casos_recuperados;
saturacion.limite_bilateral = limite_bilateral;
saturacion.limite_ideal = region12.limite_bilateral;
saturacion.tabla_resumen = tabla_resumen;

fprintf('=============================================\n');

end

function dx = dinamica_saturada(~,x,f_num,K,tau_lim)
tau_cmd = -K*x;
tau_real = min(max(tau_cmd,-tau_lim),tau_lim);
dx = f_num(x,tau_real);
dx = dx(:);
end

function recupera = evaluar_recuperacion( ...
    t,x,t_final,ventana_final, ...
    theta_tol,omega_tol,v_tol,dr_tol)

if t(end) < t_final-1e-4 || any(~isfinite(x(:)))
    recupera = false;
    return;
end

idx = t >= (t_final-ventana_final);

if ~any(idx)
    recupera = false;
    return;
end

dr = x(idx,1);
v = x(idx,2);
theta = x(idx,3);
omega_p = x(idx,4);

recupera = ...
    all(abs(theta) <= theta_tol) && ...
    all(abs(omega_p) <= omega_tol) && ...
    all(abs(v) <= v_tol) && ...
    all(abs(dr) <= dr_tol);

end

function limite = obtener_limite_bilateral(theta_vec,recupera)

paso = abs(theta_vec(2)-theta_vec(1));
magnitudes = 0:paso:max(abs(theta_vec));
limite = 0;

for mag = magnitudes

    idx_p = find(abs(theta_vec-mag) < 1e-10,1);
    idx_n = find(abs(theta_vec+mag) < 1e-10,1);

    if isempty(idx_p) || isempty(idx_n)
        break;
    end

    if recupera(idx_p) && recupera(idx_n)
        limite = mag;
    else
        break;
    end
end

end

function [value,isterminal,direction] = evento_caida(~,x)

theta = x(3);
value = pi/2-abs(theta);
isterminal = 1;
direction = -1;

end
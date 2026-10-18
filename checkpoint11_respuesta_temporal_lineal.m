function respuesta = checkpoint11_respuesta_temporal_lineal(num,control)
%CHECKPOINT 11 - RESPUESTA TEMPORAL LINEAL EN LAZO CERRADO
%Evalua el comportamiento temporal del controlador:
%   PD de equilibrio + PI de velocidad
%utilizando el modelo lineal:
%   x_dot = A*x + B*tau
%con:
%   x = [d_r; v; theta; omega_p]
%La ley equivalente de control es:
%   tau = -K*x
%   K = [Ki_v Kp_v Kp_theta Kd_theta]
%Se analiza una perturbacion angular inicial y se calculan:
%   - polos de lazo cerrado
%   - desplazamiento
%   - velocidad
%   - angulo
%   - velocidad angular
%   - par de control
%   - par aproximado por motor
%   - tiempo de establecimiento angular
%IMPORTANTE:
%Se toma d_r(0)=0, considerando el desplazamiento respecto
%al punto donde se activa el controlador.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 11 - RESPUESTA TEMPORAL LINEAL\n');
fprintf('========================================\n\n');

% ================================================================
% 1. RECUPERAR MODELO
A = num.A;
B = num.B;

if ~control.candidato_encontrado
    error(['Checkpoint 10 no contiene un candidato estable. ' ...
           'No se puede ejecutar Checkpoint 11.']);
end


% ================================================================
% 2. RECUPERAR GANANCIAS
% Checkpoint 10 guarda:
% control.ganancias = [Kp_theta Kd_theta Kp_v Ki_v]

Kp_theta = control.ganancias(1);
Kd_theta = control.ganancias(2);
Kp_v     = control.ganancias(3);
Ki_v     = control.ganancias(4);

fprintf('Ganancias utilizadas:\n\n');
fprintf('Kp_theta = %.6f\n',Kp_theta);
fprintf('Kd_theta = %.6f\n',Kd_theta);
fprintf('Kp_v     = %.6f\n',Kp_v);
fprintf('Ki_v     = %.6f\n\n',Ki_v);


% ================================================================
% 3. VECTOR EQUIVALENTE DE REALIMENTACION
% tau =
% -Ki_v*d_r
% -Kp_v*v
% -Kp_theta*theta
% -Kd_theta*omega_p
%Por tanto:
% K = [Ki_v Kp_v Kp_theta Kd_theta]
K = [ Ki_v Kp_v Kp_theta Kd_theta];

fprintf('Vector equivalente K:\n');
disp(K);

% ================================================================
% 4. MATRIZ DE LAZO CERRADO

Acl = A - B*K;
fprintf('Matriz A de lazo cerrado:\n');
disp(Acl);

% ================================================================
% 5. POLOS DE LAZO CERRADO
polos_cl = eig(Acl);

fprintf('----------------------------------------\n');
fprintf('POLOS DE LAZO CERRADO\n');
fprintf('----------------------------------------\n\n');
disp(polos_cl);

if all(real(polos_cl) < 0)
    fprintf('RESULTADO:\n');
    fprintf('Todos los polos tienen parte real negativa.\n');
    fprintf('El modelo lineal de lazo cerrado es estable.\n\n');
    estable = true;
else
    fprintf('ADVERTENCIA:\n');
    fprintf('Existe al menos un polo no estable.\n\n');
    estable = false;
end


% ================================================================
% 6. COMPARACION CON CHECKPOINT 10

fprintf('Comparacion con los polos del Checkpoint 10:\n\n');

disp(control.polos);
diferencia_polos = max(abs(sort_complex(polos_cl) - ...
            sort_complex(control.polos)));

fprintf('Diferencia maxima entre ambos calculos: %.6e\n\n', ...
    diferencia_polos);

% ================================================================
% 7. CONDICION INICIAL
%Perturbacion angular inicial:
% d_r(0)     = 0 m
% v(0)       = 0 m/s
% theta(0)   = 0.1 rad
% omega_p(0) = 0 rad/s

x0 = [ 0; 0; 0.1; 0];

theta0_deg = x0(3)*180/pi;

fprintf('----------------------------------------\n');
fprintf('CONDICION INICIAL\n');
fprintf('----------------------------------------\n\n');
fprintf('d_r(0)     = %.6f m\n',x0(1));
fprintf('v(0)       = %.6f m/s\n',x0(2));
fprintf('theta(0)   = %.6f rad\n',x0(3));
fprintf('theta(0)   = %.6f grados\n',theta0_deg);
fprintf('omega_p(0) = %.6f rad/s\n\n',x0(4));


% ================================================================
% 8. TIEMPO DE SIMULACION

t = 0:0.001:5;

% ================================================================
% 9. SISTEMA AUTONOMO DE LAZO CERRADO
C = eye(4);
Bcl = zeros(4,1);
Dcl = zeros(4,1);
sys_cl = ss(Acl,Bcl,C,Dcl);

% ================================================================
% 10. RESPUESTA A CONDICIONES INICIALES
[~,t,x] = initial(sys_cl,x0,t);

% Separar estados
d_r     = x(:,1);
v       = x(:,2);
theta   = x(:,3);
omega_p = x(:,4);
theta_deg = theta*180/pi;


% ================================================================
% 11. CALCULAR PAR DE CONTROL
% tau(t) = -K*x(t)

tau = -(x*K');
% Par ideal por motor
%Se supone que los dos motores comparten el par total
%de forma aproximadamente uniforme.
tau_motor = tau/2;

% ================================================================
% 12. METRICAS DE RESPUESTA

max_dr = max(abs(d_r));
max_v = max(abs(v));
max_theta = max(abs(theta));
max_theta_deg = max(abs(theta_deg));
max_omega = max(abs(omega_p));
max_tau = max(abs(tau));
max_tau_motor = max(abs(tau_motor));


% ================================================================
% 13. TIEMPO DE ESTABLECIMIENTO ANGULAR
%Se considera establecida la respuesta cuando: |theta| <= 0.02 rad

banda_theta = 0.02;
indices_fuera = find(abs(theta) > banda_theta);
if isempty(indices_fuera)
    ts_theta = 0;
elseif indices_fuera(end) < length(t)
    ts_theta = t(indices_fuera(end)+1);
else
    ts_theta = NaN;
end

% ================================================================
% 14. RESULTADOS NUMERICOS

fprintf('----------------------------------------\n');
fprintf('RESULTADOS TEMPORALES\n');
fprintf('----------------------------------------\n\n');

fprintf('Maximo |d_r|:\n');
fprintf('%.6f m\n\n',max_dr);

fprintf('Maximo |v|:\n');
fprintf('%.6f m/s\n\n',max_v);

fprintf('Maximo |theta|:\n');
fprintf('%.6f rad\n',max_theta);
fprintf('%.6f grados\n\n',max_theta_deg);

fprintf('Maximo |omega_p|:\n');
fprintf('%.6f rad/s\n\n',max_omega);

fprintf('Maximo |tau total|:\n');
fprintf('%.6f N*m\n\n',max_tau);

fprintf('Maximo |tau por motor|:\n');
fprintf('%.6f N*m\n\n',max_tau_motor);

if ~isnan(ts_theta)

    fprintf(['Tiempo de establecimiento angular ' '(|theta| <= %.3f rad):\n'], ...
             banda_theta);
    fprintf('%.6f s\n\n',ts_theta);
else
    fprintf(['La respuesta angular no permanece dentro de ' ...
             'la banda definida durante el intervalo simulado.\n\n']);
end


% ================================================================
% 15. VALORES FINALES

fprintf('Valores finales a t = %.2f s:\n\n',t(end));
fprintf('d_r     = %.8f m\n',d_r(end));
fprintf('v       = %.8f m/s\n',v(end));
fprintf('theta   = %.8f rad\n',theta(end));
fprintf('omega_p = %.8f rad/s\n',omega_p(end));
fprintf('tau     = %.8f N*m\n\n',tau(end));


% ================================================================
% 16. GRAFICAS

figure;
plot(t,d_r,'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('d_r (m)');
title('Desplazamiento lineal');

figure;
plot(t,v,'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('v (m/s)');
title('Velocidad lineal');

figure;
plot(t,theta_deg,'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('\theta (grados)');
title('Angulo de inclinacion');

figure;
plot(t,omega_p,'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('\omega_p (rad/s)');
title('Velocidad angular del cuerpo');

figure;
plot(t,tau,'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('\tau (N*m)');
title('Par total solicitado por el controlador');

% ================================================================
% 17. GUARDAR RESULTADOS

respuesta.Acl = Acl;
respuesta.K = K;
respuesta.polos = polos_cl;
respuesta.estable = estable;
respuesta.diferencia_polos = diferencia_polos;
respuesta.x0 = x0;
respuesta.t = t;
respuesta.x = x;
respuesta.d_r = d_r;
respuesta.v = v;
respuesta.theta = theta;
respuesta.omega_p = omega_p;
respuesta.tau = tau;
respuesta.tau_motor = tau_motor;
respuesta.max_dr = max_dr;
respuesta.max_v = max_v;
respuesta.max_theta = max_theta;
respuesta.max_theta_deg = max_theta_deg;
respuesta.max_omega = max_omega;
respuesta.max_tau = max_tau;
respuesta.max_tau_motor = max_tau_motor;
respuesta.ts_theta = ts_theta;
fprintf('========================================\n');
end

function p = sort_complex(p)
%Ordena polos para facilitar comparacion numerica.

p = p(:);
[~,idx] = sortrows([real(p) imag(p)],[1 2]);
p = p(idx);

end
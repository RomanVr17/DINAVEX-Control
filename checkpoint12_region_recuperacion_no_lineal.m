function region = checkpoint12_region_recuperacion_no_lineal( ...
    estados, p, control)
%CHECKPOINT 12 - REGION DE RECUPERACION NO LINEAL
%Utiliza directamente el modelo no lineal obtenido en:
%checkpoint3_estados mediante:estados.f = f(x,tau)
%Objetivos:
%1. Evaluar el controlador PD + PI de velocidad sobre
%   el modelo no lineal.
%2. Comparar cualitativamente el caso theta0 = 0.1 rad
%   con la respuesta lineal del Checkpoint 11.
%3. Realizar un barrido de theta0.
%4. Determinar la region de recuperacion angular ideal.
%5. Registrar para cada theta0:
%   - recuperacion del equilibrio
%   - par maximo requerido
%   - desplazamiento maximo
%   - velocidad maxima
%   - velocidad angular maxima
%Para el PI de velocidad se supone:
%
%   v_ref = 0
%   d_r(0) = 0
%de manera que la integral de velocidad corresponde
%al desplazamiento respecto al origen de control.

fprintf('\n');
fprintf('=============================================\n');
fprintf('CHECKPOINT 12 - REGION DE RECUPERACION\n');
fprintf('MODELO NO LINEAL\n');
fprintf('=============================================\n\n');

% 1. VERIFICAR CONTROLADOR
if ~control.candidato_encontrado

    error(['Checkpoint 10 no contiene un candidato estable. ' ...
           'No es posible ejecutar Checkpoint 12.']);

end

% 2. RECUPERAR GANANCIAS
Kp_theta = control.ganancias(1);
Kd_theta = control.ganancias(2);
Kp_v     = control.ganancias(3);
Ki_v     = control.ganancias(4);
fprintf('Ganancias utilizadas:\n\n');
fprintf('Kp_theta = %.6f\n',Kp_theta);
fprintf('Kd_theta = %.6f\n',Kd_theta);
fprintf('Kp_v     = %.6f\n',Kp_v);
fprintf('Ki_v     = %.6f\n\n',Ki_v);

% Vector equivalente
% tau =-Ki_v*d_r-Kp_v*v-Kp_theta*theta-Kd_theta*omega_p

K = [ Ki_v Kp_v Kp_theta Kd_theta];
fprintf('Vector equivalente K:\n');
disp(K);

% 3. RECUPERAR MODELO NO LINEAL DEL CHECKPOINT 3
f_simbolica = estados.f;
fprintf('Modelo no lineal recuperado desde estados.f\n');
fprintf('No se vuelven a derivar las ecuaciones dinamicas.\n\n');

% 4. SUSTITUIR PARAMETROS NUMERICOS
%Se conservan como variables:
%   d_r, v, theta, omega_p, tau
%Los demas simbolos presentes en estados.f se buscan
%automaticamente dentro de la estructura p.
f_parametros = sustituir_parametros( ...
    f_simbolica, estados, p);

% 5. CONVERTIR MODELO SIMBOLICO A FUNCION NUMERICA
%Se genera: f_num(x,tau) que devuelve: dx/dt
f_num = matlabFunction( f_parametros, 'Vars',{estados.x, estados.tau});
fprintf('Modelo no lineal convertido a funcion numerica.\n\n');

% 6. CONFIGURACION GENERAL DE SIMULACION
t_final = 5;
opciones = odeset( 'RelTol',1e-7, 'AbsTol',1e-9, 'Events',@evento_caida);


% 7. CRITERIOS DE RECUPERACION
%El sistema se considera recuperado si durante los
%ultimos 0.5 s permanece aproximadamente en equilibrio.
theta_tol_deg = 1.0;
theta_tol = theta_tol_deg*pi/180;
omega_tol = 0.05;      % rad/s
v_tol     = 0.02;      % m/s
dr_tol    = 0.01;      % m
ventana_final = 0.5;
fprintf('Criterios de recuperacion:\n\n');
fprintf('|theta|   <= %.3f grados\n',theta_tol_deg);
fprintf('|omega_p| <= %.3f rad/s\n',omega_tol);
fprintf('|v|       <= %.3f m/s\n',v_tol);
fprintf('|d_r|     <= %.3f m\n',dr_tol);
fprintf('durante los ultimos %.2f s.\n\n', ventana_final);

% 8. CASO DE REFERENCIA
%Misma perturbacion utilizada en Checkpoint 11:
% theta0 = 0.1 rad = 5.7296 grados
theta_ref = 0.1;
x0_ref = [ 0; 0; theta_ref; 0];

fprintf('---------------------------------------------\n');
fprintf('CASO DE REFERENCIA NO LINEAL\n');
fprintf('---------------------------------------------\n\n');
fprintf('theta0 = %.6f rad\n',theta_ref);
fprintf('theta0 = %.6f grados\n\n', theta_ref*180/pi);

[t_ref,x_ref] = ode45( @(t,x) dinamica_controlada(t,x,f_num,K), ...
    [0 t_final], ...
    x0_ref, ...
    opciones);

% Separar estados
dr_ref    = x_ref(:,1);
v_ref     = x_ref(:,2);
theta_ref_t = x_ref(:,3);
omega_ref = x_ref(:,4);

% Calcular par
tau_ref = -(x_ref*K');

% Metricas
max_dr_ref = max(abs(dr_ref));
max_v_ref = max(abs(v_ref));
max_theta_ref = max(abs(theta_ref_t));
max_omega_ref = max(abs(omega_ref));
max_tau_ref = max(abs(tau_ref));

% Determinar recuperacion

recupera_ref = evaluar_recuperacion( ...
    t_ref, x_ref, t_final, ventana_final, ...
    theta_tol, omega_tol, v_tol, dr_tol);
fprintf('Resultado del caso de referencia:\n\n');

if recupera_ref
    fprintf('RECUPERA EL EQUILIBRIO.\n\n');
else
    fprintf('NO RECUPERA EL EQUILIBRIO.\n\n');
end

fprintf('Maximo |d_r| = %.6f m\n', max_dr_ref);
fprintf('Maximo |v| = %.6f m/s\n', max_v_ref);
fprintf('Maximo |theta| = %.6f grados\n', max_theta_ref*180/pi);
fprintf('Maximo |omega_p| = %.6f rad/s\n', max_omega_ref);
fprintf('Maximo |tau|     = %.6f N*m\n\n', max_tau_ref);

% 9. GRAFICAS DEL CASO DE REFERENCIA
figure;
plot( t_ref, theta_ref_t*180/pi, 'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('\theta (grados)');
title('Modelo no lineal - Angulo de inclinacion');
figure;
plot( t_ref, tau_ref, 'LineWidth',1.3);
grid on;
xlabel('Tiempo (s)');
ylabel('\tau (N*m)');
title('Modelo no lineal - Par solicitado');

% 10. BARRIDO DE ANGULOS INICIALES
%Se exploran inicialmente:-60 grados <= theta0 <= 60 grados
%con pasos de 1 grado.
%Si todos los casos recuperan, el resultado NO significa que
%60 grados sea el limite. Significa que debe ampliarse el rango.

theta0_deg_vec = -60:1:60;
n_pruebas = length(theta0_deg_vec);

% Reservar memoria
recupera = false(n_pruebas,1);
max_tau = NaN(n_pruebas,1);
max_dr = NaN(n_pruebas,1);
max_v = NaN(n_pruebas,1);
max_omega = NaN(n_pruebas,1);
max_theta = NaN(n_pruebas,1);
ts_theta = NaN(n_pruebas,1);
tiempo_final_real = NaN(n_pruebas,1);
fprintf('---------------------------------------------\n');
fprintf('BARRIDO DE ANGULOS INICIALES\n');
fprintf('---------------------------------------------\n\n');
fprintf('Rango: %.1f a %.1f grados\n', ...
    theta0_deg_vec(1), ...
    theta0_deg_vec(end));
fprintf('Paso: %.1f grado\n\n', ...
    theta0_deg_vec(2)-theta0_deg_vec(1));


% 11. SIMULAR CADA ANGULO
for i = 1:n_pruebas
    theta0_deg = theta0_deg_vec(i);
    theta0 = theta0_deg*pi/180;
    % Condicion inicial
    x0 = [ 0; 0; theta0; 0];
    try
        [t,x] = ode45( ...
            @(t,x) dinamica_controlada( ...
                t,x,f_num,K), ...
            [0 t_final], ...
            x0, ...
            opciones);

        % Tiempo alcanzado
        tiempo_final_real(i) = t(end);


        % Estados
        dr = x(:,1);
        v = x(:,2);
        theta = x(:,3);
        omega_p = x(:,4);
        
        % Par solicitado
        tau = -(x*K');
        
        % Metricas
        max_tau(i) = max(abs(tau));
        max_dr(i) = max(abs(dr));
        max_v(i) = max(abs(v));
        max_omega(i) = max(abs(omega_p));
        max_theta(i) = ...
            max(abs(theta))*180/pi;
        % Comprobar recuperacion
        recupera(i) = evaluar_recuperacion( ...
            t, ...
            x, ...
            t_final, ...
            ventana_final, ...
            theta_tol, ...
            omega_tol, ...
            v_tol, ...
            dr_tol);
        % Tiempo de establecimiento angular
        ts_theta(i) = tiempo_establecimiento( ...
            t, ...
            theta, ...
            theta_tol, ...
            t_final);
    catch
        recupera(i) = false;
    end
end


% 12. CREAR TABLA DE RESULTADOS
tabla = table( ...
    theta0_deg_vec(:), ...
    recupera, ...
    max_tau, ...
    max_dr, ...
    max_v, ...
    max_omega, ...
    max_theta, ...
    ts_theta, ...
    tiempo_final_real, ...
    'VariableNames',{ ...
    'Theta0_deg', ...
    'Recupera', ...
    'TauMax_Nm', ...
    'DrMax_m', ...
    'VMax_m_s', ...
    'OmegaMax_rad_s', ...
    'ThetaMax_deg', ...
    'TsTheta_s', ...
    'TiempoSimulado_s'});

% 13. DETERMINAR LIMITES POSITIVO Y NEGATIVO

idx_pos = theta0_deg_vec >= 0;
theta_pos = theta0_deg_vec(idx_pos);
rec_pos = recupera(idx_pos);
idx_neg = theta0_deg_vec <= 0;
theta_neg = theta0_deg_vec(idx_neg);
rec_neg = recupera(idx_neg);

% Maximo positivo recuperado

if any(rec_pos)
    theta_max_pos = max(theta_pos(rec_pos));
else
    theta_max_pos = NaN;
end

% Maximo negativo recuperado en magnitud
if any(rec_neg)
theta_min_neg = min(theta_neg(rec_neg));
else
    theta_min_neg = NaN;
end

% 14. RANGO BILATERAL CONTIGUA
%Se busca el mayor valor A tal que:
%   -A <= theta0 <= A
%se recupere para todos los angulos probados.
paso_deg = abs( theta0_deg_vec(2)-theta0_deg_vec(1));
magnitudes = 0:paso_deg: max(abs(theta0_deg_vec));

limite_bilateral = 0;

for mag = magnitudes
    idx_p = find( ...
        abs(theta0_deg_vec-mag) < 1e-10, 1);
    idx_n = find( abs(theta0_deg_vec+mag) < 1e-10, 1);
    if isempty(idx_p) || isempty(idx_n)
        break;
    end
    if recupera(idx_p) && recupera(idx_n)
        limite_bilateral = mag;
    else
        break;
    end
end


% 15. MOSTRAR RESULTADOS

fprintf('---------------------------------------------\n');
fprintf('RESULTADO DEL BARRIDO\n');
fprintf('---------------------------------------------\n\n');

fprintf('Pruebas realizadas: %d\n', ...
    n_pruebas);

fprintf('Casos recuperados: %d\n\n', ...
    sum(recupera));

fprintf('Mayor angulo positivo recuperado:\n');
fprintf('%.2f grados\n\n', ...
    theta_max_pos);

fprintf('Mayor angulo negativo recuperado:\n');
fprintf('%.2f grados\n\n', ...
    theta_min_neg);

fprintf('Region bilateral continua verificada:\n');
fprintf('[-%.2f , +%.2f] grados\n\n', ...
    limite_bilateral, ...
    limite_bilateral);

% Avisar si se alcanzo el borde del barrido
if limite_bilateral >= ...
        max(abs(theta0_deg_vec))
    fprintf(['ADVERTENCIA:\n' ...
        'Todos los angulos del rango explorado fueron ' ...
        'recuperados.\n']);
    fprintf(['El limite real no fue encontrado; ' ...
        'se debe ampliar el barrido.\n\n']);
end


% 16. GRAFICA DE REGION DE RECUPERACION
figure;
plot( theta0_deg_vec, double(recupera), 'o-');
grid on;
xlabel('\theta_0 (grados)');
ylabel('Recuperacion');
ylim([-0.1 1.1]);
title('Region de recuperacion - Modelo no lineal');


% 17. PAR MAXIMO REQUERIDO VS ANGULO INICIAL
figure;
plot( theta0_deg_vec, max_tau, 'LineWidth',1.3);
grid on;
xlabel('\theta_0 (grados)');
ylabel('|\tau|_{max} (N*m)');
title('Par requerido vs. inclinacion inicial');

% 18. DESPLAZAMIENTO MAXIMO VS ANGULO INICIAL
figure;
plot( theta0_deg_vec, max_dr, 'LineWidth',1.3);
grid on;
xlabel('\theta_0 (grados)');
ylabel('|d_r|_{max} (m)');
title('Desplazamiento maximo vs. inclinacion inicial');

% 19. GUARDAR RESULTADOS

region.K = K;
region.Kp_theta = Kp_theta;
region.Kd_theta = Kd_theta;
region.Kp_v = Kp_v;
region.Ki_v = Ki_v;
region.f_parametros = f_parametros;


% Caso de referencia
region.referencia.theta0 = 0.1;
region.referencia.t = t_ref;
region.referencia.x = x_ref;
region.referencia.tau = tau_ref;
region.referencia.recupera = recupera_ref;
region.referencia.max_dr = max_dr_ref;
region.referencia.max_v = max_v_ref;
region.referencia.max_theta = max_theta_ref;
region.referencia.max_omega = max_omega_ref;
region.referencia.max_tau = max_tau_ref;

% Barrido
region.theta0_deg = theta0_deg_vec;
region.recupera = recupera;
region.max_tau = max_tau;
region.max_dr = max_dr;
region.max_v = max_v;
region.max_omega = max_omega;
region.max_theta = max_theta;
region.ts_theta = ts_theta;
region.tabla = tabla;


% Limites
region.theta_max_pos = theta_max_pos;
region.theta_min_neg = theta_min_neg;
region.limite_bilateral = limite_bilateral;

% Criterios
region.theta_tol = theta_tol;
region.omega_tol = omega_tol;
region.v_tol = v_tol;
region.dr_tol = dr_tol;
region.ventana_final = ventana_final;
fprintf('=============================================\n');

end

function dx = dinamica_controlada(~,x,f_num,K)
%Dinamica no lineal en lazo cerrado.
%Control PD + PI de velocidad equivalente
tau = -K*x;
%Modelo no lineal procedente del Checkpoint 3
dx = f_num(x,tau);
%Garantizar vector columna
dx = dx(:);
end

function recupera = evaluar_recuperacion( ...
    t,x,t_final,ventana_final, ...
    theta_tol,omega_tol,v_tol,dr_tol)
%Determina si el sistema recupero el equilibrio.
% Si la simulacion termino antes, se considera fallo
if t(end) < t_final-1e-4
    recupera = false;
    return;
end

% Revisar valores numericos
if any(~isfinite(x(:)))
    recupera = false;
    return;
end


% Ultima ventana temporal
idx = t >= (t_final-ventana_final);
if ~any(idx)
    recupera = false;
    return;
end

dr = x(idx,1);
v = x(idx,2);
theta = x(idx,3);
omega_p = x(idx,4);

% Criterios simultaneos
condicion_theta = all(abs(theta) <= theta_tol);

condicion_omega = all(abs(omega_p) <= omega_tol);

condicion_v = all(abs(v) <= v_tol);

condicion_dr = all(abs(dr) <= dr_tol);

recupera = condicion_theta && condicion_omega && ...
    condicion_v && condicion_dr;

end

function ts = tiempo_establecimiento( ...
    t,theta,theta_tol,t_final)
%Tiempo a partir del cual theta permanece dentro
%de la banda establecida.
if t(end) < t_final-1e-4
    ts = NaN;
    return;
end

idx_fuera = find(abs(theta) > theta_tol);
if isempty(idx_fuera)
    ts = 0;
elseif idx_fuera(end) < length(t)
    ts = t(idx_fuera(end)+1);
else
    ts = NaN;
end
end


% ================================================================
function [value,isterminal,direction] = evento_caida(~,x)
%Detiene la simulacion si el cuerpo alcanza +/-90 grados.
%Esto evita continuar simulando un estado que ya no
%representa una recuperacion desde la posicion vertical.
theta = x(3);
value = pi/2 - abs(theta);
isterminal = 1;
direction = -1;
end

function expr_num = sustituir_parametros(expr,estados,p)
%Sustituye automaticamente los parametros simbolicos
%presentes en estados.f utilizando la estructura p.
%
%No sustituye:
%
% dr, v, theta, omega_p, tau

vars = symvar(expr);

% Nombres de estados y entrada
excluir = cell(1,length(estados.x)+1);
for i = 1:length(estados.x)
    excluir{i} = char(estados.x(i));
end
excluir{end} = char(estados.tau);

% Preparar sustituciones
simbolos = sym([]);
valores = [];

for i = 1:length(vars)
    nombre = char(vars(i));
    % Ignorar estados y tau
    if any(strcmp(nombre,excluir))
        continue;
    end
    % Buscar valor en estructura de parametros
    [encontrado,valor] = ...
        buscar_parametro(p,nombre);
    if ~encontrado
        error(['No se encontro el parametro simbolico "' ...
               nombre ...
               '" dentro de la estructura proporcionada.']);
    end
    simbolos(end+1) = vars(i);
    valores(end+1) = valor;

end

% Sustituir
expr_num = subs( expr, simbolos, valores);
% Verificar que no queden parametros sin sustituir
vars_restantes = symvar(expr_num);
permitidas = excluir;


for i = 1:length(vars_restantes)
    nombre = char(vars_restantes(i));
    if ~any(strcmp(nombre,permitidas))

        error(['Quedo sin sustituir el simbolo "' nombre '".']);
    end
end
end


function [encontrado,valor] = buscar_parametro(p,nombre)
%Busca un parametro permitiendo pequeñas diferencias
%de nomenclatura, por ejemplo:
% m_h  <-> mh
% I_p  <-> Ip
% c_v  <-> cv
% c_omega <-> cw

campos = fieldnames(p);
objetivo = normalizar(nombre);
encontrado = false;
valor = NaN;

% Busqueda normal
for i = 1:length(campos)
    campo_normal = normalizar(campos{i});
    if strcmp(campo_normal,objetivo)
        valor = p.(campos{i});
        encontrado = true;
        return;
    end
end

% Alias comunes
if strcmp(objetivo,'comega')

    alias = { ...
        'cw', 'comega', 'comega'};

elseif strcmp(objetivo,'cw')
    alias = { 'cw', 'comega', 'comega'};
else
    alias = {};
end

for j = 1:length(alias)
    for i = 1:length(campos)
        campo_normal = normalizar(campos{i});
        if strcmp(campo_normal, alias{j})
            valor = p.(campos{i});
            encontrado = true;
            return;
        end
    end
end

end

function nombre = normalizar(nombre)
%Elimina guiones bajos y diferencias entre mayusculas
%y minusculas.

nombre = lower(nombre);
nombre = regexprep(nombre, '[^a-zA-Z0-9]','');

end
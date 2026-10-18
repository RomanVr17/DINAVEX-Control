function analisis = checkpoint7_analisis_lazo_abierto(ft, numerico)
%CHECKPOINT 7 - ANALISIS EN LAZO ABIERTO
%Analiza polos, ceros y estabilidad de las cuatro FT 
%Posteriormente compara los resultados con los valores
%propios de A y simula los cuatro estados ante una condicion inicial.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 7 - ANALISIS EN LAZO ABIERTO\n');
fprintf('========================================\n\n');

% Recuperar funciones de transferencia
G_dr      = ft.G_dr;
G_v       = ft.G_v;
G_theta   = ft.G_theta;
G_omega_p = ft.G_omega_p;

% Simplificar funciones de transferencia
G_dr_min      = minreal(G_dr);
G_v_min       = minreal(G_v);
G_theta_min   = minreal(G_theta);
G_omega_p_min = minreal(G_omega_p);

% Organizar funciones de transferencia
G = {G_dr_min, G_v_min, G_theta_min, G_omega_p_min};
nombres = {
    'tau -> d_r'
    'tau -> v'
    'tau -> theta'
    'tau -> omega_p'
};

analisis_ft = cell(4,1);

% Analisis de polos, ceros y estabilidad
tol = 1e-8;

for i = 1:4

    fprintf('----------------------------------------\n');
    fprintf('FT %d: %s\n', i, nombres{i});
    fprintf('----------------------------------------\n');

    polos = pole(G{i});
    ceros = zero(G{i});

    disp('Polos:');
    disp(polos);

    disp('Ceros:');
    disp(ceros);

    if any(real(polos) > tol)
        estado = 'INESTABLE';
    elseif all(real(polos) < -tol)
        estado = 'ESTABLE';
    else
        estado = 'NO ASINTOTICAMENTE ESTABLE';
    end

    fprintf('Clasificacion: %s\n\n', estado);

    analisis_ft{i}.G = G{i};
    analisis_ft{i}.polos = polos;
    analisis_ft{i}.ceros = ceros;
    analisis_ft{i}.estado = estado;
end

% Verificacion con el modelo en espacio de estados
A = numerico.A;
B = numerico.B;

valores_propios = eig(A);

fprintf('\n========================================\n');
fprintf('VERIFICACION CON EL MODELO DE ESTADOS\n');
fprintf('========================================\n\n');

disp('Valores propios de la matriz A:');
disp(valores_propios);

% Modelo de estados para simulacion
C = eye(4);
D = zeros(4,1);
sys_ss = ss(A, B, C, D);

% Condicion inicial
% x = [d_r; v; theta; omega_p]
x0 = [0; 0; 0.1; 0];

% Tiempo de simulacion
t = 0:0.001:1;

% Simulacion en lazo abierto
[~, t, x] = initial(sys_ss, x0, t);

% Graficas de los cuatro estados
figure;

subplot(2,2,1);
plot(t, x(:,1), 'LineWidth', 1.2);
grid on;
xlabel('Tiempo (s)');
ylabel('d_r (m)');
title('Desplazamiento');

subplot(2,2,2);
plot(t, x(:,2), 'LineWidth', 1.2);
grid on;
xlabel('Tiempo (s)');
ylabel('v (m/s)');
title('Velocidad lineal');

subplot(2,2,3);
plot(t, x(:,3), 'LineWidth', 1.2);
grid on;
xlabel('Tiempo (s)');
ylabel('\theta (rad)');
title('Angulo de inclinacion');

subplot(2,2,4);
plot(t, x(:,4), 'LineWidth', 1.2);
grid on;
xlabel('Tiempo (s)');
ylabel('\omega_p (rad/s)');
title('Velocidad angular');


% Guardar resultados
analisis.ft = analisis_ft;
analisis.valores_propios = valores_propios;
analisis.sys_ss = sys_ss;
analisis.x0 = x0;
analisis.t = t;
analisis.x = x;
fprintf('\n========================================\n');

end
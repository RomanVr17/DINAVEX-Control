function conjunto = checkpoint9_PD_PI_conjunto(ft)
%CHECKPOINT 9 - ANALISIS CONJUNTO PD-PI
%Analiza la arquitectura paralela propuesta:
%
%   u = u_theta + u_dr
%
%   u_theta = Kp_theta*e_theta + Kd_theta*e_omega
%   u_dr    = Kp_dr*e_dr + Ki_dr*integral(e_dr)
%
%Para referencias nulas:
%   e_theta = -theta
%   e_omega = -omega_p
%   e_dr    = -d_r
%
%El analisis se realiza a partir de las funciones de transferencia.
%Se construye el polinomio caracteristico completo, se revisan
%condiciones necesarias de estabilidad y se realiza una busqueda
%exploratoria de combinaciones de ganancias.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 9 - ANALISIS CONJUNTO PD-PI\n');
fprintf('========================================\n\n');

% Recuperar funciones de transferencia
G_dr      = minreal(ft.G_dr);
G_theta   = minreal(ft.G_theta);
G_omega_p = minreal(ft.G_omega_p);

fprintf('Arquitectura analizada:\n');
fprintf('u = u_theta + u_dr\n');
fprintf('PD para equilibrio + PI para desplazamiento\n\n');

% Tolerancia numerica
tol = 1e-7;

% Verificar G_omega_p(s) = s*G_theta(s)
s = tf([1 0],1);

error_omega = minreal(G_omega_p - s*G_theta);

[num_error_omega, ~] = tfdata(error_omega,'v');
num_error_omega = limpiar_inicio(num_error_omega,tol);

if all(abs(num_error_omega) < tol)

    relacion_angular = true;

    fprintf('VERIFICACION 1: G_omega_p(s) = s*G_theta(s)\n');

else

    relacion_angular = false;

    fprintf(['ADVERTENCIA: G_omega_p(s) no coincide con ' ...
             's*G_theta(s).\n']);

end

% Extraer polinomios de G_theta y G_dr
[Ntheta,Dtheta] = tfdata(G_theta,'v');
[Ndr,Ddr]       = tfdata(G_dr,'v');

Ntheta = limpiar_inicio(Ntheta,tol);
Dtheta = limpiar_inicio(Dtheta,tol);

Ndr = limpiar_inicio(Ndr,tol);
Ddr = limpiar_inicio(Ddr,tol);

%% Normalizar denominadores
Ntheta = Ntheta/Dtheta(1);
Dtheta = Dtheta/Dtheta(1);

Ndr = Ndr/Ddr(1);
Ddr = Ddr/Ddr(1);

% Verificar que D_dr(s) = s*D_theta(s)
Ddr_esperado = conv([1 0],Dtheta);

Ddr_esperado = Ddr_esperado/Ddr_esperado(1);

if length(Ddr) == length(Ddr_esperado) && ...
        max(abs(Ddr-Ddr_esperado)) < 1e-5

    relacion_desplazamiento = true;

    fprintf('VERIFICACION 2: D_dr(s) = s*D_theta(s)\n\n');

else

    relacion_desplazamiento = false;

    fprintf(['ADVERTENCIA: El denominador de G_dr(s) no coincide ' ...
             'con s*D_theta(s).\n']);

    fprintf(['El desarrollo algebraico automatico requiere revisar ' ...
             'la estructura de las FT.\n\n']);

end

if ~relacion_angular || ~relacion_desplazamiento

    error(['Las relaciones estructurales esperadas entre ' ...
           'las FT no se cumplen.']);

end

% Construccion del polinomio caracteristico
%
% Ecuacion caracteristica:
%
% 1 + Kp_theta*G_theta(s)
%   + Kd_theta*G_omega_p(s)
%   + (Kp_dr + Ki_dr/s)*G_dr(s) = 0
%
% Usando:
%
% G_omega_p = s*Ntheta/Dtheta
%
% G_dr = Ndr/(s*Dtheta)
%
% Al multiplicar por s^2*Dtheta:
%
% P(s) =
% s^2*Dtheta
% + Kp_theta*s^2*Ntheta
% + Kd_theta*s^3*Ntheta
% + Kp_dr*s*Ndr
% + Ki_dr*Ndr
%

P0       = conv([1 0 0],Dtheta);
P_Kpt    = conv([1 0 0],Ntheta);
P_Kdt    = conv([1 0 0 0],Ntheta);
P_Kpdr   = conv([1 0],Ndr);
P_Kidr   = Ndr;

% Igualar longitud de todos los polinomios
ncoef = max([ ...
    length(P0), ...
    length(P_Kpt), ...
    length(P_Kdt), ...
    length(P_Kpdr), ...
    length(P_Kidr)]);

P0     = pad_left(P0,ncoef);
P_Kpt  = pad_left(P_Kpt,ncoef);
P_Kdt  = pad_left(P_Kdt,ncoef);
P_Kpdr = pad_left(P_Kpdr,ncoef);
P_Kidr = pad_left(P_Kidr,ncoef);

% Mostrar polinomio general
fprintf('----------------------------------------\n');
fprintf('POLINOMIO CARACTERISTICO GENERAL\n');
fprintf('----------------------------------------\n\n');

fprintf(['P(s) = s^2*D(s) + Kp_theta*s^2*N_theta(s) ' ...
         '+ Kd_theta*s^3*N_theta(s)\n']);

fprintf(['       + Kp_dr*s*N_dr(s) ' ...
         '+ Ki_dr*N_dr(s)\n\n']);

fprintf('Coeficientes desde s^%d hasta s^0:\n\n',ncoef-1);

for i = 1:ncoef

    potencia = ncoef-i;

    fprintf(['a%d = %.6g + (%.6g)Kp_theta ' ...
             '+ (%.6g)Kd_theta ' ...
             '+ (%.6g)Kp_dr ' ...
             '+ (%.6g)Ki_dr\n'], ...
             potencia, ...
             P0(i), ...
             P_Kpt(i), ...
             P_Kdt(i), ...
             P_Kpdr(i), ...
             P_Kidr(i));

end

% Condicion necesaria de signo
%
% Para un polinomio monico Hurwitz todos sus coeficientes
% deben ser positivos.
%
% Esto NO es suficiente para orden 5, pero permite descartar
% regiones de ganancias inmediatamente.

fprintf('\n----------------------------------------\n');
fprintf('REVISION DE GANANCIAS POSITIVAS\n');
fprintf('----------------------------------------\n\n');

coef_constante_Ki = P_Kidr(end);

otros_constantes = [ ...
    P0(end) ...
    P_Kpt(end) ...
    P_Kdt(end) ...
    P_Kpdr(end)];

if max(abs(otros_constantes)) < tol && ...
        coef_constante_Ki < 0

    ganancias_positivas_viables = false;

    fprintf(['Con la convencion de signos actual, el termino ' ...
             'constante queda proporcional a Ki_dr con ' ...
             'coeficiente negativo.\n']);

    fprintf(['Para Ki_dr > 0, el termino constante es ' ...
             'negativo.\n']);

    fprintf(['Por tanto, la arquitectura con todas las ganancias ' ...
             'positivas NO puede ser Hurwitz.\n\n']);

else

    ganancias_positivas_viables = true;

    fprintf(['El signo del termino constante no descarta por si ' ...
             'solo las ganancias positivas.\n\n']);

end

% Busqueda exploratoria con ganancias con signo
%
% Esta busqueda NO sustituye una demostracion algebraica.
%
% Se permiten ganancias positivas y negativas porque el signo
% requerido depende tambien de la convencion utilizada para
% desplazamiento, inclinacion y par.

fprintf('----------------------------------------\n');
fprintf('BUSQUEDA EXPLORATORIA DE ESTABILIDAD\n');
fprintf('----------------------------------------\n\n');

%% Rangos iniciales
Kp_theta_vec = linspace(-10,10,9);

Kd_theta_vec = linspace(-1,1,9);

Kp_dr_vec = linspace(-10,10,9);

Ki_dr_vec = linspace(-10,10,9);

% Variables para guardar el mejor resultado
mejor_margen = Inf;

mejores_ganancias = [NaN NaN NaN NaN];

mejores_polos = [];

numero_estables = 0;

numero_evaluadas = 0;

tol_estabilidad = 1e-6;

% Barrido de ganancias
for i1 = 1:length(Kp_theta_vec)

    Kp_theta = Kp_theta_vec(i1);

    for i2 = 1:length(Kd_theta_vec)

        Kd_theta = Kd_theta_vec(i2);

        for i3 = 1:length(Kp_dr_vec)

            Kp_dr = Kp_dr_vec(i3);

            for i4 = 1:length(Ki_dr_vec)

                Ki_dr = Ki_dr_vec(i4);

                %% Construir polinomio para estas ganancias
                P = P0 ...
                    + Kp_theta*P_Kpt ...
                    + Kd_theta*P_Kdt ...
                    + Kp_dr*P_Kpdr ...
                    + Ki_dr*P_Kidr;

                %% Obtener polos
                polos = roots(P);

                %% Polo mas a la derecha
                max_real = max(real(polos));

                numero_evaluadas = numero_evaluadas + 1;

                %% Comprobar estabilidad
                if max_real < -tol_estabilidad

                    numero_estables = numero_estables + 1;

                    %% Guardar mejor candidato
                    if max_real < mejor_margen

                        mejor_margen = max_real;

                        mejores_ganancias = [ ...
                            Kp_theta ...
                            Kd_theta ...
                            Kp_dr ...
                            Ki_dr];

                        mejores_polos = polos;

                    end

                end

            end

        end

    end

end

% Mostrar resultados de la busqueda
fprintf('Combinaciones evaluadas: %d\n', ...
    numero_evaluadas);

fprintf('Combinaciones estables encontradas: %d\n\n', ...
    numero_estables);

if numero_estables > 0

    fprintf('Mejor candidato del barrido inicial:\n');

    fprintf('Kp_theta = %.6f\n', ...
        mejores_ganancias(1));

    fprintf('Kd_theta = %.6f\n', ...
        mejores_ganancias(2));

    fprintf('Kp_dr    = %.6f\n', ...
        mejores_ganancias(3));

    fprintf('Ki_dr    = %.6f\n\n', ...
        mejores_ganancias(4));

    disp('Polos del candidato:');
    disp(mejores_polos);

    fprintf('Polo mas a la derecha: %.6f\n\n', ...
        mejor_margen);

    % Grafica de polos del candidato
    figure;

    plot( ...
        real(mejores_polos), ...
        imag(mejores_polos), ...
        'x', ...
        'MarkerSize',9);

    grid on;
    hold on;

    yl = ylim;

    plot([0 0],yl,'--');

    xlabel('Parte real');
    ylabel('Parte imaginaria');

    title('Polos de lazo cerrado - candidato PD-PI');

    hold off;

    candidato_encontrado = true;

else

    fprintf(['No se encontro una combinacion estable dentro del ' ...
             'barrido inicial definido.\n']);

    fprintf(['Esto no constituye por si solo una demostracion de ' ...
             'imposibilidad para ganancias con signo arbitrario.\n']);

    fprintf(['Sin embargo, junto con la condicion necesaria de los ' ...
             'coeficientes, indica que conviene revisar la ' ...
             'arquitectura del lazo de desplazamiento.\n\n']);

    candidato_encontrado = false;

end

% Guardar resultados
conjunto.G_dr = G_dr;

conjunto.G_theta = G_theta;

conjunto.G_omega_p = G_omega_p;

conjunto.relacion_angular = relacion_angular;

conjunto.relacion_desplazamiento = ...
    relacion_desplazamiento;

conjunto.P0 = P0;

conjunto.P_Kp_theta = P_Kpt;

conjunto.P_Kd_theta = P_Kdt;

conjunto.P_Kp_dr = P_Kpdr;

conjunto.P_Ki_dr = P_Kidr;

conjunto.ganancias_positivas_viables = ...
    ganancias_positivas_viables;

conjunto.numero_evaluadas = ...
    numero_evaluadas;

conjunto.numero_estables = ...
    numero_estables;

conjunto.candidato_encontrado = ...
    candidato_encontrado;

conjunto.ganancias = ...
    mejores_ganancias;

conjunto.polos = ...
    mejores_polos;

conjunto.max_real = ...
    mejor_margen;

fprintf('========================================\n');

end


% ================================================================
function v = limpiar_inicio(v,tol)
%Elimina coeficientes iniciales numericamente nulos.

while length(v) > 1 && abs(v(1)) < tol

    v(1) = [];

end

end


% ================================================================
function v = pad_left(v,n)
%Completa con ceros a la izquierda hasta alcanzar longitud n.

if length(v) < n

    v = [zeros(1,n-length(v)) v];

end

end
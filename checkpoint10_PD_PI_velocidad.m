function control = checkpoint10_PD_PI_velocidad(ft)
%CHECKPOINT 10 - PD DE EQUILIBRIO + PI DE VELOCIDAD
%Analiza la arquitectura:
%   u = u_theta + u_v
%   u_theta = Kp_theta*e_theta + Kd_theta*e_omega
%   u_v     = Kp_v*e_v + Ki_v*integral(e_v)
%Se obtiene el polinomio caracteristico completo y se realiza una busqueda exploratoria de ganancias que estabilicen ...
% los cuatro polos del sistema.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 10 - PD + PI DE VELOCIDAD\n');
fprintf('========================================\n\n');

% Recuperar funciones de transferencia
G_dr      = minreal(ft.G_dr);
G_v       = minreal(ft.G_v);
G_theta   = minreal(ft.G_theta);
G_omega_p = minreal(ft.G_omega_p);

tol = 1e-7;

%1.VERIFICAR RELACIONES ENTRE LAS FUNCIONES DE TRANSFERENCIA

s = tf([1 0],1);

% Relacion angular
% omega_p = d(theta)/dt
% Por tanto:
% G_omega_p(s) = s*G_theta(s)

error_angular = minreal(G_omega_p - s*G_theta);
[num_error,~] = tfdata(error_angular,'v');       %extrae los coeficientes del numerador.
num_error = limpiar_inicio(num_error,tol);       %elimina ceros numéricos insignificantes al comienzo del vector.

if all(abs(num_error) < tol) %pregunta si todos los coeficientes son prácticamente cero. Si lo son, confirma la relación.
    relacion_angular = true;
    fprintf('VERIFICACION 1:\n');
    fprintf('G_omega_p(s) = s*G_theta(s)\n\n');
else
    relacion_angular = false;
    fprintf('ADVERTENCIA:\n');
    fprintf('G_omega_p(s) no coincide con s*G_theta(s).\n\n');
end

% Relacion lineal
% v = d(d_r)/dt
% Por tanto:
% G_v(s) = s*G_dr(s)

error_velocidad = minreal(G_v - s*G_dr);
[num_error_v,~] = tfdata(error_velocidad,'v');
num_error_v = limpiar_inicio(num_error_v,tol);

if all(abs(num_error_v) < tol)
    relacion_velocidad = true;
    fprintf('VERIFICACION 2:\n');
    fprintf('G_v(s) = s*G_dr(s)\n\n');
else
    relacion_velocidad = false;
    fprintf('ADVERTENCIA:\n');
    fprintf('G_v(s) no coincide con s*G_dr(s).\n\n');
end

%2.OBTENER POLINOMIOS DE G_theta Y G_v

[Ntheta,Dtheta] = tfdata(G_theta,'v');   %tfdata separa los coeficientes del numerador y denominador.
[Nv,Dv]         = tfdata(G_v,'v');

Ntheta = limpiar_inicio(Ntheta,tol);
Dtheta = limpiar_inicio(Dtheta,tol);

Nv = limpiar_inicio(Nv,tol);
Dv = limpiar_inicio(Dv,tol);

% Normalizar denominadores
Ntheta = Ntheta/Dtheta(1);
Dtheta = Dtheta/Dtheta(1);

Nv = Nv/Dv(1);
Dv = Dv/Dv(1);

% Verificar denominador comun
if length(Dtheta) == length(Dv) && max(abs(Dtheta-Dv)) < 1e-5
    denominador_comun = true;
    fprintf('VERIFICACION 3:\n');
    fprintf(['G_theta(s) y G_v(s) comparten el mismo ' 'denominador dinamico.\n\n']);
else
    denominador_comun = false;
    error(['G_theta(s) y G_v(s) no presentan el denominador ' 'comun esperado.']);
end

%3.CONSTRUIR POLINOMIO CARACTERISTICO
% Ecuacion caracteristica:
% 1 + Kp_theta*G_theta + Kd_theta*G_omega_p + (Kp_v + Ki_v/s)*G_v = 0
% Como:
% G_omega_p = s*Ntheta/D
% G_v = Nv/D
% multiplicando por s*D:
% P(s) = s*D + Kp_theta*s*Ntheta + Kd_theta*s^2*Ntheta + Kp_v*s*Nv + Ki_v*Nv
% El resultado es un polinomio de cuarto orden.

%La convolución de vectores de coeficientes equivale a multiplicar polinomios.
P0        = conv([1 0],Dtheta);
P_Kptheta = conv([1 0],Ntheta);
P_Kdtheta = conv([1 0 0],Ntheta);
P_Kpv     = conv([1 0],Nv);
P_Kiv     = Nv;

% Igualar longitud de todos los vectores, alinear correctamente las potencias de s antes de sumar los polinomios
ncoef = max([ length(P0), length(P_Kptheta), length(P_Kdtheta), length(P_Kpv), length(P_Kiv)]);

P0        = pad_left(P0,ncoef);
P_Kptheta = pad_left(P_Kptheta,ncoef);
P_Kdtheta = pad_left(P_Kdtheta,ncoef);
P_Kpv     = pad_left(P_Kpv,ncoef);
P_Kiv     = pad_left(P_Kiv,ncoef);

%4.MOSTRAR POLINOMIO GENERAL
fprintf('----------------------------------------\n');
fprintf('POLINOMIO CARACTERISTICO GENERAL\n');
fprintf('----------------------------------------\n\n');
fprintf('P(s) = s*D(s) + Kp_theta*s*N_theta(s) + Kd_theta*s^2*N_theta(s)+ Kp_v*s*N_v(s) + Ki_v*N_v(s)\n\n');
fprintf('Orden del polinomio: %d\n\n',ncoef-1);
fprintf('Coeficientes desde s^%d hasta s^0:\n\n',ncoef-1);

for i = 1:ncoef
    potencia = ncoef-i;
    fprintf(['a%d = %.6g + (%.6g)Kp_theta ' '+ (%.6g)Kd_theta ' ...
             '+ (%.6g)Kp_v ' '+ (%.6g)Ki_v\n'], potencia, P0(i), P_Kptheta(i), ...
             P_Kdtheta(i), P_Kpv(i), P_Kiv(i));
end

%5.REVISION DEL TERMINO CONSTANTE
fprintf('\n----------------------------------------\n');
fprintf('REVISION PRELIMINAR DE SIGNOS\n');
fprintf('----------------------------------------\n\n');

coef_Kiv_constante = P_Kiv(end); 
fprintf('Coeficiente de Ki_v en a0: %.6f\n\n', coef_Kiv_constante);
if coef_Kiv_constante < 0
    fprintf(['Con la convencion actual de signos, Ki_v positivo ' ...
             'produce un termino constante negativo.\n']);
    fprintf(['Por ello se deben permitir ganancias con signo ' ...
             'durante la busqueda.\n\n']);

else
    fprintf(['El signo del termino constante no descarta ' ...
             'Ki_v positivo.\n\n']);
end


%6.BUSQUEDA EXPLORATORIA DE GANANCIAS
% Los signos de las ganancias dependen de la convencion utilizada
% para theta, velocidad y par.
% Por ello inicialmente se permite explorar valores positivos
% y negativos.
% Estos rangos son solamente exploratorios.

fprintf('----------------------------------------\n');
fprintf('BUSQUEDA EXPLORATORIA\n');
fprintf('----------------------------------------\n\n');

% Kp_theta_vec = linspace();
% Kd_theta_vec = linspace();
% Kp_v_vec = linspace();
% Ki_v_vec = linspace();

%búsqueda paramétrica discretizada o grid search.
Kp_theta_vec = -120:5:-10;
Kd_theta_vec = -25:2.5:-2.5;
Kp_v_vec = -15:1:-1;
Ki_v_vec = -10:0.5:-0.5;

numero_evaluadas = 0;
numero_estables = 0;
tol_estabilidad = 1e-6;  %evita considerar como estable un polo demasiado cercano a cero.

%Guardar candidatos estables
candidatos = [];

%Escalas para comparar magnitud de ganancias
esc_Kpt = max(abs(Kp_theta_vec));
esc_Kdt = max(abs(Kd_theta_vec));
esc_Kpv = max(abs(Kp_v_vec));
esc_Kiv = max(abs(Ki_v_vec));

%7.BARRIDO DE GANANCIAS  prueba todas las combinaciones posibles de las cuatro ganancias
for i1 = 1:length(Kp_theta_vec)
    Kp_theta = Kp_theta_vec(i1);
    for i2 = 1:length(Kd_theta_vec)
        Kd_theta = Kd_theta_vec(i2);
        for i3 = 1:length(Kp_v_vec)
            Kp_v = Kp_v_vec(i3);
            for i4 = 1:length(Ki_v_vec)
                Ki_v = Ki_v_vec(i4);
                P = P0 ...
                    + Kp_theta*P_Kptheta ...
                    + Kd_theta*P_Kdtheta ...
                    + Kp_v*P_Kpv ...
                    + Ki_v*P_Kiv;
                polos = roots(P);                        %Calcula las raíces del polinomio
                numero_evaluadas = numero_evaluadas + 1; 
                max_real = max(real(polos));             %obtiene la parte real de todos los polos y después escoge la mayor
                if max_real < -tol_estabilidad           %aceptar esta combinación únicamente si todos los polos tienen parte real negativa
                    numero_estables = numero_estables + 1;
                    [~,idx_dom] = max(real(polos));      %idx_dom indica la posición del polo cuya parte real es mayor
                    polo_dom = polos(idx_dom);           %recupera ese polo.
                    wn_c = abs(polo_dom);
                    zeta_c = -real(polo_dom)/wn_c;      
                    ts_c = 4/abs(real(polo_dom));        %criterio del 2%
                    norma_K = sqrt( ...     %Calcula una medida de qué tan grandes son las cuatro ganancias normalizadas por sus respectivos rangos.
                        (Kp_theta/esc_Kpt)^2 + ...
                        (Kd_theta/esc_Kdt)^2 + ...
                        (Kp_v/esc_Kpv)^2 + ...
                        (Ki_v/esc_Kiv)^2);

                    candidatos(end+1,:) = [ ...    %Guardar cada candidato estable
                        Kp_theta ...
                        Kd_theta ...
                        Kp_v ...
                        Ki_v ...
                        max_real ...
                        wn_c ...
                        zeta_c ...
                        ts_c ...
                        norma_K];
                end
            end
        end
    end
end

% 8. RESULTADOS
fprintf('Combinaciones evaluadas: %d\n', numero_evaluadas);
fprintf('Combinaciones estables encontradas: %d\n\n', numero_estables);
if numero_estables > 0
    candidato_encontrado = true;
    % Tiempo mas rapido encontrado
    ts_min = min(candidatos(:,8));    %todas las filas de la columna 8

    % Permitir hasta 25% mas tiempo que el mas rapido
    %Porque el candidato absolutamente más rápido podría requerir ganancias demasiado elevadas.
    factor_tiempo = 1.25;

    idx_buenos = candidatos(:,8) <= factor_tiempo*ts_min;  %crea un vector de candidatos aptos

    candidatos_buenos = candidatos(idx_buenos,:);

    % ordenar por menor magnitud de ganancias y despues menor tiempo 
    %busca ganancias relativamente pequeñas y, en caso de candidatos parecidos, favorece el menor tiempo
    candidatos_buenos = sortrows(candidatos_buenos,[9 8]);

    mejor = candidatos_buenos(1,:);  %selecciona la primera fila, el candidato balanceado

    Kp_theta = mejor(1);
    Kd_theta = mejor(2);
    Kp_v = mejor(3);
    Ki_v = mejor(4);

    mejores_ganancias = mejor(1:4);
    
    %Recalcular los polos del candidato seleccionado
    P_mejor = P0 + Kp_theta*P_Kptheta + Kd_theta*P_Kdtheta ...
        + Kp_v*P_Kpv + Ki_v*P_Kiv;

    mejores_polos = roots(P_mejor);
    
    %verifica cuál es el más cercano al eje imaginario
    mejor_max_real = max(real(mejores_polos));

    fprintf('RESULTADO:\n');
    fprintf(['Se encontraron combinaciones capaces de estabilizar ' ...
        'el sistema.\n\n']);

    fprintf('Candidato balanceado seleccionado:\n\n');

    fprintf('Kp_theta = %.6f\n',Kp_theta);
    fprintf('Kd_theta = %.6f\n',Kd_theta);
    fprintf('Kp_v     = %.6f\n',Kp_v);
    fprintf('Ki_v     = %.6f\n\n',Ki_v);

    disp('Polos de lazo cerrado:');
    disp(mejores_polos);

    fprintf('Polo mas a la derecha: %.6f\n\n', ...
        mejor_max_real);

    % Mostrar mejores candidatos
    n_mostrar = min(10,size(candidatos_buenos,1));

    tabla_candidatos = array2table( ... %convierte la matriz numérica en una tabla con nombres de columnas
        candidatos_buenos(1:n_mostrar,:), ...
        'VariableNames',{ ...
        'Kp_theta', ...
        'Kd_theta', ...
        'Kp_v', ...
        'Ki_v', ...
        'MaxReal', ...
        'WnDom', ...
        'ZetaDom', ...
        'TsAprox_s', ...
        'NormaK'});

    fprintf('----------------------------------------\n');
    fprintf('CANDIDATOS RECOMENDADOS\n');
    fprintf('----------------------------------------\n\n');
    disp(tabla_candidatos);

    %9.ANALISIS DEL POLO DOMINANTE
    %vuelve a identificar el polo dominante, pero ahora únicamente del controlador seleccionado.
    [~,indice_dominante] = max(real(mejores_polos));
    polo_dominante = mejores_polos(indice_dominante);
    wn = abs(polo_dominante);

    if wn > tol
        zeta = -real(polo_dominante)/wn;
    else
        zeta = NaN;
    end

    if real(polo_dominante) < 0
        ts_aprox = 4/abs(real(polo_dominante));
    else
        ts_aprox = Inf;
    end

    fprintf('Polo dominante:\n');
    disp(polo_dominante);

    fprintf('Frecuencia natural aproximada: %.6f rad/s\n', wn);
    fprintf('Factor de amortiguamiento aproximado: %.6f\n', zeta);
    fprintf('Tiempo de establecimiento aproximado: %.6f s\n\n', ts_aprox);

    %10.GRAFICA DE POLOS
    figure;
    plot( real(mejores_polos), imag(mejores_polos), 'x', 'MarkerSize',10);
    grid on;
    hold on;
    limites_y = ylim;
    plot([0 0],limites_y,'--');
    xlabel('Parte real');
    ylabel('Parte imaginaria');
    title('Polos de lazo cerrado - PD + PI de velocidad');
    hold off;
else
    candidato_encontrado = false;
    mejores_ganancias = [NaN NaN NaN NaN];    %NaN = not a number
    mejores_polos = [];
    mejor_max_real = Inf;                     %Inf = infinito
    tabla_candidatos = table();
    polo_dominante = NaN;
    wn = NaN;
    zeta = NaN;
    ts_aprox = NaN;
    fprintf('RESULTADO:\n');
    fprintf(['No se encontro una combinacion estable dentro de ' 'los rangos explorados.\n']);
end

% 11. GUARDAR RESULTADOS
control.G_dr = G_dr;
control.G_v = G_v;
control.G_theta = G_theta;
control.G_omega_p = G_omega_p;
control.relacion_angular = relacion_angular;
control.relacion_velocidad = relacion_velocidad;
control.denominador_comun = denominador_comun;
control.P0 = P0;
control.P_Kp_theta = P_Kptheta;
control.P_Kd_theta = P_Kdtheta;
control.P_Kp_v = P_Kpv;
control.P_Ki_v = P_Kiv;
control.numero_evaluadas = numero_evaluadas;
control.numero_estables = numero_estables;
control.candidatos_estables = candidatos;
control.candidatos_recomendados = tabla_candidatos;
control.candidato_encontrado = candidato_encontrado;
control.ganancias = mejores_ganancias;
control.polos = mejores_polos;
control.max_real = mejor_max_real;
control.polo_dominante = polo_dominante;
control.wn = wn;
control.zeta = zeta;
control.ts_aprox = ts_aprox;

fprintf('========================================\n');
end

function v = limpiar_inicio(v,tol)
%Elimina coeficientes iniciales numericamente nulos.

while length(v) > 1 && abs(v(1)) < tol
    v(1) = [];
end
end

function v = pad_left(v,n)
%Completa con ceros a la izquierda hasta longitud n.

if length(v) < n
    v = [zeros(1,n-length(v)) v];
end
end
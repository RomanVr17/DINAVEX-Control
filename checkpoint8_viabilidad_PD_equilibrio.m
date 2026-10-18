function pd = checkpoint8_viabilidad_PD_equilibrio(ft)
%CHECKPOINT 8 - VIABILIDAD DEL CONTROLADOR PD DE EQUILIBRIO
%Analiza la relacion entre G_theta y G_omega_p, muestra polos y ceros,
%construye la ecuacion caracteristica con un PD ideal y verifica,
%mediante condiciones de Routh-Hurwitz, si el PD angular aislado
%puede estabilizar completamente el modelo de referencia.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 8 - VIABILIDAD DEL PD\n');
fprintf('========================================\n\n');

% Recuperar funciones de transferencia angulares
G_theta   = minreal(ft.G_theta);
G_omega_p = minreal(ft.G_omega_p);

disp('Funcion de transferencia angular G_theta(s):');
disp(G_theta);

disp('Funcion de transferencia de velocidad angular G_omega_p(s):');
disp(G_omega_p)

% Verificar relacion G_omega_p(s) = s*G_theta(s)
s = tf([1 0], 1);

error_relacion = minreal(G_omega_p - s*G_theta);

[num_error, ~] = tfdata(error_relacion, 'v');

tol = 1e-6;

if all(abs(num_error) < tol)
    fprintf('VERIFICACION: G_omega_p(s) = s*G_theta(s)\n\n');
    relacion_correcta = true;
else
    fprintf('ADVERTENCIA: G_omega_p(s) no coincide exactamente con s*G_theta(s).\n');
    fprintf('Revisar simplificacion o redondeo numerico.\n\n');
    relacion_correcta = false;
end

% Polos y ceros del lazo angular
polos_theta = pole(G_theta);
ceros_theta = zero(G_theta);

disp('Polos de G_theta(s):');
disp(polos_theta);

disp('Ceros de G_theta(s):');
disp(ceros_theta);

% Mapa polos-ceros
figure;
pzmap(G_theta);
grid on;
title('Mapa polos-ceros de G_\theta(s)');

% Lugar geometrico de las raices de la planta angular sin compensar
figure;
rlocus(G_theta);
grid on;
title('Lugar geometrico de las raices de G_\theta(s)');

% Obtener coeficientes de G_theta(s)
[num, den] = tfdata(G_theta, 'v');

% Eliminar ceros numericos iniciales, si existen
while length(num) > 1 && abs(num(1)) < tol
    num(1) = [];
end

while length(den) > 1 && abs(den(1)) < tol
    den(1) = [];
end

% Normalizar denominador
num = num / den(1);
den = den / den(1);

fprintf('\n----------------------------------------\n');
fprintf('ANALISIS ALGEBRAICO DEL PD\n');
fprintf('----------------------------------------\n');

if length(den) == 4 && length(num) == 2

    % G_theta(s) = (n1*s + n0)/(s^3 + d2*s^2 + d1*s + d0)
    n1 = num(1);
    n0 = num(2);

    d2 = den(2);
    d1 = den(3);
    d0 = den(4);

    fprintf('G_theta(s) = (%.6f s %+ .6f) / ', n1, n0);
    fprintf('(s^3 %+ .6f s^2 %+ .6f s %+ .6f)\n\n', d2, d1, d0);

    fprintf('Controlador propuesto:\n');
    fprintf('C_PD(s) = Kp + Kd*s\n\n');

    fprintf('Ecuacion caracteristica:\n');
    fprintf('D(s) + C_PD(s)*N(s) = 0\n\n');

    fprintf('Coeficientes del polinomio de lazo cerrado:\n');
    fprintf('a2 = %.6f + (%.6f)*Kd\n', d2, n1);
    fprintf('a1 = %.6f + (%.6f)*Kd + (%.6f)*Kp\n', d1, n0, n1);
    fprintf('a0 = %.6f + (%.6f)*Kp\n\n', d0, n0);

    % Verificacion de viabilidad mediante Routh-Hurwitz
    es_viable_pd = NaN;
    Kd_limite = NaN;
    Kp_limite = NaN;
    a1_max = NaN;

    % Caso correspondiente al modelo de referencia:
    % n1 < 0, n0 > 0, d2 > 0 y d0 < 0
    if n1 < 0 && n0 > 0 && d2 > 0 && d0 < 0

        % Para a2 > 0:
        Kd_limite = -d2/n1;

        % Para a0 > 0:
        Kp_limite = -d0/n0;

        % El mayor valor posible de a1 manteniendo a2 > 0 y a0 > 0
        a1_max = d1 + n0*Kd_limite + n1*Kp_limite;

        fprintf('Condicion a2 > 0:\n');
        fprintf('Kd < %.8f\n\n', Kd_limite);

        fprintf('Condicion a0 > 0:\n');
        fprintf('Kp > %.8f\n\n', Kp_limite);

        fprintf('Maximo teorico de a1 bajo esas restricciones:\n');
        fprintf('a1_max = %.8f\n\n', a1_max);

        if a1_max <= 0
            es_viable_pd = false;

            fprintf('RESULTADO:\n');
            fprintf(['El PD angular aislado NO puede satisfacer simultaneamente\n' ...
                     'las condiciones de Routh-Hurwitz para este modelo.\n']);
            fprintf(['Por lo tanto, no se seleccionan Kp y Kd definitivos en\n' ...
                     'este checkpoint.\n\n']);

        else
            es_viable_pd = true;

            fprintf('RESULTADO:\n');
            fprintf(['Existe una region potencial de ganancias Kp y Kd que\n' ...
                     'debe analizarse con mayor detalle mediante LGR.\n\n']);
        end

    else
        fprintf(['La estructura de signos de esta FT no coincide con el caso\n' ...
                 'de referencia. Se requiere revisar Routh-Hurwitz de forma\n' ...
                 'general para este conjunto de parametros.\n\n']);
    end

else
    fprintf(['La FT angular no tiene la forma esperada de primer orden en\n' ...
             'el numerador y tercer orden en el denominador.\n']);
    fprintf('Se omite la prueba algebraica automatica.\n\n');

    n1 = NaN;
    n0 = NaN;
    d2 = NaN;
    d1 = NaN;
    d0 = NaN;

    es_viable_pd = NaN;
    Kd_limite = NaN;
    Kp_limite = NaN;
    a1_max = NaN;
end

% Guardar resultados
pd.G_theta = G_theta;
pd.G_omega_p = G_omega_p;

pd.relacion_correcta = relacion_correcta;

pd.polos_theta = polos_theta;
pd.ceros_theta = ceros_theta;

pd.num_theta = num;
pd.den_theta = den;

pd.n1 = n1;
pd.n0 = n0;
pd.d2 = d2;
pd.d1 = d1;
pd.d0 = d0;

pd.es_viable_pd_aislado = es_viable_pd;
pd.Kd_limite = Kd_limite;
pd.Kp_limite = Kp_limite;
pd.a1_max = a1_max;

fprintf('========================================\n');

end
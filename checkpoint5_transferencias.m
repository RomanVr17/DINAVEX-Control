function transferencias = checkpoint5_transferencias(numerico)

% CHECKPOINT 5 - FUNCIONES DE TRANSFERENCIA
% Obtiene las funciones de transferencia desde el par
% de entrada hacia los cuatro estados del robot a partir
% de las matrices numericas A y B.

% RECUPERAR MATRICES DEL MODELO LINEAL

A = numerico.A;
B = numerico.B;

% DEFINIR LAS SALIDAS
% x1 = dr       -> posicion lineal
% x2 = v        -> velocidad lineal
% x3 = theta    -> angulo de inclinacion
% x4 = omega_p  -> velocidad angular

C = eye(4);

% No existe paso directo de la entrada hacia las salidas
D = zeros(4,1);

% CREAR MODELO EN ESPACIO DE ESTADOS
sistema_ss = ss(A, B, C, D);

% CONVERTIR A FUNCIONES DE TRANSFERENCIA
G = tf(sistema_ss);

% SEPARAR CADA FUNCION DE TRANSFERENCIA
G_dr      = G(1,1);   % Dr(s) / Tau(s)
G_v       = G(2,1);   % V(s) / Tau(s)
G_theta   = G(3,1);   % Theta(s) / Tau(s)
G_omega_p = G(4,1);   % Omega_p(s) / Tau(s)

% GUARDAR RESULTADOS
transferencias.sistema_ss = sistema_ss;
transferencias.G_dr      = G_dr;
transferencias.G_v       = G_v;
transferencias.G_theta   = G_theta;
transferencias.G_omega_p = G_omega_p;

end
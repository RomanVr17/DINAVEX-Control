%MODELO SIMBOLICO DINAVEX
%ARCHIVO PRINCIPAL

clear
clc

%CHECKPOINT 1 CINEMATICA Y ENERGIAS

modelo = checkpoint1_energias();

%MOSTRAR RESULTADOS DEL CHECKPOINT 1

disp('======================================')
disp('CHECKPOINT 1: CINEMATICA Y ENERGIAS')
disp('======================================')

disp(' ')
disp('Vector de posicion del centro de masa rH:')
pretty(modelo.rH)

disp(' ')
disp('Vector de velocidad del centro de masa vH:')
pretty(modelo.vH)

disp(' ')
disp('Energia cinetica de las ruedas Kr:')
pretty(modelo.Kr)

disp(' ')
disp('Energia cinetica del cuerpo Kp:')
pretty(modelo.Kp)

disp(' ')
disp('Energia cinetica total K:')
pretty(modelo.K)

disp(' ')
disp('Energia potencial V:')
pretty(modelo.V)

disp(' ')
disp('Lagrangiano L:')
pretty(modelo.L)

disp(' ')
disp('Disipacion de Rayleigh Rp:')
pretty(modelo.Rp)

%CHECKPOINT 2 ECUACIONES DE EULER-LAGRANGE

ecuaciones = checkpoint2_euler_lagrange(modelo);

%MOSTRAR RESULTADOS DEL CHECKPOINT 2

disp(' ')
disp('======================================')
disp('CHECKPOINT 2: EULER-LAGRANGE')
disp('======================================')

%PRIMERA ECUACION

disp(' ')
disp('---------- ECUACION 1 ----------')

disp(' ')
disp('d/dt(dL/dv):')
pretty(ecuaciones.d_dt_dL_dv)

disp(' ')
disp('dL/ddr:')
pretty(ecuaciones.dL_ddr)

disp(' ')
disp('dRp/dv:')
pretty(ecuaciones.dRp_dv)

disp(' ')
disp('Tau_1:')
pretty(ecuaciones.tau1)

%SEGUNDA ECUACION

disp(' ')
disp('---------- ECUACION 2 ----------')

disp(' ')
disp('d/dt(dL/domega_p):')
pretty(ecuaciones.d_dt_dL_domega)

disp(' ')
disp('dL/dtheta:')
pretty(ecuaciones.dL_dtheta)

disp(' ')
disp('dRp/domega_p:')
pretty(ecuaciones.dRp_domega)

disp(' ')
disp('Tau_2:')
pretty(ecuaciones.tau2)


%CHECKPOINT 3
estados = checkpoint3_estados(modelo, ecuaciones);
% MOSTRAR RESULTADOS DEL CHECKPOINT 3
disp(' ')
disp('======================================')
disp('CHECKPOINT 3: MODELO EN ESPACIO DE ESTADOS')
disp('======================================')

disp(' ')
disp('Matriz dinamica M:')
pretty(estados.M)

disp(' ')
disp('Vector b:')
pretty(estados.b)

disp(' ')
disp('x1_dot = f1:')
pretty(estados.f1)

disp(' ')
disp('x2_dot = f2:')
pretty(estados.f2)

disp(' ')
disp('x3_dot = f3:')
pretty(estados.f3)

disp(' ')
disp('x4_dot = f4:')
pretty(estados.f4)

%CHECKPOINT 4
lineal = checkpoint4_linealizacion(estados);

% MOSTRAR RESULTADOS DEL CHECKPOINT 4
disp(' ')
disp('======================================')
disp('CHECKPOINT 4: LINEALIZACION')
disp('======================================')

disp(' ')
disp('Modelo evaluado en el punto de equilibrio:')
pretty(lineal.f_eq)

disp(' ')
disp('Matriz A linealizada:')
pretty(lineal.A)

disp(' ')
disp('Matriz B linealizada:')
pretty(lineal.B)

% =========================================================
% EVALUACION NUMERICA CON LOS PARAMETROS DE ARTURO
% =========================================================

p_arturo = parametros_arturo();
num_arturo = evaluar_modelo_lineal(modelo, lineal, p_arturo);

disp(' ')
disp('======================================')
disp('MODELO NUMERICO CON PARAMETROS DE ARTURO')
disp('======================================')

disp(' ')
disp('Matriz A numerica:')
disp(num_arturo.A)

disp(' ')
disp('Matriz B numerica:')
disp(num_arturo.B)

% =========================================================
% CHECKPOINT 5 - FUNCIONES DE TRANSFERENCIA
% =========================================================

ft_arturo = checkpoint5_transferencias(num_arturo);

disp(' ')
disp('======================================')
disp('FUNCIONES DE TRANSFERENCIA ARTURO')
disp('======================================')

disp(' ')
disp('Posicion / Par:')
ft_arturo.G_dr

disp(' ')
disp('Velocidad / Par:')
ft_arturo.G_v

disp(' ')
disp('Angulo / Par:')
ft_arturo.G_theta

disp(' ')
disp('Velocidad angular / Par:')
ft_arturo.G_omega_p

% =========================================================
% CHECKPOINT 6 - CONTROLABILIDAD
% =========================================================

control_arturo = checkpoint6_controlabilidad(num_arturo.A, num_arturo.B);

% =========================================================
% CHECKPOINT 7 - ANALISIS EN LAZO ABIERTO
% =========================================================

analisis_arturo = checkpoint7_analisis_lazo_abierto(ft_arturo, num_arturo);
                
% =========================================================
% CHECKPOINT 8 - VIABILIDAD DEL PD DE EQUILIBRIO
% =========================================================
pd_arturo = checkpoint8_viabilidad_PD_equilibrio(ft_arturo);

% =========================================================
% CHECKPOINT 9 - ANALISIS CONJUNTO PD-PI
conjunto_arturo = checkpoint9_PD_PI_conjunto(ft_arturo);
% =========================================================

% =========================================================
% CHECKPOINT 10 - PD DE EQUILIBRIO + PI DE VELOCIDAD
control_vel_arturo = checkpoint10_PD_PI_velocidad(ft_arturo);
% =========================================================

% =========================================================
% CHECKPOINT 11 - RESPUESTA TEMPORAL LINEAL
respuesta_arturo = checkpoint11_respuesta_temporal_lineal( ...
    num_arturo, control_vel_arturo);
% =========================================================

% =========================================================
% CHECKPOINT 12 - REGION DE RECUPERACION NO LINEAL
region_arturo = checkpoint12_region_recuperacion_no_lineal( ...
    estados, p_arturo, control_vel_arturo);
% =========================================================

% =========================================================
%CHECKPOINT 13
tau_limites = [0.5 1 2 4 6 8 10 12]; % N*m, par total
saturacion_arturo = checkpoint13_saturacion_actuador(estados, control_vel_arturo, region_arturo, tau_limites);
% =========================================================

% =========================================================
% EVALUACION NUMERICA CON LOS PARAMETROS DE DINAVEX
% =========================================================

% p_dinavex = parametros_DINAVEX();
% num_dinavex = evaluar_modelo_lineal( modelo, lineal, p_dinavex);
% 
% disp(' ')
% disp('======================================')
% disp('MODELO NUMERICO DE DINAVEX')
% disp('======================================')
% 
% disp(' ')
% disp('Matriz A de DINAVEX:')
% disp(num_dinavex.A)
% 
% disp(' ')
% disp('Matriz B de DINAVEX:')
% disp(num_dinavex.B)

% =========================================================
% CHECKPOINT 5 - FUNCIONES DE TRANSFERENCIA DINAVEX
% =========================================================

% ft_dinavex = checkpoint5_transferencias(num_dinavex);
% 
% disp(' ')
% disp('======================================')
% disp('FUNCIONES DE TRANSFERENCIA - DINAVEX')
% disp('======================================')
% 
% disp(' ')
% disp('Posicion / Par:')
% ft_dinavex.G_dr
% 
% disp(' ')
% disp('Velocidad / Par:')
% ft_dinavex.G_v
% 
% disp(' ')
% disp('Angulo / Par:')
% ft_dinavex.G_theta
% 
% disp(' ')
% disp('Velocidad angular / Par:')
% ft_dinavex.G_omega_p

% =========================================================
% CHECKPOINT 6 - CONTROLABILIDAD DINAVEX
% =========================================================

% control_dinavex = checkpoint6_controlabilidad(num_dinavex.A, num_dinavex.B);
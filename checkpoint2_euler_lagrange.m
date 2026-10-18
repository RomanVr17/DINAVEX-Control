function ecuaciones = checkpoint2_euler_lagrange(modelo)

%Aplica Euler-Lagrange al modelo simbolico para obtener
%las dos ecuaciones dinamicas generales del robot.

% CHECKPOINT 2
% ECUACIONES DE EULER-LAGRANGE

% RECUPERAR VARIABLES DEL CHECKPOINT 1
L  = modelo.L;
Rp = modelo.Rp;
dr      = modelo.dr;
theta   = modelo.theta;
v       = modelo.v;
omega_p = modelo.omega_p;

%ACELERACIONES SIMBOLICAS

syms dv domega_p real

%dv       = derivada de v respecto al tiempo
%domega_p = derivada de omega_p respecto al tiempo

% VARIABLES DEL SISTEMA Y SUS DERIVADAS
z = [dr;
     theta;
     v;
     omega_p];

zdot = [v;
        omega_p;
        dv;
        domega_p];


% ECUACION 1 DE EULER-LAGRANGE
% d/dt(dL/dv) - dL/ddr + dRp/dv = Tau_1

dL_dv  = diff(L,v);
dL_ddr = diff(L,dr);
dRp_dv = diff(Rp,v);

% Derivada total respecto al tiempo
d_dt_dL_dv = simplify( jacobian(dL_dv,z)*zdot );

% Tau_1
tau1 = simplify(d_dt_dL_dv - dL_ddr + dRp_dv );

% ECUACION 2 DE EULER-LAGRANGE
% d/dt(dL/domega_p) - dL/dtheta
% + dRp/domega_p = Tau_2

dL_domega  = diff(L,omega_p);
dL_dtheta  = diff(L,theta);
dRp_domega = diff(Rp,omega_p);

% Derivada total respecto al tiempo
d_dt_dL_domega = simplify(jacobian(dL_domega,z)*zdot );

% Tau_2
tau2 = simplify(d_dt_dL_domega - dL_dtheta + dRp_domega );

% GUARDAR RESULTADOS
% Aceleraciones simbolicas
ecuaciones.dv       = dv;
ecuaciones.domega_p = domega_p;

% Terminos de la ecuacion 1
ecuaciones.dL_dv      = dL_dv;
ecuaciones.dL_ddr     = dL_ddr;
ecuaciones.dRp_dv     = dRp_dv;
ecuaciones.d_dt_dL_dv = d_dt_dL_dv;

% Terminos de la ecuacion 2
ecuaciones.dL_domega      = dL_domega;
ecuaciones.dL_dtheta      = dL_dtheta;
ecuaciones.dRp_domega     = dRp_domega;
ecuaciones.d_dt_dL_domega = d_dt_dL_domega;

% Ecuaciones finales
ecuaciones.tau1 = tau1;
ecuaciones.tau2 = tau2;

end
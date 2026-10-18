function estados = checkpoint3_estados(modelo, ecuaciones)

%Resuelve simultaneamente las ecuaciones dinamicas para
%obtener x_dot = f(x,tau) del modelo no lineal.

%CHECKPOINT 3
%MODELO NO LINEAL EN ESPACIO DE ESTADOS

%RECUPERAR VARIABLES DEL CHECKPOINT 1
dr      = modelo.dr;
v       = modelo.v;
theta   = modelo.theta;
omega_p = modelo.omega_p;
R       = modelo.R;

%RECUPERAR ACELERACIONES DEL CHECKPOINT 2
dv       = ecuaciones.dv;
domega_p = ecuaciones.domega_p;

%ENTRADA DEL SISTEMA
syms tau real

%Para robot autoequilibrado:
%tau_1 = tau/R
%tau_2 = 0

%ECUACIONES DINAMICAS
eq1 = ecuaciones.tau1 == tau/R;
eq2 = ecuaciones.tau2 == 0;

%SEPARAR LAS ACELERACIONES EN FORMA MATRICIAL
%Las ecuaciones tienen la forma:
%M(q) * [dv; domega_p] = b(q, dq, tau)
%donde las incógnitas son dv y domega_p.

[M, b] = equationsToMatrix([eq1; eq2], [dv; domega_p]);
M = simplify(M);
b = simplify(b);

%RESOLVER EL SISTEMA SIMULTANEO

aceleraciones = M \ b;
dv_sol = simplify(aceleraciones(1));
domega_p_sol = simplify(aceleraciones(2));

%VECTOR DE ESTADOS
%x1 = dr
%x2 = v
%x3 = theta
%x4 = omega_p

x = [dr; v; theta; omega_p];

%ECUACIONES DE ESTADO
%x1_dot = x2
%x2_dot = dv
%x3_dot = x4
%x4_dot = domega_p

f1 = v;
f2 = dv_sol;
f3 = omega_p;
f4 = domega_p_sol;

f = [f1; f2; f3; f4];

%GUARDAR RESULTADOS

estados.tau = tau;
estados.eq1 = eq1;
estados.eq2 = eq2;
estados.M = M;
estados.b = b;
estados.dv = dv_sol;
estados.domega_p = domega_p_sol;
estados.x = x;
estados.f1 = f1;
estados.f2 = f2;
estados.f3 = f3;
estados.f4 = f4;
estados.f = f;

end
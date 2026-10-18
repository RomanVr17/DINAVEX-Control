function lineal = checkpoint4_linealizacion(estados)

%Calcula las matrices Jacobianas A y B del modelo no lineal
%y las evalua alrededor del punto de equilibrio vertical.

%CHECKPOINT 4
%LINEALIZACION MEDIANTE MATRICES JACOBIANAS
%RECUPERAR VECTOR DE ESTADOS, MODELO NO LINEAL Y ENTRADA

x   = estados.x;
f   = estados.f;
tau = estados.tau;

%MATRICES JACOBIANAS GENERALES
%Modelo no lineal:
%x_dot = f(x,tau)
%La linealizacion de primer orden se obtiene mediante:
%A = df/dx
%B = df/dtau

A_sym = jacobian(f, x);
B_sym = jacobian(f, tau);

%PUNTO DE EQUILIBRIO
%Para el equilibrio vertical:
%dr = 0
%v = 0
%theta = 0
%omega_p = 0
%tau = 0

xe   = sym(zeros(4,1));
taue = sym(0);

%VARIABLES Y VALORES A SUSTITUIR
variables_eq = [x; tau];
valores_eq = [xe; taue];

%VERIFICAR EL PUNTO DE EQUILIBRIO
f_eq = simplify( subs(f, variables_eq, valores_eq),'Steps', 100);

%EVALUAR LAS JACOBIANAS EN EL PUNTO DE EQUILIBRIO
A_lin = subs(A_sym, variables_eq, valores_eq);
B_lin = subs(B_sym, variables_eq, valores_eq);

%SIMPLIFICAR DESPUES DE EVALUAR EL EQUILIBRIO
A_lin = simplify(A_lin, 'Steps', 100);
B_lin = simplify(B_lin, 'Steps', 100);

%GUARDAR RESULTADOS
lineal.xe   = xe;
lineal.taue = taue;
lineal.f_eq = f_eq;
lineal.A_sym = A_sym;
lineal.B_sym = B_sym;
lineal.A = A_lin;
lineal.B = B_lin;

end
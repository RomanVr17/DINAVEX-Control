function numerico = evaluar_modelo_lineal(modelo, lineal, p)

%Sustituye un conjunto de parametros fisicos en las matrices
%simbolicas A y B para obtener su forma numerica.

%EVALUACION NUMERICA DEL MODELO LINEAL
%PARAMETROS SIMBOLICOS

parametros_sym = [modelo.mh, modelo.mr, modelo.Ip, modelo.Ir, modelo.l, ...
    modelo.R, modelo.c_v, modelo.c_w, modelo.g ];

%VALORES NUMERICOS
parametros_num = [ p.mh, p.mr, p.Ip, p.Ir, p.l, p.R, p.c_v, ...
    p.c_w, p.g ];

%SUSTITUIR EN LAS MATRICES SIMBOLICAS

A_num = subs(lineal.A, parametros_sym, parametros_num);
B_num = subs(lineal.B, parametros_sym, parametros_num);

% CONVERTIR DE SIMBOLICO A DOUBLE

A_num = double(A_num);
B_num = double(B_num);

% GUARDAR RESULTADOS

numerico.A = A_num;
numerico.B = B_num;

end
function modelo = checkpoint1_energias()
%Calcula simbolicamente la cinematica, energias, Lagrangiano
%y funcion de disipacion de Rayleigh del robot.

%CINEMATICA, ENERGIAS, LAGRANGIANO Y RAYLEIGH
syms mh mr Ip Ir l R c_v c_w g real

%VARIABLES GENERALIZADAS Y VELOCIDADES
syms dr theta v omega_p real

%POSICION DEL CENTRO DE MASA DEL CUERPO
rH = [dr + l*sin(theta);
      l*cos(theta);
      0];

%VELOCIDAD DEL CENTRO DE MASA
vH = [v + l*omega_p*cos(theta);
     -l*omega_p*sin(theta);
      0];

%VELOCIDAD ANGULAR DE LAS RUEDAS
omega_r = v/R;

%ENERGIA CINETICA DE LAS RUEDAS

Kr = (1/2)*mr*v^2 + (1/2)*Ir*omega_r^2;
Kr = simplify(Kr);

%ENERGIA CINETICA DEL CUERPO
Kp = (1/2)*mh*(vH.'*vH) + (1/2)*Ip*omega_p^2;
Kp = simplify(expand(Kp));

%ENERGIA CINETICA TOTAL
K = simplify(Kp + Kr);

%ENERGIA POTENCIAL
Vr = sym(0);
Vp = mh*g*l*cos(theta);
V = simplify(Vp + Vr);

%LAGRANGIANO
L = simplify(K - V);

%FUNCION DE DISIPACION DE RAYLEIGH
Rp = (1/2)*( c_v*v^2 + c_w*(v/R - omega_p)^2 );
Rp = simplify(Rp);

%GUARDAR RESULTADOS EN UNA ESTRUCTURA
%Parametros
modelo.mh = mh;
modelo.mr = mr;
modelo.Ip = Ip;
modelo.Ir = Ir;
modelo.l  = l;
modelo.R  = R;
modelo.c_v = c_v;
modelo.c_w = c_w;
modelo.g  = g;

%Variables
modelo.dr      = dr;
modelo.theta   = theta;
modelo.v       = v;
modelo.omega_p = omega_p;

%Cinematica
modelo.rH      = rH;
modelo.vH      = vH;
modelo.omega_r = omega_r;

%Energias
modelo.Kr = Kr;
modelo.Kp = Kp;
modelo.K  = K;
modelo.Vr = Vr;
modelo.Vp = Vp;
modelo.V  = V;

%Lagrangiano y Rayleigh
modelo.L  = L;
modelo.Rp = Rp;

end
function p = parametros_DINAVEX()

%PARAMETROS DEL ROBOT DINAVEX
%Distancia del eje de las ruedas al centro de masa [m]
p.l = 0;         

%Masa del cuerpo del robot [kg]
p.mh = 0;        

%Momento de inercia del cuerpo respecto al eje de balance [kg*m^2]
p.Ip = 0;          

%Masa equivalente de las ruedas [kg]
p.mr = 0;         

%Momento de inercia de las ruedas [kg*m^2]
p.Ir = 0;        

%Radio de las ruedas [m]
p.R = 0;          

%Coeficiente de friccion rotacional [N*m*s/rad]
p.c_w = 0;       

%Coeficiente de friccion lineal [N*s/m]
p.c_v = 0;       

%Gravedad  [m/s^2]
p.g = 9.81;         

end
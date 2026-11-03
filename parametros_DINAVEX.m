function p = parametros_DINAVEX()
%PARAMETROS DEL ROBOT DINAVEX

%Distancia del eje de las ruedas al centro de masa [m]
p.l = 0.3535;         

%Masa del cuerpo del robot [kg]
p.mh = 22.937;        

%Momento de inercia del cuerpo respecto al eje de pitch que pasa por el CM [kg*m^2]
p.Ip = 1.1178;

%Masa equivalente de las ruedas [kg]
p.mr = 7.4;         

%Momento de inercia equivalente de ambas ruedas [kg*m^2]
p.Ir = 0.02185;        

%Radio de las ruedas [m]
p.R = 0.0825;          

%Coeficiente de friccion rotacional [N*m*s/rad]
p.c_w = 0.0067;       

%Coeficiente de friccion lineal [N*s/m]
p.c_v = 3.5;       

%Gravedad  [m/s^2]
p.g = 9.81;         
end
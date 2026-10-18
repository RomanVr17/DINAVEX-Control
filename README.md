# DINAVEX-Control

Modelado matemático, análisis y control del robot autoequilibrado DINAVEX

# DINAVEX



Modelado matemático y diseño del sistema de control de un robot didáctico autoequilibrado desarrollado como parte del Proyecto Integrador.



\## Objetivo del repositorio



Este repositorio contiene los scripts utilizados para:



\- derivación del modelo dinámico;

\- obtención del modelo no lineal;

\- linealización alrededor del equilibrio;

\- obtención de funciones de transferencia;

\- análisis de controlabilidad;

\- análisis de estabilidad en lazo abierto;

\- diseño y evaluación del controlador;

\- simulaciones lineales y no lineales;

\- análisis de región de recuperación;

\- análisis de saturación del actuador;

\- discretización del controlador.



\## Estructura principal



\- `Main.m`: programa principal.

\- `checkpoint1\_energias.m`: cálculo de energías.

\- `checkpoint2\_euler\_lagrange.m`: formulación mediante Euler-Lagrange.

\- `checkpoint3\_estados.m`: modelo no lineal en espacio de estados.

\- `checkpoint4\_linealizacion.m`: linealización del modelo.

\- `checkpoint5\_transferencias.m`: funciones de transferencia.

\- `checkpoint6\_controlabilidad.m`: análisis de controlabilidad.

\- `checkpoint7\_analisis\_lazo\_abierto.m`: análisis de polos y ceros.

\- `checkpoint8\_viabilidad\_PD\_equilibrio.m`: análisis del controlador PD.

\- `checkpoint9\_PD\_PI\_conjunto.m`: análisis PD + PI.

\- `checkpoint10\_PD\_PI\_velocidad.m`: controlador PD + PI de velocidad.

\- `checkpoint11\_respuesta\_temporal\_lineal.m`: respuesta temporal lineal.

\- `checkpoint12\_region\_recuperacion\_no\_lineal.m`: región de recuperación no lineal.

\- `checkpoint13\_saturacion\_actuador.m`: análisis de saturación.

\- `parametros\_arturo.m`: parámetros de referencia.

\- `parametros\_DINAVEX.m`: parámetros del prototipo DINAVEX.



\## Variables de estado



El modelo utiliza:



x = \[d\_r, v, theta, omega\_p]^T



donde:



\- `d\_r`: desplazamiento longitudinal.

\- `v`: velocidad lineal.

\- `theta`: ángulo de inclinación.

\- `omega\_p`: velocidad angular del cuerpo.



\## Control



La arquitectura de control considerada utiliza:



\- controlador PD para equilibrio angular;

\- controlador PI de velocidad;

\- estimación de theta y omega\_p mediante IMU y filtro de Kalman;

\- estimación de velocidad lineal mediante encoders.



La acción total de control se expresa como:



u = u\_theta + u\_v



\## Parámetros



Se consideran dos conjuntos de parámetros:



\- `parametros\_arturo.m`: utilizados para validación del modelo.

\- `parametros\_DINAVEX.m`: parámetros propios del prototipo obtenidos a partir del CAD, mediciones físicas y estimaciones experimentales.



\## Ejecución



Ejecutar:



`Main.m`



El programa llama secuencialmente a los distintos checkpoints de modelado, análisis y control.



\## Estado actual



Actualmente se encuentran implementados los Checkpoints 1 al 13.



El siguiente desarrollo corresponde a la selección de frecuencia de muestreo y discretización del controlador.



\## Proyecto



Proyecto Integrador  

Ingeniería Mecatrónica  

Universidad Autónoma de San Luis Potosí


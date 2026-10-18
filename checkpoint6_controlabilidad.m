function control = checkpoint6_controlabilidad(A, B) 
%CHECKPOINT 6 - CONTROLABILIDAD
%Calcula la matriz de controlabilidad para un sistema de 4 estados:
%Mc = [B  AB  A^2B  A^3B]
%verifica si su rango es igual al numero de estados.

fprintf('\n========================================\n');
fprintf('CHECKPOINT 6 - CONTROLABILIDAD\n');
fprintf('========================================\n\n');

% Numero de estados
n = size(A,1);

% Matriz de controlabilidad
Mc = [B A*B A^2*B A^3*B];

% Mostrar matriz
disp('Matriz de controlabilidad Mc = ');
disp(Mc);

% Calcular rango
rango_Mc = rank(Mc);

fprintf('Numero de estados: %d\n', n);
fprintf('Rango de Mc:       %d\n\n', rango_Mc);

% Determinar controlabilidad
if rango_Mc == n
    fprintf('RESULTADO: El sistema es completamente controlable.\n');
    es_controlable = true;
else
    fprintf('RESULTADO: El sistema NO es completamente controlable.\n');
    es_controlable = false;
end

% Guardar resultados
control.Mc = Mc;
control.rango = rango_Mc;
control.n_estados = n;
control.es_controlable = es_controlable;

fprintf('\n========================================\n');

end
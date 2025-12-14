function [p, m] = pr_hmmMINE2(o, a, b, pi)
    % INPUTS:
    % o = Secuencia de observaciones
    % a = Matriz de transición ampliada (TRANS_HAT)
    % b = Matriz de emisión ampliada (EMIS_HAT)
    % pi = Probabilidades iniciales (ya incluida en a)

    % OUTPUTS:
    % p = Probabilidad total de la secuencia
    % m = Matriz de probabilidades Forward (T x n)

    T = length(o);  % Longitud de la secuencia de observaciones
    n = size(a, 1); % Número de estados (incluyendo el ficticio)

    % Inicialización
    m = zeros(T, n); % Matriz de probabilidades intermedias
    for i = 2:n % Excluye el estado ficticio
        m(1, i) = b(i, o(1)) * a(1, i); % Usa la fila 1 de a (PI_HAT)
    end

    % Recursión
    for t = 1:(T-1)
        for j = 2:n % Excluye el estado ficticio
            z = 0;
            for i = 2:n % Excluye el estado ficticio
                z = z + a(i, j) * m(t, i);
            end
            m(t+1, j) = z * b(j, o(t+1));
        end
    end

    % Terminación
    p = sum(m(T, 2:end)); % Sumar solo los estados válidos
end
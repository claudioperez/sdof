function x = GammaSDOF(scheme, model, f, u0, v0, maxiter)
    if nargin < 5, u0 = 0; end
    if nargin < 6, v0 = 0; end
    if nargin < 7, maxiter = 10; end

    nt = length(f);
    x = zeros(nt, 3);
    
    % Given initial conditions
    x(1, 1) = u0;
    x(1, 2) = v0;
    % Consistent initial acceleration
    [p, k, c, m] = evaluate(model, [u0, v0, 0]);
    x(1, 3) = -p / m;

    ip = scheme.iy;

    for i = 1:nt-1
        fa = increment(scheme, i, f);

        % Prediction of augmented state xa
        xa = predict(scheme, x(i, :), scheme.iy);

        % Iterative solve for equilibriated xa
        for iter = 1:maxiter
            [p, k, c, m] = evaluate(model, xa);
            ga = p - fa;
            Dg = scheme.ga(3)*m + scheme.ga(2)*c + scheme.ga(1)*k;
            dy = -ga / Dg;

            % Equation (8)
            for r = 1:3
                xa(r) = xa(r) + scheme.ga(r) * dy;
            end
        end

        % Advance from augmented state xa to true state x(i+1,:) at time t+dt
        x(i+1, :) = advance(scheme, x(i, :), xa);
    end
    
    x = x.';
end

function fa = increment(scheme, i, f)
    fa = (1 - scheme.lam) * f(i) + scheme.lam * f(i + 1);
end

function [p, k, c, m] = evaluate(model, xa)
    k = model(1);
    c = model(2);
    m = model(3);
    p = m * xa(3) + c * xa(2) + k * xa(1);
end

function xa = predict(scheme, xo, ip)
    if nargin < 3, ip = scheme.iy; end

    xa = zeros(1, 3);
    for r = 1:3
        sum_term = 0;
        for s = 1:3
            fact = scheme.Ga(r, s) + scheme.ga(r) / scheme.ga(ip) * (double(ip == s) - scheme.Ga(ip, s));
            sum_term = sum_term + fact * xo(s);
        end
        xa(r) = sum_term;
    end
end

function x = advance(scheme, xo, xa)
    iy = scheme.iy;
    x = zeros(1, 3);

    for s = 1:3
        for r = 1:3
            fact = scheme.Gn(s, r) - scheme.gn(r) / scheme.ga(iy) * scheme.Ga(iy, r);
            x(s) = x(s) + fact * xo(r);
        end
        x(s) = x(s) + xa(iy) * scheme.gn(s) / scheme.ga(iy);
    end
end


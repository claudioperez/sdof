main

function main
    % Define parameters
    f = [0.0000, 5.0000, 8.6603, 10.0000, 8.6603, 5.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000];
    m = 1;
    k = 2;
    c = 0;
    T = 2 * pi * sqrt(m / k);
    dt = 10 * T; % Size of time step
    nt = 40;     % Number of time steps
    f = zeros(1, nt);
    v0 = 0;
    u0 = 1;
    Ey = "a";
    
    % Set parameters
    beta = 0.3025;
    gamma = 0.6;
    model = [k, c, m];
    
    % Run Newmark gamma method
    [Gn, gn]  = NewmarkScheme(Ey, dt, beta, gamma);
    scheme    = GammaScheme(Ey, Gn, gn, Gn, gn, Ey);
    [u, v, a] = GammaSDOF(scheme, model, f, v0, u0);
    plot(u, 'o', 'DisplayName', 'newmark');
    hold on;

    % Modified average gamma
    alpha   = -1.0
    alpha_u = 1 + alpha;
    alpha_a = 1;
    gamma   = 1 / 2 + alpha_a - alpha_u;
    beta    = (1 + alpha_a - alpha_u)^2 / 4;
    [Ga, ga] = AverageScheme([alpha_u, 1, alpha_a], {Gn, gn});
    scheme = GammaScheme(Ey, Gn, gn, Ga, ga, alpha_u, Ey);
    [u, v, a] = GammaSDOF(scheme, model, f, v0, u0);
    plot(u, 'x', 'DisplayName', sprintf('\\alpha = %.2f', alpha));

    % Collocation method
    beta = 0.3025;
    gamma = 0.6;
    theta = 1.021712;
    [Gn, gn] = NewmarkScheme(Ey, dt, beta, gamma);
    [Gc, gc] = CollocateScheme(Ey, dt, theta, beta, gamma);
    scheme   = GammaScheme(Ey, Gn, gn, Gc, gc, theta, Ey);
    [u, v, a] = GammaSDOF(scheme, model, f, v0, u0);
    plot(u, '-', 'DisplayName', sprintf('\\theta = %.2f', theta));

    % Show plot
    legend show;
    hold off;
end


function scheme = GammaScheme(Ey, Gn, gn, Ga, ga, lam, predictor)
    scheme.iy = find('uva' == Ey);
    scheme.Gn = Gn;
    scheme.gn = gn;
    scheme.Ga = Ga;
    scheme.ga = ga;
    scheme.lam = lam;
end


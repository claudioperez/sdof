function [Ga, ga] = AverageScheme(alpha, scheme_a)
    G = scheme_a{1};
    g = scheme_a{2};
    A = diag(alpha);
    Ga = A * (G - eye(3)) + eye(3);
    ga = A * g;
end


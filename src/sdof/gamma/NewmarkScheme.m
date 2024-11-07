function [Ga, ga] = NewmarkScheme(Ey, dt, beta, gamma, betas, betass, gammas)
    if nargin < 5 || isempty(gammas)
        gammas = 1 - gamma;
    end
    if nargin < 6 || isempty(betass)
        betass = 1;
    end
    if nargin < 7 || isempty(betas)
        betas = (1 - 2 * beta) / 2;
    end

    % First form for Ey == "a"
    Ga = [1, betass * dt, betas * dt^2;
          0,           1,  gammas * dt;
          0,           0,            0];

    ga = [beta * dt^2; gamma * dt; 1];

    if ~strcmp(Ey, "a")
        [Ga, ga] = convert({Ga, ga}, "a", Ey);
    end
end

function [Gz, gz] = convert(scheme, Ey, Ez)
    if strcmp(Ey, Ez)
        Gz = scheme{1};
        gz = scheme{2};
        return;
    end

    Gy = scheme{1};
    gy = scheme{2};
    iy = find('uva' == Ey);
    iz = find('uva' == Ez);
    iw = find(~ismember(1:3, [iy, iz]));

    mu = 1 / gy(iz);
    gz = gy * mu;

    Gz = zeros(3, 3);
    for i = 1:3
        if i ~= iz
            for j = 1:3
                Gz(i, j) = Gy(i, j) - mu * gy(i) * Gy(iz, j);
                if i == iy
                    Gz(i, j) = 0;
                end
            end
        else
            Gz(i, :) = [0, 0, 0];
        end
    end
end


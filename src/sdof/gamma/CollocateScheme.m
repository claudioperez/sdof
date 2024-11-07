function [G, g] = CollocateScheme(Ey, dt, theta, varargin)

    [Gc, gc] = NewmarkScheme("a", dt * theta, varargin{:});
    [Gy, gy] = NewmarkScheme("a", dt, varargin{:});

    iy = find('uva' == Ey);
    ia = 3;

    G = Gc;
    g = zeros(1, 3);
    
    for i = 1:3
        for j = 1:3
            G(i, j) = Gc(i, j) + gc(i) * ((1 - theta) * (j == ia) - theta / gy(iy) * Gy(iy, j));
        end
        g(i) = gc(i) * theta / gy(iy);
    end
end


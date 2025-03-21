import sympy as sp
from gamma import newmark_gamma, convert, average_gamma, collocate
from integrate import advance, predict, GammaScheme


def print_method(method):
    Ey = "a"
    iy = "uva".index(Ey)

    dt = sp.Symbol(r"\Delta t")

    if method == "GS4":
        lamda = sp.symbols(" ".join(f"lambda_{i}" for i in range(1,6)))


        betas = lamda[1] - lamda[2]
        betass = lamda[0]
        beta = lamda[2]
        gamma = lamda[4]
        gammas = lamda[3] - lamda[4]

        params = beta, gamma, betas, betass, gammas

    else:
        params = sp.symbols("beta gamma")
        alpha = sp.symbols(" ".join(f"alpha_{i}" for i in "uva"))
        theta = sp.symbols("vartheta")


    xc = sp.symbols(r"\tilde{u} \tilde{v} \tilde{a}")
    xo = sp.symbols(r"\bm{u}_{n} \bm{v}_{n} \bm{a}_{n}")
    xn = sp.symbols(r"\bm{u}_{n+1} \bm{v}_{n+1} \bm{a}_{n+1}")

    Gn, gn = newmark_gamma(Ey, dt, *params)
    Gc, gc = collocate(Ey, dt, theta, *params)
    scheme = GammaScheme(Ey, Gn, gn, Gc, gc, theta)

    xa = predict(scheme, xo)
    for i in xa:
        print("$$")
        print(sp.latex(i))
        print("$$")

    xa[iy] = xc[iy]
    for i in advance(scheme, xo, xa):
        print("$$")
        print(sp.latex(i))
        print("$$")

def print_table(method):
    dt = sp.Symbol(r"\Delta t")

    alpha = sp.symbols(" ".join(f"alpha_{i}" for i in "uva"))
    theta = sp.symbols("vartheta")

    if method=="GS4":
        lamda = sp.symbols(" ".join(f"lambda_{i}" for i in range(1,6)))


        betas = lamda[1] - lamda[2]
        betass = lamda[0]
        beta = lamda[2]
        gamma = lamda[4]
        gammas = lamda[3] - lamda[4]

        params = beta, gamma, betas, betass, gammas

    else:
        params = sp.symbols("beta gamma")

    def create(Ey):
        if method == "newmark":
            return newmark_gamma(Ey, dt, *params)

        elif method == "alpha":
            scheme = newmark_gamma(Ey, dt, *params)
            return average_gamma(alpha, scheme)
        else:
            return collocate(Ey, dt, theta, *params)


    print("$$")
    print(r"\begin{array}{l|ccc}")
    for i,Ei in enumerate("uva"):

        for j, Ej in enumerate("uva"):

            print(rf"\tilde{{\Gamma}}_{{ {Ei} {Ej} }} & ", end="")

            for iy, Ey in enumerate("uva"):
                if iy:
                    print("     & ", end="")

                G, g = create(Ey)
                Gij = sp.simplify(G[i][j])
#               Gij = sp.collect(Gij, alpha[i]*dt)
#               if c:= Gij.coeff(alpha[i]*dt):
#                   Gij = alpha[i]*dt*sp.expand(c)
                print(sp.latex(
                    sp.nsimplify(Gij)
                ).replace("1.0", "1")
                )
            print(r"\\[0.3cm]")

    print("\\hline ")

    for i,Ei in enumerate("uva"):
        print(rf"\tilde{{\gamma}}_{{ {Ei} }} & ", end="")
        for iy, Ey in enumerate("uva"):
            if iy:
                print("     & ", end="")
            _, g = create(Ey)
            gi = sp.simplify(g[i])
            gi = sp.collect(gi, alpha[i]*dt)
#           if c:= gi.coeff(alpha[i]*dt):
#               gi = alpha[i]*dt*sp.expand(c)
            print(sp.latex(gi).replace("1.0", "1"))
        print("\n" + r"\\[0.3cm]")

    print(r"\end{array}")
    print("$$")

    return


    for Ey in "uv":
        Gy, gy = newmark_gamma(Ey, dt, beta, gamma, betas, betass, gammas)
        for Ez in "uva":
            if Ez != Ey:
                print(Ey, "->", Ez)
                Gaz, gaz = newmark_gamma(Ez, dt, beta, gamma, betas, betass, gammas)
                Gaz = [[sp.simplify(ij) for ij in ii] for ii in Gaz]

                Gz, gz = convert((Gy, gy), Ey, Ez)


                for i in range(3):
                    for j in range(3):
                        Gzij = sp.simplify(Gz[i][j])
                        assert (Gzij == Gaz[i][j])

def test_newmark():
    beta, gamma = sp.symbols("beta gamma")
    dt = sp.Symbol(r"\Delta t")

#    for Ey in "uva":
#        print(Ey)
#        Gy, gy = newmark_gamma(Ey, dt, beta, gamma)
#        print([sp.simplify(j) for j in gy])
#        for i in Gy:
#            print([sp.simplify(j) for j in i])


    for Ey in "uv":
        Gy, gy = newmark_gamma(Ey, dt, beta, gamma)
        for Ez in "uva":
            if Ez != Ey:
                Gaz, gaz = newmark_gamma(Ez, dt, beta, gamma)
                Gaz = [[sp.simplify(ij) for ij in ii] for ii in Gaz]

                Gz, gz = convert((Gy, gy), Ey, Ez)



                for i in range(3):
                    for j in range(3):
                        Gzij = sp.simplify(Gz[i][j])
                        assert (Gzij == Gaz[i][j])

#       print("$$")
#       print(sp.latex(sp.Matrix(newmark(E, dt, beta, gamma))))
#       print("$$")

if __name__ == "__main__":
    test_newmark()

    for table in "Newmark", "Wilson", "Alpha":
        print(f"# {table}\n")
        print_table(table.lower())
#       print_method()
        print()


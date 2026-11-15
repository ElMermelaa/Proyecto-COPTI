"""    mean_field_linear_blocks(p, lambda)

Bloques del jacobiano de la red todos-con-todos alrededor de un equilibrio
homogéneo. J_parallel=J_individual y J_perp=J_individual+diag(lambda,0).
Su validez NO se extiende sin cambios a una topología arbitraria.
"""
function mean_field_linear_blocks(p::FHNParameters, lambda::Real)
    isfinite(lambda) || throw(ArgumentError("lambda debe ser real finito."))
    Jparallel = jacobian_fhn(equilibrium(p)[1], p)
    Jperp = copy(Jparallel)
    Jperp[1,1] += lambda
    eigenparallel = eigvals(Jparallel)
    eigenperp = eigvals(Jperp)
    transverse_stable = all(real(z) < 0 for z in eigenperp)
    return (parallel=Jparallel, transverse=Jperp,
            eigen_parallel=eigenparallel, eigen_transverse=eigenperp,
            transverse_stable=transverse_stable,
            full_equilibrium_stable=transverse_stable && all(real(z)<0 for z in eigenparallel))
end

"""    transverse_floquet(orbit, p, lambda; reltol=1e-10, abstol=1e-10)

Calcula los multiplicadores transversales del ciclo sincronizado para la red
todos-con-todos (lambda constante). Integra SIMULTÁNEAMENTE la órbita y su
matriz fundamental, evitando el error de interpolación de orbit.solution(t).

La órbita se inicializa desde el estado refinado por shooting de la Etapa B:

    z' = F(z),     Phi' = A_perp(z,lambda) Phi,     Phi(0)=I_2
    A_perp = [1-v^2+lambda  -1; 1/tau  -b/tau].

Para lambda = 0 debe existir un multiplicador aproximadamente igual a 1.
Para lambda != 0, ambos multiplicadores son transversales a la variedad de
sincronía y no tienen por qué incluir un multiplicador igual a 1.
"""
function transverse_floquet(orbit, p::FHNParameters, lambda::Real;
                             reltol::Real=1e-10, abstol::Real=1e-10)
    orbit.status == :converged || throw(ArgumentError("Ciclo no refinado."))
    isfinite(lambda) || throw(ArgumentError("lambda debe ser finito."))
    reltol > 0 && abstol > 0 ||
        throw(ArgumentError("Las tolerancias deben ser positivas."))
    isapprox(orbit.current, p.I; atol=1e-8, rtol=0) ||
        throw(ArgumentError("La corriente del ciclo y de p no coincide."))

    # y = [v, w, Phi_11, Phi_21, Phi_12, Phi_22]; orden por columnas.
    # Al integrar la órbita junto a Phi se evita evaluar A_perp sobre una
    # interpolación del ciclo guardado con saveat relativamente grueso.
    function orbit_and_variational!(dy, y, _, t)
        v, w = y[1], y[2]
        dy[1] = fhn_nonlinearity(v) - w + p.I
        dy[2] = (v + p.a - p.b*w) / p.tau

        a11 = 1 - v^2 + lambda
        a21 = 1 / p.tau
        a22 = -p.b / p.tau

        dy[3] = a11*y[3] - y[4]
        dy[4] = a21*y[3] + a22*y[4]
        dy[5] = a11*y[5] - y[6]
        dy[6] = a21*y[5] + a22*y[6]
        return nothing
    end

    y0 = [Float64(orbit.initial_state[1]), Float64(orbit.initial_state[2]),
          1.0, 0.0, 0.0, 1.0]
    prob = ODEProblem(orbit_and_variational!, y0, (0.0, Float64(orbit.period)))
    sol = solve(prob, Tsit5(); reltol=reltol, abstol=abstol,
                save_everystep=false, save_start=false, save_end=true)

    yend = sol.u[end]
    M = [yend[3] yend[5]; yend[4] yend[6]]
    mus = eigvals(M)
    radius = maximum(abs.(mus))
    closure_error = norm(yend[1:2] - y0[1:2])
    return (monodromy=M, multipliers=mus,
            spectral_radius=radius,
            stability=radius < 1-1e-4 ? :attractive :
                      radius > 1+1e-4 ? :unstable : :near_critical,
            orbit_closure_error=closure_error)
end

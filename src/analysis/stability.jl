"""    jacobian_fhn(v, p)

Jacobiano del sistema individual en `(v,w)` (independiente de `w`):

    [ 1-v²      -1  ]
    [ 1/tau    -b/tau ]
"""
function jacobian_fhn(v::Real, p::FHNParameters)
    J = Matrix{Float64}(undef, 2, 2)
    J[1, 1] = 1 - v^2
    J[1, 2] = -1.0
    J[2, 1] = 1 / p.tau
    J[2, 2] = -p.b / p.tau
    return J
end

"""    stability_info(p; tol=1e-9)

Información de estabilidad local del equilibrio único.
`classification` puede ser `:stable`, `:unstable` o `:nonhyperbolic`.
Cuando la traza es casi nula, NO se afirma estabilidad no lineal.
"""
function stability_info(p::FHNParameters; tol::Real=1e-9)
    tol >= 0 || throw(ArgumentError("tol debe ser no negativo."))
    z = equilibrium(p)
    J = jacobian_fhn(z[1], p)
    trJ = tr(J)
    detJ = det(J)
    lambda = eigvals(J)
    classification = trJ < -tol ? :stable : trJ > tol ? :unstable : :nonhyperbolic
    discriminant = trJ^2 - 4 * detJ
    geometry = discriminant < -tol ? :focus :
               discriminant > tol ? :node : :degenerate
    return (
        equilibrium=z,
        jacobian=J,
        trace=trJ,
        determinant=detJ,
        eigenvalues=lambda,
        classification=classification,
        geometry=geometry,
    )
end

"""    hopf_currents(p)

Corrientes de traza cero cuando `0 < b < 1` y `tau > b`, con cruce transversal.
Retorna `nothing` cuando estas hipótesis no se cumplen. Este cálculo identifica
**candidatos teóricos**: la continuación numérica de ramas y ciclos corresponde
a la etapa B.
"""
function hopf_currents(p::FHNParameters)
    (p.b > 0 && p.tau > p.b) || return nothing
    s = sqrt(1 - p.b / p.tau)
    i_of_v(v) = v^3 / 3 - v + (v + p.a) / p.b
    return (
        v_minus=-s,
        v_plus=s,
        I_minus=i_of_v(-s),
        I_plus=i_of_v(s),
        omega=sqrt(p.tau - p.b^2) / p.tau,
    )
end

"""    equilibrium_sweep(p, currents)

Evalúa equilibrio y estabilidad al variar únicamente `I`. No realiza
continuación de órbitas periódicas.
"""
function equilibrium_sweep(p::FHNParameters, currents)
    Ivals = Float64.(collect(currents))
    isempty(Ivals) && throw(ArgumentError("El barrido no puede estar vacío."))
    all(isfinite, Ivals) || throw(ArgumentError("Las corrientes deben ser finitas."))

    vstars = Float64[]
    wstars = Float64[]
    traces = Float64[]
    determinants = Float64[]
    classes = Symbol[]
    for I in Ivals
        info = stability_info(with_current(p, I))
        push!(vstars, info.equilibrium[1])
        push!(wstars, info.equilibrium[2])
        push!(traces, info.trace)
        push!(determinants, info.determinant)
        push!(classes, info.classification)
    end
    return (I=Ivals, v=vstars, w=wstars, trace=traces,
            determinant=determinants, classification=classes)
end

"""    parameter_sweep(p, parameter, values)

Barrido del equilibrio y de la traza al modificar solamente uno de
`:a`, `:b`, `:tau` o `:I`. No integra las trayectorias.
"""
function parameter_sweep(p::FHNParameters, parameter::Symbol, values)
    vals = Float64.(collect(values))
    isempty(vals) && throw(ArgumentError("El barrido no puede estar vacío."))
    params = [with_parameter(p, parameter, val) for val in vals]
    infos = stability_info.(params)
    return (
        parameter=parameter,
        values=vals,
        v=[info.equilibrium[1] for info in infos],
        w=[info.equilibrium[2] for info in infos],
        trace=[info.trace for info in infos],
        determinant=[info.determinant for info in infos],
        classification=[info.classification for info in infos],
    )
end

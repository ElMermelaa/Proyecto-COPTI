"""
    FHNParameters(a, b, tau, I)
    FHNParameters(; a=0.7, b=0.5, tau=12.5, I=1.4)

Parámetros de una neurona FitzHugh–Nagumo del enunciado MA4703:

    v̇ = v - v³/3 - w + I
    τ ẇ = v + a - b*w

Hipótesis: `a ≥ 0`, `0 ≤ b < 1`, `tau > 0`, `I ∈ ℝ`.
Los valores predeterminados son **solo un ejemplo numérico**, no están fijados
por el enunciado ni por el cuaderno teórico.
"""
struct FHNParameters
    a::Float64
    b::Float64
    tau::Float64
    I::Float64

    function FHNParameters(a::Real, b::Real, tau::Real, I::Real)
        all(isfinite, (a, b, tau, I)) ||
            throw(ArgumentError("Todos los parámetros deben ser finitos."))
        a >= 0 || throw(ArgumentError("Se requiere a ≥ 0."))
        0 <= b < 1 || throw(ArgumentError("Se requiere 0 ≤ b < 1."))
        tau > 0 || throw(ArgumentError("Se requiere tau > 0."))
        return new(Float64(a), Float64(b), Float64(tau), Float64(I))
    end
end

FHNParameters(; a::Real=0.7, b::Real=0.5, tau::Real=12.5, I::Real=1.4) =
    FHNParameters(a, b, tau, I)

"""    with_current(p, I)

Devuelve parámetros nuevos cambiando únicamente la corriente externa.
"""
with_current(p::FHNParameters, I::Real) = FHNParameters(p.a, p.b, p.tau, I)

"""    with_parameter(p, name, value)

Copia parámetros modificando solamente `:a`, `:b`, `:tau` o `:I`.
Las restricciones siguen verificándose en `FHNParameters`.
"""
function with_parameter(p::FHNParameters, name::Symbol, value::Real)
    name in (:a, :b, :tau, :I) ||
        throw(ArgumentError("Parámetro desconocido: $name. Usa :a, :b, :tau o :I."))
    return FHNParameters(
        name === :a ? value : p.a,
        name === :b ? value : p.b,
        name === :tau ? value : p.tau,
        name === :I ? value : p.I,
    )
end

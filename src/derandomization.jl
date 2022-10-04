using Random: AbstractRNG, GLOBAL_RNG

"""
	hits(obs, op)

Return wether operator `op` hits observable `obs`.

# Examples

```jldoctest
julia> hits("XIZI", "XYZX")
true
```
"""
function hits(obs::String, op::String)::Bool
    if length(obs) * length(op) == 0
        throw(DomainError((obs, op), "zero length"))
    end
    if length(obs) != length(op)
        throw(DomainError((obs, op), "not the same length"))
    end
    f(x::Tuple{Char,Char}) = x[1] ∉ [x[2], 'I']
    filter(f, collect(zip(obs, op))) |> length == 0
end

"""
	weight(o)

Returns the weight of the given PauliId string.
"""
function weight(o::String)::Int
    filter(x -> x ∈ ('X', 'Y', 'Z'), o) |> length
end

"""
    partial_confidence_bound_expectation(obs, measurements, M, ε)

Computes the confidence bound expectation value over the remaining unspecified
measurements for the Pauli observables `obs` and the already chosen Pauli
measurements `measurements` over a total budget `M`, for an additive error `ε`.

See the paper to choose `M` and `ε` correctly.

Example:
```
    partial_confidence_bound_expectation(Set(["XIZ"]), ["XXX", "YYZ", "X"], 4, 0.5)
```
"""
function partial_confidence_bound_expectation(
    obs::Set{String},
    measurements::Vector{String},
    M::Int,
    ε::Float64,
)::Float64
    n = collect(obs)[1] |> length # number of qubits inferred from the observable list
    nu = 1 - exp(-(ε^2) / 2)
    result = 0

    m = length(measurements)
    k = last(measurements) |> length

    for o in obs
        # 1. Fixed measurements
        exp1 = sum([hits(o, meas) for meas in measurements[1:m-1]])
        contrib1 = exp(-ε^2 / 2 * exp1)

        # 2. Partially fixed measurement
        contrib2 = 1 - nu * hits(o[1:k], last(measurements)) * 3.0^(-weight(o[k+1:n]))

        # 3. Unspecified measurements
        contrib3 = (1 - nu * 3.0^(-weight(o)))^(M - m)

        result += contrib1 * contrib2 * contrib3
    end
    # TODO(lvignoli): write it as a function and broadcast-sum it on all the observables (nicer to read I think)
    return result
end

"""
    derandomization(obs, M, ε)

Compute the `M` Pauli measurements to operate to get estimators with additive
error `ε` for the expectation value of `obs` according to the derandomization
procedure.
"""
function derandomization(obs::Set{String}, M::Int, ε::Float64)::Vector{String}
    n = collect(obs)[1] |> length
    # TODO(lvignoli): the exepected shape of arguments of the confidence bound is
    # very inconvenient and force this implementation. Refactor it (custom
    # struct?)
    function _extend(P::Vector{String}, W::String)
        if length(P) == 0
            return [W]
        end
        if length(P[end]) < n
            return [P[1:end-1]; P[end] * W]
        else
            return [P; W]
        end
    end
    P♯::Vector{String} = []
    while length(foldl(*, P♯)) < n * M
        W = argmin(
            W -> partial_confidence_bound_expectation(obs, _extend(P♯, W), M, ε),
            ["X", "Y", "Z"],
        )
        P♯ = _extend(P♯, W)
    end
    return P♯
end

function check_confidence_bound(
    obs::Set{String},
    measurements::Vector{String},
    ε::Float64,
)::Tuple{Bool,Float64}
    a = sum(exp(-(ε^2 / 2) * sum(hits.(o, measurements))) for o in obs)
    δ = 2 * a
    return 0 < δ < 1, δ
end

function expect_derand(
    reg::AbstractRegister,
    obs::String,
    ops::Vector{String};
    rng::AbstractRNG = GLOBAL_RNG,
)::Float64
    res = 0
    n = length(obs)
    for op in ops
        s = apply(reg, pauli_basis_change(op))
        bs = measure(s; rng)[1]

        # Yao bitstrings are in little endian. Iterating on them correctly take
        # bits from right to left.
        if hits(obs, op)
            q = [b == 0 ? 1 : -1 for b in bs]  # Map bits to corresponding eigenvalue: 0 => 1, 1 => -1
            res += prod([q[i] for i in 1:n if obs[i] != 'I'])
        end
    end
    htot = hits.(obs, ops) |> sum
    res / htot
end

"""
Estimate the expectation values of the Pauli observables in the state `reg` using the Pauli
measurements `measurements`.
"""
function expect_derand(
    reg::AbstractRegister,
    observables::Set{String},
    measurements::Vector{String};
    rng::AbstractRNG = GLOBAL_RNG,
)::Dict{String,Float64}
    Dict(o => expect_derand(reg, o, measurements; rng) for o in observables)
end

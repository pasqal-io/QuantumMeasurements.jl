"""
    number_samples(k, precision, probability, M, n; robust)

Compute the number of samples recquired to estimate the expected values of `M`
that are `k` local with the desired `precision` with the given `probability`.
"""
function number_samples(
    k::Int, precision::Float64, probability::Float64, M::Int, n=1; robust=false
)::NTuple{4,Int64}
    δ = 1 - probability
    Ke = log(2 * M / δ)

    Ne = exp(log(4) * k) * 34 #See formal version of theorem 1 of arXiv:2002.08953v2
    Ne /= precision^2

    Nc = 0
    Kc = 0

    if robust
        Ne *= 1 #TODO(rmartin): Take noise into account and add the multiplying factor

        Kc = log(2 / δ) + k * log(n)

        Nc = Ne
        # Nc *= 2 * exp(log(3) * k)  It is too much, makes the program too slow without a substantial gain in precision
        Nc *= 1 #TODO(rmartin): Take noise into account and add the multiplying factor
    end

    Ne = trunc(Int, Ne)
    Ke = trunc(Int, Ke)
    Nc = trunc(Int, Nc)
    Kc = trunc(Int, Kc)
    return (Ne, Ke, Nc, Kc)
end

"""
    median_of_mean(estimations, N, K)

Split the array `estimations` in `K` equally-sized part. Compute the `K` means of these parts. Return the median of the `K` means.
"""
function median_of_mean(estimations::Array{Float64}, N::Int, K::Int)::Float64
    median([mean([estimations[r] for r = ((k - 1) * N + 1):(k * N)]) for k = 1:K])
end

"""
    single_estimate_obs(shadow, obs)

Based on one `shadows`, give an estimation of the observable `obs`'s expected value
which equals the expected value in average.
Uses the kronecker factorization.
"""
function single_estimate_obs(shadow::Shadow, obs::PauliObs)::Float64
    U = shadow.gate
    b = shadow.measure
    prod([
        real(tr(Matrix(mat(pauli_from_char[obs.pauli[j]]))' * chain_inverse(Matrix(mat(U[j])), b[j]))) for
        j in keys(obs.pauli)
    ])
end

"""
    estimate_obs(shadows, obs, N, K)

Based on `shadows` of a state, return an estimation of the observable `obs`'s expected value.
Uses the median of mean protocol to have a good estimation.
"""
function estimate_obs(shadows::Array{Shadow}, obs::PauliObs, Ne::Int, Ke::Int)::Float64
    R = Ne * Ke
    median_of_mean([single_estimate_obs(shadows[r], obs) for r = 1:R], Ne, Ke)
end

"""
    estimate_set_obs(shadows, observables, N, K)

Return a dictionary which assciociates to each observable of the set `observables`
an estimation of its expected value.
"""
function estimate_set_obs(
    shadows::Array{Shadow}, observables::Set{PauliObs}, Ne::Int, Ke::Int
)::Dict{String,Float64}
    Dict(PauliObs_to_string(o) => estimate_obs(shadows, o, Ne, Ke) for o in observables)
end

"""
    robust_estimate_set_obs(shadows, f, observables, N, K)

Return a dictionary which assciociates to each observable of the set `observables`
an estimation of its expected value.

NOTE(rmartin): Calculating each observable's weight (instead of taking the max) and compute only the usefull bit strings may be faster in some cases.
"""
function robust_estimate_set_obs(
    shadows::Array{Shadow}, calibration::Calibration, observables::Set{PauliObs}, Ne::Int, Ke::Int
)::Dict{String,Float64}
    Dict(PauliObs_to_string(o) => robust_estimate_obs(calibration, shadows, o, Ne, Ke) for o in observables)
end

"""
    robust_estimate_obs(f, shadows, observable, N, K, bit_string)

Return an estimation of the `observable`'s expected value taking noise into account.
It is computed as the median of `K` means, each one being a mean of `N` single estimation.
"""
function robust_estimate_obs(
    calibration::Calibration, shadows::Array{Shadow}, observable::PauliObs, Ne::Int, Ke::Int
)::Float64
    R = Ne * Ke
    median_of_mean(
        [robust_single_estimate(calibration, shadows[r], observable) for r = 1:R], Ne, Ke
    )
end

"""
    robust_single_estimate(f, bit_string, shadow, observable)

Based on a single `shadow` and noise correlation factors `f`, compute a single estimation
which equals the `observable`'s expected value in average.
"""
function robust_single_estimate(
    calibration::Calibration, shadow::Shadow, observable::PauliObs
)::Float64
    n = shadow.n_qubits
    f = calibration.coefficients
    bit_string = zeros(Int, n)
    terms = observable.pauli
    product = 1
    for j in keys(terms)
        bit_string[j] = 1
        U = Matrix(mat(shadow.gate[j]))
        b = shadow.measure[j] == 0 ? [1, 0] : [0, 1]
        Oj = Matrix(mat(pauli_from_char[terms[j]]))
        product *= b' * U * Oj * U' * b 
    end

    product /= f[bit_string]
    
    return real(product)
end

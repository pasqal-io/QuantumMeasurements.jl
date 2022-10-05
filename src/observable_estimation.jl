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
function single_estimate_obs(shadow::Shadow, obs::String)::Float64
    U = shadow.gate
    b = shadow.measure
    prod([
        real(tr(Matrix(mat(pauli_from_char[obs[j]]))' * chain_inverse(Matrix(mat(U[j])), b[j]))) for
        j = 1:(shadow.n_qubits)
    ])
end

"""
    estimate_obs(shadows, obs, N, K)

Based on `shadows` of a state, return an estimation of the observable `obs`'s expected value.
Uses the median of mean protocol to have a good estimation.
"""
function estimate_obs(shadows::Array{Shadow}, obs::String, Ne::Int, Ke::Int)::Float64
    R = Ne * Ke
    median_of_mean([single_estimate_obs(shadows[r], obs) for r = 1:R], Ne, Ke)
end

"""
    estimate_set_obs(shadows, observables, N, K)

Return a dictionary which assciociates to each observable of the set `observables`
an estimation of its expected value.
"""
function estimate_set_obs(
    shadows::Array{Shadow}, observables::Set{String}, Ne::Int, Ke::Int
)::Dict{String,Float64}
    Dict(o => estimate_obs(shadows, o, Ne, Ke) for o in observables)
end

"""
    robust_estimate_set_obs(shadows, f, observables, N, K)

Return a dictionary which assciociates to each observable of the set `observables`
an estimation of its expected value.

NOTE(rmartin): Calculating each observable's weight (instead of taking the max) and compute only the usefull bit strings may be faster in some cases.
"""
function robust_estimate_set_obs(
    shadows::Array{Shadow}, calibration::Calibration, observables::Set{String}, Ne::Int, Ke::Int
)::Dict{String,Float64}
    Dict(o => robust_estimate_obs(calibration, shadows, o, Ne, Ke) for o in observables)
end

"""
    robust_estimate_obs(f, shadows, observable, N, K, bit_string)

Return an estimation of the `observable`'s expected value taking noise into account.
It is computed as the median of `K` means, each one being a mean of `N` single estimation.
"""
function robust_estimate_obs(
    calibration::Calibration, shadows::Array{Shadow}, observable::String, Ne::Int, Ke::Int
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

The inverse of the quantum channel representing the measurement protocol is a sum of terms.
Each term is associated to a bitstring of {0, 1}ⁿ and must be multiplied by the corresponding correlation factor of `f`.
The only terms that contribute to the estimation are those for which the corresponding bitstring
does not have a "1" in the same position as an "I" in the `observable`.
This justifies why only bitstring with less than k "1" have been considered (where k is the maximum weight of observables)

NOTE(rmartin): Computing robust_single_estimate for each shadow is the bottleneck of the
implementation because its complexity is proportional to nᵏ.
"""
function robust_single_estimate(
    calibration::Calibration, shadow::Shadow, observable::String
)::Float64
    n = shadow.n_qubits
    estimations = 0
    f = calibration.coefficients
    for bit_string in keys(f)
        product = 1 / f[bit_string]
        for j = 1:n
            if observable[j] == 'I'
                if bit_string[j] == 1 #The corresponding term does not contribute.
                    product = 0
                    break
                end #No "else" beacause in the other case we only have to multiply product by one (I expected value)
            else
                if bit_string[j] == 0
                    product *= tr(Matrix(mat(pauli_from_char[observable[j]])))
                else
                    U = Matrix(mat(shadow.gate[j]))
                    b = shadow.measure[j] == 0 ? [1, 0] : [0, 1]
                    Oj = Matrix(mat(pauli_from_char[observable[j]]))
                    product *= b' * U * Oj * U' * b - tr(Oj) / 2
                end
            end
        end
        estimations += product
    end
    return real(estimations)
end

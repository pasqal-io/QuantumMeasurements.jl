# TODO(lvignoli): allow to specify the measurements, or cache the derandomized measurements,
# so that we do not recompute the derandomized measurement everytime, it is so slow.

# TODO(lvignoli,mdagrada): unify the exposed API for dispatching between tomography and
# derandomization. Ideas:
# - use Singleton type to specify the method. But we need additional parameters for config.
# - use a strategy pattern, while staying on idiomatic julia.

"""
    estimate_expect(state, observables, precision; rng)

# Arguments

  - `state::AbstractRegister`: a quantum state given as a Yao register. Usually constructed
    using a Yao block like so.
    
    ```julia
    s = zero_state(2)
    block = chain(3.1 * kron(X, Y), put(1 => H), -0.3 * kron(Z, Z))
    apply!(block, s)
    ```

  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values in `state`.
  - `precision::Float64`: the approximate additive error of the predicitons with respect to
    the true expectation values. It is not a hard upper bound. See Theorem 1 of
    arXiv:2103.07510 for details.

# Keywords

  - `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
    number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

  - `Dict{String,Float64}`: the expectation values in
    `state` for each observable in `observables`.
"""
function estimate_expect(
    state::AbstractRegister,
    observables::Set{String},
    precision::Float64;
    rng::AbstractRNG=GLOBAL_RNG,
)::Dict{String,Float64}
    L = length(observables)

    M = log(L) * max(exp.(weight.(observables))...)
    M /= precision^2
    M = trunc(Int, M)
    println("Using $M derandomized Pauli measurements")

    measurements = derandomization(observables, M, precision)

    ms = unique(measurements)
    measurements_count = Dict(m => count(==(m), measurements) for m in ms)

    println("Derandomized measurements:")
    println(measurements_count)

    expect_derand(state, observables, measurements; rng)
end

"""
    estimate_expect_shadows(state, observables, precision, probability; noise_model, robust, rng)

For each observable in `observables`, give an estimation of the expected value with `state`
that is ϵ accurate with a probability of at least δ.

# Arguments

  - `state::AbstractRegister`: a quantum state given as a Yao register. Usually constructed
    using a Yao block like so.
    
    ```julia
    s = zero_state(2)
    block = chain(3.1 * kron(X, Y), put(1 => H), -0.3 * kron(Z, Z))
    apply!(block, s)
    ```

  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values in `state`.
  - `precision::Float64`: the approximate additive error of the predicitons with respect to
    the true expectation values. It is not a hard upper bound. See formal version of Theorem 1 of
    arXiv:2002.08953v2 for details.
  - `probability::Float64` : the desired probability of having the estimations of expection values
    with the `precision` given as an input. The default value is 0.95 if unspecified.

# Keywords

  - `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
    measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
  - `robust::Bool`: If this boolean is true then the estimation is done with the method
    robust to noise taken as an argument (`noise_model`). It defaults on `false` if unspecified.
  - `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
    number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

  - `Dict{String,Float64}`: the expectation values in
    `state` for each observable in `observables`.
"""
function estimate_expect_shadows(
    state::AbstractRegister,
    observables::Set{String},
    precision::Float64,
    probability::Float64=0.95;
    noise_model::Noise=Noise(identity, identity),
    robust::Bool=false,
    rng::AbstractRNG=GLOBAL_RNG,
)::Dict{String,Float64}
    pauli_observables = Set([PauliObs(o) for o in observables])
    M = length(pauli_observables)
    k = max(weight.(pauli_observables)...)
    n = nqubits(state)

    Ne, Ke, Nc, Kc = number_samples(k, precision, probability, M, n; robust)

    println("Predictions have an accuracy of $precision with probability at least $probability")

    println("Using $(Ne*Ke) classical shadows split into $Ke equally-sized parts")
    shadows_list = classical_shadows(state, Ne * Ke; noise_model, rng)
    println("Classical shadows generated")

    if robust
        println("Generating $(Nc*Kc) noise shadows split into $Kc equally-sized parts")
        noise_shadows_list = classical_shadows(zero_state(n), Nc * Kc; noise_model, rng)
        println("Noise shadows generated")
        println("Computing calibration coefficients")
        bit_strings = observables_bit_strings(pauli_observables)
        f = Calibration(noise_shadows_list, bit_strings, Nc, Kc)
        println("Calibration done")

        println("Estimating expectation values")
        return robust_estimate_set_obs(shadows_list, f, pauli_observables, Ne, Ke)
    else
        println("Estimating expectation values")
        return estimate_set_obs(shadows_list, pauli_observables, Ne, Ke)
    end
end

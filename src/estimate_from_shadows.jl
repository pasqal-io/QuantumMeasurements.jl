"""
    estimate_from_shadows(shadows, observables, Re, Ke)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. Uses a median of mean estimation using `Ke` groups, each of `Ne` shadows.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `Re::Int`: The number of shadows to use to compute the average of the expected values in the
    median of mean estimation protocol.
  - `Ke::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow}, observables::Set{String}, Re::Int, Ke::Int=1
)::Dict{String,Float64}
    nb_shadows = length(shadows)
    Re > nb_shadows ? Re = nb_shadows : Re
    Ne = Re ÷ Ke
    pauli_observables = Set([PauliObs(o) for o in observables])
    println("Using $(Ne*Ke) classical shadows split into $Ke equally-sized parts")
    estimate_set_obs(shadows, pauli_observables, Ne, Ke)
end

"""
    estimate_from_shadows(shadows, observables, precision, probability)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. Uses a median of mean estimation.
The estimations are ϵ accurate with a probability of at least δ if enough shadows are given in input.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `precision::Float64`: the approximate additive error of the predicitons with respect to
    the true expectation values. It is not a hard upper bound.
  - `probability::Float64` : the desired probability of having the estimations of expection values
    with the precision given as an input. The default value is 0.95 if unspecified.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow}, observables::Set{String}, precision::Float64, probability::Float64=0.95
)::Dict{String,Float64}
    k = k = max(weight.(observables)...)
    M = length(observables)
    n = shadows[1].n_qubits
    Ne, Ke = number_samples(k, precision, probability, M, n; robust=false)

    if Ne * Ke <= length(shadows)
        println("Predictions have an accuracy of $precision with probability at least $probability")
    else
        println(
            "Not enough shadows to return an estimation with the desired precision and probability"
        )
    end
    estimate_from_shadows(shadows, observables, Ne * Ke, Ke)
end

"""
    estimate_from_shadows(shadows, calibration, observables, Re, Ke)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. Uses a median of mean estimation using `Ke` groups, each of `Ne` shadows.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `calibration::Calibration`: The calibration to have an estimation protocol that is robust to noise.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `Re::Int`: The number of shadows to use to compute the average of the expected values in the
    median of mean estimation protocol.
  - `Ke::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow}, calibration::Calibration, observables::Set{String}, Re::Int, Ke::Int=1
)::Dict{String,Float64}
    pauli_observables = Set([PauliObs(o) for o in observables])
    nb_shadows = length(shadows)
    Re > nb_shadows ? Re = nb_shadows : Re
    Ne = Re ÷ Ke

    println("Using $(Ne*Ke) classical shadows split into $Ke equally-sized parts")
    robust_estimate_set_obs(shadows, calibration, pauli_observables, Ne, Ke)
end

"""
    estimate_from_shadows(shadows, calibration, observables, precision, probability)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. Uses a median of mean estimation.
The estimations are ϵ accurate with a probability of at least δ if enough shadows are given in input.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `calibration::Calibration`: The calibration to have an estimation protocol that is robust to noise.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `precision::Float64`: the approximate additive error of the predicitons with respect to
    the true expectation values. It is not a hard upper bound.
  - `probability::Float64` : the desired probability of having the estimations of expection values
    with the precision given as an input. The default value is 0.95 if unspecified.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow},
    calibration::Calibration,
    observables::Set{String},
    precision::Float64,
    probability::Float64=0.95,
)::Dict{String,Float64}
    k = k = max(weight.(observables)...)
    M = length(observables)
    n = shadows[1].n_qubits
    Ne, Ke = number_samples(k, precision, probability, M, n; robust=false)

    if Ne * Ke <= length(shadows)
        println("Predictions have an accuracy of $precision with probability at least $probability")
    else
        println(
            "Not enough shadows to return an estimation with the desired precision and probability"
        )
    end
    estimate_from_shadows(shadows, calibration, observables, Ne * Ke, Ke)
end

"""
    estimate_from_shadows(shadows, noise_shadows, observables, Re, Ke, Rc, Kc)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. Uses a median of mean estimation using `Ke` groups, each of `Ne` shadows.
The estimation is robust to noise. It computes a Calibration with the shadows from `noise_shadows` using
`Nc`*`Kc` shadows split into `Kc` equally sized parts.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `noise_shadows::Array{Shadow}`: Shadows of the zero state to gain information of the noise and
    compute the Calibration. Can be generated with the function `generate_noise_shadows`.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `Re::Int`: The number of shadows to use to compute the average of the expected values in the
    median of mean estimation protocol.
  - `Ke::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the expected values' estimations.
  - 'Rc::Int': The number of noise shadows to use to compute the Calibration for a
    median of mean estimation protocol.
  - 'Kc::Int': The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow},
    noise_shadows::Array{Shadow},
    observables::Set{String},
    Re::Int,
    Ke::Int=1,
    Rc::Int=Re,
    Kc::Int=Ke,
)::Dict{String,Float64}
    nb_shadows = length(shadows)
    nb_noise_shadows = length(noise_shadows)
    Re > nb_shadows ? Re = nb_shadows : Re
    Rc > nb_noise_shadows ? Rc = nb_noise_shadows : Rc
    Ne = Re ÷ Ke
    Nc = Rc ÷ Kc

    println("Using $(Nc*Kc) classical shadows split into $Kc equally-sized parts for calibration")
    calibration = Calibration(noise_shadows, observables_bit_strings(observables), Nc, Kc)

    estimate_from_shadows(shadows, calibration, observables, Ne * Ke, Ke)
end

"""
    estimate_from_shadows(shadows, noise_shadows, observables, precision, probability)

For each observable in `observables`, give an estimation of the expected value for the state
used to generate the array `shadows`. The estimation is robust to noise. It computes
a Calibration with the shadows from `noise_shadows`.
The estimations are ϵ accurate with a probability of at least δ if enough shadows are given in input.

# Arguments

  - `shadows::Array{Shadow}`: The array of Shadow generated from the state
    in which observables' expected values must be estimated.
  - `noise_shadows::Array{Shadow}`: Shadows of the zero state to gain information of the noise and
    compute the Calibration. Can be generated with the function `generate_noise_shadows`.
  - `observables::Set{String}`: the set of Pauli string observables for which we estimate the
    expectation values.
  - `precision::Float64`: the approximate additive error of the predicitons with respect to
    the true expectation values. It is not a hard upper bound.
  - `probability::Float64` : the desired probability of having the estimations of expection values
    with the precision given as an input. The default value is 0.95 if unspecified.

# Returns

  - `Dict{String,Float64}`: the expectation values for each observable in `observables` in
    the state used to generate the shadows.
"""
function estimate_from_shadows(
    shadows::Vector{Shadow},
    noise_shadows::Vector{Shadow},
    observables::Set{String},
    precision::Float64,
    probability::Float64=0.95,
)::Dict{String,Float64}
    k = k = max(weight.(observables)...)
    M = length(observables)
    n = shadows[1].n_qubits

    Ne, Ke, Nc, Kc = number_samples(k, precision, probability, M, n; robust=true)
    if Ne * Ke <= length(shadows) && Nc * Kc <= length(noise_shadows)
        println("Predictions have an accuracy of $precision with probability at least $probability")
    else
        println(
            "Not enough shadows to return an estimation with the desired precision and probability"
        )
    end
    estimate_from_shadows(shadows, noise_shadows, observables, Ne * Ke, Ke, Nc * Kc, Kc)
end

"""
    calibrate(bit_strings, noise_shadows, Rc, Kc)

Return a Calibration using the noise shadows of the array `noise_shadows`. Uses a median of mean estimation
with `Nc`*`Kc` shadows split into `Kc equally` sized parts.

# Arguments

  - `bit_strings::Set{Array{Int}}`: The array of bitstrings for which the corresponding coeffcient will be  computed.
  - `noise_shadows::Array{Shadow}`: Shadows of the zero state to gain information of the noise and
    compute the Calibration. Can be generated with the function `generate_noise_shadows`.
  - `Rc::Int`: The number of noise shadows to use to compute the Calibration for
    median of mean estimation protocol.
  - `Kc::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Returns

  - `Dict{Array{Int},Float64}`: The coeffcient of the Calibration for each bitstring of the argument
    array `bit_strings`
"""
function calibrate(
    bit_strings::Set{Array{Int}}, noise_shadows::Array{Shadow}, Nc::Int, Kc::Int
)::Calibration
    nb_shadows = length(noise_shadows)
    Nc * Kc > nb_shadows ? Nc = nb_shadows ÷ Kc : Nc
    Calibration(noise_shadows, bit_strings, Nc, Kc)
end

"""
    calibrate(observables, noise_shadows, Rc, Kc)

Return a Calibration using the noise shadows of the array `noise_shadows`. Uses a median of mean estimation
with `Nc`*`Kc` shadows split into `Kc equally` sized parts. Only Calibration's coeffcients that contribute to
the estimation of observables' expected values are computed.

# Arguments

  - `observables::Set{String}`: the set of Pauli string observables for which the expected values will
    be estimated using this Calibration.
  - `noise_shadows::Array{Shadow}`: Shadows of the zero state to gain information of the noise and
    compute the Calibration. Can be generated with the function `generate_noise_shadows`.
  - `Rc::Int`: The number of noise shadows to use to compute the Calibration for
    median of mean estimation protocol.
  - `Kc::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Returns

  - `Dict{Array{Int},Float64}`: The coeffcient of the Calibration for each bitstring that contribute to
    the estimation of observables' expected values.
"""
function calibrate(
    observables::Set{String}, noise_shadows::Array{Shadow}, Nc::Int, Kc::Int
)::Calibration
    pauli_observables = Set([PauliObs(o) for o in observables])
    bit_strings = observables_bit_strings(pauli_observables)
    calibrate(bit_strings, noise_shadows, Nc, Kc)
end

"""
    calibrate(observables, Rc, Kc; noise_model, rng)

Return a Calibration for the noise model given as an argument. Uses a median of mean estimation
with `Nc`*`Kc` shadows split into `Kc equally` sized parts. Only Calibration's coeffcients that contribute to
the estimation of observables' expected values are computed.

# Arguments

  - `observables::Set{String}`: the set of Pauli string observables for which the expected values will
    be estimated using this Calibration.
  - `Rc::Int`: The number of noise shadows to use to compute the Calibration for
    median of mean estimation protocol.
  - `Kc::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Keywords

  - `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
    measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
  - `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
    number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

  - `Dict{Array{Int},Float64}`: The coeffcient of the Calibration for each bitstring that contribute to
    the estimation of observables' expected values.
"""
function calibrate(
    observables::Set{String},
    Nc::Int,
    Kc::Int;
    noise_model::Noise=Noise(identity, identity),
    rng::AbstractRNG=GLOBAL_RNG,
)::Calibration
    n = length(first(observables))
    noise_shadows = classical_shadows(zero_state(n), Nc * Kc; noise_model, rng)
    calibrate(observables, noise_shadows, Nc, Kc)
end

"""
    calibrate(observables, noise_shadows, Rc, Kc)

Return a Calibration using the noise shadows of the array `noise_shadows`. Uses a median of mean estimation
with `Nc`*`Kc` shadows split into `Kc equally` sized parts. Only Calibration's coeffcients that contribute to
the estimation of expected values for `k`-local or less observables are computed.

# Arguments

  - `k::Int`: The maximum weight of observables which expected values will be estimated using the Calibration.
  - `noise_shadows::Array{Shadow}`: Shadows of the zero state to gain information of the noise and
    compute the Calibration. Can be generated with the function `generate_noise_shadows`.
  - `Rc::Int`: The number of noise shadows to use to compute the Calibration for
    median of mean estimation protocol.
  - `Kc::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Returns

  - `Dict{Array{Int},Float64}`: The coeffcient of the Calibration for each bitstring that contribute to
    the estimation of the expected value of an observable which locality is `k` or less.
"""
function calibrate(k::Int, noise_shadows::Array{Shadow}, Nc::Int, Kc::Int)::Calibration
    n = first(noise_shadows).n_qubits
    bit_strings = bit_string_k(n, k)
    calibrate(bit_strings, noise_shadows, Nc, Kc)
end

"""
    calibrate(observables, Rc, Kc; noise_model, rng)

Return a Calibration for the noise model given as an argument. Uses a median of mean estimation
with `Nc`*`Kc` shadows split into `Kc equally` sized parts. Only Calibration's coeffcients that contribute to
the estimation of expected values for `k`-local or less observables are computed.

# Arguments

  - `k::Int`: The maximum weight of observables which expected values will be estimated using the Calibration.
  - `Rc::Int`: The number of noise shadows to use to compute the Calibration for
    median of mean estimation protocol.
  - `Kc::Int`: The number of means to compute to then take the median for the median of mean estimation
    protocol for the Calibration.

# Keywords

  - `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
    measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
  - `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
    number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

  - `Dict{Array{Int},Float64}`: The coeffcient of the Calibration for each bitstring that contribute to
    the estimation of the expected value of an observable which locality is `k` or less.
"""
function calibrate(
    n::Int,
    k::Int,
    Nc::Int,
    Kc::Int;
    noise_model::Noise=Noise(identity, identity),
    rng::AbstractRNG=GLOBAL_RNG,
)::Calibration
    noise_shadows = classical_shadows(zero_state(n), Nc * Kc; noise_model, rng)
    calibrate(k, noise_shadows, Nc, Kc)
end

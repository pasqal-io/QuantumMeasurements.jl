"""
    generate_classical_shadows(state, R; noise_model, rng)

Return an array of `R` classical shadows of `state` that can be used to estimate the expected value of any observable.

# Arguments

- `state::AbstractRegister`: a quantum state given as a Yao register. Usually constructed
  using a Yao block like so.
  ```julia
  s = zero_state(2)
  block = chain(3.1*kron(X,Y), put(1=>H), -0.3*kron(Z,Z))
  apply!(block, s)
  ```
- `R::Int` : the desired number of classical shadows.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of `R` classical shadows of the `state`. 
"""
function generate_classical_shadows(
    state::AbstractRegister,
    R::Int;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    classical_shadows(state, R; noise_model, rng)
end

"""
    generate_classical_shadows(state, k, M, precision, probability; noise_model, rng)

Return an array of classical shadows of `state` that can be used to estimate the expected value of some observables.
Estimations based on these shadows will have the desired `precision` with the given `probability` 
for a set of `M` observables that are `k` local.
See formal version of Theorem 1 of arXiv:2002.08953v2 for details.

# Arguments

- `state::AbstractRegister`: a quantum state given as a Yao register. Usually constructed
  using a Yao block like so.
  ```julia
  s = zero_state(2)
  block = chain(3.1*kron(X,Y), put(1=>H), -0.3*kron(Z,Z))
  apply!(block, s)
  ```
- `k::Int`: The locality of the observables which expected values will be estimated using these classical shadows.
  "IIXIZ" is a 2 local observable.
- `M::Int`: The number of observables which expected values will be estimated using these classical shadows.
- `precision::Float64`: the approximate additive error of the predicitons with respect to
  the true expectation values. It is not a hard upper bound.
- `probability::Float64` : the desired probability of having the estimations of expection values
  with the precision given as an input. The default value is 0.95 if unspecified.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of `R` classical shadows of the `state`. 
"""
function generate_classical_shadows(
    state::AbstractRegister,
    k::Int,
    M::Int,
    precision::Float64,
    probability::Float64 = 0.95;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    Ne, Ke = number_samples(k, precision, probability, M)

    generate_classical_shadows(state, Ne * Ke; noise_model, rng)
end

"""
    generate_classical_shadows(state, observables, precision, probability; noise_model, rng)

Return an array of classical shadows of `state` that can be used to estimate the expected value of some observables.
Estimations based on these shadows will have the desired `precision` with the given `probability` 
for observables in the set `observables`.
See formal version of Theorem 1 of arXiv:2002.08953v2 for details.

# Arguments

- `state::AbstractRegister`: a quantum state given as a Yao register. Usually constructed
  using a Yao block like so.
  ```julia
  s = zero_state(2)
  block = chain(3.1*kron(X,Y), put(1=>H), -0.3*kron(Z,Z))
  apply!(block, s)
  ```
- `observables::Set{String}`: the set of Pauli string observables for which we estimate the
  expectation values in `state`.
- `precision::Float64`: the approximate additive error of the predicitons with respect to
  the true expectation values. It is not a hard upper bound. 
- `probability::Float64` : the desired probability of having the estimations of expection values
  with the precision given as an input. The default value is 0.95 if unspecified.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of `R` classical shadows of the `state`. 
"""
function generate_classical_shadows(
    state::AbstractRegister,
    observables::Set{String},
    precision::Float64,
    probability::Float64 = 0.95;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    M = length(observables)
    k = max(weight.(observables)...)

    generate_classical_shadows(state, k, M, precision, probability; noise_model, rng)

end

"""
    generate_noise_shadows(n, R; noise_model, rng)

Return an array of `R` noise shadows that are classical shadows of the zero state. They can be used 
to compute coefficients to obtain estimations of any obsrvable's expected value that are robust to the noise.

# Arguments

- `n::Int`: the number of qu-bits of the system.
- `R::Int` : the desired number of noise shadows.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of `R` noise shadows (which are classical shadows of the
   zero state) that give information on the noise.
"""
function generate_noise_shadows(
    n::Int,
    R::Int;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    generate_classical_shadows(zero_state(n), R; noise_model, rng)

end

"""
    generate_noise_shadows(n, k, M, precision, probability; noise_model, rng)

Return an array of noise shadows that are classical shadows of the zero state. They can be used 
to compute coefficients to obtain estimations of any obsrvable's expected value that are robust to the noise.
Estimations based on coefficients based on these noise shadows will have the desired `precision` 
with the given `probability` for a set of `M` observables that are `k` local.
See Theorem 10 of arXiv:2011.09636v2 for details.

# Arguments

- `n::Int`: the number of qu-bits of the system.
- `k::Int`: The locality of the observables which expected values will be estimated using these classical shadows.
  "IIXIZ" is a 2 local observable.
- `M::Int`: The number of observables which expected values will be estimated using these classical shadows.
- `precision::Float64`: the approximate additive error of the predicitons with respect to
  the true expectation values. It is not a hard upper bound. 
- `probability::Float64` : the desired probability of having the estimations of expection values
  with the precision given as an input. The default value is 0.95 if unspecified.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of noise shadows (which are classical shadows of the
   zero state) that give information on the noise.
"""
function generate_noise_shadows(
    n::Int,
    k::Int,
    M::Int,
    precision::Float64,
    probability::Float64 = 0.95;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    sample_number = number_samples(k, precision, probability, M, n; robust = true)
    Nc, Kc = sample_number[3], sample_number[4]

    generate_noise_shadows(n, Nc * Kc; noise_model, rng)
end

"""
    generate_noise_shadows(n, observables, precision, probability; noise_model, rng)

Return an array of noise shadows that are classical shadows of the zero state. They can be used 
to compute coefficients to obtain estimations of any obsrvable's expected value that are robust to the noise.
Estimations based on coefficients based on these noise shadows will have the desired `precision`
with the given `probability` for observables in the set `observables`.
See Theorem 10 of arXiv:2011.09636v2 for details.

# Arguments

- `n::Int`: the number of qu-bits of the system.
- `observables::Set{String}`: the set of Pauli string observables for which we estimate the
  expectation values in `state`.
- `precision::Float64`: the approximate additive error of the predicitons with respect to
  the true expectation values. It is not a hard upper bound. 
- `probability::Float64` : the desired probability of having the estimations of expection values
  with the precision given as an input. The default value is 0.95 if unspecified.

# Keywords

- `noise_model::Noise`: use this keyword to make measurements noisy. No noise is applied to
  measurements if unspecified. Functions from the file noise.jl generate common noise_model functions.
- `rng::AbstractRNG`: use this keyword to see the measurement process with a given random
  number generator for repeatable results. It defaults on `GLOBAL_RNG` if unspecified.

# Returns

-  `Array{Shadow}`: An array of noise shadows (which are classical shadows of the
   zero state) that give information on the noise.
"""
function generate_noise_shadows(
    n::Int,
    observables::Set{String},
    precision::Float64,
    probability::Float64 = 0.95;
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Array{Shadow}

    M = length(observables)
    k = max(weight.(observables)...)
    generate_noise_shadows(n, k, M, precision, probability; noise_model, rng)
end

using Random: AbstractRNG, GLOBAL_RNG
using LinearAlgebra, Statistics

"""
    Shadow

Represent a single classical shadow of a state which consists in a measurement and
a unitary gate that has been applied before the measurement.
`gate` is the unitary gate that had been applied before the measurement which result is `measure`.
`gate` is represented by the n factors of its kronecker factorization.
"""
struct Shadow
    n_qubits::Int
    gate::Array{Yao.AbstractBlock}
    measure::Array{Int8}

    function Shadow(
        state::AbstractRegister;
        noise_model::Noise=Noise(identity, identity),
        rng::AbstractRNG=GLOBAL_RNG,
    )
        n = nqubits(state)
        gate = random_clifford_gate(n; rng)
        b = pauli_measure(gate, state; noise_model)
        new(n, gate, b)
    end
end

"""
    Calibration

Store the coeffcients necessary to compute the inverse of the noisy protocol measurement (seen as a quantum channel).
To eachbit string of {0,1}ⁿ corresponds a coeffcient.
See equations (15) and (16) of arXiv:2002.08953v2 for more details.

The inner constructor returns a `Calibration` where every bitstring of the argument `bit_strings`
is added to the dictionary but with a corresponding value set to 0.
"""
struct Calibration
    n_qubits::Int
    coefficients::Dict{Array{Int},Float64}

    function Calibration(bit_strings::Set{Array{Int}})::Calibration
        n = length(first(bit_strings))
        coefficients = Dict{Array{Int},Float64}([z => 0.0 for z in bit_strings])
        new(n, coefficients)
    end

    function Calibration(n::Int)
        new(n, Dict{Array{Int},Float64}())
    end
end

"""
    add_single_calibration!(noise_shadow, calibration)

For each bitstring in the dictionary of `calibration` compute the coeffcient corresponding to the `shadow`
and add it to the existing value. In average this computation equals the factors of the inverse channel.
"""
function add_single_calibration!(noise_shadow::Shadow, calibration::Calibration)
    b = noise_shadow.measure #Noisy measurement from a known state (usually the zero state)
    gate = noise_shadow.gate
    n = noise_shadow.n_qubits

    #Results of noisy measurements allow to compute noise calibration coefficients from these factors
    corelation_factors = [
        if b[j] == 0
            [1, 0]' * Matrix(mat(gate[j])) * [1 0; 0 -1] * Matrix(mat(gate[j]))' * [1, 0]
        else
            [0, 1]' * Matrix(mat(gate[j])) * [1 0; 0 -1] * Matrix(mat(gate[j]))' * [0, 1]
        end for j = 1:n
    ]

    for bit_string in keys(calibration.coefficients)
        coeff = 1.0
        for j = 1:n
            if bit_string[j] == 1
                coeff *= corelation_factors[j]
            end
        end
        calibration.coefficients[bit_string] += coeff
    end
end

"""
    mean_calibration!(noise_shadows, calibration)

For each bitstring in the dictionary of `calibration` compute the coeffcient corresponding to the `shadow`
by calculating the mean of coeffcients computed for each shadow of the array `noise_shadows`.
"""
function mean_calibration!(noise_shadows::Array{Shadow}, calibration::Calibration)
    R = length(noise_shadows)
    for shadow in noise_shadows
        add_single_calibration!(shadow, calibration)
    end
    for bit_string in keys(calibration.coefficients)
        calibration.coefficients[bit_string] /= R
    end
end

"""
    Calibration(noise_shadows, bit_strings, N, K)

Return a Calibration computing each coeffcient corresponding to a bitstring of `bit_string` threw
a median of mean estimation splitting the shadows of the array `noise_shadows` into `K` equally sized parts.
"""
function Calibration(
    noise_shadows::Array{Shadow}, bit_strings::Set{Array{Int}}, Nc::Int, Kc::Int
)::Calibration
    n = first(noise_shadows).n_qubits
    calibrate!(Calibration(n), bit_strings, noise_shadows, Nc, Kc)
end

"""
    calibrate!(calibration, bit_strings, noise_shadows, N, K)

Modify and return the argument `calibration` computing each coeffcient corresponding to a bitstring of `bit_string` threw
a median of mean estimation splitting the shadows of the array `noise_shadows` into `K` equally sized parts.
"""
function calibrate!(
    calibration::Calibration,
    bit_strings::Set{Array{Int}},
    noise_shadows::Array{Shadow},
    Nc::Int,
    Kc::Int,
)
    mean_calibrations = Array{Calibration}(undef, Kc)
    for t = 1:Kc
        mean_calibrations[t] = Calibration(bit_strings)
        mean_calibration!(noise_shadows[((t - 1) * Nc + 1):(t * Nc)], mean_calibrations[t])
    end
    for bit_string in bit_strings
        calibration.coefficients[bit_string] = median([
            mean_calibrations[t].coefficients[bit_string] for t = 1:Kc
        ])
    end
    calibration
end

"""
    random_clifford_gate(n)

Generate an array of `n` random unitary Clifford gates.
"""
function random_clifford_gate(n::Int; rng::AbstractRNG=GLOBAL_RNG)::Array{Yao.AbstractBlock}
    gate = Array{Yao.AbstractBlock}(undef, n)
    for j = 1:n
        U1 = pauli_from_char[rand(rng, ('I', 'X', 'Y', 'Z'))]
        U2 = pauli_from_char[rand(rng, ('I', 'V', 'W', 'H', 'K', 'L'))]
        gate[j] = chain(U2, U1)
    end
    return gate
end

"""
    pauli_measure(gates, s; noise_model)

Return a measure of the state `s` after applying the unitary gate represented by the array `gates` of its kronecker factors.
It applies the noise functions to the state after the unitary gate and to the measurement result.
If no `noise_model` is specified, no noise is applied to the measurement.

# Examples

bit_flip_noise(0.1) returns a noise_model that flips each qubit with propability O.1

'''jldoctest
julia> pauli_measure(random_clifford_gate(2), ArrayReg(bit"00"), noise_model=bit_flip_noise(0.1))
2-element Vector{Int8}:
1
0
'''
"""
function pauli_measure(
    gates::Array{Yao.AbstractBlock},
    s::AbstractRegister;
    noise_model::Noise=Noise(identity, identity),
)::Array{Int8}
    state = copy(s)
    n = nqubits(state)

    #Applying the random unitary gate
    U = chain(n, repeat(gates[i], [i]) for i = 1:n)
    apply!(state, U)

    #Applying the noise model on the state
    noise_model.noise_circuit(state)

    m = measure(state, (1:n))[1][1:n]

    #Applying the noise model on the measurement
    noise_model.noise_measure(m)

    return m
end

"""
    classical_shadows(state, R; noise_model, rng)

Generate an array of `R` classical snapshots of an unknown `state`.
It gives a classical representation of `state` that can be used
to predict expection values of any set of observables.
"""
function classical_shadows(
    state::AbstractRegister,
    R::Int;
    noise_model::Noise=Noise(identity, identity),
    rng::AbstractRNG=GLOBAL_RNG,
)::Array{Shadow}
    [Shadow(state; noise_model, rng) for _ = 1:R]
end

"""
    observables_bit_strings(observables)   

Return a set containing only the bitstrings of {0,1}ⁿ that contribute to the estimation of the expected value of
at least one observable in the set `observables`. A bitstring contribute to the estimation of an observable's expected value
if there is no 1 in the bitstring in the same position as a I in the observable.

# Examples

```jldoctest
julia> QuantumMeasurements.observables_bit_strings(Set(["XIII", "IIZY"]))
Set{Array{Int64}} with 5 elements:
  [0, 0, 0, 1]
  [0, 0, 0, 0]
  [0, 0, 1, 1]
  [0, 0, 1, 0]
  [1, 0, 0, 0]
```
"""
function observables_bit_strings(observables::Set{String})::Set{Array{Int}}
    s = Set{Array{Int}}()
    n = length(first(observables))
    for o in observables #For each observable in observables we add to the set s the bitstrings that contribute to the estimation of its expected value.
        k = weight(o)
        bit_string = [zeros(Int, n) for _ = 1:(2^k)]

        non_I = Array{Int}(undef, k) #Index of qubits where the observable o acts non trivally.
        i = 1
        for j = 1:n
            if o[j] != 'I'
                non_I[i] = j
                i += 1
            end
        end
        # At an index where the observable's factor is not I, the bitstring
        # contribute to the inverse channel whether it is 0 or 1, so, after
        # setting the bits in the same poistion as I to 0, we generate all the
        # possibilities for the other positions.
        bit = bitarray(collect(0:(2^k - 1)), k)
        for z = 1:(2^k)
            for j = 1:k
                if bit[j, z]
                    bit_string[z][non_I[j]] = 1
                end
            end
        end

        push!(s, bit_string...)
    end
    return s
end

"""
    bit_string_k(n, k)

Return an Set that contains every bitstring of {0, 1}ⁿ which number of 1 is less or equal to `k`.

# Examples

```jldoctest
julia> QuantumMeasurements.bit_string_k(3, 2)
Set{Array{Int64}} with 7 elements:
  [0, 1, 0]
  [1, 0, 1]
  [1, 0, 0]
  [0, 0, 1]
  [0, 0, 0]
  [1, 1, 0]
  [0, 1, 1]
```
"""
function bit_string_k(n::Int, k::Int)::Set{Array{Int}}
    O_string = [0 for _ = 1:n]
    bit_strings = Set{Vector{Int}}([O_string])
    recursive_bit_string_k(O_string, 0, 0, k, bit_strings)
    return bit_strings
end

"""
    recursive_bit_string_k(string, lastOne, nbOne, k, bit_strings)

Add to the array `bit_strings` all the bit strings obtained from `string` replacing one 0 after the last 1 of `string`
"""
function recursive_bit_string_k(
    string::Array{Int}, lastOne::Int, NbOne::Int, k::Int, bit_strings::Set{Vector{Int}}
)
    n = length(string)
    if !(lastOne >= n || NbOne >= k)
        for i = (lastOne + 1):n
            new_string = copy(string)
            new_string[i] = 1
            push!(bit_strings, new_string)
            recursive_bit_string_k(new_string, i, NbOne + 1, k, bit_strings)
        end
    end
end

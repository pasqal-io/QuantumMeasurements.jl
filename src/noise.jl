using Random: AbstractRNG, GLOBAL_RNG

"""
    Noise

Simulate a noisy operation and a noisy measurement.
`noise_circuit` is a function that takes a state as unique argument and modfies it to
simulate a noisy gate.
`noise_measure` is a function that takes a Array of 0 and 1 and modifies it to simulate
a noisy measurement. It can also simulate a noisy gate followed by a noisy measurement as a whole.

"""
struct Noise
    noise_circuit::Function
    noise_measure::Function
end

"""
	noise_composition(noise_models)

Return a Noise struct which is a composition of each Noise of the array noise_models.

"""
function noise_composition(noise_models::Array{Noise})::Noise
    function noise_circuit!(state::AbstractRegister)
        for Λ in noise_models
            Λ.noise_circuit(state)
        end
    end
    function noise_measure!(measrue::Array{Int})
        for Λ in noise_models
            Λ.noise_measure(measrue)
        end
    end
    return Noise(noise_circuit!, noise_measure!)
end

"""
	bit_flip_noise(rate, subsystem)

Return a Noise that randomly flip each qubit which index is in the array `subsystem` with probability `rate`.
If subsystem is unspecified, the Noise acts on every single qubit of the state.

"""
function bit_flip_noise(rate::T, subsystem::Array{Int})::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        bit_flip = [rand() for _ in subsystem]
        Λ = chain(n)
        for i in eachindex(subsystem)
            #Adding bit-flip to the noise channel
            if bit_flip[i] < rate
                Λ = chain(Λ, put(n, subsystem[i] => X))
            end
        end
        apply!(state, Λ)
    end
    return Noise(noise!, identity)
end

function bit_flip_noise(rate::T)::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        uniform_noise = bit_flip_noise(rate, [1:n...])
        uniform_noise.noise_circuit(state)
    end
    return Noise(noise!, identity)
end

"""
    depolarization_noise(rate)

Return a Noise that randomly depolarize each qubit which index is in the array `subsystem` with probability `rate`.
If subsystem is unspecified, the Noise acts on every single qubit of the state.

"""
function depolarization_noise(rate::T, subsystem::Array{Int})::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        Λ = chain(n)
        depolarizing = [rand() for _ in subsystem]
        for i in eachindex(subsystem)
            #Adding depolarizing to the noise channel
            if depolarizing[i] < rate / 4
                Λ = chain(Λ, put(n, subsystem[i] => X))
            elseif depolarizing[i] < rate / 2
                Λ = chain(Λ, put(n, subsystem[i] => Y))
            elseif depolarizing[i] < 3 * rate / 4
                Λ = chain(Λ, put(n, subsystem[i] => Z))
            end
        end
        apply!(state, Λ)
    end
    return Noise(noise!, identity)
end

function depolarization_noise(rate::T)::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        uniform_noise = depolarization_noise(rate, [1:n...])
        uniform_noise.noise_circuit(state)
    end
    return Noise(noise!, identity)
end

"""
	amplitude_damping_noise(rate)

Return a Noise that change the measurement so that a result of 1 is changed to 0 
with probability `rate` and result of 0 are unchanged.

"""
function amplitude_damping_noise(rate::T)::Noise where {T<:Real}
    function noise!(measure::Array{Int})
        n = length(measure)
        random_draw = [rand() for _ = 1:n]
        for j = 1:n
            if measure[j] == 1 && random_draw[j] < rate
                measure[j] = 0
            end
        end
    end
    return Noise(identity, noise!)
end

"""
	amplitude_damping_noise_gate(rate)

Return a function that change a qubit from the state 1 to 0 with probability `rate` and does not modifiy a qubit in state 0.

"""
function amplitude_damping_noise_gate(rate::T)::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        Λ = chain(2 * n)
        θ = 2 * asin(sqrt(rate))
        append_qubits!(state, n)

        #Adding Amplitude damping noise
        A(i, j) = chain(2 * n, control(i, j => rot(Y, θ)), control(j, i => X))
        B = chain(A(i, i + n) for i = 1:n)
        Λ = chain(Λ, B)
        Λ = chain(Λ, Measure(2 * n, locs = n+1:2*n))

        apply!(state, Λ)
    end
    return Noise(noise!, identity)
end

"""
	X_error(θ)

Return a Noise that apply a rotation of `θ` in the X axis to each qubit.

"""
function X_error(θ::T)::Noise where {T<:Real}
    function noise!(state::AbstractRegister)
        n = nqubits(state)
        Λ = chain(n, put(i => rot(X, θ)) for i = 1:n)
        apply!(state, Λ)
    end
    return Noise(noise!, identity)
end

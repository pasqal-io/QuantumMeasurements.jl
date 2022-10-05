"""
	chain_inverse(Pauli, measure)

The measurement protocol acts as a quantum channel on the state.
This function computes its inverse : it takes the measurement gate and outcome to return an
estimation of a kronecker factor of the density matrix.

"""
function chain_inverse(gate::Matrix, measure::Int8)::Matrix
    measure == 0 ? M = [1 0; 0 0] : M = [0 0; 0 1]
    return 3 * gate' * M * gate - [1 0; 0 1]
end

"""
	single_estimate(shadow)

Generate single classical snapshot of the unknown state from a single `shadow`.
This snapshot exactly reproduces the underlying state in expectation.

"""
function single_estimate(shadow::Shadow)::Matrix{ComplexF64}
    n = shadow.n_qubits
    gate = shadow.gate
    b = shadow.measure
    if n > 1
        ρ = [chain_inverse(Matrix(mat(gate[j])), b[j]) for j = 1:n]
        ρt = kron([ρ[n-j+1] for j = 1:n]...)
    else
        ρt = chain_inverse(Matrix(mat(gate[1])), b[1])
    end
    return ρt
end

"""
    robust_single_estimate_state(shadow, f, bit_string)

Generate single classical snapshot of the unknown state from a single shadow based on the noise calibration `f`.
This snapshot exactly reproduces the underlying state in expectation.

"""
function robust_single_estimate_state(
    shadow::Shadow,
    calibration::Calibration,
)::Matrix{ComplexF64}
    n = shadow.n_qubits
    gate = shadow.gate
    b = shadow.measure

    f = calibration.coefficients

    ρ = Set{Matrix}()
    for z in keys(f)
        ρj = Array{Matrix}(undef, n)
        for j = 1:n
            vec = b[j] == 0 ? [1, 0] : [0, 1]

            #Compute the inverse of the channel measurement to return an estimation of the state's density matrix
            Π(A::Matrix) = tr(A) * [1 0; 0 1] / 2
            U = Matrix(mat(gate[j]))

            if z[j] == 0
                P = Π(U' * vec * vec' * U)
            else
                P = U' * vec * vec' * U - Π(U' * vec * vec' * U)
            end
            ρj[j] = P
        end
        push!(ρ, kron([ρj[n-j+1] for j = 1:n]...) / f[z])
    end
    sum(ρ)
end

"""
state_reconstruction(state, R; robust, noise_model)

Return an estimation of `state`'s density matrix using `R` classical shadows of that `state`.
If `robust` == true, uses the robust algorithm.

"""
function state_reconstruction(
    state::AbstractRegister,
    R::Int;
    robust::Bool = false,
    noise_model::Noise = Noise(identity, identity),
    rng::AbstractRNG = GLOBAL_RNG,
)::Matrix{ComplexF64}
    shadows = classical_shadows(state, R; noise_model, rng)
    println("Classical shadows generated")
    n = nqubits(state)
    if robust
        noise_shadow_list = classical_shadows(zero_state(n), R; noise_model, rng)
        state_reconstruction(shadows, noise_shadow_list, R)
    else
        state_reconstruction(shadows, R)
    end
end

"""
state_reconstruction(shadows, R)

Return an estimation of the state's density matrix used to generate the classical shadows of the array `shadows` 
Uses `R` classical shadows of the array `shadows` and a non robust estimation method.

"""
function state_reconstruction(shadows::Array{Shadow}, R::Int)::Matrix{ComplexF64}
    nb_shadows = length(shadows)
    R = R > nb_shadows ? nb_shadows : R
    mean([single_estimate(shadows[r]) for r = 1:R])
end

"""
state_reconstruction(shadows, calibration, R)

Return an estimation of the state's density matrix used to generate the classical shadows of the array `shadows` 
Uses `R` classical shadows of the array `shadows` and the `calibration` to use a robust estimation method.

"""
function state_reconstruction(
    shadows::Array{Shadow},
    calibration::Calibration,
    R::Int,
)::Matrix{ComplexF64}
    nb_shadows = length(shadows)
    R = R > nb_shadows ? nb_shadows : R
    mean([robust_single_estimate_state(shadows[r], calibration) for r = 1:R])
end

"""
state_reconstruction(shadows, calibration, R)

Return an estimation of the state's density matrix used to generate the classical shadows of the array `shadows` 
Uses `R` classical shadows of the array `shadows` and the `R` noise shadows to compute a calibration and
use a robust method.

"""
function state_reconstruction(
    shadows::Array{Shadow},
    noise_shadows::Array{Shadow},
    R::Int,
)::Matrix{ComplexF64}
    n = first(shadows).n_qubits
    nb_shadows = length(noise_shadows)
    Rc = R > nb_shadows ? nb_shadows : R
    calibration = Calibration(noise_shadows, bit_string_k(n, n), Rc, 1)
    println("Calibration done")
    state_reconstruction(shadows, calibration, R)
end

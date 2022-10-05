"""
    filter_pauli_operations(obs::KronBlock)

Given an observable represented as a Kroenecker block of operators, this function
filter the qubit locations of the blocks containing a Pauli operation (X, Y or Z)
and not the trivial identity operator. The output of this function is a Vector{Int64}
with the qubit indices where Pauli operations are applied in the input observable
"""
function filter_pauli_operations(obs::KronBlock, op_to_exclude::Vector{T}) where {T<:PauliGate}
    locs::Vector{Int} = []
    for (l, b) in zip(obs.locs, obs.blocks)
        if b ∉ op_to_exclude
            push!(locs, l[1])
        end
    end
    return locs
end

filter_pauli_operations(obs) = filter_pauli_operations(obs, [I2])

"""
    rotate!(circuit, locs::Vector{Int}, op::PauliGate)

Rotate the input circuit inplace in order to measure the corresponding non-local observable
given by the `op` argument
"""
function rotate(circuit, locs::Vector{Int}, op::T) where {T<:PauliGate}
    n_qubits = nqubits(circuit)
    for i in locs
        if op == X
            circuit = chain(n_qubits, circuit, put(i => Ry(-π / 2)))
        elseif op == Y
            circuit = chain(n_qubits, circuit, put(i => Rx(π / 2)))
        end
    end

    return circuit
end

"""
    rotate_to_measurement_basis(term::KronBlock, circuit::AbstractBlock)

Add a final rotation to a quantum circuit for measuring non-local observables in the Z basis which
is usually the only one accessible to real devices.

Return a circuit with the corresponding rotations already applied.
"""
function rotate_to_measurement_basis(term::KronBlock, circuit::AbstractBlock)
    rot_circuit = copy(circuit)

    # perform rotations for terms containing the X Pauli operation
    rot_circuit = rotate(rot_circuit, filter_pauli_operations(term, [I2, Z, Y]), X)

    # perform rotations for terms containing the Y Pauli operation
    rot_circuit = rotate(rot_circuit, filter_pauli_operations(term, [I2, Z, X]), Y)

    return rot_circuit
end

"""
    compute_expectation_value(meas::Vector{T}, observable::KronBlock) where{T<:DitStr}

Compute the expectation value of a given observable built as a tensor product of Pauli
operators starting from a set of Z-basis measurements
"""
function compute_expectation_value(meas::Vector{T}, observable::KronBlock) where {T<:DitStr}
    budget = length(meas)
    locs = filter_pauli_operations(observable)

    val = map(meas) do m
        (-1)^(sum(m[locs]))
    end |> sum

    return val / budget
end

"""
    state_tomography(observables::Vector{AbstractBlock}, circuit::AbstractBlock; budget::Int = 1000)

Perform a naive tomography procedure which allows to reconstruct the state of a quantum system starting
from its classical (bitstring) representations measured in the Z basis. This function takes as input a
(parametrized) quantum circuit U(θ) and a list of observables and uses the tomography procedure to
evaluate the quantum mechanical expectation value of each observable ̂O defined as:

    ψ = U(θ)|ψ0>
    E[̂O] = <ψ|O|ψ>
"""
function expect_tomography(
    observables::Vector{T}, circuit::AbstractBlock; budget::Int=1000
) where {T<:AbstractBlock}
    values = []

    for obs in observables

        # collect the various term on the observable
        pauli_obs = PauliObservable(obs)

        # compute the expectation value for each term and sum them together
        n = 0
        value = map(pauli_obs.terms) do term
            n += 1
            coeff = pauli_obs.coefficients[n]

            # first find non Z operators and perform the necessary circuit rotations
            rot_circuit = rotate_to_measurement_basis(term, circuit)

            # perform the measurement with the given budget
            ψ0 = zero_state(nqubits(term))
            ψ = apply(ψ0, rot_circuit)
            meas = measure(ψ, nshots=budget)

            # compute the expectation value
            val = compute_expectation_value(meas, term)
            val * coeff
        end |> sum
        push!(values, value)
    end

    return values
end

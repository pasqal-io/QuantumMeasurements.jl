using QuantumMeasurements
using Test

@testset "YaoConverter" begin
    n_qubits = 2

    # create a random circuit
    circuit = chain(n_qubits, [])

    # add two layers of single-qubit rotations
    layer_1 = chain(n_qubits, put(1 => Rx(1.0)), put(2 => Rx(1.0)))
    layer_2 = chain(n_qubits, put(1 => Ry(2.0)), put(2 => Ry(2.0)))
    layer_3 = chain(n_qubits, put(1 => Rz(3.0)), put(2 => Rz(3.0)))

    # add controlledX (CNOT)
    control_gate = chain(n_qubits, control(1, 2 => X))

    # chain all together
    push!(circuit, chain(n_qubits, layer_1, layer_2, layer_3, control_gate))
    plot(circuit)
end

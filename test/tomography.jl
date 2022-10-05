@testset "Tomography with Z basis only" begin
    budget = 1000000
    n_qubits = 2

    bs = []
    for i = 1:n_qubits
        b = chain(n_qubits, put(i => Rz(2.0)), put(i => Rx(2.0)))
        push!(bs, b)
    end
    circuit = chain(n_qubits, bs...)

    O1 = pauli_string(2, 2 => Z)
    O2 = pauli_string(Z, Z)

    observables = [O1, O2]
    values_tomo = expect_tomography(observables, circuit, budget=budget)

    ψ0 = zero_state(n_qubits)
    ψ = apply(ψ0, circuit)
    values_exact = [expect(O1, ψ), expect(O2, ψ)]

    @test isapprox(values_tomo, values_exact, rtol=1e-2, atol=1e-2)
end

@testset "Tomography complete operators" begin
    budget = 1000000
    n_qubits = 5

    bs = []
    for i = 1:n_qubits
        b = chain(n_qubits, put(i => Ry(2.0)), put(i => Ry(2.0)))
        push!(bs, b)
    end
    circuit = chain(n_qubits, bs...)

    O1 = pauli_string(n_qubits, 2 => X)
    O2 = pauli_string(Z, X, X, Z, Z)
    O3 = 2.0 * pauli_string(n_qubits, 1 => Y, 3 => Z) + 3.0 * pauli_string(X, X, Z, I2, I2)

    observables = [O1, O2, O3]
    values_tomo = expect_tomography(observables, circuit, budget=budget)

    ψ0 = zero_state(n_qubits)
    ψ = apply(ψ0, circuit)
    values_exact = [expect(O1, ψ), expect(O2, ψ), expect(O3, ψ)]

    @test isapprox(values_tomo, values_exact, rtol=1e-2, atol=1e-2)
end

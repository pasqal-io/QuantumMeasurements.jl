@testset "Correct promotion rule and type hierarchy" begin
    @test promote_type(XGate, YGate) == PauliGate
    @test [I2, X, Y, Z] isa Vector{PauliGate}
    @test PauliGate <: ConstantGate
end

@testset "No gate on invalid qubits (inherited from Yao.kron)" begin
    @test_throws ErrorException pauli_string(0, 1 => X)
    @test_throws ErrorException pauli_string(1, 2 => X)
    @test_throws LocationConflictError pauli_string(1, 1 => X, 1 => X)
end

@testset "No operators other than Pauli matrices and identity" begin
    @test_throws MethodError pauli_string(1, 1 => H)
    @test_throws MethodError pauli_string(1, 1 => Rx(1.0))
    @test_throws MethodError pauli_string(H)
    @test_throws MethodError pauli_string(H, X)
    @test_throws MethodError pauli_string(X, Rx(1.0))
end

@testset "Good number of qubits" begin
    @test pauli_string(I2, X, Y, Z).n == 4
    @test pauli_string(4, 2 => X).n == 4
end

@testset "Can create Hamiltonians as sum of scaled Pauli strings" begin
    @test 2 * pauli_string(X) isa Any
    @test pauli_string(X) + pauli_string(Y) isa Any
    @test 3 * pauli_string(X, I2, Y) + 4.5 * pauli_string(3, 2 => Z) -
          pauli_string(Y, I2, I2) isa Any
end

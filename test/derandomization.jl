@testset "Hits" begin
    @test hits("XX", "XX") == true
    @test hits("XI", "XX") == true
    @test hits("IX", "XX") == true
    @test hits("II", "XX") == true
    @test hits("ZZ", "XX") == false
    @test hits("IZ", "XX") == false
    @test hits("ZI", "XX") == false
    @test hits("ZY", "XX") == false

    @test hits("XI", "XX") == true
    @test hits("XI", "XY") == true
    @test hits("XI", "XZ") == true
    @test hits("IIIIII", "IIXIZY") == true
    @test hits("IIIIZY", "IIXIZY") == true
    @test hits("IIIIZX", "IIXIZY") == false
    @test hits("XIXIZY", "IIXIZY") == false

    @test_throws DomainError hits("", "X")
    @test_throws DomainError hits("", "")
end

@testset "Weight" begin
    @test weight("I") == 0
    @test weight("X") == 1
    @test weight("XYZ") == 3
    @test weight("IXYZIXYZIXYZ") == 9
end

@testset "Confidence bound expectation" begin
    ε0 = √2
    ε1 = 0.987654321
    nu(x::Float64) = 1 - exp(-(x^2) / 2)

    cases = []

    # 1 qubit
    push!(
        cases,
        (Set(["I"]), ["X"], 2, ε0, exp(-2)),
        (Set(["Z"]), ["X"], 2, ε0, 1 - nu(ε0) / 3),
        (Set(["X", "Z"]), ["X"], 2, ε0, ((3 - nu(ε0)) * (2 - nu(ε0))) / 3),
        (Set(["X", "Z"]), ["X", "Y"], 3, ε0, (1 + exp(-1)) * (1 - nu(ε0) / 3)),
    )

    # 2 qubits
    push!(
        cases,
        # IX
        (Set(["IX"]), ["X"], 1, ε0, 1 - nu(ε0) / 3),
        (Set(["IX"]), ["Y"], 1, ε0, 1 - nu(ε0) / 3),
        (Set(["IX"]), ["Z"], 1, ε0, 1 - nu(ε0) / 3),
        (Set(["IX"]), ["X"], 2, ε0, (1 - nu(ε0) / 3)^2),
        (Set(["IX"]), ["XY"], 2, ε0, 1 - nu(ε0) / 3),
        (Set(["IX"]), ["XX"], 2, ε0, (1 - nu(ε0)) * (1 - nu(ε0) / 3)),
        (Set(["IX"]), ["XY", "Z"], 2, ε0, 1 - nu(ε0) / 3),
        (Set(["IX"]), ["XX", "X"], 2, ε0, (1 / exp(1)) * (1 - nu(ε0) / 3)),
        (Set(["IX"]), ["ZX", "YX", "X"], 4, ε0, (exp(-2) * (1 - nu(ε0) / 3)^2)),
        # XX
        (Set(["XX"]), ["X"], 1, ε0, 1 - nu(ε0) / 3),
        (Set(["XX"]), ["Y"], 1, ε0, 1),
        (Set(["XX"]), ["Z"], 1, ε0, 1),
    )

    # 2 2-qubit observables
    push!(
        cases,
        (
            Set(["XY", "ZI"]),
            ["XY", "ZY", "X"],
            4,
            ε0,
            exp(-1) * (1 - nu(ε0) / 3) * (1 + (1 - nu(ε0) / 3^2)),
        ),
    )

    # Test dependency in budget M
    push!(
        cases,
        (Set(["IX"]), ["XX", "X"], 12, ε0, (1 / exp(1)) * (1 - nu(ε0) / 3)^11),
        (Set(["I"]), ["X"], 100, ε0, exp(-100)),
        (Set(["Z"]), ["X"], 100, ε0, (1 - nu(ε0) / 3)^99),
    )

    # Test dependency in precision ε
    push!(
        cases,
        (Set(["IX"]), ["XX"], 2, 0.1234, (1 - nu(0.1234)) * (1 - nu(0.1234) / 3)),
        (Set(["IX"]), ["XX", "Z"], 2, 0.1234, exp(-(0.1234^2) / 2) * (1 - nu(0.1234) / 3)),
    )

    # Big tests: multiple observables with large budget and many measurements
    push!(
        cases,
        (
            Set(["XIY", "IIZ", "XXI"]),
            ["XXY", "YXZ", "ZY"],
            4,
            ε0,
            exp(-1) * (2 * (1 - nu(ε0) / 3^2) + (1 - nu(ε0) / 3)^2),
        ),
        (
            Set(["XIY", "IIZ", "XXI"]),
            ["XXY", "YXZ", "ZY"],
            4,
            ε1,
            exp(-(ε1^2) / 2) * (2 * (1 - nu(ε1) / 3^2) + (1 - nu(ε1) / 3)^2),
        ),
    )

    for (observables, assigned_measurements, M, ε, want) in cases
        got = partial_confidence_bound_expectation(observables, assigned_measurements, M, ε)
        @test got ≈ want
    end
end

@testset "Derandomization" begin
    cases = [
        (Set(["X"]), 2, √2, ["X", "X"]),
        (Set(["Y"]), 2, √2, ["Y", "Y"]),
        (Set(["Z"]), 2, √2, ["Z", "Z"]),
        (Set(["XZ", "YI", "IX"]), 2, √2, ["YX", "YX"]),
        (Set(["XZ", "YI", "IX"]), 3, √2, ["YX", "XZ", "YX"]),
    ]

    for (observables, M, ε, want) in cases
        got = derandomization(observables, M, ε)
        @test got == want
    end
end

@testset "Statistics on derandomization" begin
    cases = [
        (Set(["X"]), 1000, 0.1, Dict("X" => 1000)),
        (Set(["X", "Y"]), 1000, 0.1, Dict("X" => 500, "Y" => 500)),
        (Set(["X", "Y", "Z"]), 999, 0.1, Dict("X" => 333, "Y" => 333, "Z" => 333)),
        (Set(["XY", "ZZ"]), 1000, 0.1, Dict("XY" => 500, "ZZ" => 500)),
        (Set(["XX", "XI"]), 1000, 0.1, Dict("XX" => 1000)),
        (Set(["XXY", "YYZ"]), 1000, 0.1, Dict("XXY" => 500, "YYZ" => 500)),
    ]

    for (observables, M, ε, expected_meas_count) in cases
        measurements = derandomization(observables, M, ε)
        meas_count = Dict([o => count(==(o), measurements) for o in unique(measurements)])
        @test meas_count == expected_meas_count
    end
end

@testset "Confidence bound" begin
    cases = [
        (Set(["X"]), ["X", "X"], 0.1, 2 * exp(-0.1^2 / 2 * 2)),  # <1
        (
            Set(["XZ", "YX", "IY"]),
            ["XZ", "YX", "XY", "XY"],
            0.1,
            2 * (2 * exp(-0.1^2 / 2) + exp(-0.1^2 / 2 * 2)),
        ),  # >1
    ]
    for (obs, measurements, ε, δ_true) in cases
        ok, δ = check_confidence_bound(obs, measurements, ε)
        @test δ ≈ δ_true
        @test ok == (δ_true < 1)
    end
end

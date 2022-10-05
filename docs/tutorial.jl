### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 725a9ee8-b419-4de7-a27a-f3ffb12c5033
begin
    using Pkg
    Pkg.activate("..")
    Pkg.instantiate()

    using QuantumMeasurements
    using Yao
end

# ╔═╡ 2b9cb58b-5162-4aba-919d-9fb6c1d38146
md"""
## Write Pauli observables

QuantumMeasurements.jl offers a simple primitive for representing observables made of Pauli operators. If the observables are written in this way, they can be fed into all available QuantumMeasurements.jl measurement procedures. 

For example, let's write an observable acting on two qubits with a Pauli $Z$ operator

$\hat{O} = Z_1Z_2$
"""

# ╔═╡ c256d1a8-a774-4826-8803-7c1982dca70c
begin
    obs_simple = pauli_string(Z, Z)
    Matrix(mat(obs_simple))
end

# ╔═╡ 01df233e-ad25-4180-9588-222461de9c94
md"""
It is also easy to add more complex operators composed by multiple factors. For example, given a 5 qubit system, let's write the following operator:


$\hat{O} = Z_1Z_2 + 3 Z_1X_2X_3 + 4 Y_4Y_5$
"""

# ╔═╡ 073eca16-7619-4715-b6bf-fa109c47c41b
begin
    n_qubits = 5
    obs_full =
        pauli_string(n_qubits, 1 => Z, 2 => Z) +
        3 * pauli_string(n_qubits, 1 => Z, 2 => X, 3 => X) +
        4 * pauli_string(n_qubits, 4 => Y, 5 => Y)
end

# ╔═╡ 849d1a8d-9629-4705-a2ca-b666c4b6e76c
md"""
## Derandomized Pauli measurements

Following reference arxiv:2103:07510, QuantumMeasurements.jl offers estimators for the expectation values of a list of observables based on derandomized Pauli measurements.
"""

# ╔═╡ 8e48294f-2344-4c91-9fde-d45102e166c7
md"""
Let us pick a few 3-qubit Pauli observables.
We would like to estimate their expectation value in some quantum state $\ket{\psi}$ with the least number of measurements.
"""

# ╔═╡ 2324c375-fca3-40c4-9c8a-45ae6e3fe83c
observables = Set(["XXI", "IYZ", "IZI"])

# ╔═╡ 8ec5215a-b6b3-406b-9205-71424e85c78d
md"""
These Pauli observables are incompatible with each other: when measuring in the XXX basis, we gain information only on the first expectation value.
Similarly, measuring in the, say, XZZ basis, we gain information only for the last expectation value.

Derandomization of Pauli measurements provides a deterministic measurement setup ensuring the estimation of all the expectation values equally with the fewest possible measurements.
"""

# ╔═╡ 1f63a4f6-4e0a-43e1-907c-4857135c6e50
md"Let's construct a random state on these 3 qubits"

# ╔═╡ 39a08f0d-9b95-4e92-9fb8-e40e75a1afa4
s = rand_state(3);

# ╔═╡ 9cab82a6-91f8-4526-854d-acaf3e244a3b
statevec(s)

# ╔═╡ 0a86771b-a71b-4791-bb3b-9b601624f2f5
md"The true expectation values of the three observables in `s` are:"

# ╔═╡ e12daa31-7e9c-4a9b-8fea-af4f409f37a9
expect_vals_true = Dict(o => expect(pauli_string(o), s) for o in observables)

# ╔═╡ 69b8ae08-611a-4bd7-9703-43129e8bfa67
md"""
We can estimate these values using derandomized Pauli measurements using the `estimate_expect()` method of `QuantumMeasurements.jl`.
The `precision` argument is the additive error we would like to achieve with respect to the true values.
It is not a hard upper-bound, but rather the measurements are expected to fall in this error range with a high probability.

The computation of the necessary Pauli measurement is made on the fly and may be a bit slow.
"""

# ╔═╡ c2f945e0-185a-4cf1-8d93-8eab6e8fd3ca
expect_vals = estimate_expect(s, observables, 0.05)

# ╔═╡ 54153dc8-9f42-4961-99c8-df0b1bd6b039
md"""
It takes about 3300 measurements to estimate the expectations values of these 3 observables with an error close to 0.05.

We can check that the estimation is close enough to the true values.
"""

# ╔═╡ 37175781-d191-40cb-bfda-b37158ef4a62
for o in observables
    println(abs(expect_vals[o] - expect_vals_true[o]))
end

# ╔═╡ c930571c-947e-4506-aaa3-3d258f6d6af6
md"""
## Quantum state tomography

QuantumMeasurements.jl offers also an estimator based on standard quite state tomography. Here, the expectation value of the Pauli observable is computed directly from the Z-basis measurements given a certain number of shots.

The interface is similar to the derandomization estimator. Let's first choose a number of qubit and a measurement budget and the construct few observables to measure.
"""

# ╔═╡ 3edb6f3f-fccf-4f7a-848e-2258080d5939
n = 5

# ╔═╡ 1f2b388d-11a8-4605-a073-f9c80d40559a
budget = 1000000

# ╔═╡ f02dda39-ccf0-4f88-bc58-e2994343948e
begin
    O1 = pauli_string(n, 2 => X)
    O2 = pauli_string(Z, X, X, Z, Z)
    O3 = 2.0 * pauli_string(n, 1 => Y, 3 => Z) + 3.0 * pauli_string(X, X, Z, I2, I2)
    nothing
end

# ╔═╡ 7c376299-1685-4ad5-9dc0-f96a280ba54a
md"""
Then construct a simple quantum circuit for with some single-qubit rotations for preparing our quantum state.
"""

# ╔═╡ 23dcbe82-0b2c-490e-bdfd-e30a9315169d
begin
    bs = []
    for i = 1:n
        b = chain(n, put(i => Ry(2.0)), put(i => Ry(2.0)))
        push!(bs, b)
    end
    circuit = chain(n, bs...)
    nothing
end

# ╔═╡ 7a67930f-a149-4282-a460-9620ddb63f44
md"""
Let's now compute the expectation values of the observables defined above using quantum state tomography and compare with the expected value computed using statevector propagation.
"""

# ╔═╡ e7e4a1a8-eebd-4304-88f8-1ce10b80a7f6
begin
    all_obs = [O1, O2, O3]
    expect_vals_tomo = expect_tomography(all_obs, circuit, budget=budget)
end

# ╔═╡ 28a439fd-ac55-4c23-bacb-68571fce1e4a
begin
    ψ0 = zero_state(n_qubits)
    ψ = apply(ψ0, circuit)
    expect_vals_tomo_true = [expect(O1, ψ), expect(O2, ψ), expect(O3, ψ)]
end

# ╔═╡ b41d0b52-72a4-4f13-8585-c09b6b66b1f3
for i = 1:length(all_obs)
    println(abs(expect_vals_tomo[i] - expect_vals_tomo_true[i]))
end

# ╔═╡ 5af67659-0651-4f6c-8bde-83e63eead6da
md"""
## Classical Shadows

Following reference arXiv:2002.08953v2, QuantumMeasurements offers estimators for the expectation values of a list of observables based on classical shadows.
Random unitary matrices that rotate the state before measurement are randomly chosen from the local Clifford group $Cl_2^{\otimes n}$. See subsection C of the section 5 for more details.
"""

# ╔═╡ b97abec7-79a3-4d9f-8a47-e45509c31cd1
md"""
Let us pick a few 4-qubit Pauli observables.
We would like to estimate their expectation value in some quantum state $\ket{\psi}$.
"""

# ╔═╡ b6d249d5-5a18-4754-b91e-278330013244
obs = Set(["XXII", "ZIZI", "YIIZ"])

# ╔═╡ e319a652-6ad6-4cae-9193-784e744338df
md"Let's construct a random state on these 4 qubits"

# ╔═╡ b933ee85-9f12-42b3-9224-64f282320982
ρ = rand_state(4)

# ╔═╡ 3d80105b-122e-440c-8ad0-8c01d9aa8e91
md"The true expectation values of the three observables in `s` are:"

# ╔═╡ 90079df8-f0e1-4de0-8622-270bceca778f
expected_vals_true = Dict(o => expect(pauli_string(o), ρ) for o in obs)

# ╔═╡ 875f2714-0d84-47a5-8dd4-82841b00841d
md"""
We can estimate these values with random classical shadows using the `estimate_expect_shadows()` method of `QuantumMeasurements`.
The `precision` argument is the additive error we would like to achieve with respect to the true values.
It is not a hard upper-bound, but rather the measurements are expected to fall in this error range with a probability given as argument. If unspecified the probability is set to 0.95
"""

# ╔═╡ 66ff89a5-d3c0-4fe8-b726-ef73a319e070
expected_vals = estimate_expect_shadows(ρ, obs, 0.1)

# ╔═╡ 3ead0cc8-dcce-4369-be88-da54e7635867
md"""
It takes about 217596 measurements to estimate the expectations values of these 3 observables. This represents quite a lot of measurements but this is the minimum required to ensure the precision desired. The estimations will often be a lot more precise.

We can check that the estimation is close enough to the true values.
"""

# ╔═╡ a6caa431-6fa3-42d2-b719-23782dd3172e
for o in obs
    println(abs(expected_vals[o] - expected_vals_true[o]))
end

# ╔═╡ 92903d46-06b0-4008-8fb4-7d1186e96b7a
md"""
A strong property of classical shadows is that the recquired number of shadows to be $\epsilon$ precise with probability $1-\delta $ for M observbales with maximum weight (number of non identity operator) of $k$ is $\mathcal{O}(log(\frac{M}{\delta})\frac{3^k}{\epsilon^2})$ which is independant of the number of qubits.

Let's estimate some observables with the same weight as before but on 6-qubits state.
"""

# ╔═╡ bd46282f-9280-4af1-9a2f-a76e5bbfe3e4
begin
    obs_6 = Set(["XXIIII", "ZIIIZI", "IIYXII"])
    ρ_6 = rand_state(6)
    estimate_expect_shadows(ρ_6, obs_6, 0.1)
end

# ╔═╡ d069a1cd-7af7-4499-9122-aeae3fd3832c
md"""
We can observe that 217596 classical shadows have been used as for previous estimations
"""

# ╔═╡ 154a4e22-9126-4034-bad3-742c1841e5b2
md"""
One other strengh of classical shadows is that they can be generated before determining which observables they will be used to compute expected values.
`generate_classical_shadows()` and `estimate_from_shadows()` methods of `QuantumMeasurements.jl` allow to manipulate shadows.
"""

# ╔═╡ 78fad7f8-7c98-487f-8205-75b3fc853e87
md"""
Let's generate 100000 classiscal shadows of the state ρ.
"""

# ╔═╡ 6886b7f3-5086-43d7-ba87-fe975cf32de0
shadows = generate_classical_shadows(ρ, 100000)

# ╔═╡ a156a87e-98fe-459a-a804-7177dd41cbdc
md"""
We can again use these shadows to estimate expected values of observables from the set `obs` using for exemple only 50000 classical shadows split into 10 equally sized parts for a median of mean estimation.
"""

# ╔═╡ 20b0ea20-dfd9-4dc2-8352-85a699c5d5d7
estimate_from_shadows(shadows, obs, 50000, 10)

# ╔═╡ 4d207e08-7c31-4f1b-874b-009d84881f08
md"""
It is possible to specify some knowledge on observables when generating classical shadows or specify what is the desired precision of the estimation so the number of required shadows is automatically determined.

Let's see some exemples :
"""

# ╔═╡ 3bd66405-bf3f-499c-82dd-9b1e7c5ba37b
md"For $4$ observables that are $2$-local with a precision of 0.1:"

# ╔═╡ 64e0f826-bb1c-4dd7-b331-e776ff1c9d75
sh = generate_classical_shadows(ρ, 2, 3, 0.1)

# ╔═╡ 57aea7f4-39f8-44b0-bd36-9e28dafb30b7
length(sh)

# ╔═╡ a48f84b1-48d7-429a-a652-ace28b0569b7
md"""
We then determine our observables and estimate their expected values:
"""

# ╔═╡ ad7bd5e0-0a9a-4158-9b8c-6a21126909cb
begin
    obs_set = Set(["IIIX", "IIXI", "IIZZ"])
    estimate_from_shadows(sh, obs_set, 0.1)
end

# ╔═╡ 80ff8bf0-1917-4971-95e5-a86add6f7be1
md"""
However if we want to estimate the expected values of $3$-local observables or if we want a more precised estimation, the number of classical shadows generated will not be enough (an estimation using all available classical shadows is then returned) :
"""

# ╔═╡ a60a84e2-c7fb-4ae5-8724-38062fab25cb
estimate_from_shadows(sh, obs_set, 0.01)

# ╔═╡ 5d63fd01-9674-4fcf-8b14-bdbbe5392040
begin
    obs_3_local = Set(["ZZIX", "YYXI", "IIZZ"])
    estimate_from_shadows(sh, obs_3_local, 0.1)
end

# ╔═╡ cbff155f-5a6b-49af-b7f5-da2604863397
md"""
Classical shadows allow to reconstruct the density matrix of the unknown state. However the time complexity is exponential in the number of qubits $n$.
"""

# ╔═╡ d667206c-49a2-4a28-a996-3261badb75bf
ρ_matrix = state_reconstruction(ρ, 20000)

# ╔═╡ d2b2e48e-6d6a-43cc-bdab-eb47fa376a35
begin
    using LinearAlgebra
    function fidelity(state, matrix)::Float64
        fidelity(density_matrix(state).state, matrix)
    end

    function fidelity(m1::Matrix, m2::Matrix)::Float64
        A = sqrt(m1)
        return real((tr(sqrt(A * m2 * A))))
    end

    fidelity(ρ, ρ_matrix)
end

# ╔═╡ 2e3cfbff-bb55-4cce-b7c2-75f5638b2d13
md"We can eavluate the fidelity between the state and its reconstruction :"

# ╔═╡ 8604baa1-5be4-4519-92f9-2a8608359f55
md"""
## Noisy gates and measurements

QuantumMeasurements.jl offers the possibility to add noise each time a gate is applied or a measurement is done.
Methods from QuantumMeasurements.jl allow to construct common noise models :
"""

# ╔═╡ 9ab1889c-6d27-49a1-9354-7ecf166bdd4e
bf = bit_flip_noise(0.1)

# ╔═╡ 6320d315-35e5-4597-9b2f-55f00c7b10fe
md"""
When applied, each quibit has a probability of $0.1$ to be flipped :
"""

# ╔═╡ 71148961-f04c-4716-a510-33a3dd6a3c51
begin
    N = 10000
    n_qubit = 3
    average_measure = zeros(Float64, n_qubit)
    for i = 1:N
        zero = zero_state(n_qubit)
        bf.noise_circuit(zero)
        average_measure += measure(zero)[1][1:n_qubit]
    end
    average_measure /= N
end

# ╔═╡ 610844f6-3593-496f-9249-45e665bf8189
md"""
For each single qubit, we measure $1$ instead of $0$ 10% of the time.

We can also have several noise models :
"""

# ╔═╡ 216ab5ff-fa3a-457f-bd19-0aebbc4c9ec5
noise_model = noise_composition([bit_flip_noise(0.1), depolarization_noise(0.2)])

# ╔═╡ 7a582eea-291f-4a4c-9403-90a677779fe8
md"""
Noise that applies only on a subpart of the system can be defined specifying which qubits it will affect :
"""

# ╔═╡ f40de983-54a0-4112-abc5-8b6fec2774b5
begin
    noise_12 = bit_flip_noise(0.1, [1, 2]) #Bit flip with probability of 0.1 on qubits 1 and 2
    noise_3 = bit_flip_noise(0.2, [3]) #Bit flip with probability of 0.2 on qubit 3
    noise_not_uniform = noise_composition([noise_12, noise_3])

    average_meas = zeros(Float64, n_qubit)

    for i = 1:N
        zero = zero_state(n_qubit)
        noise_not_uniform.noise_circuit(zero)
        average_meas += measure(zero)[1][1:n_qubit]
    end
    average_meas /= N
end

# ╔═╡ c751189e-05e8-4e6e-aa7d-41d1ba104dfe
md"""
Estimations using classical shadows can be performed with noisy gates and measurements specifying `noise_model` as a keyword argument:
"""

# ╔═╡ 48b2092b-d965-4001-b0b7-d111e2ca4fc1
noisy_expected_vals = estimate_expect_shadows(s, observables, 0.1; noise_model)

# ╔═╡ 6e98e5bf-2bdc-4076-9c0d-86a4dfe8f348
md"""
Due to the noise, the estimation is not accurate anymore. The gap with theorical expected values being :
"""

# ╔═╡ d2a9c1cb-f4fa-4245-be3b-1e1782e2e703
begin
    noisy_expected_vals_true = Dict(o => expect(pauli_string(o), s) for o in observables)
    for o in observables
        println(abs(noisy_expected_vals[o] - noisy_expected_vals_true[o]))
    end
end

# ╔═╡ 1ddd531c-8754-45b8-b8ae-bd63a10fda88
md"""
To tackle this issue QuantumMeasurements.jl offers a robust classical shadows algorithm.
"""

# ╔═╡ c4706f67-b7a1-4dfd-b24d-b3eb25140be9
md"""
## Robust Classical Shadows

Following reference arXiv:2011.09636v2, QuantumMeasurements offers a robust estimation protocol based on classical shadows.

This protocol can be easily used using the keyword argument `robust` set to true :
"""

# ╔═╡ 59e89559-1df4-4d7a-bfdd-e445b634a78d
begin
    robust_expected_vals = estimate_expect_shadows(s, observables, 0.1; noise_model, robust=true)
    println("\n Gap between estimations and true values")
    for o in observables
        println(abs(robust_expected_vals[o] - noisy_expected_vals_true[o]))
    end
end

# ╔═╡ 4fb21ef8-fe2a-4d18-a7c6-77d4ae7d55ab
md"""
Robust classical shadows protocol is based on a calibration which consists in coeffcients that modify the expression of the function that maps an classical shadow to an estimation. This calibration is implemented by a Julia struct. It can be generated using the `calibrate()` method with different set of arguments depending on what is known by the user.

"""

# ╔═╡ 1647ae66-0667-4427-9118-fc1838f23276
calibrate(observables, 10000, 10; noise_model)

# ╔═╡ 22e8fb1a-3b90-4c57-950f-f4aa8ba81b0b
md"""
For a system of n qubits, every coeffecient is associated to a bit string of $\left\{ 0,1\right\}^n$. For a given observable, if a bit string has a $1$ at the same position as an $I$ Pauli observable then the coefficient associated to that bitstring does not contribute to the estimation of this observable.

As a consequence only coefficients that contribute to the estimation of observables given as arguments are computed when calling the method `calibrate()`. It is also possible to compute coeffcients for every $k$-local observable of a $n$-qubits system.
"""

# ╔═╡ 165d71b6-20b2-4707-8512-44feaad6cdc1
md"For $n=4$ and $k=2$ :"

# ╔═╡ aa716221-69d8-4c9c-9dca-a46b83e7eb20
calibrate(4, 2, 10000, 10; noise_model).coefficients

# ╔═╡ 4d0313c7-abf5-4e95-a369-7c32c181cf22
md"""
These coefficients are computed by generating classical shadows of a known state (the zero state). As a consequence, noise shadows (classical shadows of the zero state) can be first generated before determining which observables' expected values will be estimated :
"""

# ╔═╡ 8f9bbcc7-798a-4ad5-b051-9ae7fe904a85
begin
    #Get a classical representation of the system
    noise_shadows = generate_noise_shadows(4, 150000; noise_model)
    classical_shadows = generate_classical_shadows(ρ, 100000; noise_model)

    #Determine some observables'expected values :
    calibration = calibrate(obs, noise_shadows, 15000, 10)
    estimate_from_shadows(classical_shadows, calibration, obs, 0.15)
end

# ╔═╡ 30ad58bb-6201-4c49-8192-5fd4495b4781
md"Or directly :"

# ╔═╡ 37b82968-ae03-4645-bbf1-3f986f7298d0
estimate_from_shadows(classical_shadows, noise_shadows, obs, 0.15)

# ╔═╡ 3d00c5b5-451c-4b83-9de6-9aaac60f348a
md"""
True values being :
"""

# ╔═╡ 4e3d0ac7-ead1-418a-bfe2-3dfabf7a067e
Dict(o => expect(pauli_string(o), ρ) for o in obs)

# ╔═╡ 036072ca-cea0-493f-86f6-be4ebe51dbfd
md"""
Let's compare the standard and robust protocol for different noise rates :
"""

# ╔═╡ b64df966-e649-4dd4-8247-ecf74a6b4525
rates = [0.0, 0.05, 0.1, 0.2]

# ╔═╡ 4ede67ca-9041-4bd6-852b-05e7545d04f4
md"We will use $100000$ classical shadows and noise shadows"

# ╔═╡ f411f819-68ed-4802-bcb1-15841ac0bd85
R = 100000

# ╔═╡ c75d5423-5ad6-445e-b35b-b73f4ab1dd0a
begin
    function fidelity_compare(
        rate::Float64, state::AbstractRegister, R::Int
    )::Tuple{Float64,Float64}
        noise_model = bit_flip_noise(rate)
        n = nqubits(state)

        shadows = generate_classical_shadows(state, R; noise_model)
        noise_shadows = generate_noise_shadows(n, R; noise_model)

        standard = state_reconstruction(shadows, R)
        robust = state_reconstruction(shadows, noise_shadows, R)

        return (fidelity(s, standard), fidelity(s, robust))
    end
end

# ╔═╡ 355d74dc-3bac-4faf-89cd-5ad469e1612f
begin
    standard = Array{Float64}(undef, length(rates))
    robust = Array{Float64}(undef, length(rates))
    for i in eachindex(rates)
        r = rates[i]
        standard[i], robust[i] = fidelity_compare(r, s, R)
    end
end

# ╔═╡ 3bd41953-b7ba-4942-a7fb-ab3d3d270917
begin
    Pkg.add("Plots")
    using Plots
    plot(
        rates,
        [standard, robust],
        lab=["Standard" "Robust"],
        markershapes=:circle,
        legend_position=(0.85, 0.8),
    )
    ylabel!("Fidelity")
    xlabel!("bit flip rate")
end

# ╔═╡ 9d9da279-1f6e-4236-898f-cdb234d3dd1e
md"Results can now be displayed :"

# ╔═╡ e0f58157-86e7-4666-b6d2-46bb7cf6c17e

# ╔═╡ Cell order:
# ╠═725a9ee8-b419-4de7-a27a-f3ffb12c5033
# ╟─2b9cb58b-5162-4aba-919d-9fb6c1d38146
# ╠═c256d1a8-a774-4826-8803-7c1982dca70c
# ╠═01df233e-ad25-4180-9588-222461de9c94
# ╠═073eca16-7619-4715-b6bf-fa109c47c41b
# ╟─849d1a8d-9629-4705-a2ca-b666c4b6e76c
# ╟─8e48294f-2344-4c91-9fde-d45102e166c7
# ╠═2324c375-fca3-40c4-9c8a-45ae6e3fe83c
# ╟─8ec5215a-b6b3-406b-9205-71424e85c78d
# ╠═1f63a4f6-4e0a-43e1-907c-4857135c6e50
# ╠═39a08f0d-9b95-4e92-9fb8-e40e75a1afa4
# ╠═9cab82a6-91f8-4526-854d-acaf3e244a3b
# ╠═0a86771b-a71b-4791-bb3b-9b601624f2f5
# ╠═e12daa31-7e9c-4a9b-8fea-af4f409f37a9
# ╟─69b8ae08-611a-4bd7-9703-43129e8bfa67
# ╠═c2f945e0-185a-4cf1-8d93-8eab6e8fd3ca
# ╟─54153dc8-9f42-4961-99c8-df0b1bd6b039
# ╠═37175781-d191-40cb-bfda-b37158ef4a62
# ╟─c930571c-947e-4506-aaa3-3d258f6d6af6
# ╟─3edb6f3f-fccf-4f7a-848e-2258080d5939
# ╟─1f2b388d-11a8-4605-a073-f9c80d40559a
# ╠═f02dda39-ccf0-4f88-bc58-e2994343948e
# ╟─7c376299-1685-4ad5-9dc0-f96a280ba54a
# ╠═23dcbe82-0b2c-490e-bdfd-e30a9315169d
# ╟─7a67930f-a149-4282-a460-9620ddb63f44
# ╠═e7e4a1a8-eebd-4304-88f8-1ce10b80a7f6
# ╠═28a439fd-ac55-4c23-bacb-68571fce1e4a
# ╠═b41d0b52-72a4-4f13-8585-c09b6b66b1f3
# ╟─5af67659-0651-4f6c-8bde-83e63eead6da
# ╟─b97abec7-79a3-4d9f-8a47-e45509c31cd1
# ╠═b6d249d5-5a18-4754-b91e-278330013244
# ╟─e319a652-6ad6-4cae-9193-784e744338df
# ╟─b933ee85-9f12-42b3-9224-64f282320982
# ╟─3d80105b-122e-440c-8ad0-8c01d9aa8e91
# ╠═90079df8-f0e1-4de0-8622-270bceca778f
# ╠═875f2714-0d84-47a5-8dd4-82841b00841d
# ╠═66ff89a5-d3c0-4fe8-b726-ef73a319e070
# ╟─3ead0cc8-dcce-4369-be88-da54e7635867
# ╠═a6caa431-6fa3-42d2-b719-23782dd3172e
# ╟─92903d46-06b0-4008-8fb4-7d1186e96b7a
# ╠═bd46282f-9280-4af1-9a2f-a76e5bbfe3e4
# ╠═d069a1cd-7af7-4499-9122-aeae3fd3832c
# ╟─154a4e22-9126-4034-bad3-742c1841e5b2
# ╟─78fad7f8-7c98-487f-8205-75b3fc853e87
# ╠═6886b7f3-5086-43d7-ba87-fe975cf32de0
# ╟─a156a87e-98fe-459a-a804-7177dd41cbdc
# ╠═20b0ea20-dfd9-4dc2-8352-85a699c5d5d7
# ╟─4d207e08-7c31-4f1b-874b-009d84881f08
# ╟─3bd66405-bf3f-499c-82dd-9b1e7c5ba37b
# ╠═64e0f826-bb1c-4dd7-b331-e776ff1c9d75
# ╠═57aea7f4-39f8-44b0-bd36-9e28dafb30b7
# ╟─a48f84b1-48d7-429a-a652-ace28b0569b7
# ╠═ad7bd5e0-0a9a-4158-9b8c-6a21126909cb
# ╠═80ff8bf0-1917-4971-95e5-a86add6f7be1
# ╠═a60a84e2-c7fb-4ae5-8724-38062fab25cb
# ╠═5d63fd01-9674-4fcf-8b14-bdbbe5392040
# ╠═cbff155f-5a6b-49af-b7f5-da2604863397
# ╠═d667206c-49a2-4a28-a996-3261badb75bf
# ╟─2e3cfbff-bb55-4cce-b7c2-75f5638b2d13
# ╠═d2b2e48e-6d6a-43cc-bdab-eb47fa376a35
# ╠═8604baa1-5be4-4519-92f9-2a8608359f55
# ╠═9ab1889c-6d27-49a1-9354-7ecf166bdd4e
# ╠═6320d315-35e5-4597-9b2f-55f00c7b10fe
# ╠═71148961-f04c-4716-a510-33a3dd6a3c51
# ╠═610844f6-3593-496f-9249-45e665bf8189
# ╠═216ab5ff-fa3a-457f-bd19-0aebbc4c9ec5
# ╟─7a582eea-291f-4a4c-9403-90a677779fe8
# ╠═f40de983-54a0-4112-abc5-8b6fec2774b5
# ╟─c751189e-05e8-4e6e-aa7d-41d1ba104dfe
# ╠═48b2092b-d965-4001-b0b7-d111e2ca4fc1
# ╠═6e98e5bf-2bdc-4076-9c0d-86a4dfe8f348
# ╠═d2a9c1cb-f4fa-4245-be3b-1e1782e2e703
# ╟─1ddd531c-8754-45b8-b8ae-bd63a10fda88
# ╟─c4706f67-b7a1-4dfd-b24d-b3eb25140be9
# ╠═59e89559-1df4-4d7a-bfdd-e445b634a78d
# ╠═4fb21ef8-fe2a-4d18-a7c6-77d4ae7d55ab
# ╠═1647ae66-0667-4427-9118-fc1838f23276
# ╟─22e8fb1a-3b90-4c57-950f-f4aa8ba81b0b
# ╟─165d71b6-20b2-4707-8512-44feaad6cdc1
# ╠═aa716221-69d8-4c9c-9dca-a46b83e7eb20
# ╟─4d0313c7-abf5-4e95-a369-7c32c181cf22
# ╠═8f9bbcc7-798a-4ad5-b051-9ae7fe904a85
# ╟─30ad58bb-6201-4c49-8192-5fd4495b4781
# ╠═37b82968-ae03-4645-bbf1-3f986f7298d0
# ╟─3d00c5b5-451c-4b83-9de6-9aaac60f348a
# ╠═4e3d0ac7-ead1-418a-bfe2-3dfabf7a067e
# ╟─036072ca-cea0-493f-86f6-be4ebe51dbfd
# ╠═b64df966-e649-4dd4-8247-ecf74a6b4525
# ╟─4ede67ca-9041-4bd6-852b-05e7545d04f4
# ╠═f411f819-68ed-4802-bcb1-15841ac0bd85
# ╠═c75d5423-5ad6-445e-b35b-b73f4ab1dd0a
# ╠═355d74dc-3bac-4faf-89cd-5ad469e1612f
# ╟─9d9da279-1f6e-4236-898f-cdb234d3dd1e
# ╠═3bd41953-b7ba-4942-a7fb-ab3d3d270917
# ╠═e0f58157-86e7-4666-b6d2-46bb7cf6c17e

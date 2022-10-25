### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 320355b8-dc20-4c5b-a99d-69827b554923
# ╠═╡ show_logs = false
begin
    using Pkg
    Pkg.activate(".")

    using QuantumMeasurements
    using Yao
    using YaoPlots
    using LinearAlgebra
    using Plots
end

# ╔═╡ 9cb9f01c-2920-11ed-331b-337a297c8f86
md"""
## Variational Quantum Eigensolver

See arXiv:1803.11173v1 for more details.

One appplication of quantum computers is to implement a variational quantum eigensolver. The objective is to find the ground state of an hamiltonian. To achieve this, the idea is to have a parametrization of quantum states. Then use an optimization algorithm to minimize the energy. At each step the gradient can be computed with  classical shadows as the energy is the expected value of the hamiltonian.

"""

# ╔═╡ 5ce95f76-dd5a-4660-b0d8-efdbdb2bad7c
md"Let's define a simple hamiltonian : $\hat{\mathcal{H}} = Z_1Z_2 + X_1$ on a $3$ qubits system"

# ╔═╡ fb361f36-b364-4e2e-a683-b71ac2a84aa8
H_set = Set(["ZZI", "XII"])

# ╔═╡ 1b7e5d77-c459-44cd-a595-96fd928772ff
md"""
We now create a parametrized space of state using $p = 6$ layers of rotation gates and entangling gates.

We first generate $p*n$ axis of rotations
"""

# ╔═╡ 9fc0c652-ee37-4eb0-bc7f-92e187308776
begin
    n = 3
    p = 6
    axis = [QuantumMeasurements.pauli_from_char[rand(('X', 'Y', 'Z'))] for _ = 1:n, _ = 1:p]
end

# ╔═╡ 2a1aff7f-56db-444c-932a-26b91522c9ba
md"""
Then we need a function which maps a matrix of $n*p$ angles to a state. To do so for each angle we apply a rotation around the direction defined by the matrix `axis` to the qubit which corresponds to the line of the angle in the matix. Between these successive rotations we apply entangling gates.
"""

# ╔═╡ 774be95b-a95d-48bd-93c2-b097d7ed7d9a
begin
    function entangling_gate()::Yao.AbstractBlock
        chain(n, control(n, 1 => Z), chain(control(i, i + 1 => Z) for i = 1:(n - 1)))
    end

    function rotation_gate(angle::Matrix{Float64}, l::Int)::Yao.AbstractBlock
        chain(n, put(i => rot(axis[i, l], angle[i, l])) for i = 1:n)
    end

    function unitary_gate(angle::Matrix{Float64})::Yao.AbstractBlock
        chain(chain(entangling_gate(), rotation_gate(angle, l)) for l = 1:p)
    end

    function state_angle(angle::Matrix{Float64})::Yao.AbstractRegister
        s = zero_state(n)
        apply!(s, unitary_gate(angle))
    end

    YaoPlots.plot(unitary_gate(π * 2 * rand(Float64, (n, p))))
end

# ╔═╡ 11e9a310-0ee9-4166-bb49-d8d5f402a9f1
md"""
We now use classical shadows to evaluate the energy of one state defined by the matrix angle. We will use $4*1200$ classical shadows :
"""

# ╔═╡ 5c2e4a5f-e019-4a34-956f-9d920a0fdd80
begin
    N = 1200
    K = 4
    function energy(angle::Matrix{Float64})::Float64
        state = state_angle(angle)
        shadows = generate_classical_shadows(state, N * K)
        estimations = estimate_from_shadows(shadows, H_set, K * N, K)
        return estimations["ZZI"] + estimations["XII"]
    end
end

# ╔═╡ ac2a60ae-2f5c-4169-b370-e50cd4a38922
md"""
To minimize the energy the gradient of the energy is required. A litteral expression can be found. See arXiv:1811.11184v1 for more details.
"""

# ╔═╡ 4f6a846b-2766-4a56-b129-e4d1cac5a652
begin
    function gradient(angle::Matrix{Float64})::Matrix{Float64}
        [gradient_element(angle, i, l) for i = 1:n, l = 1:p]
    end

    function gradient_element(angle::Matrix{Float64}, i::Int, l::Int)::Float64
        angle[i, l] += π / 2
        Ep = energy(angle)
        angle[i, l] -= π
        Em = energy(angle)
        angle[i, l] += π / 2
        (Ep - Em) / 2
    end
end

# ╔═╡ 3c37e8fa-2740-4d7c-af65-a01945c1f70d
md"""
We now can implement the optimization algorithm with a learning rate η.
"""

# ╔═╡ 14a7d8ef-8258-4a2a-ae93-27f7e349cd3f
begin
    function iteration!(angle::Matrix{Float64}, η::Float64)
        grad = gradient(angle)
        for i = 1:n, l = 1:p
            angle[i, l] -= η * grad[i, l]
        end
    end
end

# ╔═╡ 394c2791-e4bc-422a-871f-124c5cf78504
md"""
To evaluate the algorithm we need to compute the theorical groundstate energy.
"""

# ╔═╡ 611917b4-2d06-4a09-926c-0e65c1d4f4b8
begin
    H_matrix = zeros(Complex, (2^n, 2^n))
    for H_term in H_set
        global H_matrix += Matrix(
            mat(kron([QuantumMeasurements.pauli_from_char[H_term[j]] for j = 1:n]...))
        )
    end

    Eth = eigvals(H_matrix)[1]
end

# ╔═╡ 97d5602f-fc4e-4408-85d3-7a3b3d147ae4
md"""
Let's plot a trajectory of the VQE starting from a random matrix with a learning rate η=0.2 and 20 iterations.
"""

# ╔═╡ a32c3bda-b88b-403d-9d6d-ada7f33f7f45
begin
    max_iter = 20
    E_trajectory = Array{Float64}(undef, max_iter)
    η = 0.2

    angle = π * 2 * rand(Float64, (n, p))
    for i = 1:max_iter
        iteration!(angle, η)
        E_trajectory[i] = energy(angle) - Eth
    end
end

# ╔═╡ 215182bf-3a50-4714-887f-cdb8d86f7d7a
begin
    Plots.plot(1:max_iter, E_trajectory, lab="E-Eth")
    xlabel!("iteration")
    ylabel!("E-Eth")
end

# ╔═╡ Cell order:
# ╠═320355b8-dc20-4c5b-a99d-69827b554923
# ╟─9cb9f01c-2920-11ed-331b-337a297c8f86
# ╟─5ce95f76-dd5a-4660-b0d8-efdbdb2bad7c
# ╠═fb361f36-b364-4e2e-a683-b71ac2a84aa8
# ╟─1b7e5d77-c459-44cd-a595-96fd928772ff
# ╠═9fc0c652-ee37-4eb0-bc7f-92e187308776
# ╟─2a1aff7f-56db-444c-932a-26b91522c9ba
# ╠═774be95b-a95d-48bd-93c2-b097d7ed7d9a
# ╠═11e9a310-0ee9-4166-bb49-d8d5f402a9f1
# ╠═5c2e4a5f-e019-4a34-956f-9d920a0fdd80
# ╟─ac2a60ae-2f5c-4169-b370-e50cd4a38922
# ╠═4f6a846b-2766-4a56-b129-e4d1cac5a652
# ╟─3c37e8fa-2740-4d7c-af65-a01945c1f70d
# ╠═14a7d8ef-8258-4a2a-ae93-27f7e349cd3f
# ╟─394c2791-e4bc-422a-871f-124c5cf78504
# ╠═611917b4-2d06-4a09-926c-0e65c1d4f4b8
# ╟─97d5602f-fc4e-4408-85d3-7a3b3d147ae4
# ╠═a32c3bda-b88b-403d-9d6d-ada7f33f7f45
# ╠═215182bf-3a50-4714-887f-cdb8d86f7d7a

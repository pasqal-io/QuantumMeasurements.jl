import Base: promote_rule

"""
PauliGate is a the union type of Pauli matrices supplemented with the identity.
"""
const PauliGate = Union{I2Gate,XGate,YGate,ZGate}

# Define the promote rules for subtypes of PauliGate. Vectors like [X,Y,Z] will
# be of type Vector{PauliGate} instead of Vecctor{ConstantGate{1,2}}. Same thing
# for other containers, like Vararg. It makes possible to define functions of
# signature f(vs::Vararg{PauliGate}).
promote_rule(::Type{S}, ::Type{T}) where {S<:PauliGate,T<:PauliGate} = PauliGate

"""
Construct a Pauli string block.
Enforces that it consists of Pauli operators acting on different qubits.
"""
pauli_string(n::Int, vs::Vararg{Pair{Int,<:PauliGate}}) = Yao.kron(n, vs...)
pauli_string(vs::Vararg{<:PauliGate}) = Yao.kron(vs...)

"""
    PauliObs

Represent a pauli observable storing only non identity terms.
`n_qubits` is the number of qubits on which the observable acts.
`pauli` is a dictionary which keys are the position of the non identity terms
and values are the character refering to the pauli matrix (namely 'X', 'Y' or 'Z').
"""
struct PauliObs
    n_qubits::Int
    pauli::Dict{Int64, Char}

    function PauliObs(obs::String)::PauliObs
        n_qubits = length(obs)
        pauli = Dict{Int64, Char}()
        for i in 1:n_qubits
            if obs[i] != 'I'
                pauli[i] = obs[i]
            end
        end
        new(n_qubits, pauli)
    end

    function Pauli(n_qubits::Int, pauli::Dict{Int64, Char})
        new(n_qubits, pauli)
    end
end

function PauliObs_to_string(obs::PauliObs)::String
    string = collect("I"^obs.n_qubits)
    for p in obs.pauli
        string[p[1]] = p[2]
    end
    join(string)
end

function weight(pauli_observable::PauliObs)::Int
    length(pauli_observable.pauli)
end

struct PauliObservable
    n_qubits::Int
    terms::Vector{AbstractBlock}
    coefficients::Vector{Float64}
end

function PauliObservable(observable::Yao.KronBlock)
    n_qubits = Yao.nqubits(observable)
    PauliObservable(n_qubits, [observable], [1.0])
end

function PauliObservable(observable::Yao.Scale)
    n_qubits = Yao.nqubits(observable)
    PauliObservable(n_qubits, [observable.content], [observable.alpha])
end

function PauliObservable(observable::Yao.Add)
    n_qubits = Yao.nqubits(observable)
    terms = []
    coefficients = []
    for b in observable.list
        obs = PauliObservable(b)
        terms = vcat(terms, obs.terms)
        coefficients = vcat(coefficients, obs.coefficients)
    end
    PauliObservable(n_qubits, terms, coefficients)
end

@const_gate V = [1 1; -im im] / sqrt(2)

@const_gate S = [1 0; 0 im]

@const_gate W = [1 im; 1 -im] / sqrt(2)

@const_gate F = [1 -im; 1 im] / sqrt(2)

@const_gate K = [1 im; im 1] / sqrt(2)

@const_gate L = [1 0; 0 im]

pauli_from_char = Dict(
    'I' => I2,
    'X' => X,
    'Y' => Y,
    'Z' => Z,
    'H' => H,
    'S' => S,
    'F' => F,
    'V' => V,
    'W' => W,
    'K' => K,
    'L' => L,
)

re_pauli_string = r"^[IXYZ]+$"

function pauli_string(s::String)
    if match(re_pauli_string, s) === nothing
        throw(DomainError)
    end

    kron([pauli_from_char[c] for c in s]...)
end

basis_change = Dict('X' => H, 'Y' => chain(shift(-π / 2), H), 'Z' => I2)

function pauli_basis_change(s::String)
    kron([basis_change[c] for c in s]...)
end

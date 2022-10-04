# TODO: This should be removed at first since we only treat Yao circuits

abstract type AbstractCircuit end

abstract type QiskitCircuit <: AbstractCircuit end
abstract type CirqCircuit <: AbstractCircuit end

const CircuitFormat = Union{ChainBlock,QiskitCircuit,CirqCircuit}

# generic structur to hold a converted circuit
# the concrete implementation of each converter will perform the
# conversion from the initial circuit platform to Yao
struct Circuit{T<:CircuitFormat}
    total_n_qubits::Int
    original_circuit::T
    yao_circuit::ChainBlock
end

# If giving as input a Yao circuit no conversion is needed
function convert(n_qubits::Int, circuit::ChainBlock{D}) where {D}
    @assert n_qubits == nqubits(circuit)
    Circuit{ChainBlock}(n_qubits, circuit, circuit)
end

function convert(circuit::ChainBlock{D}) where {D}
    Circuit{ChainBlock}(nqubits(circuit), circuit, circuit)
end

function convert(total_n_qubits::Int, circuit::QiskitCircuit)
    # TODO: Implement conversion from Qiskit serialized format (CPY)
    yao_circuit = nothing
    Circuit{QiskitCircuit}(total_n_qubits, yao_circuit, circuit)
end

function convert(total_n_qubits::Int, circuit::CirqCircuit)
    # TODO: Implement conversion from Cirq serialized format (protobuf)
    yao_circuit = nothing
    Circuit{CirqCircuit}(total_n_qubits, yao_circuit, circuit)
end

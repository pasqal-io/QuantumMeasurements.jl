```@meta
CurrentModule = QuantumMeasurements
```

# QuantumMeasurements.jl

## Presentation

`QuantumMeasurements.jl` is a Julia package for designing and tuning Quantum Measurement protocols.
This package can help:

- Generating and processing randomized measurement protocols, including the *Classical Shadows*
protocol.
- Preprocessing a derandomized sequence of measurements from a list of observables.
- Tuning the measurement protocol to a specific quantum circuit.

It also includes tools for evaluating the impact of *noise* in the measurement protocol.

The project repository can be found in the GitHub repository [QuantumMeasurements](https://github.com/pasqal-io/QuantumMeasurements.jl).

## Documentation

### Introduction

Standard quantum mechanics tells us that the state of a quantum system can be described by a density matrix $\rho$. To find out the expectation value of a given observable $\hat {\mathcal O}$, we calculate $$\langle \hat {\mathcal O} \rangle = \mathrm{Tr}(\hat {\mathcal O} \rho).$$

In practice, this means operating a quantum processing device that produces $\rho$, and performs a measurement of the qubits in some basis that is suitable for $\hat {\mathcal O}$. A (large) number of measurement outcomes is then used to estimate $\langle \hat{\mathcal O} \rangle$.

The same type of argument is used to express a full quantum state classically, by expanding the density matrix in a basis of observables. This is known as a *tomography* of the state. Likewise, one can attempt to measure powers of the density matrix, $\mathrm{Tr}(\rho^m)$.

Full tomography of an $n$-qubit quantum state requires an exponential (in $n$) number of measurements, as well as a large amount of classical postprocessing. 

Moreover, many observables are expressed as a linear combination of others, $\hat{\mathcal O} = \sum_{i=1}^L \hat{\mathcal O_i}$, where $L$ can be a very large number. Thus, the feasible number of measurements becomes a scarce resource in practical quantum computing.
### Randomized Measurements

A random measurement consists of random unitary rotation followed by a fixed measurement on each copy produced by the Quantum processing device. Appropriately averaging over these measurements, produces an efficient estimator of the desired expectation value. Thus, this protocol creates a robust classical representation of the quantum state, in the sense that the measurement information can be used for different purposes and is moreover available to postprocess in order to mitigate noise and errors.

`QuantumMeasurements.jl` provides methods for evaluating  $\varepsilon$-accurate confidence bounds, given an empirical average of randomized measurements. The protocol (if accepted) will perform a prediction of $L$ *local Pauli* expectation values $\mathrm{Tr}(\hat O_1\rho), \ldots, \mathrm{Tr}(\hat O_L \rho)$.  Given a set of Pauli observables `obs_list` , and a trial set of $M$ random Pauli measurements `meas_list`, one can simply evaluate:

```julia
conf_b = confidence_bound(obs_list, meas_list, epsilon)  # Check for example that it is less than some desired delta/2
empirical_avg = empirical_average(obs_list, meas_list) # Will be below epsilon with probability 1-delta
```

### Derandomization

By optimizing the resulting expected confidence bound for a given list of measurements (originally random), one can obtain a derandomized measurement protocol. `QuantumMeasurements.jl` performs this by using the function `derandomize`:

```julia    
derandomized_meas_list = derandomize(obs_list, meas_list, epsilon, delta)
```

This will produce an optimized *derandomized* list of measurements which we can now pass to our Quantum processing device (or even simulate). 





```@index
```

```@autodocs
Modules = [QuantumMeasurements]
```

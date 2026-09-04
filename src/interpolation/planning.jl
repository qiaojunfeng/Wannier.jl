export interpolate

abstract type Observable end

function result_name end
function requirements end
function _allocate_observable_result end
function _assemble_observable! end

struct _ObservableRequirements{P <: Tuple, I <: Tuple}
    primitive_operators::P
    intermediates::I
end

struct _InterpolationPlan{M, O <: Tuple}
    model::M
    observables::O
    intermediates::Set{Symbol}
    batch_size::Int
end

struct _InterpolationWorkspace{P, DP, H, DH, EV, U, B, OP, E}
    phase::P
    derivative_phase::DP
    hamiltonian::H
    hamiltonian_gradient::DH
    eigenvalues::EV
    eigenvectors::U
    derivative_values::B
    operator_product::OP
    diagonalization::E
end

function _normalize_observables(observables::Tuple)
    isempty(observables) && throw(ArgumentError("at least one observable is required"))
    all(observable -> observable isa Observable, observables) ||
        throw(ArgumentError("every interpolation request must be an observable recipe"))
    names = map(result_name, observables)
    length(unique(names)) == length(names) ||
        throw(ArgumentError("an observable can be requested only once"))
    return observables
end

function _plan_interpolation(
        model::InterpolationModel,
        observables::Tuple,
        number_kpoints::Integer,
    )
    recipes = _normalize_observables(observables)

    required = Set{Symbol}()
    for observable in recipes
        observable_requirements = requirements(observable)
        for operator_name in observable_requirements.primitive_operators
            haskey(model.operators, operator_name) || throw(
                ArgumentError(
                    "$(nameof(typeof(observable))) requires the missing primitive " *
                        "operator :$operator_name",
                ),
            )
        end
        union!(required, observable_requirements.intermediates)
    end
    supported = Set((:eigensystem, :hamiltonian_gradient))
    unsupported = setdiff(required, supported)
    isempty(unsupported) ||
        error("unsupported interpolation requirements: $(sort!(collect(unsupported)))")
    derivative_order = :hamiltonian_gradient in required ? 1 : 0
    batch_size = _interpolation_batch_size(
        model, number_kpoints; hamiltonian_derivative_order = derivative_order
    )
    return _InterpolationPlan(model, recipes, required, batch_size)
end

function _allocate_interpolation_workspace(plan::_InterpolationPlan)
    model = plan.model
    number_wannier = n_wannier(model)
    number_vectors = n_Rvectors(model)
    batch_size = plan.batch_size
    T = eltype(model.operators.hamiltonian.coefficients)
    RT = typeof(real(zero(T)))

    phase = Matrix{T}(undef, number_vectors, batch_size)
    hamiltonian = Array{T}(undef, number_wannier, number_wannier, batch_size)
    eigenvalues = Matrix{RT}(undef, number_wannier, batch_size)
    eigenvectors = similar(hamiltonian)
    diagonalization = _HermitianDiagonalizationWorkspace(T, number_wannier)
    if :hamiltonian_gradient in plan.intermediates
        derivative_phase = similar(phase)
        hamiltonian_gradient = Array{T}(
            undef, number_wannier, number_wannier, 3, batch_size
        )
        derivative_values = Matrix{T}(
            undef, number_wannier^2, batch_size
        )
        operator_product = Matrix{T}(undef, number_wannier, number_wannier)
    else
        derivative_phase = nothing
        hamiltonian_gradient = nothing
        derivative_values = nothing
        operator_product = nothing
    end
    return _InterpolationWorkspace(
        phase,
        derivative_phase,
        hamiltonian,
        hamiltonian_gradient,
        eigenvalues,
        eigenvectors,
        derivative_values,
        operator_product,
        diagonalization,
    )
end

function _allocate_interpolation_result(
        model::InterpolationModel,
        observables::Tuple,
        number_kpoints::Integer,
    )
    arrays = map(observables) do observable
        _allocate_observable_result(observable, model, number_kpoints)
    end
    return NamedTuple{map(result_name, observables)}(arrays)
end

function _last_axis_view(array::AbstractArray, indices)
    selectors = (ntuple(_ -> Colon(), ndims(array) - 1)..., indices)
    return view(array, selectors...)
end

function _result_batch_view(result::NamedTuple, indices)
    arrays = map(array -> _last_axis_view(array, indices), values(result))
    return NamedTuple{keys(result)}(arrays)
end

function _interpolation_intermediates!(
        plan::_InterpolationPlan,
        kpoint_batch::AbstractVector,
        workspace::_InterpolationWorkspace,
    )
    number_kpoints = length(kpoint_batch)
    number_kpoints <= plan.batch_size ||
        throw(DimensionMismatch("k-point batch exceeds the planned batch size"))

    phase = view(workspace.phase, :, 1:number_kpoints)
    _fourier_phase_block!(phase, plan.model.real_space, kpoint_batch)
    hamiltonian = view(workspace.hamiltonian, :, :, 1:number_kpoints)
    _evaluate_real_space_operator!(
        hamiltonian, plan.model.operators.hamiltonian, phase
    )

    hamiltonian_gradient = nothing
    if :hamiltonian_gradient in plan.intermediates
        derivative_phase = view(workspace.derivative_phase, :, 1:number_kpoints)
        derivative_values = view(workspace.derivative_values, :, 1:number_kpoints)
        hamiltonian_gradient = view(
            workspace.hamiltonian_gradient, :, :, :, 1:number_kpoints
        )
        _evaluate_hamiltonian_gradient!(
            hamiltonian_gradient,
            plan.model.operators.hamiltonian,
            plan.model.real_space,
            phase,
            derivative_phase,
            derivative_values,
        )
    end

    eigenvalues = view(workspace.eigenvalues, :, 1:number_kpoints)
    eigenvectors = view(workspace.eigenvectors, :, :, 1:number_kpoints)
    _diagonalize_hermitian_batch!(
        eigenvalues, eigenvectors, hamiltonian, workspace.diagonalization
    )
    return (;
        hamiltonian,
        hamiltonian_gradient,
        eigenvalues,
        eigenvectors,
    )
end

function interpolate!(
        result_view::NamedTuple,
        plan::_InterpolationPlan,
        kpoint_batch::AbstractVector,
        workspace::_InterpolationWorkspace,
    )
    length(kpoint_batch) > 0 || throw(ArgumentError("k-point batches cannot be empty"))
    intermediates = _interpolation_intermediates!(plan, kpoint_batch, workspace)
    for observable in plan.observables
        destination = getproperty(result_view, result_name(observable))
        _assemble_observable!(destination, observable, intermediates, workspace)
    end
    return result_view
end

"""
    interpolate(model, kpoints, observable)
    interpolate(model, kpoints, observables)

Interpolate one observable recipe or a tuple of recipes at fractional-coordinate
`kpoints`. The result is always a named tuple, with the k-point index on the
final axis of every result array. A combined request shares Fourier phases,
Hamiltonian derivatives, and eigendecompositions.
"""
function interpolate(
        model::InterpolationModel,
        kpoints::AbstractVector,
        observable::Observable,
    )
    return interpolate(model, kpoints, (observable,))
end

function interpolate(
        model::InterpolationModel,
        kpoints::AbstractVector,
        observables::Tuple,
    )
    number_kpoints = length(kpoints)
    number_kpoints > 0 || throw(ArgumentError("kpoints cannot be empty"))
    plan = _plan_interpolation(model, observables, number_kpoints)
    result = _allocate_interpolation_result(model, plan.observables, number_kpoints)
    workspace = _allocate_interpolation_workspace(plan)

    for first_index in 1:plan.batch_size:number_kpoints
        indices = first_index:min(first_index + plan.batch_size - 1, number_kpoints)
        interpolate!(
            _result_batch_view(result, indices),
            plan,
            view(kpoints, indices),
            workspace,
        )
    end
    return result
end

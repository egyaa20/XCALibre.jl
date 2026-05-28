export AbstractFluid, AbstractIncompressible, AbstractCompressible
export Fluid
export Incompressible, WeaklyCompressible, Compressible
export Phase, Fluid, Multiphase
export AbstractModel, AbstractEosModel, AbstractViscosityModel
export AbstractMultiphaseModel, VOF, Mixture

abstract type AbstractFluid end
abstract type AbstractIncompressible <: AbstractFluid end
abstract type AbstractCompressible <: AbstractFluid end
abstract type AbstractMultiphase <: AbstractFluid end

abstract type AbstractPhase <: AbstractMultiphase end
abstract type AbstractModel end
abstract type AbstractEosModel <: AbstractModel end
abstract type AbstractViscosityModel <: AbstractModel end


Base.show(io::IO, fluid::AbstractFluid) = print(io, typeof(fluid).name.wrapper)


"""
    Fluid <: AbstractFluid

Abstract fluid model type for constructing new fluid models.

### Fields
- 'args' -- Model arguments.

"""
struct Fluid{T,ARG}
    args::ARG
end

"""
    Incompressible <: AbstractIncompressible

Incompressible fluid model containing fluid field parameters for incompressible flows.

### Fields
- 'nu'   -- Fluid kinematic viscosity.
- 'rho'  -- Fluid density.

### Examples
- `Fluid{Incompressible}(nu=0.001, rho=1.0)` - Constructor with default values.
"""
@kwdef struct Incompressible{S1, S2, F1, F2} <: AbstractIncompressible
    nu::S1
    rho::S2
    nuf::F1
    rhof::F2
end
Adapt.@adapt_structure Incompressible

Fluid{Incompressible}(; nu, rho=1.0) = begin
    coeffs = (nu=nu, rho=rho)
    ARG = typeof(coeffs)
    Fluid{Incompressible,ARG}(coeffs)
end

(fluid::Fluid{Incompressible, ARG})(mesh) where ARG = begin
    coeffs = fluid.args
    (; rho, nu) = coeffs
    nu = ConstantScalar(nu)
    nuf = nu
    rho = ConstantScalar(rho)
    rhof = rho
    Incompressible(nu, rho, nuf, rhof)
end

"""
    WeaklyCompressible <: AbstractCompressible

Weakly compressible fluid model containing fluid field parameters for weakly compressible 
    flows with constant parameters - ideal gas with constant viscosity.

### Fields
- 'nu'   -- Fluid kinematic viscosity.
- 'cp'   -- Fluid specific heat capacity.
- `gamma` -- Ratio of specific heats.
- `Pr`   -- Fluid Prandtl number.

### Examples
- `Fluid{WeaklyCompressible}(; nu=1E-5, cp=1005.0, gamma=1.4, Pr=0.7)` - Constructor with 
default values.
"""
struct WeaklyCompressible{S1, S2, F1, F2, T} <: AbstractCompressible
    nu::S1
    rho::S2
    nuf::F1
    rhof::F2
    cp::T
    gamma::T
    Pr::T
    R::T
end
Adapt.@adapt_structure WeaklyCompressible

Fluid{WeaklyCompressible}(; nu, cp, gamma, Pr) = begin
    coeffs = (nu=nu, cp=cp, gamma=gamma, Pr=Pr)
    ARG = typeof(coeffs)
    Fluid{WeaklyCompressible,ARG}(coeffs)
end

(fluid::Fluid{WeaklyCompressible, ARG})(mesh) where ARG = begin
    coeffs = fluid.args
    (; nu, cp, gamma, Pr) = coeffs
    cp = ConstantScalar(cp)
    gamma = ConstantScalar(gamma)
    Pr = ConstantScalar(Pr)
    R = ConstantScalar(cp.values*(1.0 - (1.0/gamma.values)))

    nu = ConstantScalar(nu)
    rho = ScalarField(mesh)
    nuf = nu
    rhof = FaceScalarField(mesh)
    WeaklyCompressible(nu, rho, nuf, rhof, cp, gamma, Pr, R)
end

"""
    Compressible <: AbstractCompressible

Compressible fluid model containing fluid field parameters for compressible flows with 
    constant parameters - ideal gas with constant viscosity.

### Fields
- 'nu'   -- Fluid kinematic viscosity.
- 'cp'   -- Fluid specific heat capacity.
- `gamma` -- Ratio of specific heats.
- `Pr`   -- Fluid Prantl number.

### Examples
- `Fluid{Compressible}(; nu=1E-5, cp=1005.0, gamma=1.4, Pr=0.7)` - Constructur with default values.
"""
@kwdef struct Compressible{S1, S2, F1, F2, T} <: AbstractCompressible
    nu::S1
    rho::S2
    nuf::F1
    rhof::F2
    cp::T
    gamma::T
    Pr::T
    R::T
end
Adapt.@adapt_structure Compressible

Fluid{Compressible}(; nu=1E-5, cp=1005.0, gamma=1.4, Pr=0.7 ) = begin
    coeffs = (nu=nu, cp=cp, gamma=gamma, Pr=Pr)
    ARG = typeof(coeffs)
    Fluid{Compressible,ARG}(coeffs)
end

(fluid::Fluid{Compressible, ARG})(mesh) where ARG = begin
    coeffs = fluid.args
    (; nu, cp, gamma, Pr) = coeffs
    cp = ConstantScalar(cp)
    gamma = ConstantScalar(gamma)
    Pr = ConstantScalar(Pr)
    R = ConstantScalar(cp.values*(1.0 - (1.0/gamma.values)))

    nu = ConstantScalar(nu)
    rho = ScalarField(mesh)
    nuf = nu
    rhof = FaceScalarField(mesh)
    Compressible(nu, rho, nuf, rhof, cp, gamma, Pr, R)
end


"""
    Phase <: AbstractPhase

Configuration structure for a single fluid phase.

### Fields
- `rho` -- Density model (Equation of State) for the phase.
- `mu`  -- Viscosity model for the phase.
- `cp`  -- Specific heat capacity at constant pressure [J/(kg·K)]. Used by the
           multiphase temperature energy model. Defaults to 0.0 (ignored when
           the energy model is `Isothermal`).
- `k`   -- Thermal conductivity [W/(m·K)]. Same role as `cp`.
"""
struct Phase{E<:AbstractEosModel, V<:AbstractViscosityModel, CP, K} <: AbstractPhase
    rho::E
    mu::V
    cp::CP
    k::K
end

function Phase(; rho, mu, cp=0.0, k=0.0) # cp/k optional — only used by MultiphaseTemperature energy model
    rho_model = rho isa AbstractFloat ? ConstEos(rho) : rho
    mu_model = mu  isa AbstractFloat ? ConstMu(mu) : mu
    return Phase(rho_model, mu_model, cp, k)
end

@kwdef struct PhaseState{E<:AbstractEosModel, V<:AbstractViscosityModel, S1,S2,S3,S4,S5} <: AbstractPhase
    rho_model::E
    mu_model::V

    rho::S1
    mu::S2
    k::S3
    cp::S4
    beta::S5
end
Adapt.@adapt_structure PhaseState


function build_phase(phase_setup::Phase, mesh)
    rho   = phase_setup.rho isa ConstEos ? ConstantScalar(phase_setup.rho.rho) : ScalarField(mesh)
    mu    = phase_setup.mu  isa ConstMu ? ConstantScalar(phase_setup.mu.mu) : ScalarField(mesh)
    cp    = ConstantScalar(phase_setup.cp)
    k     = ConstantScalar(phase_setup.k)
    beta  = ScalarField(mesh)

    return PhaseState(
        rho_model = phase_setup.rho,
        mu_model = phase_setup.mu,
        rho=rho,
        mu=mu,
        k=k,
        cp=cp,
        beta=beta
    )
end

"""
    build_phase_table_mode(phase_setup, mesh)

Variant of `build_phase` used when `fluid_properties::HelmholtzTable` is
present on `Fluid{Multiphase}`. All four T-dependent properties (rho,
mu, cp, k) are allocated as `ScalarField`s rather than `ConstantScalar`s
so the live update routine can mutate them per cell at every outer
iteration. Initial values are populated from the supplied `phase_setup`
constants (which themselves came from the snapshot at `T_snapshot`).
"""
function build_phase_table_mode(phase_setup::Phase, mesh)
    rho   = ScalarField(mesh); initialise!(rho, phase_setup.rho.rho)
    mu    = ScalarField(mesh); initialise!(mu,  phase_setup.mu.mu)
    cp    = ScalarField(mesh); initialise!(cp,  phase_setup.cp)
    k     = ScalarField(mesh); initialise!(k,   phase_setup.k)
    beta  = ScalarField(mesh)

    return PhaseState(
        rho_model = phase_setup.rho,
        mu_model = phase_setup.mu,
        rho=rho,
        mu=mu,
        k=k,
        cp=cp,
        beta=beta
    )
end


"""
    AbstractMultiphaseModel

Tag type for multiphase interface/coupling models. Concrete subtypes (`VOF`,
`Mixture`) carry the solver-specific knobs and are used to dispatch inside the
multiphase solver (e.g. `model.fluid.model isa VOF`).
"""
abstract type AbstractMultiphaseModel end

"""
    VOF(; sigma=0.0, cAlpha=1.0, smooth=0, cycles=1) <: AbstractMultiphaseModel

Volume-of-Fluid interface-capturing settings.

### Fields
- `sigma`  -- Surface tension coefficient [N/m].
- `cAlpha` -- Interface compression coefficient (MULES).
- `smooth` -- Number of Laplacian passes to smooth α before computing the
              compression-flux normal `n̂f`. `0` (default) keeps the raw
              gradient. `> 0` reduces per-cell sawtooth artefacts on thin
              interfaces at low/zero σ.
- `cycles` -- Number of MULES sub-cycles per outer time step. `1` (default)
              disables sub-cycling. `> 1` splits the outer dt into N MULES
              α-updates with `dt_sub = dt/N`, increasing the effective α-CFL
              headroom for adaptive time-stepping by a factor of N.
"""
@kwdef struct VOF{T1,T2,T3,T4} <: AbstractMultiphaseModel
    sigma::T1        = 0.0
    cAlpha::T2       = 1.0
    smooth::T3       = 0
    cycles::T4       = 1
end
Adapt.@adapt_structure VOF

"""
    Mixture(; diameter=1.0e-6) <: AbstractMultiphaseModel

Manninen drift-flux mixture-model settings.

### Fields
- `diameter` -- Dispersed-phase particle/bubble diameter [m].
"""
@kwdef struct Mixture{T1} <: AbstractMultiphaseModel
    diameter::T1 = 1.0e-6
end
Adapt.@adapt_structure Mixture


"""
    Multiphase <: AbstractMultiphase

Multiphase fluid model containing multiple phases and their interaction properties.

### Fields
- 'phases'             -- Tuple of PhaseState structures.
- 'physics_properties' -- NamedTuple of physical models (drag, surface tension, etc.).
- 'volume_fraction'    -- Index of the phase tracked by the volume fraction field.
- 'alpha'              -- Volume fraction ScalarField.
- 'alphaf'             -- Volume fraction FaceScalarField.
- 'rho'                -- Mixture density ScalarField.
- 'rhof'               -- Mixture density FaceScalarField.
- 'nu'                 -- Mixture kinematic viscosity ScalarField.
- 'nuf'                -- Mixture kinematic viscosity FaceScalarField.
- 'p_rgh'              -- Dynamic pressure ScalarField.
- 'p_rghf'             -- Dynamic pressure FaceScalarField.
"""
@kwdef struct Multiphase{M,P1,P2,S1,F1,S2,F2,S3,F3,S4,F4} <: AbstractMultiphase
    model::M
    phases::P1
    physics_properties::P2
    volume_fraction::Int
    alpha::S1
    alphaf::F1
    rho::S2
    rhof::F2
    nu::S3
    nuf::F3
    p_rgh::S4
    p_rghf::F4
end
Adapt.@adapt_structure Multiphase


Fluid{Multiphase}(; phases::NTuple{2, Phase}, model=nothing, kwargs...) = begin
    @assert model isa AbstractMultiphaseModel "Expected `model = VOF(...)` or `model = Mixture(...)`, got: $(typeof(model))"
    coeffs = (; phases, model, kwargs...)
    ARG = typeof(coeffs)
    Fluid{Multiphase, ARG}(coeffs)
end

(fluid::Fluid{Multiphase, ARG})(mesh) where {ARG} = begin
    coeffs = fluid.args
    physics_properties = Base.structdiff(coeffs, (phases = nothing, model = nothing))

    phase_setups = coeffs.phases
    @assert phase_setups isa Tuple{Phase, Phase} "Phases must be a plain Tuple of exactly two Phase objects, e.g. (Phase(...), Phase(...))"

    volume_fraction = 1  # First phase is always the tracked phase

    build_multiphase(coeffs.model, phase_setups, physics_properties, mesh, volume_fraction)
end

build_property(property, mesh) = property
build_property(setup::Gravity, mesh) = build_gravityModel(setup, mesh)

function build_multiphase(model::AbstractMultiphaseModel, phase_setups::Tuple{<:AbstractPhase, <:AbstractPhase}, physics_properties_setup::NamedTuple, mesh, volume_fraction::Int)
    # Phase-1 / Phase-2 high-fidelity properties: if the user supplied a
    # `fluid_properties = HelmholtzTable(...)` keyword:
    #   - the phase setups are first rewritten to use the table snapshot
    #     at `T_snapshot` (so initial values are physical);
    #   - then phases are allocated with ScalarField backings (Phase-2
    #     mode) so that the live per-cell update routine has somewhere to
    #     write into each iteration.
    phase_setups = maybe_snapshot_from_fluid_properties(phase_setups, physics_properties_setup)

    if haskey(physics_properties_setup, :fluid_properties) && physics_properties_setup.fluid_properties isa HelmholtzTable
        phases = map(setup -> build_phase_table_mode(setup, mesh), phase_setups)
    else
        phases = map(setup -> build_phase(setup, mesh), phase_setups)
    end

    built_properties = map(prop_setup -> build_property(prop_setup, mesh), physics_properties_setup)

    alpha  = ScalarField(mesh)
    alphaf = FaceScalarField(mesh)

    rho  = ScalarField(mesh)
    rhof = FaceScalarField(mesh)

    nu  = ScalarField(mesh)
    nuf = FaceScalarField(mesh)

    p_rgh  = ScalarField(mesh)
    p_rghf = FaceScalarField(mesh)

    Multiphase(model=model, phases=phases, physics_properties=built_properties, volume_fraction=volume_fraction, alpha=alpha, alphaf=alphaf, rho=rho, rhof=rhof, nu=nu, nuf=nuf, p_rgh=p_rgh, p_rghf=p_rghf)
end

# Default: no fluid_properties kwarg, return phase setups unchanged.
maybe_snapshot_from_fluid_properties(phase_setups::Tuple, physics_properties::NamedTuple) =
    haskey(physics_properties, :fluid_properties) ?
        _snapshot_from_fluid_properties(phase_setups, physics_properties.fluid_properties) :
        phase_setups

# Generic catch-all so a user-supplied `fluid_properties = SomethingElse()`
# is left alone. Specialised on `HelmholtzTable` in the FluidProperties
# module (HelmholtzTable.jl is loaded later in the module include order).
_snapshot_from_fluid_properties(phase_setups, _other) = phase_setups
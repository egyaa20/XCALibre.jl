export RPIWallState, init_rpi_wall_state, update_rpi_wall_state!

"""
    RPIWallState

Per-iteration RPI wall-state container. Holds the per-face partition values
and the per-cell *correction* sources that the multiphase solver folds back
into the T and α equations.

How it integrates (no new BC type needed):
- The user's `NeumannFunction(:Wall_*_Heated, heated_flux)` keeps applying
  the full Joule flux `q_w` to the wall via the standard BC.
- This state's `cell_source_T` deposits `-q_e · A_face / V_cell` at the
  *wall-adjacent cell* for each heated face. Net heat into the cell becomes
  `(q_w − q_e) · A_face = (q_c + q_q) · A_face` — RPI's partition exactly.
- Similarly, `cell_source_alpha` deposits the vapour-generation mass
  source at wall-adjacent cells.

Fields:
- `rpi`              : RPI configuration
- `sat`              : SatProps (T_sat, h_fg, ρ_l, ρ_v, …)
- `q_w_func`         : closure (coords, time, fID) → q_w  (the Joule flux)
- `h_c`              : convective HTC estimate fed to rpi_partition
- `q_to_fluid_face`  : FaceScalarField  (q_c + q_q) — diagnostic only
- `mdot_e_face`      : FaceScalarField  evaporation rate — diagnostic only
- `cell_source_T`    : ScalarField  −q_e · A / V at wall cells [W/m³]  (folded into T_source)
- `cell_source_alpha`: ScalarField  −ṁ_e · A / (ρ_l V) at wall cells [1/s] (α correction)
- `heated_patches`   : Vector{Symbol} — patches where partition fires
"""
struct RPIWallState{F<:AbstractFloat,QF,FF,SF}
    rpi::RPI{F}
    sat::SatProps{F}
    q_w_func::QF
    h_c::F
    q_to_fluid_face::FF
    mdot_e_face::FF
    cell_source_T::SF
    cell_source_alpha::SF
    heated_patches::Vector{Symbol}
end


"""
    init_rpi_wall_state(rpi, sat, q_w_func, h_c, mesh, heated_patches)

Allocate the field containers backing an `RPIWallState`. Pass to the
multiphase fluid via `wall_boiling = state`; the solver picks it up
automatically each outer corrector.
"""
function init_rpi_wall_state(rpi::RPI, sat::SatProps, q_w_func,
                              h_c::Real, mesh, heated_patches::Vector{Symbol})
    F = _get_float(mesh)
    q_to_fluid_face   = FaceScalarField(mesh)
    mdot_e_face       = FaceScalarField(mesh)
    cell_source_T     = ScalarField(mesh)
    cell_source_alpha = ScalarField(mesh)
    return RPIWallState{F,typeof(q_w_func),typeof(q_to_fluid_face),typeof(cell_source_T)}(
        rpi, sat, q_w_func, F(h_c),
        q_to_fluid_face, mdot_e_face,
        cell_source_T, cell_source_alpha, heated_patches,
    )
end


"""
    update_rpi_wall_state!(state, model, time, config)

Once per outer corrector: for every heated wall face, read T at the
adjacent boundary cell (lagged proxy for T_w), Newton-solve the
partition, then deposit:
  - q_e · A_face into `cell_source_T[cID]`     as a NEGATIVE source  (so that
    net q into cell = q_w − q_e = q_c + q_q after BC + this correction)
  - ṁ_e · A_face / (ρ_l V) into `cell_source_alpha[cID]` as a NEGATIVE source
    (liquid disappears at wall)

A small CPU-side loop over the heated boundary faces — fine for typical
mesh sizes. If profiling shows this hot, kernelise per patch later.
"""
function update_rpi_wall_state!(state::RPIWallState, model, time, config)
    mesh        = model.domain
    rpi         = state.rpi
    sat         = state.sat
    q_w_func    = state.q_w_func
    h_c         = state.h_c

    # Everything we touch from `mesh` (cells, faces, boundaries, boundary_cellsID)
    # lives on the device when `mesh = model.domain` is the adapted mesh.
    # The function runs as a CPU loop, so pull each one to host once.
    T_host          = Array(model.energy.T.values)
    cells_host      = Array(mesh.cells)
    Vol_host        = [c.volume for c in cells_host]
    faces_host      = Array(mesh.faces)
    bnd_cells       = Array(mesh.boundary_cellsID)
    boundaries_host = Array(mesh.boundaries)

    fluid_nt = (
        ρ_l  = sat.rho_l,
        ρ_v  = sat.rho_v,
        μ_l  = sat.mu_l,
        cp_l = sat.cp_l,
        k_l  = sat.k_l,
        h_fg = sat.h_fg,
        σ    = sat.sigma,
    )

    qfacebuf    = Array(state.q_to_fluid_face.values);   fill!(qfacebuf, 0)
    mdotbuf     = Array(state.mdot_e_face.values);       fill!(mdotbuf, 0)
    srcT_buf    = Array(state.cell_source_T.values);     fill!(srcT_buf, 0)
    srcA_buf    = Array(state.cell_source_alpha.values); fill!(srcA_buf, 0)

    for patch in state.heated_patches
        ID = boundary_index(boundaries_host, patch)
        boundary = boundaries_host[ID]
        for fID in boundary.IDs_range
            face   = faces_host[fID]
            cID    = bnd_cells[fID]
            T_l    = T_host[cID]
            q_w    = q_w_func(face.centre, time, fID)

            try
                sol = rpi_solve_wall_temperature(
                    q_w, T_l, sat.T_sat, h_c, fluid_nt, rpi;
                    T_w_guess = max(T_l + 1.0, sat.T_sat + 1.0),
                    tol = 1.0e-3, max_iter = 30)
                part = sol.partition

                qfacebuf[fID] = part.q_c + part.q_q
                mdotbuf[fID]  = part.mdot_e

                # Negative S_T correction: subtract the evaporation portion
                # from the wall-adjacent cell's heat budget (BC already
                # applied full q_w; we remove q_e).
                #   units: q_e [W/m²] · A_face [m²] / V_cell [m³] = [W/m³]
                srcT_buf[cID] += -(part.q_e * face.area) / Vol_host[cID]

                # Negative S_α correction: liquid disappearing at wall.
                #   units: ṁ_e [kg/(m²·s)] · A_face / (ρ_l · V_cell) = [1/s]
                srcA_buf[cID] += -(part.mdot_e * face.area) / (sat.rho_l * Vol_host[cID])
            catch
                # Newton diverged (typically q_w too small for boiling).
                # Leave face at zero correction → no partition; BC's q_w
                # all goes to fluid via the regular convection path.
                qfacebuf[fID] = q_w
                mdotbuf[fID]  = 0.0
            end
        end
    end

    copyto!(state.q_to_fluid_face.values,   qfacebuf)
    copyto!(state.mdot_e_face.values,       mdotbuf)
    copyto!(state.cell_source_T.values,     srcT_buf)
    copyto!(state.cell_source_alpha.values, srcA_buf)

    return nothing
end

update_rpi_wall_state!(::Nothing, _model, _time, _config) = nothing
build_property(s::RPIWallState, mesh) = s

# Boiling submodels.
#
# Currently includes the RPI wall heat-flux partitioning model
# (Kurul & Podowski, 1990) with these closures:
#   - Lemmert-Chawla     : nucleation-site density
#   - Tolubinsky-Kostanchuk: bubble departure diameter
#   - Cole               : bubble departure frequency
#
# The contents are pure functions plus a configuration struct. Solver
# integration (alpha-source from `mdot_e`, latent-heat sink for T,
# iterative wall-T solve when `q_w` is prescribed) is not done here —
# this module produces the per-face quantities that those couplings
# would consume.

include("RPI.jl")
include("sat_props.jl")
include("Lee.jl")
include("RPI_wall_update.jl")

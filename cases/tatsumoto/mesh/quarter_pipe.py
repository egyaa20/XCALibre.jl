#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
SALOME 9.14 — parameterised 90 degree quarter-pipe structured hex mesh.

Driven by a JSON parameter file rather than edited constants, so a case TOML
can define the geometry and the cluster can build it unattended:

    SALOME_MESH_PARAMS=params.json SALOME_MESH_OUT=mesh.unv \
        salome -t quarter_pipe.py

Parameters are read from the JSON file (see PARAM_DEFAULTS below). The mesh is
exported as UNV, which XCALibre reads via UNV3D_mesh(file, scale=0.001) — this
script works in millimetres, the scale factor is applied at read time.

Geometry: three axial zones stacked along +Z (the flow direction the tatsumoto
stages assume), inlet at Z=0:

    Z = 0                       Inlet
      |  L_entrance             unheated development length
    Z = L_entrance
      |  L_heated               HEATED section
    Z = L_entrance + L_heated
      |  L_exit                 unheated outlet extension (may be 0)
    Z = total                   Outlet

Cross-section is a quarter disc meshed with an O-grid (butterfly) so there is
no singular cell on the axis and the curved wall carries a clean boundary
layer. In the X-Y plane (apex O on the pipe axis):

      Y=R  Wy
           |\
           | \  (arc)
      Y=c  Py--Pd--Wm      Wm on the arc at 45 degrees
           |   |    \
           | 1 | 3 /
      O----Px--Wx----  X
           X=c  X=R

    block 1 = core square (O,Px,Pd,Py)
    block 2 = lower outer (Px,Wx,Wm,Pd)
    block 3 = upper outer (Py,Pd,Wm,Wy)
    c = core_ratio * R

Emitted face groups — these names are a contract with the tatsumoto stage
scripts (cases/tatsumoto/stages/*.jl) and must not be renamed:

    Inlet, Outlet, Symmetry1 (Y=0), Symmetry2 (X=0),
    Wall_Heated, Wall_Unheated  (Wall_Unheated spans BOTH the entrance and
                                 exit extensions)
"""

import json
import math
import os
import sys

import salome
salome.salome_init()

import GEOM
from salome.geom import geomBuilder
geompy = geomBuilder.New()

import SMESH
from salome.smesh import smeshBuilder
smesh = smeshBuilder.New()

# ============================================================================
# 1.  PARAMETERS
# ============================================================================
PARAM_DEFAULTS = {
    # -- geometry [mm] -------------------------------------------------------
    "radius":     2.7,     # pipe radius; Tatsumoto tube ID is 5.4 mm
    "L_entrance": 200.0,   # unheated development length (paper: 200 mm = 37d)
    "L_heated":   100.0,   # heated section (paper: 100 mm)
    "L_exit":     0.0,     # unheated extension past the heated section
    "core_ratio": 0.4,     # O-grid core half-size as a fraction of R

    # -- resolution ----------------------------------------------------------
    "n_quarter":  10,      # segments on core + circumferential edges
    "n_radial":   20,      # segments on radial (wall-normal) edges; carries BL
    "n_axial":    300,     # total axial segments across the FULL pipe length;
                           # split across entrance/heated/exit proportionally
                           # to zone length (uniform axial cell size). Ignored
                           # when the explicit per-zone counts below are set.
    "n_ax_entrance": None, # explicit per-zone axial counts; all three must be
    "n_ax_heated":   None, # set together (entrance/exit may be 0 when the
    "n_ax_exit":     None, # zone length is 0)
    "axial_ratio": 1.0,    # largest/smallest axial cell in the entrance and
                           # exit zones, fine end towards the heated section;
                           # heated zone is always uniform

    # -- boundary layer ------------------------------------------------------
    # Preferred: state the target first wall-normal cell height directly and
    # let the grading be solved for it. Set first_cell to null to specify the
    # raw SALOME scale factor (last/first cell length) via bl_growth instead.
    "first_cell": 0.01,
    "bl_growth":  20.0,
    "bl_flip":    False,   # True if the BL ends up clustered at the axis
}


def load_params():
    path = os.environ.get("SALOME_MESH_PARAMS")
    p = dict(PARAM_DEFAULTS)
    if path:
        with open(path) as fh:
            user = json.load(fh)
        unknown = set(user) - set(PARAM_DEFAULTS)
        if unknown:
            raise SystemExit("unknown mesh parameter(s): %s" % sorted(unknown))
        p.update(user)
    return p


def solve_scale(length, n, first_cell):
    """SALOME scale factor (last/first) giving `first_cell` on a graded edge.

    Cells follow a geometric progression with per-cell ratio q, so
    first*(q^n - 1)/(q - 1) = length and scale = q^(n-1). Inverted by
    bisection: the first cell shrinks monotonically as q grows.
    """
    if n < 2:
        return 1.0
    uniform = length / n
    if first_cell >= uniform:
        print("  WARN first_cell %.4g >= uniform spacing %.4g mm; using uniform"
              % (first_cell, uniform))
        return 1.0

    def first_of(q):
        if abs(q - 1.0) < 1e-12:
            return length / n
        return length * (q - 1.0) / (q ** n - 1.0)

    lo, hi = 1.0 + 1e-12, 5.0
    if first_of(hi) > first_cell:
        print("  WARN cannot reach first_cell=%.4g mm with n_radial=%d "
              "(min achievable %.4g); clamping" % (first_cell, n, first_of(hi)))
        return hi ** (n - 1)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if first_of(mid) > first_cell:
            lo = mid
        else:
            hi = mid
    return (0.5 * (lo + hi)) ** (n - 1)


def first_cell_of(length, n, scale):
    """Length of the first cell on an edge graded by `scale` (last/first)."""
    if n < 2 or abs(scale - 1.0) < 1e-9:
        return length / max(n, 1)
    q = scale ** (1.0 / (n - 1))
    return length * (q - 1.0) / (q ** n - 1.0)


P = load_params()

R          = float(P["radius"])
L_ENT      = float(P["L_entrance"])
L_HEAT     = float(P["L_heated"])
L_EXIT     = float(P["L_exit"])
CORE_RATIO = float(P["core_ratio"])
N_QUARTER  = int(P["n_quarter"])
N_RADIAL   = int(P["n_radial"])
N_AXIAL    = int(P["n_axial"])
BL_FLIP    = bool(P["bl_flip"])

if L_HEAT <= 0:
    raise SystemExit("L_heated must be > 0")
if not (0.0 < CORE_RATIO < 1.0 / math.sqrt(2.0)):
    raise SystemExit("core_ratio must be in (0, 0.707) so the core stays inside the arc")

Z_HEAT_START = L_ENT
Z_HEAT_END   = L_ENT + L_HEAT
Z_TOTAL      = Z_HEAT_END + L_EXIT

c = CORE_RATIO * R                      # core square half-size
s = R / math.sqrt(2.0)                  # arc point at 45 degrees

# Radial edge lengths differ between the axis-aligned and diagonal edges. One
# SALOME hypothesis covers both, so the grading is solved for the axis edges
# and the diagonal's resulting first cell is reported for information.
L_RAD_AXIS = R - c
L_RAD_DIAG = R - c * math.sqrt(2.0)

if P["first_cell"] is None:
    BL_SCALE = float(P["bl_growth"])
else:
    BL_SCALE = solve_scale(L_RAD_AXIS, N_RADIAL, float(P["first_cell"]))

fc_axis = first_cell_of(L_RAD_AXIS, N_RADIAL, BL_SCALE)
fc_diag = first_cell_of(L_RAD_DIAG, N_RADIAL, BL_SCALE)

if N_AXIAL < 1:
    raise SystemExit("n_axial must be >= 1")

# Per-zone axial counts: explicit when n_ax_* are set, otherwise N_AXIAL split
# proportionally to zone length (uniform axial cell size). Each present zone
# keeps at least 1 segment even if its share rounds to 0.
explicit_ax = [P["n_ax_entrance"], P["n_ax_heated"], P["n_ax_exit"]]
if any(v is not None for v in explicit_ax):
    if any(v is None for v in explicit_ax):
        raise SystemExit("n_ax_entrance/n_ax_heated/n_ax_exit must be set together")
    n_ax_ent, n_ax_heat, n_ax_exit = (int(v) for v in explicit_ax)
    if (L_ENT > 0) != (n_ax_ent > 0) or (L_EXIT > 0) != (n_ax_exit > 0):
        raise SystemExit("per-zone axial counts must be > 0 exactly for zones with length > 0")
    if n_ax_heat < 1:
        raise SystemExit("n_ax_heated must be >= 1")
else:
    AXIAL_CELL = Z_TOTAL / N_AXIAL
    n_ax_ent  = max(1, int(round(L_ENT  / AXIAL_CELL))) if L_ENT  > 0 else 0
    n_ax_heat = max(1, int(round(L_HEAT / AXIAL_CELL)))
    n_ax_exit = max(1, int(round(L_EXIT / AXIAL_CELL))) if L_EXIT > 0 else 0
n_ax_total = n_ax_ent + n_ax_heat + n_ax_exit

AXIAL_RATIO = float(P["axial_ratio"])
if AXIAL_RATIO < 1.0:
    raise SystemExit("axial_ratio must be >= 1 (largest/smallest)")

cross_cells = N_QUARTER ** 2 + 2 * N_QUARTER * N_RADIAL

print("=" * 68)
print("  Quarter pipe : R = %g mm, total L = %g mm" % (R, Z_TOTAL))
print("    entrance %g mm (%d cells) | heated %g mm (%d) | exit %g mm (%d)"
      % (L_ENT, n_ax_ent, L_HEAT, n_ax_heat, L_EXIT, n_ax_exit))
print("  core_ratio = %g (core half-size %.4g mm)" % (CORE_RATIO, c))
print("  n_quarter=%d, n_radial=%d, n_axial=%d (heated cell %.4g mm, axial_ratio %g)"
      % (N_QUARTER, N_RADIAL, n_ax_total, L_HEAT / n_ax_heat, AXIAL_RATIO))
print("  BL scale factor = %.4g (flip=%s)" % (BL_SCALE, BL_FLIP))
print("  first wall cell ~ %.5g mm (axis edges), %.5g mm (diagonal)"
      % (fc_axis, fc_diag))
print("  cross-section cells = %d, total ~ %d"
      % (cross_cells, cross_cells * n_ax_total))
print("=" * 68)

# ============================================================================
# 2.  GEOMETRY
# ============================================================================
O  = geompy.MakeVertex(0, 0, 0)
Px = geompy.MakeVertex(c, 0, 0)
Py = geompy.MakeVertex(0, c, 0)
Pd = geompy.MakeVertex(c, c, 0)
Wx = geompy.MakeVertex(R, 0, 0)
Wy = geompy.MakeVertex(0, R, 0)
Wm = geompy.MakeVertex(s, s, 0)

e_OWx = geompy.MakeEdge(O, Wx)
arc   = geompy.MakeArc(Wx, Wm, Wy)
e_WyO = geompy.MakeEdge(Wy, O)
sector_face = geompy.MakeFace(geompy.MakeWire([e_OWx, arc, e_WyO]), True)

# Butterfly cutters meeting at Pd; their endpoints split the sector boundary.
cut1 = geompy.MakeEdge(Px, Pd)
cut2 = geompy.MakeEdge(Py, Pd)
cut3 = geompy.MakeEdge(Pd, Wm)

block_2d = geompy.MakePartition(
    [sector_face], [cut1, cut2, cut3], [], [],
    geompy.ShapeType["FACE"], 0, [], 0)

vz    = geompy.MakeVectorDXDYDZ(0, 0, 1)
solid = geompy.MakePrismVecH(block_2d, vz, Z_TOTAL)

# Split the extrusion at each heated-section boundary so the wall can be
# grouped into heated / unheated patches.
cut_planes = []
for z in (Z_HEAT_START, Z_HEAT_END):
    if 1e-9 < z < Z_TOTAL - 1e-9:
        cut_planes.append(geompy.MakePlane(
            geompy.MakeVertex(0, 0, z),
            geompy.MakeVectorDXDYDZ(0, 0, 1), 4 * R))

domain = geompy.MakePartition(
    [solid], cut_planes, [], [], geompy.ShapeType["SOLID"], 0, [], 0) \
    if cut_planes else solid
geompy.addToStudy(domain, "QuarterPipe")

# ============================================================================
# 3.  FACE GROUPS
# Order matters: Z-planes first, then the curved wall by radius, then the two
# flat symmetry cuts by centroid.
# ============================================================================
tol    = 1.0e-3
r_wall = 0.90 * R

inlet_f, outlet_f, sym1_f, sym2_f, wall_uh_f, wall_h_f = [], [], [], [], [], []

for f in geompy.SubShapeAll(domain, geompy.ShapeType["FACE"]):
    x, y, z = geompy.PointCoordinates(geompy.MakeCDG(f))
    r = math.sqrt(x * x + y * y)
    if abs(z) < tol:
        inlet_f.append(f)
    elif abs(z - Z_TOTAL) < tol:
        outlet_f.append(f)
    elif r > r_wall:
        if Z_HEAT_START + tol < z < Z_HEAT_END - tol:
            wall_h_f.append(f)
        else:
            wall_uh_f.append(f)      # entrance AND exit extensions
    elif abs(y) < tol:
        sym1_f.append(f)
    elif abs(x) < tol:
        sym2_f.append(f)


def make_group(items, name, shape_type):
    g = geompy.CreateGroup(domain, geompy.ShapeType[shape_type])
    geompy.UnionList(g, items)
    geompy.addToStudyInFather(domain, g, name)
    return g


if not wall_h_f:
    raise SystemExit("no heated wall faces found — check L_entrance/L_heated")

g_inlet   = make_group(inlet_f,   "Inlet",         "FACE")
g_outlet  = make_group(outlet_f,  "Outlet",        "FACE")
g_sym1    = make_group(sym1_f,    "Symmetry1",     "FACE")
g_sym2    = make_group(sym2_f,    "Symmetry2",     "FACE")
g_wall_uh = make_group(wall_uh_f, "Wall_Unheated", "FACE")
g_wall_h  = make_group(wall_h_f,  "Wall_Heated",   "FACE")

print("  faces: Inlet=%d Outlet=%d Sym1=%d Sym2=%d Wall_H=%d Wall_UH=%d"
      % (len(inlet_f), len(outlet_f), len(sym1_f), len(sym2_f),
         len(wall_h_f), len(wall_uh_f)))

# ============================================================================
# 4.  EDGE GROUPS — axial (per zone) and radial-to-wall (for the BL)
# ============================================================================
axial_heat = []
axial_ent_lo, axial_ent_hi = [], []    # graded zones, split by edge direction:
axial_exit_lo, axial_exit_hi = [], []  # *_lo = edge start vertex at the lower Z
bl_v0, bl_v1 = [], []          # radial-to-wall, split by which vertex is the wall
r_tol = 0.05 * R

for e in geompy.SubShapeAll(domain, geompy.ShapeType["EDGE"]):
    verts = geompy.SubShapeAll(e, geompy.ShapeType["VERTEX"])
    if len(verts) != 2:
        continue
    p0 = geompy.PointCoordinates(verts[0])
    p1 = geompy.PointCoordinates(verts[1])
    dx, dy, dz = abs(p1[0] - p0[0]), abs(p1[1] - p0[1]), abs(p1[2] - p0[2])

    # Axial edge (along Z) -> bucket by the zone it spans; the entrance and
    # exit zones also split by direction so the grading can put the fine end
    # towards the heated section on every edge.
    if dz > tol and dx < tol and dy < tol:
        zmid = 0.5 * (p0[2] + p1[2])
        lo_first = p0[2] < p1[2]
        if zmid < Z_HEAT_START - tol:
            (axial_ent_lo if lo_first else axial_ent_hi).append(e)
        elif zmid > Z_HEAT_END + tol:
            (axial_exit_lo if lo_first else axial_exit_hi).append(e)
        else:
            axial_heat.append(e)
        continue

    # In-plane edge: keep only the radial-to-wall ones (exactly one endpoint
    # on the arc). Arcs and inner edges are left to the global hypothesis.
    if dz > tol:
        continue
    on0 = abs(math.sqrt(p0[0] ** 2 + p0[1] ** 2) - R) < r_tol
    on1 = abs(math.sqrt(p1[0] ** 2 + p1[1] ** 2) - R) < r_tol
    if on0 == on1:
        continue
    (bl_v0 if on0 else bl_v1).append(e)

print("  axial edges: entrance=%d heated=%d exit=%d"
      % (len(axial_ent_lo) + len(axial_ent_hi), len(axial_heat),
         len(axial_exit_lo) + len(axial_exit_hi)))
print("  BL radial edges: wall-at-v0=%d wall-at-v1=%d" % (len(bl_v0), len(bl_v1)))

# ============================================================================
# 5.  MESH
# ============================================================================
mesh = smesh.Mesh(domain, "QuarterPipe_Mesh")
mesh.Segment().NumberOfSegments(N_QUARTER)      # global default
mesh.Quadrangle()
mesh.Hexahedron()

if axial_heat and n_ax_heat > 0:
    mesh.Segment(geom=make_group(axial_heat, "Axial_heated", "EDGE")) \
        .NumberOfSegments(n_ax_heat)

# BL grading. The wall must be the FINE end: the v0 group (wall at the edge
# start) grows start->end so it takes scale > 1; the v1 group (wall at the
# end) must shrink towards the end, so it takes the reciprocal.
factor_v0 = (1.0 / BL_SCALE) if BL_FLIP else BL_SCALE
factor_v1 = BL_SCALE if BL_FLIP else (1.0 / BL_SCALE)


def apply_graded(edges, name, n, scale):
    if not edges:
        return
    sub = mesh.Segment(geom=make_group(edges, name, "EDGE"))
    hyp = sub.NumberOfSegments(n, scale)
    # DistrType 1 = scale (geometric). Set explicitly: it also defeats
    # SALOME's hypothesis reuse, which would otherwise share one object
    # between the two groups and give them the same direction.
    hyp.SetDistrType(1)
    hyp.SetScaleFactor(scale)
    hyp.SetNumberOfSegments(n)
    print("  [%s] n=%d scale=%.4g" % (name, hyp.GetNumberOfSegments(),
                                      hyp.GetScaleFactor()))


# Entrance/exit axial grading: fine end towards the heated section. The scale
# factor is last/first along the edge, so it depends on each edge's direction.
if abs(AXIAL_RATIO - 1.0) < 1e-12:
    for edges, n, label in ((axial_ent_lo + axial_ent_hi, n_ax_ent, "entrance"),
                            (axial_exit_lo + axial_exit_hi, n_ax_exit, "exit")):
        if edges and n > 0:
            mesh.Segment(geom=make_group(edges, "Axial_" + label, "EDGE")) \
                .NumberOfSegments(n)
else:
    apply_graded(axial_ent_lo,  "Axial_entrance_lo", n_ax_ent, 1.0 / AXIAL_RATIO)
    apply_graded(axial_ent_hi,  "Axial_entrance_hi", n_ax_ent, AXIAL_RATIO)
    apply_graded(axial_exit_lo, "Axial_exit_lo",     n_ax_exit, AXIAL_RATIO)
    apply_graded(axial_exit_hi, "Axial_exit_hi",     n_ax_exit, 1.0 / AXIAL_RATIO)

apply_graded(bl_v0, "BL_WallAtV0", N_RADIAL, factor_v0)
apply_graded(bl_v1, "BL_WallAtV1", N_RADIAL, factor_v1)

print("[MESH] computing ...")
if not mesh.Compute():
    raise SystemExit("mesh computation FAILED — see the SALOME message log")

for grp, name in ((g_inlet, "Inlet"), (g_outlet, "Outlet"),
                  (g_sym1, "Symmetry1"), (g_sym2, "Symmetry2"),
                  (g_wall_uh, "Wall_Unheated"), (g_wall_h, "Wall_Heated")):
    mesh.GroupOnGeom(grp, name, SMESH.FACE)

# ============================================================================
# 6.  EXPORT + STATS
# ============================================================================
out = os.environ.get("SALOME_MESH_OUT")
if not out:
    raise SystemExit("SALOME_MESH_OUT is not set — nowhere to write the mesh")
out = os.path.abspath(out)
outdir = os.path.dirname(out)
if outdir and not os.path.isdir(outdir):
    os.makedirs(outdir)

mesh.ExportUNV(out)

stats = {
    "nodes":       mesh.NbNodes(),
    "hexahedra":   mesh.NbHexas(),
    "volumes":     mesh.NbVolumes(),
    "faces":       mesh.NbFaces(),
    "first_cell_axis_mm": fc_axis,
    "first_cell_diag_mm": fc_diag,
    "bl_scale":    BL_SCALE,
    "axial_cells": {"entrance": n_ax_ent, "heated": n_ax_heat, "exit": n_ax_exit},
    "axial_ratio": AXIAL_RATIO,
    "z_heated":    [Z_HEAT_START, Z_HEAT_END],
    "z_total":     Z_TOTAL,
    "params":      P,
}
with open(os.path.splitext(out)[0] + ".mesh.json", "w") as fh:
    json.dump(stats, fh, indent=2, sort_keys=True)

print("=" * 68)
print("  nodes        : %d" % stats["nodes"])
print("  hexahedra    : %d" % stats["hexahedra"])
print("  volume cells : %d" % stats["volumes"])
print("  face cells   : %d" % stats["faces"])
print("  wrote %s" % out)
print("=" * 68)

if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()

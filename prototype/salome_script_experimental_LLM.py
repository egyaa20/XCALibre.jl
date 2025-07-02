#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mesh_structured_hex.py – 1×0.2×0.2 block with a structured, undistorted hexahedral mesh.
Segments each axis independently (nx, ny, nz) and uses the Hexahedron() algorithm.
Tested with SALOME 9.14 (Windows, Python 3.9).
"""

import salome
# Initialize SALOME (GUI or headless)
try:
    salome.salome_init()
except Exception:
    salome.salome_init_without_session()

from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
import SMESH

# User parameters
Lx, Ly, Lz = 1.0, 0.2, 0.2    # block dimensions (m)
nx, ny, nz = 50, 10, 10      # subdivisions along X, Y, Z
OUT_BASENAME = "structured_hex"

# 1) Build geometry
geom = geomBuilder.New()
box  = geom.MakeBoxDXDYDZ(Lx, Ly, Lz)
geom.addToStudy(box, "HexBlock")

# 2) Create mesh object
smesh = smeshBuilder.New()
mesh  = smesh.Mesh(box, "HexMesh")

# 3) Identify edges by orientation
edges = geom.SubShapeAll(box, geom.ShapeType["EDGE"])
x_edges, y_edges, z_edges = [], [], []
for e in edges:
    # Get the two end vertices of edge e
    verts = geom.SubShapeAll(e, geom.ShapeType["VERTEX"])
    coords = [geom.PointCoordinates(geom.MakeCDG(v)) for v in verts]
    dx = abs(coords[0][0] - coords[1][0])
    dy = abs(coords[0][1] - coords[1][1])
    dz = abs(coords[0][2] - coords[1][2])
    if dx > 1e-6 and dy < 1e-6 and dz < 1e-6:
        x_edges.append(e)
    elif dy > 1e-6 and dx < 1e-6 and dz < 1e-6:
        y_edges.append(e)
    elif dz > 1e-6 and dx < 1e-6 and dy < 1e-6:
        z_edges.append(e)

# 4) 1D segmentation: set number of segments per axis
segx = mesh.Segment()
segx.NumberOfSegments(nx, x_edges)
segx.Propagation()

segy = mesh.Segment()
segy.NumberOfSegments(ny, y_edges)
segy.Propagation()

segz = mesh.Segment()
segz.NumberOfSegments(nz, z_edges)
segz.Propagation()

# 5) 2D quadrangle on all faces
mesh.Quadrangle()

# 6) 3D hexahedral algorithm (no parameters)
mesh.Hexahedron()

# 7) Compute and export
if not mesh.Compute():
    raise RuntimeError("Hexahedral meshing failed – check segmentation and geometry")

mesh.ExportUNV(f"{OUT_BASENAME}.unv", renumber=True)
print(f"Structured hexahedral mesh ({nx}×{ny}×{nz}) written to {OUT_BASENAME}.unv/.med")







# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# mesh_boundary_layer_numbered.py – 1×0.2×0.2 block with a 5-layer viscous prism layer only at the x=0 face,
# with proper element numbering for UNV export.
# Tested with SALOME 9.14 (Windows, Python 3.9).
# """

# import salome
# # Initialize SALOME (GUI or headless)
# try:
#     salome.salome_init()
# except Exception:
#     salome.salome_init_without_session()

# import SMESH
# from salome.geom import geomBuilder
# from salome.smesh import smeshBuilder

# # User parameters
# Lx, Ly, Lz = 1.0, 0.2, 0.2      # block dimensions (m)
# thickness = 0.01               # total prism layer thickness (m)
# n_layers  = 5                  # number of prism layers
# stretch   = 1.3                # growth factor between layers
# OUT_BASENAME = "mesh_bl_only_x0"

# # 1) Geometry: create block
# geom = geomBuilder.New()
# box  = geom.MakeBoxDXDYDZ(Lx, Ly, Lz)
# geom.addToStudy(box, "Block")

# # 2) Identify face at x=0 via centroid
# faces = geom.SubShapeAll(box, geom.ShapeType["FACE"])
# face0 = None
# for f in faces:
#     cdg = geom.MakeCDG(f)
#     cx, cy, cz = geom.PointCoordinates(cdg)
#     if abs(cx) < 1e-6:
#         face0 = f
#         break
# if face0 is None:
#     raise RuntimeError("Could not find face at x=0")

# # 3) Map faces to their IDs
# face_ids    = geom.SubShapeAllIDs(box, geom.ShapeType["FACE"])
# shape_to_id = {f: fid for f, fid in zip(faces, face_ids)}
# ignore_ids  = [shape_to_id[f] for f in faces if f is not face0]

# # 4) Create base mesh and viscous-layer builder
# smesh = smeshBuilder.New()
# mesh  = smesh.Mesh(box, "BL_Only_x0")
# vl    = mesh.ViscousLayerBuilder()
# vl.setBuilderParameters(thickness, n_layers, stretch, ignore_ids, groupName="BL_x0")

# # 5) Mesh the shrunk geometry (tetrahedra)
# shrink = vl.GetShrinkGeometry()
# # Identify corresponding shrunk face for quadrangle
# shrink_faces = geom.SubShapeAll(shrink, geom.ShapeType["FACE"])
# face_shrink0 = None
# for f in shrink_faces:
#     cdg = geom.MakeCDG(f)
#     x0, y0, z0 = geom.PointCoordinates(cdg)
#     if abs(x0 - thickness) < 1e-6:
#         face_shrink0 = f
#         break
# if face_shrink0 is None:
#     raise RuntimeError("Could not find shrunk face for prism layers")

# shrinkMesh = smeshBuilder.New().Mesh(shrink, "Shrink_x0")
# shrinkMesh.Segment().NumberOfSegments(8)
# shrinkMesh.Triangle()
# shrinkMesh.Quadrangle(face_shrink0)
# shrinkMesh.Tetrahedron()
# if not shrinkMesh.Compute():
#     raise RuntimeError("Shrink mesh failed")

# # 6) Build final mesh with prism layers
# final = vl.AddLayers(shrinkMesh)

# # 7) Define a global 1D algorithm to avoid warnings
# seg_final = final.Segment()
# seg_final.LocalLength(thickness / n_layers)
# seg_final.Propagation()

# # 8) Compute final mesh
# if not final.Compute():
#     raise RuntimeError("Final mesh compute failed")

# # 9) Export results
# # not necessary

# print(f"Boundary-layer mesh written to {OUT_BASENAME}.med")
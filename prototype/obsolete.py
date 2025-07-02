# #!/usr/bin/env python3
# # bfs_hexmesh_fixed.py – BFS hex mesh with viscous layers, no GetID() calls

# import salome
# try:
#     salome.salome_init()
# except:
#     salome.salome_init_without_session()

# import GEOM, SMESH
# from salome.geom import geomBuilder
# from salome.smesh import smeshBuilder

# # -----------------------------------------------------------------------------
# # 0) PARAMETERS (SI units)
# # -----------------------------------------------------------------------------
# x_inlet   = -0.1016      # m
# x_step    =   0.0000     # m
# x_outlet  =   0.2813     # m
# H_up      =   0.045008   # m
# h_step    =   0.011252   # m
# H_dn      =   H_up - h_step
# z_thick   =   0.0100     # m

# tol       =   1e-6       # geometric tolerance

# # mesh resolution (~160 k cells)
# n_x_up    =   80
# n_x_dn    =  160
# n_y_up    =   30
# n_y_dn    =   24
# n_z       =    3

# # viscous prism layers (Δy₁ = 1 µm)
# first_layer = 1.0e-6
# n_layers    = 10
# growth      = 1.20

# # -----------------------------------------------------------------------------
# # 1) BUILD & STUDY GEOMETRY
# # -----------------------------------------------------------------------------
# geom = geomBuilder.New()

# # 2D profile → 3D extrusion (BFS: tall→short)
# P = {
#   'A': geom.MakeVertex(x_inlet , 0.0   , 0.0),
#   'B': geom.MakeVertex(x_step  , 0.0   , 0.0),
#   'C': geom.MakeVertex(x_step  , H_dn  , 0.0),
#   'D': geom.MakeVertex(x_outlet, H_dn  , 0.0),
#   'E': geom.MakeVertex(x_outlet, H_up  , 0.0),
#   'F': geom.MakeVertex(x_inlet , H_up  , 0.0),
# }
# edges = [geom.MakeEdge(P[i],P[j]) for i,j in
#          [('A','B'),('B','C'),('C','D'),
#           ('D','E'),('E','F'),('F','A')]]
# wire2d = geom.MakeWire(edges)
# face2d = geom.MakeFace(wire2d, True)
# solid3d= geom.MakePrismVecH(face2d,
#            geom.MakeVectorDXDYDZ(0,0,z_thick), 1)
# geom.addToStudy(solid3d, "BFS_solid")

# # -----------------------------------------------------------------------------
# # 2) PICK & TAG FACES (via BoundingBox)
# # -----------------------------------------------------------------------------
# def pick_face(test):
#     for f in geom.SubShapeAll(solid3d, geom.ShapeType["FACE"]):
#         xi,xa,yi,ya,zi,za = geom.BoundingBox(f)
#         if test(xi,xa,yi,ya,zi,za):
#             return f
#     raise RuntimeError("face not found")

# face_inlet  = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(xi-x_inlet)<tol and abs(xa-x_inlet)<tol)
# face_outlet = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(xi-x_outlet)<tol and abs(xa-x_outlet)<tol)
# face_step   = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(xi-x_step)<tol and abs(xa-x_step)<tol and ya>tol)
# face_base   = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(yi)<tol and xa<x_step+tol)
# face_bottom = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(yi)<tol and xi>x_step-tol)
# face_top    = pick_face(lambda xi,xa,yi,ya,zi,za:
#                   abs(ya-H_up)<tol)

# # tag geometry faces in study (optional)
# for name,f in [
#     ("inlet",face_inlet), ("outlet",face_outlet),
#     ("stepWall",face_step), ("stepBase",face_base),
#     ("bottomWall",face_bottom), ("topWall",face_top)
# ]:
#     geom.addToStudyInFather(solid3d, f, name)

# # -----------------------------------------------------------------------------
# # 3) MAP GEOMETRY FACES → IDs (for viscous‐layer ignore list)
# # -----------------------------------------------------------------------------
# geo_faces   = geom.SubShapeAll(solid3d, geom.ShapeType["FACE"])
# geo_face_ids= geom.SubShapeAllIDs(solid3d, geom.ShapeType["FACE"])

# # Helper: find the ID of a picked face by comparing its bounding box
# def find_geo_id(target_face):
#     bbox_t = geom.BoundingBox(target_face)
#     for fid, face in zip(geo_face_ids, geo_faces):
#         if all(abs(a - b) < tol for a, b in zip(bbox_t, geom.BoundingBox(face))):
#             return fid
#     raise RuntimeError("Could not map face to ID")

# # Build the ignore list by matching bounding boxes
# ignore_ids = [
#     find_geo_id(face_inlet),
#     find_geo_id(face_outlet)
# ]
# # -----------------------------------------------------------------------------
# # 4) BUILD STRUCTURED HEX MESH
# # -----------------------------------------------------------------------------
# smesh = smeshBuilder.New()
# mesh  = smesh.Mesh(solid3d, "BFS_hex")

# seg = mesh.Segment()
# seg.Propagation()

# def set_length(hyp, L):
#     if L > 1e-12:
#         hyp.LocalLength(L)

# for e in geom.SubShapeAll(solid3d, geom.ShapeType["EDGE"]):
#     v1,v2 = geom.SubShapeAll(e, geom.ShapeType["VERTEX"])
#     x1,y1,_= geom.PointCoordinates(v1)
#     x2,y2,_= geom.PointCoordinates(v2)

#     # z-edges
#     if abs(x1-x2)>tol and abs(y1-y2)>tol:
#         seg.NumberOfSegments(n_z); continue

#     # horizontal
#     if abs(y1-y2)<tol:
#         L=abs(x2-x1)
#         if max(x1,x2)<=x_step+tol:
#             set_length(seg, L/n_x_up)
#         else:
#             set_length(seg, L/n_x_dn)
#         continue

#     # vertical
#     if abs(x1-x2)<tol:
#         L=abs(y2-y1)
#         if x1 < x_step+tol:
#             set_length(seg, L/n_y_up)
#         else:
#             set_length(seg, L/n_y_dn)

# assert mesh.Compute(), "base mesh failed"

# # -----------------------------------------------------------------------------
# # 5) ADD VISCOUS PRISM LAYERS (ignore inlet/outlet)
# # -----------------------------------------------------------------------------
# vl = mesh.ViscousLayerBuilder()
# totalTh = first_layer*(growth**n_layers - 1)/(growth - 1)
# vl.setBuilderParameters(totalTh, n_layers, growth,
#                         ignore_ids, groupName="prism")
# mesh_with_prisms = vl.AddLayers(mesh)

# # -----------------------------------------------------------------------------
# # 6) TAG MESH GROUPS FOR OPENFOAM BCs
# # -----------------------------------------------------------------------------
# for name, f in [
#     ("inlet",face_inlet), ("outlet",face_outlet),
#     ("stepWall",face_step), ("stepBase",face_base),
#     ("bottomWall",face_bottom), ("topWall",face_top)
# ]:
#     mesh_with_prisms.GroupOnGeom(f, name, SMESH.FACE)

# print("Hex cells   :", mesh_with_prisms.NbHexas())
# print("Prism cells :", mesh_with_prisms.NbPrisms())
# print("Total cells :", mesh_with_prisms.NbElements())


#!/usr/bin/env python3
# bfs_hexmesh_fixed.py – BFS hex mesh with viscous layers, no GetID() calls
# Corrected version (v2) to fix face identification logic

import salome
try:
    salome.salome_init()
except:
    # Fallback for environments where salome_init() might fail without a session
    salome.salome_init_without_session()

import GEOM
import SMESH
from salome.geom import geomBuilder
from salome.smesh import smeshBuilder
import math # Import math for isnan check if needed, although tol comparison is better

# -----------------------------------------------------------------------------
# 0) PARAMETERS (SI units)
# -----------------------------------------------------------------------------
x_inlet   = -0.1016      # m
x_step    =  0.0000      # m
x_outlet  =  0.2813      # m
H_up      =  0.045008    # m
h_step    =  0.011252    # m
H_dn      =  H_up - h_step # Calculated height downstream (Y-coordinate of bottom wall after step)
z_thick   =  0.0100      # m

tol       =  1e-6        # geometric tolerance (important for comparisons)

# mesh resolution (~160 k cells)
n_x_up    =  80
n_x_dn    = 160
n_y_up    =  30
n_y_dn    =  24
n_z       =   3

# viscous prism layers (Δy₁ = 1 µm)
first_layer = 1.0e-6
n_layers    = 10
growth      = 1.20

# -----------------------------------------------------------------------------
# 1) BUILD GEOMETRY
# -----------------------------------------------------------------------------
geom = geomBuilder.New()

# 2D profile → 3D extrusion (BFS: tall→short)
P = {
    'A': geom.MakeVertex(x_inlet , 0.0   , 0.0),
    'B': geom.MakeVertex(x_step  , 0.0   , 0.0),
    'C': geom.MakeVertex(x_step  , H_dn  , 0.0), # Point at bottom corner after step
    'D': geom.MakeVertex(x_outlet, H_dn  , 0.0), # Point at bottom corner, outlet
    'E': geom.MakeVertex(x_outlet, H_up  , 0.0),
    'F': geom.MakeVertex(x_inlet , H_up  , 0.0),
}
# Ensure points are created before making edges
geom.addToStudy(P['A'], 'P_A')
geom.addToStudy(P['B'], 'P_B')
geom.addToStudy(P['C'], 'P_C')
geom.addToStudy(P['D'], 'P_D')
geom.addToStudy(P['E'], 'P_E')
geom.addToStudy(P['F'], 'P_F')

edges = [geom.MakeEdge(P[i],P[j]) for i,j in
         [('A','B'),('B','C'),('C','D'), # Edges defining the profile
          ('D','E'),('E','F'),('F','A')]]
wire2d = geom.MakeWire(edges)
face2d = geom.MakeFace(wire2d, True) # Use planar face creation
solid3d= geom.MakePrismVecH(face2d,
           geom.MakeVectorDXDYDZ(0,0,z_thick), 1) # Extrude along Z
geom.addToStudy(solid3d, "BFS_solid")

# -----------------------------------------------------------------------------
# 2) GET ALL FACES AND IDs
# -----------------------------------------------------------------------------
# Retrieve all face objects and their corresponding IDs from the solid
geo_faces    = geom.SubShapeAll(solid3d, geom.ShapeType["FACE"])
geo_face_ids = geom.SubShapeAllIDs(solid3d, geom.ShapeType["FACE"])

# Create a dictionary mapping the retrieved face objects to their IDs
# This map will be used later to get IDs for the ignore list
shape_to_id  = {f: fid for f, fid in zip(geo_faces, geo_face_ids)}

# -----------------------------------------------------------------------------
# 3) IDENTIFY & TAG SPECIFIC FACES (using BoundingBox on retrieved faces)
# -----------------------------------------------------------------------------
# Initialize variables to store the identified face objects
face_inlet = None
face_outlet = None
face_step = None
face_base = None
face_bottom = None
face_top = None
face_front = None # Initialize front/back faces
face_back = None

# Iterate through the list of faces retrieved in step 2
for f in geo_faces:
    # Get the bounding box coordinates for the current face
    try:
        xi, xa, yi, ya, zi, za = geom.BoundingBox(f)
    except Exception as e:
        print(f"Warning: Could not get BoundingBox for a face: {e}")
        continue # Skip this face if BBox fails

    # Check bounding box coordinates against known dimensions using the tolerance
    # Use absolute difference comparisons with tolerance for floating point numbers

    # Inlet face check (constant minimum and maximum X, matching x_inlet)
    if abs(xi - x_inlet) < tol and abs(xa - x_inlet) < tol:
        face_inlet = f
        # print(f"Identified Inlet: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue # Move to next face once identified

    # Outlet face check (constant minimum and maximum X, matching x_outlet)
    elif abs(xi - x_outlet) < tol and abs(xa - x_outlet) < tol:
        face_outlet = f
        # print(f"Identified Outlet: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

    # Step wall face check (constant min/max X at x_step, Y between 0 and H_dn)
    # Ensure it's the vertical face connecting y=0 and y=H_dn at x=x_step
    elif abs(xi - x_step) < tol and abs(xa - x_step) < tol and abs(yi - 0.0) < tol and abs(ya - H_dn) < tol:
         face_step = f
         # print(f"Identified Step Wall: {f} (ID: {shape_to_id.get(f)})") # Debug print
         continue

    # Step base face check (constant min/max Y at 0, X < x_step)
    # This is the bottom face before the step
    elif abs(yi - 0.0) < tol and abs(ya - 0.0) < tol and xa < x_step + tol:
        face_base = f
        # print(f"Identified Step Base: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

    # *** CORRECTED *** Bottom wall face check (constant min/max Y at H_dn, X > x_step)
    # This is the bottom face *after* the step
    elif abs(yi - H_dn) < tol and abs(ya - H_dn) < tol and xi > x_step - tol:
        face_bottom = f
        # print(f"Identified Bottom Wall: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

    # Top wall face check (constant min/max Y at H_up)
    elif abs(yi - H_up) < tol and abs(ya - H_up) < tol:
        face_top = f
        # print(f"Identified Top Wall: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

    # Front face check (constant min/max Z at 0)
    elif abs(zi - 0.0) < tol and abs(za - 0.0) < tol:
        face_front = f
        # print(f"Identified Front Face: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

    # Back face check (constant min/max Z at z_thick)
    elif abs(zi - z_thick) < tol and abs(za - z_thick) < tol:
        face_back = f
        # print(f"Identified Back Face: {f} (ID: {shape_to_id.get(f)})") # Debug print
        continue

# --- Verification ---
# Check if all required faces were successfully identified
identified_faces = {
    "inlet": face_inlet, "outlet": face_outlet, "stepWall": face_step,
    "stepBase": face_base, "bottomWall": face_bottom, "topWall": face_top,
    "front": face_front, "back": face_back
}
missing = [name for name, face in identified_faces.items() if face is None]

if missing:
    # Raise an error if any face wasn't found
    raise RuntimeError(f"Could not find faces based on BBox criteria: {', '.join(missing)}")
else:
    print("All required faces identified successfully.")


# --- Tagging ---
# Add the identified faces to the study tree with descriptive names (optional but good practice)
geom.addToStudyInFather(solid3d, face_inlet, "inlet")
geom.addToStudyInFather(solid3d, face_outlet, "outlet")
geom.addToStudyInFather(solid3d, face_step, "stepWall")
geom.addToStudyInFather(solid3d, face_base, "stepBase")
geom.addToStudyInFather(solid3d, face_bottom, "bottomWall")
geom.addToStudyInFather(solid3d, face_top, "topWall")
geom.addToStudyInFather(solid3d, face_front, "front")
geom.addToStudyInFather(solid3d, face_back, "back")


# --- Prepare Ignore List ---
# Get the IDs for the inlet and outlet faces using the shape_to_id map.
ignore_ids = [
    shape_to_id[face_inlet],
    shape_to_id[face_outlet]
]
print(f"Ignoring faces with IDs: {ignore_ids} for viscous layers (Inlet, Outlet)")

# -----------------------------------------------------------------------------
# 4) BUILD STRUCTURED HEX MESH (Base Mesh)
# -----------------------------------------------------------------------------
smesh = smeshBuilder.New()
mesh  = smesh.Mesh(solid3d, "BFS_hex_base") # Name the base mesh

# --- Define Meshing Algorithms ---
# Use Hexahedron algorithm for the 3D solid
hex_algo = mesh.Hexahedron(algo=smeshBuilder.Hexa) # Specify algo type if needed

# Use Quadrangle mapping for the 2D faces (helps structure)
quad_algo = mesh.Quadrangle(algo=smeshBuilder.Quadrangle) # Specify algo type if needed

# Use Wire Discretization for edge segmentation
# This is now done via hypotheses added directly to edges below

# --- Define Hypotheses (Number of Segments) ---
# Create hypotheses for edge divisions
# Note: These are created but applied directly to edges in the loop
nx_up_hyp = smesh.CreateHypothesis('NumberOfSegments')
nx_up_hyp.SetNumberOfSegments(n_x_up)
nx_up_hyp.SetDistrType(smeshBuilder.UniformDistribution) # Ensure uniform distribution

nx_dn_hyp = smesh.CreateHypothesis('NumberOfSegments')
nx_dn_hyp.SetNumberOfSegments(n_x_dn)
nx_dn_hyp.SetDistrType(smeshBuilder.UniformDistribution)

ny_up_hyp = smesh.CreateHypothesis('NumberOfSegments')
ny_up_hyp.SetNumberOfSegments(n_y_up)
ny_up_hyp.SetDistrType(smeshBuilder.UniformDistribution)

ny_dn_hyp = smesh.CreateHypothesis('NumberOfSegments')
ny_dn_hyp.SetNumberOfSegments(n_y_dn)
ny_dn_hyp.SetDistrType(smeshBuilder.UniformDistribution)

nz_hyp    = smesh.CreateHypothesis('NumberOfSegments')
nz_hyp.SetNumberOfSegments(n_z)
nz_hyp.SetDistrType(smeshBuilder.UniformDistribution)

# --- Apply Hypotheses to Edges ---
# Iterate through all edges of the solid geometry
assigned_edges = 0
for edge in geom.SubShapeAll(solid3d, geom.ShapeType["EDGE"]):
    # Get the vertices of the edge
    vertices = geom.SubShapeAll(edge, geom.ShapeType["VERTEX"])
    if len(vertices) != 2: continue # Skip if not a simple edge

    # Get coordinates of the vertices
    p1 = geom.PointCoordinates(vertices[0])
    p2 = geom.PointCoordinates(vertices[1])
    x1, y1, z1 = p1
    x2, y2, z2 = p2

    # Determine edge orientation and apply the correct hypothesis
    # Use tolerance for floating point comparisons

    hypothesis_to_add = None
    # Z-direction edges (constant X and Y)
    if abs(x1 - x2) < tol and abs(y1 - y2) < tol:
        hypothesis_to_add = nz_hyp

    # X-direction edges (constant Y and Z)
    elif abs(y1 - y2) < tol and abs(z1 - z2) < tol:
        # Check if the edge is in the upstream or downstream region based on X-coordinate
        # Use midpoint or max X to determine region
        if max(x1, x2) <= x_step + tol: # Upstream part (includes edge ending exactly at x_step)
            hypothesis_to_add = nx_up_hyp
        else: # Downstream part
            hypothesis_to_add = nx_dn_hyp

    # Y-direction edges (constant X and Z)
    elif abs(x1 - x2) < tol and abs(z1 - z2) < tol:
        # Check if the edge is associated with the upstream or downstream height
        if x1 < x_step + tol: # Vertical edges in the upstream part (x < x_step)
            hypothesis_to_add = ny_up_hyp
        elif abs(x1 - x_step) < tol: # The vertical step edge itself (x = x_step)
             # This edge connects y=0 to y=H_dn. Use n_y_dn? Or n_y_up?
             # Let's check the original logic - it used n_y_up here. Let's stick to that for consistency.
             # If this causes issues, consider n_y_dn or (n_y_up - n_y_dn).
             hypothesis_to_add = ny_up_hyp # Match original script's implied logic
        else: # Vertical edges in the downstream part (x > x_step)
            hypothesis_to_add = ny_dn_hyp

    # Add the determined hypothesis to the edge
    if hypothesis_to_add:
        mesh.AddHypotheses(edge, [hypothesis_to_add])
        assigned_edges += 1
    else:
        print(f"Warning: Edge not classified for hypothesis assignment. Coords: ({x1},{y1})-({x2},{y2})")

print(f"Assigned hypotheses to {assigned_edges} edges.")

# --- Compute the Base Mesh ---
try:
    is_computed = mesh.Compute()
    if not is_computed:
        raise RuntimeError("Base mesh computation failed.")
    print("Base mesh computed successfully.")
    print(f"Base mesh cells: {mesh.NbElements()}")

except Exception as e:
    print(f"Error during base mesh computation: {e}")
    # Consider adding more detailed error handling or inspection here if needed
    raise # Re-raise the exception to stop the script

# -----------------------------------------------------------------------------
# 5) ADD VISCOUS PRISM LAYERS (ignore inlet/outlet)
# -----------------------------------------------------------------------------
# Create a Viscous Layers hypothesis object
vl_hyp = smesh.CreateHypothesis('ViscousLayers')

# Set parameters: total thickness, number of layers, stretching factor
# Calculate total thickness based on first layer, growth, and number of layers
total_thickness = first_layer * (growth**n_layers - 1) / (growth - 1) if abs(growth - 1.0) > tol else first_layer * n_layers
vl_hyp.SetTotalThickness(total_thickness)
vl_hyp.SetNumberLayers(n_layers)
vl_hyp.SetStretchingFactor(growth)
vl_hyp.SetExtrusionMethod(smeshBuilder.SurfaceOffset) # Specify extrusion method if needed

# Specify the face IDs to *exclude* from layer generation (Inlet and Outlet)
# Use the 'ignore_ids' list derived earlier
vl_hyp.SetFacesToExclude(ignore_ids)

# Add the viscous layer hypothesis to the mesh's 3D algorithm
# This tells the Hexahedron algorithm to incorporate these layers
algo_3d = mesh.GetAlgorithm(smeshBuilder.VOLUME, 0) # Get the main 3D algorithm
if algo_3d:
     algo_3d.AddHypothesis(vl_hyp)
     print("Viscous layer hypothesis added to the 3D algorithm.")
else:
     print("Warning: Could not retrieve 3D algorithm to add viscous layers hypothesis.")
     # As a fallback, try adding to the mesh directly, though less standard
     # mesh.AddHypothesis(vl_hyp)


# --- Compute Mesh with Viscous Layers ---
# Recompute the mesh to apply the viscous layers
try:
    is_computed_vl = mesh.Compute()
    if not is_computed_vl:
        raise RuntimeError("Mesh computation with viscous layers failed.")
    print("Mesh with viscous layers computed successfully.")
    # Rename the final mesh for clarity
    smesh.SetName(mesh.GetMesh(), "BFS_mesh_with_prisms")
    mesh_with_prisms = mesh # Use the same mesh object updated with layers

except Exception as e:
    print(f"Error during viscous layer computation: {e}")
    raise

# -----------------------------------------------------------------------------
# 6) CREATE MESH GROUPS FOR OPENFOAM BCs
# -----------------------------------------------------------------------------
# Create mesh groups based on the original geometry faces for boundary conditions
# Use the face objects identified in step 3

# Function to safely create group and print warning on failure
def create_group(mesh_obj, geom_face, name):
    group = mesh_obj.GroupOnGeom(geom_face, name, SMESH.FACE)
    if not group:
        print(f"Warning: Could not create mesh group '{name}'.")
    # else: # Debug print
    #     print(f"Created mesh group '{name}' with {group.NbElements()} faces.")
    return group

# Create groups for all identified faces
inlet_group = create_group(mesh_with_prisms, face_inlet, "inlet")
outlet_group = create_group(mesh_with_prisms, face_outlet, "outlet")
stepWall_group = create_group(mesh_with_prisms, face_step, "stepWall")
stepBase_group = create_group(mesh_with_prisms, face_base, "stepBase")
bottomWall_group = create_group(mesh_with_prisms, face_bottom, "bottomWall")
topWall_group = create_group(mesh_with_prisms, face_top, "topWall")
front_group = create_group(mesh_with_prisms, face_front, "front")
back_group = create_group(mesh_with_prisms, face_back, "back")


# --- Final Output ---
print("\n--- Mesh Statistics ---")
print(f"Hexahedra cells : {mesh_with_prisms.NbHexas()}")
print(f"Prism cells     : {mesh_with_prisms.NbPrisms()}")
# Note: Viscous layers often generate prisms.
print(f"Wedge cells     : {mesh_with_prisms.NbWedges()}")
print(f"Pyramid cells   : {mesh_with_prisms.NbPyramids()}")
print(f"Tetrahedra cells: {mesh_with_prisms.NbTetras()}")
print(f"Total cells     : {mesh_with_prisms.NbElements()}")
print(f"Total nodes     : {mesh_with_prisms.NbNodes()}")

# Optional: Update the study viewer if running interactively
if salome.sg.hasDesktop():
    salome.sg.updateObjBrowser()

print("\nScript finished.")


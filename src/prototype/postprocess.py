import pyvista as pv
import pandas as pd
import os
import re

# --- Configuration ---
field_name = "p"
start_point = [0.5, 0.0, 0.0]
end_point   = [0.5, 1.0, 0.0]
resolution  = 100
output_csv  = "probe_results.csv"

def get_step_number(filepath):
    """Extracts the number from 'time_#.vtk' for correct sorting."""
    match = re.search(r'time_(\d+)', filepath)
    return int(match.group(1)) if match else 0

# --- File Discovery ---
files = [f for f in os.listdir('.') if f.endswith('.vtk') and 'time_' in f]
files.sort(key=get_step_number)

if not files:
    print("No .vtk files found with 'time_' in the name.")
    exit()

print(f"🚀 Found {len(files)} files. Starting extraction...")

# Results dictionary: {arc_length: [...], time_1: [...], time_2: [...]}
results = {}

for i, file in enumerate(files):
    # Load the mesh
    mesh = pv.read(file)
    
    # Check if data is in Points or Cells
    # Interpolation along a line requires Point Data
    if field_name not in mesh.point_data:
        if field_name in mesh.cell_data:
            mesh = mesh.cell_data_to_point_data()
        else:
            print(f"⚠️ Warning: Field '{field_name}' not found in {file}. Skipping.")
            continue

    # Sample the line
    # 
    sampled = mesh.sample_over_line(start_point, end_point, resolution=resolution)
    
    # On first file, save the distance (x-axis for your plots)
    if i == 0:
        results['arc_length'] = sampled["Distance"]
    
    # Save the field data
    col_name = file.replace(".vtk", "")
    results[col_name] = sampled.point_data[field_name]
    
    print(f"✅ Processed {file}")

# --- Save to CSV ---
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False)
print(f"\n🎉 Success! Data saved to {output_csv}")


# import matplotlib.pyplot as plt
# df.set_index('arc_length').plot()
# plt.ylabel(field_name)
# plt.show()
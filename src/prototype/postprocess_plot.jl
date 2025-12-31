using Plots
using CSV
using DataFrames

# Load the data we just saved
df = CSV.read("probe_results.csv", DataFrame)

# Create a plot
p = plot(xlabel="Distance along line", ylabel="Pressure (p)", title="Pressure over Line over Time")

# Plot every time step (columns 2 to end)
for name in names(df)[2:end]
    plot!(df.arc_length, df[!, name], label=name, lw=2, alpha=0.8)
end

display(p)
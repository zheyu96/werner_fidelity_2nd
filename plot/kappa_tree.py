import matplotlib.pyplot as plt
import pandas as pd

# Data from the user
# data = {
#     "x": [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3],
#     "line1": [0.759751, 0.767366, 0.772817, 0.7768, 0.779758, 0.781983, 0.783669, 0.784952, 0.785925, 0.786656,
#               0.787196, 0.787582, 0.787842, 0.787998, 0.788067, 0.788062, 0.787996, 0.787876, 0.78771, 0.787503, 0.787261],
#     "line2": [0.759751, 0.767137, 0.77235, 0.776102, 0.778845, 0.780871, 0.782373, 0.783483, 0.784294, 0.784872,
#               0.785266, 0.785512, 0.785638, 0.785664, 0.785606, 0.785479, 0.785292, 0.785054, 0.784771, 0.78445, 0.784095]
# }

data = {
    "x": [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.1, 2.2, 2.3],
    "line1": [0.759751, 0.767366, 0.772817, 0.7768, 0.779758, 0.781983, 0.783669, 0.784952, 0.785925, 0.786656,
              0.787196, 0.787582, 0.787842, 0.787998],
    "line2": [0.759751, 0.767137, 0.77235, 0.776102, 0.778845, 0.780871, 0.782373, 0.783483, 0.784294, 0.784872,
              0.785266, 0.785512, 0.785638, 0.785664]
}

# data = {
#     "x": [1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2],
#     "line1": [0.759751, 0.767366, 0.772817, 0.7768, 0.779758, 0.781983, 0.783669, 0.784952, 0.785925, 0.786656,
#               0.787196],
#     "line2": [0.759751, 0.767137, 0.77235, 0.776102, 0.778845, 0.780871, 0.782373, 0.783483, 0.784294, 0.784872,
#               0.785266]
# }


# Create a DataFrame
df = pd.DataFrame(data)

plt.rcParams.update({'font.size': 20})
plt.rc('text', usetex=True)
plt.rc('mathtext', fontset='stix')
plt.rc('font', family='serif')

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(df["x"], df["line1"], label=r"Skewed strategy tree", marker='o', markersize=9)
plt.plot(df["x"], df["line2"], label=r"Complete strategy tree", marker='x', markersize=9)

# Customize the plot
# plt.title("Plot of Line 1 and Line 2", fontsize=14)
plt.xlabel(r"$\kappa$", fontsize=24)
plt.ylabel(r"Fidelity", fontsize=24, labelpad=10)
plt.xticks([1, 1.25, 1.5, 1.75, 2, 2.25])
plt.legend(fontsize=16)
plt.grid(True)
plt.savefig('comF1F2.png', dpi=1300, bbox_inches='tight')
plt.savefig('comF1F2.eps', dpi=1000, bbox_inches='tight')

# Show the plot
plt.show()
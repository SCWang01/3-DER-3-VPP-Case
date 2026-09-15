import matplotlib.pyplot as plt


plt.rcParams["font.family"] = "Calibri"
plt.rcParams["axes.unicode_minus"] = False

periods = [1, 2, 3, 4, 5]
time_edges = [0, 1, 2, 3, 4, 5]
outputs = {
    "1": [15, 15, 5, 25, 15],
    "2": [-10, 5, 5, -15, 5],
    "3": [-5, -20, -10, -10, -20],
}


fig, axes = plt.subplots(3, 1, figsize=(8, 7), sharex=True, sharey=True)

for ax, (name, values) in zip(axes, outputs.items()):
    for left, right, value in zip(time_edges[:-1], time_edges[1:], values):
        color = "#5B9BD5" if value >= 0 else "#C00000"
        ax.hlines(value, left, right, color=color, linewidth=3)
        ax.vlines(left, 0, value, color=color, linewidth=1.2, alpha=0.55)
        ax.vlines(right, 0, value, color=color, linewidth=1.2, alpha=0.55)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel(f"VPP(DER) {name}\noutput (kW)")
    ax.set_yticks([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25])
    ax.set_ylim(-27, 27)
    ax.set_xticks(time_edges)
    ax.set_xlim(time_edges[0], time_edges[-1])
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_title(f"VPP(DER) {name}")

axes[-1].set_xlabel("One period (5 days)")
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig("prosumers_output_3_subplots.png", dpi=600)
print("Plot saved as 'prosumers_output_3_subplots.png'")

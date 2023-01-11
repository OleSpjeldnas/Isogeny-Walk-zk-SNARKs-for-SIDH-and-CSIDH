import matplotlib.pyplot as plt

# Define the array A
A = [
    14170,
    19736,
    21912,
    27356,
    33211,
    36087,
    42974,
    51108,
    55199,
    63010,
    71110,
    76863,
    84327,
    92906,
    99174,
    109709,
]
A = [a / 1000 for a in A]
# Create a figure and axis object
fig, ax = plt.subplots()

# Set the x-values to be 2^(4+i) for i in range(15)
x = [2 ** (4 + i) for i in range(16)]
# Plot the points (x, A[i])
ax.plot(x, A, "o-", markersize=10)

# Set the x-axis label and y-axis label
ax.set_xlabel("Witness Size")
ax.set_ylabel("Proof Size (kB)")
# Set the title of the plot
ax.set_title("Proof Size")
ax.set_xscale("log", base=2)
# ax.set_yscale('log', base=2)
# ax.set_xticklabels([f"2^{i}" for i in range(4,20)])
# ax.set_yticklabels([f"2^{i}" for i in range(4,20)])


# Show the plot
plt.show()

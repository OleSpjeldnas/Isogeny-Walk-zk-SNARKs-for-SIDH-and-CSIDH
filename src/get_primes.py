import matplotlib.pyplot as plt

# Define the array A
A = [
20508,
24031,
29604,
37771,
41319,
46825,
57426,
61586,
70200,
79322,
86556,
94525,
103009,
110175,
122696,
132168
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
ax.set_ylabel("Size (kB)")
# Set the title of the plot
ax.set_title("Proof Size")
ax.set_xscale("log", base=2)
#ax.set_yscale('log', base=2)
# ax.set_xticklabels([f"2^{i}" for i in range(4,20)])
# ax.set_yticklabels([f"2^{i}" for i in range(4,20)])


# Show the plot
plt.show()

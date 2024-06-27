import numpy as np
import matplotlib.pyplot as plt

# Updated coefficients for the 20 cubic polynomials based on the latest provided data
coefficients_updated_3 = [
    [-0.254626, -6.87491, -60.8332, -176.126],
    [0.205731, 4.17367, 27.5555, 59.5768],
    [0.245201, 5.00253, 33.3575, 73.1148],
    [-0.216734, -3.3123, -16.5315, -26.6632],
    [-0.00676361, -0.162739, -0.783675, -0.41684],
    [0.298589, 3.50149, 13.8732, 19.1257],
    [-1.73299, -14.7827, -40.9794, -35.7269],
    [6.22188, 32.9465, 54.479, 27.912],
    [-11.6901, -20.7895, 0.743018, 10],
    [11.0037, -20.7895, 0.743018, 10],
    [-4.8684, 26.8267, -46.8732, 25.8721],
    [0.685617, -6.49736, 19.775, -18.5601],
    [0.365929, -3.62017, 11.1434, -9.92847],
    [-0.237734, 3.62379, -17.8325, 28.706],
    [-0.139292, 2.14716, -10.4493, 16.4007],
    [0.0870017, -1.92613, 13.9904, -32.4787],
    [0.161685, -3.49448, 24.9689, -58.0952],
    [-0.319443, 8.05258, -67.4076, 188.242],
    [0.354185, -10.1354, 96.2838, -302.832],
    [-0.163397, 5.3921, -58.9907, 214.75]
]

x_points = np.arange(-9, 12)  # from -9 to 11
y_points = [
    0.1270, 0.9134, 0.6324, 0.0975, 0.2785, 0.5469, 0.9575, 0.9649, 0.1576, 10,
    0.9572, 0.4854, 0.8003, 0.1419, 0.4218, 0.9157, 0.9157, 0.7922, 0.9595, 0.6557, 0.8147
]

# Create plot with updated coefficients and the given points
plt.figure(figsize=(14, 8))

# Plot each polynomial in its respective interval
for i, coeffs in enumerate(coefficients_updated_3):
    # Define the interval
    x = np.linspace(-9 + i, -8 + i, 100)
    # Calculate polynomial values
    p = np.poly1d(coeffs)
    y = p(x)
    # Plot the polynomial
    plt.plot(x, y, label=f'Interval {-9+i} to {-8+i}')

# Add points
plt.scatter(x_points, y_points, color='red', label='Given Points')

plt.title('Graphs of 20 Cubic Polynomials with Given Points (Latest Coefficients)')
plt.xlabel('x')
plt.ylabel('y')
plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
plt.grid(True)
plt.tight_layout()
plt.show()


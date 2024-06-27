import matplotlib.pyplot as plt

def read_points_from_file(filename):
    """读取文件中的点坐标并返回它们的列表"""
    points = []
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                x, y = map(float, line.split(','))
                points.append((x, y))
    return points

def plot_points(points, color, label):
    """将点列表绘制到图上"""
    x_vals, y_vals = zip(*points)
    plt.scatter(x_vals, y_vals, color=color, label=label)

def main():
    filenames = ['0.txt', '1.txt', '2.txt']
    colors = ['red', 'green', 'blue']
    labels = ['0', '1', '2']

    plt.figure(figsize=(8, 6))

    for filename, color, label in zip(filenames, colors, labels):
        points = read_points_from_file(filename)
        plot_points(points, color, label)

    plt.title('Projected points')
    plt.xlabel('Coordinate 1')
    plt.ylabel('Coordinate 2')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == '__main__':
    main()

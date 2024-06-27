import numpy as np
import time


np.set_printoptions(precision=20, suppress=True)


def angle(vector1: np.array([]), vector2: np.array([])):
    """
    Calculate the angle between two vectors
    :param vector1: np.array([])
    :param vector2: np.array([])
    :return:
    """
    # Calculate the dot product
    dot_product = np.dot(vector1, vector2)

    # Calculate L2 norm of two vectors
    norm_1 = np.linalg.norm(vector1)
    norm_2 = np.linalg.norm(vector2)

    # Calculate the cosine value of the angle
    cos_angle = dot_product/(norm_1*norm_2)

    return np.arccos(cos_angle)


def angle2point(angle_ra):
    """
    Given an angle, calculate the corresponding point on unit circle. Suppose the angle is the angle between the vector
    from the origin to the point and the negative x-axis. The angle is in radian.
    :param angle_ra:
    :return:
    """
    x = -np.cos(angle_ra)
    y = np.sin(angle_ra)
    dot = np.array([x, y])
    return dot


def point2angle(point: np.array([])):
    return abs(np.arctan2(point[1], -point[0]))


def find_t(floor, ceiling, error, depth=0):
    """
    Recursively find the reflection point on the unit circle using binary search.
    :param floor:
    :param ceiling:
    :param error:
    :param depth:
    :return:
    """
    # Calculate median value
    mid_angle = (floor + ceiling)/2
    mid_point = angle2point(mid_angle)

    vector1 = P - mid_point
    vector2 = Q - mid_point

    angle1 = angle(vector1, mid_point)
    angle2 = angle(vector2, mid_point)
    delta = angle1 - angle2

    if abs(delta) < error:
        print("Recursion depth:", depth)
        return mid_point
    if delta > 0:
        return find_t(floor, mid_angle, error, depth+1)
    if delta < 0:
        return find_t(mid_angle, ceiling, error, depth+1)


def findimage():
    d1 = np.linalg.norm(P-T)
    d2 = np.linalg.norm(Q-T)
    vector = (T-P)*d2/d1
    return vector + T


# Read 8 testcases
for i in range(0, 8):

    P_inputList = input().strip().split(' ')
    Q_inputList = input().strip().split(' ')

    start_time = time.time()

    P = np.array([float(P_inputList[0]), float(P_inputList[1])])
    Q = np.array([float(Q_inputList[0]), float(Q_inputList[1])])
    print("----------Testcase", i + 1, "----------")
    T = find_t(0, point2angle(Q), 0.00000001)
    print(T)
    print(findimage())

    end_time = time.time()

    print('Time:', end_time - start_time, 's')
    print()

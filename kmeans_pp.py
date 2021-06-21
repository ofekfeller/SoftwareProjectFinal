import sys

import numpy as np
import pandas as pd
import typing
import csv
import mykmeanssp


def read_csv_files_to_numpy(filename1, filename2):
    file = pd.read_csv(filename1, header=None)
    other_file = pd.read_csv(filename2, header=None)
    merged_data = pd.merge(file, other_file, on=0)
    merged_data.set_index(0, inplace=True)

    return merged_data.sort_index().to_numpy()


def find_min_d(vec, centers):
    min_val = np.infty
    for i in range(len(centers)):
        norm = np.linalg.norm(np.subtract(centers[i], vec))
        norm = norm**2
        if norm < min_val:
            min_val = norm

    return min_val


# is there better approach we can handle this with?
def is_int(num):
    int(num)


def get_args(args):
    if len(args) == 4:
        is_int(args[1])
        return int(args[1]), 300, args[2], args[3]
    elif len(args) == 5:
        is_int(args[1])
        is_int(args[2])
        return int(args[1]), int(args[2]), args[3], args[4]
    else:
        raise Exception("Arguments are corrupted")


def main():
    np.random.seed(0)
    args = sys.argv
    K, max_iter, filename1, filename2 = get_args(sys.argv)
    points = read_csv_files_to_numpy(filename1, filename2)
    z = 1
    points_num = points.shape[0]
    assert (K < points_num)
    rnd = np.random.choice(range(points_num))
    indexes = str(rnd)
    centers = [points[rnd, :]]
    for i in range(K-1):
        prob_vec = [find_min_d(points[j, :], centers) for j in range(points_num)]
        total = sum(prob_vec)
        prob_vec = [val/total for val in prob_vec]
        p = np.random.choice([i for i in range(points_num)], p=prob_vec)
        indexes += f",{p}"
        centers.append(points[p, :])

    print(indexes)

    centers = [center.tolist() for center in centers]

    res = mykmeanssp.fit(points.tolist(), centers, [1 for i in range(max_iter)])

    res = np.round(np.array(res),decimals=4)

    for row in res:
        print(",".join(map(str, row)))


if __name__ == '__main__':
    main()

from sklearn.cluster import KMeans
from sklearn.datasets import load_iris
from matplotlib import pyplot


def main():
    data = load_iris().data
    res = list()
    for i in range(1,11):
        kmeans = KMeans(init='k-means++', random_state=0, n_clusters=i).fit(data)
        res.append(kmeans.inertia_)

    max_dif = 0
    max_ind = 0
    for i in range(1,10):
        if res[i]-res[i-1] > max_dif:
            max_ind = i

    fig, ax = pyplot.subplots()
    pyplot.plot([i for i in range(1, 11)], res)
    circ = pyplot.Circle(((2-0.47)/10.45 ,152/711), color='0', radius=0.018, fill=False,transform=ax.transAxes)
    ax.add_patch(circ)
    pyplot.savefig('elbow.png', format='png')


if __name__ == '__main__':
    main()
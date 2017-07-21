import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from sklearn.decomposition import PCA


def parse(filepath):
    sample_names = []
    data = []
    label = []
    hdr_line = ""
    with open(filepath, "r") as f:
        hdr_line = f.readline() # header
        for line in f.readlines():
            tokens = line.strip().split('\t')
            sample_names.append(tokens[0])
            data.append([float(i) for i in tokens[1:-1]])
            label.append(int(tokens[-1]))
### END - with open

    hdr_tokens = hdr_line.strip().split('\t')[1:]
    header = np.asarray(hdr_tokens)
    data = np.asarray([np.asarray(i) for i in data])
    label = np.asarray(label)
    return sample_names, header, data, label

def plot2d(X, y, names, title, figName):
    fig = plt.figure(figsize=(400, 400))
    plt.scatter(X[:, 0], X[:, 1])
    #plt.scatter(X[:, 0], X[:, 1], c=y, cmap=plt.cm.spectral, s=20, label=y)
    print "ok2"
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(title)
    print "ok3"
    #plt.show()
    plt.savefig(figName + '.png')

    

def main():
    filepath = "deep_input_c3_pert_v2LogN_top_cs_3x"
    sample_names, header, X, y = parse(filepath)
    print X.shape, y.shape
    pca = PCA(n_components=2)
    X_2d = pca.fit_transform(X)
    names = ["True", "False"]
    print "ok"
    print X_2d.shape
    plot2d(X_2d[:10], y, names, "PCA Analysis top_cs_3x", "pca_top_cs_3x")
    
    
main()
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score

from node2vec import Node2Vec


def read_labels(csv_path):
    df = pd.read_csv(csv_path)

    label_col = df.columns[-1]
    node_ids = df.index.astype(int)
    labels = df[label_col].astype(str).values
    
    y = pd.Series(labels, index=node_ids)
    return y


def run_node2vec(G, dimensions, p, q, workers, seed):
    model = Node2Vec(
        G,
        dimensions=dimensions,
        walk_length=80,
        num_walks=10,
        p=p,
        q=q,
        workers=workers,
        seed=seed,
    ).fit(window=10, min_count=1, batch_words=4)

    keys = list(model.wv.index_to_key)
    X = np.vstack([model.wv[k] for k in keys])

    node_ids = []
    for k in keys:
        node_ids.append(int(k))
    emb = pd.DataFrame(X, index=node_ids)
    return emb


def cluster_and_compute_nmi(emb, y_true, seed):
    common = emb.index.intersection(y_true.index)

    X = emb.loc[common].values
    y = y_true.loc[common].values

    K = len(pd.unique(y))

    km = KMeans(n_clusters=K, n_init="auto", random_state=seed)
    pred = km.fit_predict(X)

    nmi = normalized_mutual_info_score(y, pred)
    return float(nmi)


def main():
    np.random.seed(42)

    G = nx.read_edgelist("cell_graph.edgelist", nodetype=int).to_undirected()
    y_true = read_labels("Levine_matrix.csv")

    # Problem 1
    emb = run_node2vec(
        G,
        dimensions=128,
        p=1.0,
        q=1.0,
        workers=4,
        seed=42,
    )
    nmi_default = cluster_and_compute_nmi(emb, y_true, seed=42)
    print(f"(1) NMI={nmi_default:.4f}")

    # Problem 2
    dim_list = [32, 64, 128, 256, 512]
    nmis_dim = []
    for d in dim_list:
        emb_d = run_node2vec(G, dimensions=d, p=1.0, q=1.0, workers=4, seed=42)
        nmi_d = cluster_and_compute_nmi(emb_d, y_true, seed=42)
        nmis_dim.append(nmi_d)
        print(f"(2) d={d} - NMI={nmi_d:.4f}")

    plt.figure()
    plt.plot(dim_list, nmis_dim, marker="o")
    plt.xlabel("Embedding dimension")
    plt.ylabel("NMI")
    plt.title("NMI vs embedding dimension")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("p3res/nmi_vs_dimensions.png")
    print("[Part 2] Saved plot: nmi_vs_dimensions.png")

    # Problem 3
    _target = "p" # or "q", I chose "p" in my experiments
    _vals = [0.25, 0.5, 1, 2, 4]

    nmis_sweep = []
    for v in _vals:
        p, q = 1.0, 1.0
        if _target == "p":
            p = v
        else:
            q = v

        emb_v = run_node2vec(G, dimensions=128, p=p, q=q, workers=4, seed=42)
        nmi_v = cluster_and_compute_nmi(emb_v, y_true, seed=42)
        nmis_sweep.append(nmi_v)
        print(f"(3) {_target}={v} (d={128}) - NMI={nmi_v:.4f}")

    plt.figure()
    plt.plot(_vals, nmis_sweep, marker="o")
    plt.xlabel(_target)
    plt.ylabel("NMI")
    plt.title(f"NMI vs {_target} (dims={128})")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"p3res/nmi_vs_{_target}.png")
    print(f"[Part 3] Saved plot: nmi_vs_{_target}.png")


if __name__ == "__main__":
    main()
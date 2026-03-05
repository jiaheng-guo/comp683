import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import normalized_mutual_info_score, roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression


def greedy_modularity_partition(G):
    comms = nx.algorithms.community.greedy_modularity_communities(G)
    labels = np.full(G.number_of_nodes(), -1, dtype=int)
    for cid, nodes in enumerate(comms):
        for v in nodes:
            labels[v] = cid
    return labels


def main():
    np.random.seed(42)

    df = pd.read_csv("Levine_matrix.csv")

    X = df.iloc[:, :-1].to_numpy(dtype=float)
    y_raw = df.iloc[:, -1]

    labeled_mask = ~pd.isna(y_raw)
    X_labeled = X[labeled_mask.values]
    y_labeled = y_raw[labeled_mask].astype(int).to_numpy()

    K = len(np.unique(y_labeled))

    scaler = StandardScaler()
    X_labeled_z = scaler.fit_transform(X_labeled)

    km = KMeans(n_clusters=K, n_init="auto", random_state=42)
    pred_X_labeled = km.fit_predict(X_labeled_z)

    nmi_X = normalized_mutual_info_score(y_labeled, pred_X_labeled)
    print(f"(1) Clustering on X: NMI={nmi_X:.4f}")

    n = X.shape[0]
    G = nx.read_edgelist("cell_graph.edgelist", nodetype=int).to_undirected()
    G.add_nodes_from(range(n))

    pred_G_all = greedy_modularity_partition(G)
    pred_G_labeled = pred_G_all[labeled_mask.values]

    nmi_G = normalized_mutual_info_score(y_labeled, pred_G_labeled)
    print(f"(2) Graph partitioning: NMI={nmi_G:.4f}")

    if nmi_X > nmi_G:
        print("(3) X clustering seems better (higher NMI).")
    else:
        print("(3) Graph partitioning seems better (higher NMI).")

    target = 21
    target_mask = (~pd.isna(y_raw)) & (y_raw.astype("Int64") == target)
    target_idx = np.where(target_mask.values)[0]

    X_all_z = StandardScaler().fit_transform(X)
    pred_X_all = KMeans(n_clusters=K, n_init="auto", random_state=42).fit_predict(X_all_z)

    if len(target_idx) == 0:
        pdc_X = 0
        pdc_G = 0
    else:
        pdc_X = len(set(pred_X_all[target_idx].tolist()))
        pdc_G = len(set(pred_G_all[target_idx].tolist()))

    print(f"(4) #clusters containing pDC(label=21): X={pdc_X}, G={pdc_G}")

    T = {11, 12, 17, 18}
    M = {1, 2, 3}

    keep_mask = (~pd.isna(y_raw)) & (y_raw.astype("Int64").isin(list(T | M)))
    keep_idx = np.where(keep_mask.values)[0]

    X_bin = X[keep_idx]
    y_keep = y_raw.iloc[keep_idx].astype(int).to_numpy()
    y_bin = np.array([0 if lab in T else 1 for lab in y_keep], dtype=int)

    X_bin_z = StandardScaler().fit_transform(X_bin)

    X_tr, X_te, y_tr, y_te = train_test_split(
        X_bin_z, y_bin, test_size=0.2, random_state=42, stratify=y_bin
    )

    clf = LogisticRegression(max_iter=5000, class_weight="balanced", random_state=42)
    clf.fit(X_tr, y_tr)

    prob = clf.predict_proba(X_te)[:, 1]
    fpr, tpr, _ = roc_curve(y_te, prob)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.3f}")
    plt.plot([0, 1], [0, 1], linestyle="--")
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC: T-cells (0) vs Monocytes (1)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig("p2res/roc_t_vs_mono.png")

    print(f"(5) ROC-AUC={roc_auc:.4f}")


if __name__ == "__main__":
    main()
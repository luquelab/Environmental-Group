---
title: Code snippets
parent: Pipeline Overview
nav_order: 1
---

# Code snippets

[Run this notebook in Colab](https://colab.research.google.com/github/luquelab/Environmental-Group/blob/main/Bioinfo_fetch_FASTA.ipynb)

## Upload FASTA File or Fetch

```python
from Bio import SeqIO, Entrez
from io import StringIO
import pandas as pd
import re

Entrez.email = "your_email@example.com"

def df_from_records(records):
    df = pd.DataFrame([{
        "ID": r.id,
        "Description": r.description,
        "Sequence": str(r.seq),
        "Length": len(r.seq)
    } for r in records])

    df = df.drop_duplicates(subset="Sequence").reset_index(drop=True)
    return df
```

## Calculate Properties
```python
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqRecord import SeqRecord
import pandas as pd

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

def calculate_properties(data):
    rows = []

    if isinstance(data, pd.DataFrame):
        source_data = data.to_dict(orient="records")
    else:
        source_data = data

    for item in source_data:
        if isinstance(item, SeqRecord):
            seq = str(item.seq).strip().upper()
            if not seq:
                continue
            rows.append({
                "ID": item.id,
                "Description": item.description,
                "Sequence": seq
            })

        elif isinstance(item, dict):
            seq = str(item.get("Sequence", "")).strip().upper()
            if not seq:
                continue
            rows.append({
                "ID": item.get("ID", "Unknown"),
                "Description": item.get("Description", ""),
                "Sequence": seq
            })

    results = []

    for row in rows:
        seq = row["Sequence"]
        invalid_chars = set(seq) - VALID_AA
        if invalid_chars:
            continue

        analysis = ProteinAnalysis(seq)

        results.append({
            "ID": row["ID"],
            "Length": len(seq),
            "Molecular_Weight": analysis.molecular_weight(),
            "Isoelectric_Point": analysis.isoelectric_point(),
            "Extinction_Coefficient_Reduced": analysis.molar_extinction_coefficient()[0],
            "Extinction_Coefficient_Oxidized": analysis.molar_extinction_coefficient()[1],
            "Aromaticity": analysis.aromaticity(),
            "Instability_Index": analysis.instability_index(),
            "Gravy": analysis.gravy()
        })

    return pd.DataFrame(results)
```

## Compare Sequences
```python
from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np

def calculate_pairwise_similarity(data):
    rows = []

    if isinstance(data, pd.DataFrame):
        source_data = data.to_dict(orient="records")
    else:
        source_data = data

    for item in source_data:
        if isinstance(item, SeqRecord):
            seq = str(item.seq).strip()
            if not seq:
                continue
            rows.append({
                "ID": item.id,
                "Sequence": seq
            })

        elif isinstance(item, dict):
            seq = str(item.get("Sequence", "")).strip()
            if not seq:
                continue
            rows.append({
                "ID": item.get("ID", "Unknown"),
                "Sequence": seq
            })

    ids = [row["ID"] for row in rows]
    sequences = [row["Sequence"] for row in rows]

    similarity_matrix = pd.DataFrame(
        np.zeros((len(ids), len(ids))),
        index=ids,
        columns=ids
    )

    for i in range(len(sequences)):
        for j in range(len(sequences)):
            seq1 = sequences[i]
            seq2 = sequences[j]

            alignment = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]
            matches = alignment.score
            max_len = max(len(seq1), len(seq2))
            similarity = (matches / max_len) * 100

            similarity_matrix.iloc[i, j] = similarity

    return similarity_matrix
```

## Heatmap

```python
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="white", font_scale=0.95)

plt.figure(figsize=(12, 10))

ax = sns.heatmap(
    similarity_df,
    annot=True,
    fmt=".1f",
    cmap="YlGnBu",
    linewidths=0.5,
    linecolor="white",
    square=True,
    cbar_kws={"label": "Similarity (%)"},
    annot_kws={"size": 9}
)

plt.title("Pairwise Sequence Similarity Heatmap", fontsize=16, weight="bold", pad=14)
plt.xlabel("Protein ID", fontsize=12, weight="bold")
plt.ylabel("Protein ID", fontsize=12, weight="bold")
plt.xticks(rotation=45, ha="right", fontsize=9)
plt.yticks(rotation=0, fontsize=9)
plt.tight_layout()
plt.show()
```

## Physicochemical Dendrogram

```python
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

physchem_data = properties_df[[
    "Length",
    "Molecular_Weight",
    "Isoelectric_Point",
    "Extinction_Coefficient_Reduced",
    "Extinction_Coefficient_Oxidized",
    "Aromaticity",
    "Instability_Index",
    "Gravy"
]]

scaler = StandardScaler()
scaled_physchem = scaler.fit_transform(physchem_data)

linked_physchem = linkage(scaled_physchem, method="ward")

fig, ax = plt.subplots(figsize=(13, 7), facecolor="white")

dendrogram(
    linked_physchem,
    labels=properties_df["ID"].tolist(),
    leaf_rotation=90,
    leaf_font_size=12,
    ax=ax
)

ax.set_title("Dendrogram Based on Physicochemical Properties")
plt.tight_layout()
plt.show()
```

## Sequence Dendrogram

```python
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt

distance_df = 100 - similarity_df.copy()
np.fill_diagonal(distance_df.values, 0)

condensed_distance = squareform(distance_df.values)
linked_seq = linkage(condensed_distance, method="average")

fig, ax = plt.subplots(figsize=(13, 7), facecolor="white")

dendrogram(
    linked_seq,
    labels=similarity_df.index.tolist(),
    leaf_rotation=90,
    leaf_font_size=12,
    ax=ax
)

ax.set_title("Dendrogram Based on Sequence Similarity")
plt.tight_layout()
plt.show()
```
## Combined Dendrogram
```python
from sklearn.preprocessing import StandardScaler
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

physchem_data = properties_df[[
    "Length",
    "Molecular_Weight",
    "Isoelectric_Point",
    "Extinction_Coefficient_Reduced",
    "Extinction_Coefficient_Oxidized",
    "Aromaticity",
    "Instability_Index",
    "Gravy"
]]

scaler = StandardScaler()
scaled_physchem = scaler.fit_transform(physchem_data)

physchem_dist = squareform(pdist(scaled_physchem, metric="euclidean"))

seq_dist = 100 - similarity_df.values
np.fill_diagonal(seq_dist, 0)

physchem_dist_norm = physchem_dist / physchem_dist.max()
seq_dist_norm = seq_dist / seq_dist.max()

combined_dist = 0.5 * physchem_dist_norm + 0.5 * seq_dist_norm

linked_combined = linkage(squareform(combined_dist), method="average")

fig, ax = plt.subplots(figsize=(13, 7), facecolor="white")

dendrogram(
    linked_combined,
    labels=properties_df["ID"].tolist(),
    leaf_rotation=90,
    leaf_font_size=12,
    ax=ax
)

ax.set_title("Dendrogram Based on Combined Information")
plt.tight_layout()
plt.show()
```

---
title: Code
nav_order: 8
---

# Code

Main tasks in the code:

- read FASTA data
- fetch sequences
- calculate properties
- compare sequences
- plot results

## Example

```python
from Bio import SeqIO
import pandas as pd

records = list(SeqIO.parse("sample.fasta", "fasta"))

protein_df = pd.DataFrame([
    {
        "ID": r.id,
        "Description": r.description,
        "Sequence": str(r.seq),
        "Length": len(r.seq)
    }
    for r in records
])

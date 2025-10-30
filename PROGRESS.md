# BIO214 Project 2 - Implementation Progress

**Repository:** https://github.com/saijaya/bio214-project2  
**Due Date:** 10/31/25 at 5pm  
**Files to Submit:** `ppi.py`, `gsea.py` (and optional `node2vec_pretrained`)

---

## 🎯 Project Overview

Analyzing endometriosis gene expression data using:
1. **Part 1 (PPI):** Node2Vec embeddings to predict disease-related genes
2. **Part 2 (GSEA):** Gene Set Enrichment Analysis to identify enriched pathways

---

## ✅ Part 1 (PPI) - COMPLETE

### **Implementation Status:**
- ✅ `load_data()` - Reads disease genes and interaction network (1,476 nodes, 142K edges)
- ✅ `calculate_embedding()` - Node2Vec with exact parameters (seed=42, dimensions=64)
- ✅ `get_close_genes()` - Batched cosine distance computation (efficient!)
- ✅ Pre-trained model saved for autograder speed

### **Key Concepts:**
- **Graph** = discrete representation of gene associations
- **Random walks** = "sentences" generated from graph structure  
- **Word2Vec** = learns embeddings from those "sentences"
- **Node2Vec** = wrapper that does walks + Word2Vec
- **Embeddings** = continuous vectors enabling fast similarity search

### **Critical Parameters:**
```python
Node2Vec: dimensions=64, walk_length=30, num_walks=100, workers=1, seed=42
Training: window=3, min_count=1, batch_words=4
```

### **Data Files:**
- `data/disease_gene_list.txt` - 8 disease genes
- `data/interaction_network_filtered.txt` - 142,406 interactions, 1,476 genes

### **Status:** ✅ Fully implemented, tested, and ready for submission!

---

## ⏳ Part 2 (GSEA) - IN PROGRESS

### **Current Progress:**
- ✅ `load_data()` - Implemented (reads 3 files)
- ✅ `get_gene_rank_order()` - Implemented (ranks by logFC)
- ⏳ `get_enrichment_score()` - **NEXT TO IMPLEMENT**
- ⏳ `get_sig_sets()` - After that

### **Data Files:**
- `data/GSE25628_filtered_expression.txt` - 10,702 genes × 13 samples (log-normalized)
- `data/GSE25628_samples.txt` - Sample labels (0 = healthy [6], 1 = endometriosis [7])
- `data/c2.cp.kegg.v6.2.symbols.filtered.gmt` - 145 KEGG pathways

---

## 📊 GSEA Algorithm Flowchart

```
STEP 1: Calculate actual enrichment scores
├─ Rank genes by logFC (disease - healthy)
│  └─ logFC = mean_disease - mean_healthy (already logged!)
└─ For each pathway:
    └─ Brownian bridge walk → ES

STEP 2: Permutation testing (100 times)
├─ For i in 1 to 100:
│   ├─ Shuffle labels (randomly assign samples to healthy/disease)
│   ├─ Recalculate logFC with shuffled labels
│   ├─ Re-rank genes
│   └─ For each pathway:
│       └─ Brownian bridge walk → permuted_ES[i]

STEP 3: Calculate p-values
└─ For each pathway:
    ├─ Count: how many permuted_ES >= actual_ES?
    └─ p_value = count / 100

STEP 4: Bonferroni correction
└─ corrected_p = p_value × 145 (number of pathways)

STEP 5: Filter significant
└─ Return pathways where corrected_p < 0.05
```

---

## 🎯 Brownian Bridge Walk Algorithm

**For implementing `get_enrichment_score(geneset)`:**

```python
# Step 1: Filter pathway genes to those in expression data
pathway_genes_in_data = set(ranked_genes) & set(pathway_genes)

# Step 2: Calculate scores
n_in_set = len(pathway_genes_in_data)
n_not_in_set = len(ranked_genes) - n_in_set

up_score = 1.0 / n_in_set        # Reward for pathway gene
down_score = 1.0 / n_not_in_set  # Penalty for non-pathway gene

# Step 3: Walk through ranked gene list
running_sum = 0
max_sum = 0

for gene in ranked_genes:
    if gene in pathway_genes_in_data:
        running_sum += up_score
    else:
        running_sum -= down_score
    
    max_sum = max(max_sum, running_sum)

# Step 4: Return ES (supremum of walk)
ES = round(max_sum, 2)  # Round to 2 decimals!
```

---

## 🔑 Implementation Details

### **Instance Variables:**
```python
self.human_gene_expression_df   # DataFrame: genes × samples
self.human_condition_dict       # Dict: {sample_id: 0 or 1}
self.geneset_definitions_dict   # Dict: {pathway_name: [genes]}
self.ranked_genes               # List of gene names (cached)
```

### **LogFC Calculation (Implemented):**
```python
# Separate samples by condition
healthy_samples = [s for s, c in self.human_condition_dict.items() if c == 0]
disease_samples = [s for s, c in self.human_condition_dict.items() if c == 1]

# Calculate means
mean_healthy = self.human_gene_expression_df[healthy_samples].mean(axis=1)
mean_disease = self.human_gene_expression_df[disease_samples].mean(axis=1)

# LogFC (NO additional log - already log-normalized!)
logFC = mean_disease - mean_healthy

# Sort descending
ranked_genes = logFC.sort_values(ascending=False).index.tolist()
```

---

## ⚠️ Critical Requirements

### **Efficiency:**
- ❌ Avoid: `for-loops`, `iterrows()`, `append()`
- ✅ Use: vectorized pandas/numpy operations
- **Target:** `get_sig_sets()` must complete in **< 3 minutes**

### **Exact Parameters:**
- **100 permutations** (not 500!)
- **Bonferroni correction** = p_value × 145 (number of pathways)
- **Round ES to 2 decimal places**

### **Return Types:**
```python
load_data(exp, samp, kegg) → None

get_gene_rank_order() → list[str]
# List of gene names sorted by logFC (highest first)

get_enrichment_score(geneset) → float
# Enrichment score rounded to 2 decimals

get_sig_sets(p) → list[str]
# List of significant pathway names, or [] if none
```

### **Command Line Usage:**
```bash
python3 ppi.py diseaseGeneFile interactionNetworkFile
python3 gsea.py expfile sampfile keggfile
```

---

## 📋 Next Steps

### **1. Implement `get_enrichment_score(geneset)` ⏳**
- Get ranked gene list (call `get_gene_rank_order()`)
- Filter pathway genes to those in expression data
- Calculate up_score and down_score
- Perform Brownian bridge walk
- Return supremum (max) of walk, rounded to 2 decimals

**Test:**
```python
es = gsea.get_enrichment_score("KEGG_CITRATE_CYCLE_TCA_CYCLE")
print(f"ES: {es}")  # Should be a float like 4.52
```

### **2. Implement `get_sig_sets(p_threshold)` ⏳**
- Calculate actual ES for all 145 pathways
- Run 100 permutations:
  - Shuffle sample labels
  - Recalculate rankings
  - Recalculate ES for all pathways
- Calculate p-values (count >= actual / 100)
- Apply Bonferroni correction (p × 145)
- Return pathways where corrected_p < p_threshold

**Test:**
```python
import time
start = time.time()
sig_sets = gsea.get_sig_sets(0.05)
elapsed = time.time() - start
print(f"Found {len(sig_sets)} significant pathways in {elapsed:.1f}s")
# Must be < 180 seconds (3 minutes)!
```

### **3. Final Testing & Submission**
- Test edge cases (no significant sets, extreme thresholds)
- Verify all return types match requirements
- Update header comments (name, date, collaborators)
- Submit to Gradescope: `ppi.py`, `gsea.py`, `node2vec_pretrained`

---

## 🚀 Environment Setup

**Activate virtual environment:**
```bash
cd /Users/krshnajaya/Desktop/bio214-project2
source bio214-project2/bin/activate
```

**Run scripts:**
```bash
# Part 1
python3 src/ppi.py data/disease_gene_list.txt data/interaction_network_filtered.txt

# Part 2
python3 src/gsea.py data/GSE25628_filtered_expression.txt data/GSE25628_samples.txt data/c2.cp.kegg.v6.2.symbols.filtered.gmt
```

---

## 📝 Git Status

**Latest Commits:**
- `521ab0d` - Implement GSEA load_data and get_gene_rank_order methods
- `060826e` - WIP: Start GSEA implementation Part 2
- `8121113` - Complete PPI implementation with Node2Vec embeddings

**Uncommitted:**
- `src/helper_functions.py` - Helper functions for reference (kept local)

---

## 🎓 Working Approach

- ✅ Understand concepts before coding
- ✅ Implement incrementally, test each method
- ✅ Use vectorized operations for efficiency
- ✅ Test locally before using autograder
- ✅ Focus on conceptual clarity over code generation

---

**Last Updated:** In progress  
**Status:** Part 1 complete, Part 2 ~50% done


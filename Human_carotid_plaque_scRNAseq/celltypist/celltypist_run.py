import scanpy as sc
import celltypist
from celltypist import models
import numpy as np
from matplotlib import pyplot as plt

RNA = sc.read('/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/carotid_plaque_without_PBMCs/Whole_dataset_annotated_human_RNA.h5ad')
integrated = sc.read('/gpfs/data/giannarellilab/Mike_G/Immune_checkpoint_inhibitors/carotid_plaque_without_PBMCs/Whole_dataset_annotated_human_integrated.h5ad')
sc.pp.pca(integrated)
# Transfer integrated PCA for clustering used by majority voting. X_pca is used by default
RNA.obsm['X_pca'] = integrated.obsm['X_pca']
del(integrated)

model = models.Model.load(model = 'Immune_All_Low.pkl')

# save normalized counts in raw slot.
#RNA = RNA.raw

# normalize to depth 10 000
sc.pp.normalize_total(RNA, target_sum=1e4)

# logaritmize
sc.pp.log1p(RNA)

# Predict
predictions = celltypist.annotate(RNA, model = 'Immune_All_Low.pkl', majority_voting = True)
predictions.predicted_labels.to_csv('celltypist_labels.csv')

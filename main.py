
import pyranges1 as pr
import pyfaidx


ann = pr.example_data.ncbi_gff
print("original data:\n", ann)
print("original columns:\n",ann.columns)

ann = ann.get_with_loc_columns(['Feature', 'Parent', 'ID'])
print("3 selected columns and genome locations:\n", ann)

# ann = ann.get(['Feature', 'Parent', 'ID'])
# print("3 selected columns:\n", ann)

# fasta file
genome_file = pr.example_data.files['ncbi.fasta']
with open(genome_file) as fh:
  for _ in range(8):
    print(fh.readline().strip())

print("column feature:\n", ann.Feature)

# Create new column
ann['midpoint'] = (ann.Start + ann.End) // 2
print("new column:\n", ann.get_with_loc_columns(['midpoint']))

# filter rows and create a new PyRanges object called cds:
selector = (ann.Feature == 'CDS')
cds = ann [selector]
print("cds1:\n", cds)
cds = cds.drop( ['Feature', 'Parent', 'midpoint'], axis=1 )
print("cds2:\n", cds)

# selecting intervals
reg = cds.loci['CAJFCJ010000097.1', '+', 50000:55000]
print(reg)

# code below derives the first codon of each CDS group; grouping is defined by their ID:
first = cds.slice_ranges(start=0, end=3, group_by='ID')
print("first codons:\n", first)

print("cds head:\n", cds.head())

# find promotor
g = cds.extend_ranges(ext_5=300, group_by='ID')
print(g.head())

prom = g.slice_ranges(0, 300, group_by='ID')
print("promotor:\n", prom.head())

cds.upstream(length=300, group_by='ID').head()

pyf=pyfaidx.Fasta(genome_file)
cor_prom = prom.clip_ranges(chromsizes=pyf)
cor_prom.head()

select_diff_end = cor_prom.End != prom.End
print(prom[select_diff_end])

print(cor_prom[select_diff_end])

# Detecting overlaps
print("overlaps:\n", cor_prom.overlap(cds))

j = cor_prom.join_overlaps(cds)
print("both objects:\n", j)

clu_cds = cds.cluster_overlaps()
clu_cds
c = clu_cds.Cluster.value_counts()
multi_clusters = c[ c > 1 ].index
multi_clu_cds = clu_cds[ clu_cds.Cluster.isin(multi_clusters) ].copy()
print(multi_clu_cds)

# SORTING
print("sorting:\n", multi_clu_cds.sort_ranges())


import pyranges1 as pr
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
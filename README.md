# gene_network
Workflow for gene co-expression and regulatory network 
## Gene Regulatory Network 
Download and make grn.tab file
```bash
# download gene regulatory network from Sonawane, A.R. et al
wget https://www.cell.com/cms/10.1016/j.celrep.2017.10.001/attachment/e7309c03-e579-4119-a95e-376ab2066cbb/mmc2.csv
# reduce file to specific tissue (in this case Uterus)
(head -n 1 mmc2.csv; grep 'Uterus' mmc2.csv) > grn.csv
cat mmc2.csv | sed 's/\"//g' | sed 's/,/\t/4;s/,/\t/3;s/,/\t/1;s/,/\t/1' > grn.tab #Convert to tab-delimited format
```
## Gene Co-Expression Network

You will build the gene coexpression file by:
1. manually by downloading the supplemental data from Hickmann A et. al. At the following link: https://gsajournals.figshare.com/ndownloader/files/31186971
2. Manually parse out and place in a tab-delimited file. Add column headers `GeneA` and `GeneB`.


## Download Biomart Mapping Table
Download geneID mapping table from [https://www.ensembl.org/](https://www.ensembl.org/) 
Make sure to include these headers:
- Gene stable ID
- Gene stable ID version
- Transcript stable ID
- Transcript stable ID version
- Gene name
##  Create Database with Sqlite
Prepare both files for database loading. 
Ensure both `gcn.tab` and `grn.tab` are tsv or tab files 

```bash
module load sqlite
sqlite3 mapping
.separator "\t"
.import grn.tab grn
.separator "\t"
.import gcn.tab gcn
.mode csv
.import mart_export.txt names
.separator "\t"
.headers on
.output GRN-REMAPPED.tsv
SELECT TF, names."Gene Name", Tissues, TargetGene
FROM grn
INNER JOIN names on names."Gene stable ID"=grn.TargetGene
WHERE Tissues LIKE '%Uterus%'; 
.quit
```

## Merge GRN and GCN

```bash
# Process GRN file - ensure tab separation
cat GRN-REMAPPED.tsv | awk -F'\t' 'BEGIN {OFS="\t"} {print $1,$2}' > temp
awk 'BEGIN {OFS="\t"} {print "GRN", $0}' temp > temp2
sed '1s/GRN/NetworkType/g' temp2 > temp3
sed '1s/"Gene/GeneTarget/g' temp3 > GRN_edges.tab
rm temp temp2 temp3
# Process GCN file - convert from comma to tab if needed
cat gcn.tab | awk -F',' 'BEGIN {OFS="\t"} {print $1,$2}' > temp
awk 'BEGIN {OFS="\t"} {print "GCN", $0}' temp > temp2
sed '1s/GCN/NetworkType/g' temp2 > GCN_edges.tab
rm temp temp2
# concate the files
(tail -n +2 GRN_edges.tab; tail -n +2 GRN_edges.tab; tail -n +2 GCN_edges.tab) > merged.gcn.grn.tab
# remove duplicates
cat merged.gcn.grn.tab | uniq | sed 's/\s/\t/g' > unique.merged.gcn.grn.tab
```

Now you have `unique.merged.gcn.grn.tab` which can be used to create networks with Cytoscape or NetworkX python library.

# overlap_network
This code will create an overlap of the top 200 co-expressed genes.This code runs on the `gcn.tab` and remove the first column. And uses the `GRN-REMAPPED.tsv` file as the grn.
Example usage:
```
chmod +x overlap-network.py
python3 overlap-network.py --grn GRN-REMAPPED.tsv --gcn v2_gcn.tab --output overlap --visualize
```
![overlap_overlap_only_visualization (1)](https://github.com/user-attachments/assets/ec3080ce-c023-46f9-9d4e-003c8ca26326)


# ActiveDriverWGS

ActiveDriverWGS is a driver discovery tool for analysis of whole genome sequencing data. It analyzes the mutational burden of SNVs and short INDELs in functionally defined regions of interest and retrieves regions which are significantly mutated compared to a background window. 

For more information, please refer to the ActiveDriverWGS publication. https://www.biorxiv.org/content/early/2017/12/19/236802

### Data Input
##### NOTE: Currently ActiveDriverWGS only works on regions/mutations from GRCh37. Stay tuned for an update.
ActiveDriverWGS requires three input files:

#### 1. Mutations (all_mutations_WGS.txt)
6 Column Tab Separated File (Chromosome, Start, End, Reference, Alternate, Patient ID)

```
chr1 1485768 1485768 G A Patient_1
chr1 3092744 3092744 T C Patient_1
chr1 4352513 4352513 A G Patient_1
chr1 4415835 4415835 C T Patient_1
chr1 4551961 4551961 T A Patient_1
```

#### 2. Regions of Interest (regions.bed)
Regions of interest can be coding or noncoding should be in a BED12 format

```
chr1 69090 70008 gc19_pc.cds::gencode::OR4F5::ENSG00000186092.4 0 + 69090 70008 0 1 918, 0,
chr1 138529 139309 gc19_pc.cds::gencode::AL627309.1::ENSG00000237683.5 0 - 138529 139309 0 1 780, 0,
chr1 367658 368597 gc19_pc.cds::gencode::OR4F29::ENSG00000235249.1 0 + 367658 368597 0 1 939, 0,
chr1 621095 622034 gc19_pc.cds::gencode::OR4F16::ENSG00000185097.2 0 - 621095 622034 0 1 939, 0,
chr1 738531 739137 gc19_pc.cds::gencode::AL669831.1::ENSG00000269831.1 0 - 738531 739137 0 3 87,25,17, 0,256,589

```

#### 3. Sites of Interest (tf_bindingsites.bed)
Sites of interest (e.g. Transcription factor binding sites) which may be active sites in regions of interest should be in a BED4 format

```
chr1 24810 24820 M0082_1.02:TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E
chr1 132482 132492 M0082_1.02:TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E
chr1 228608 228618 M0082_1.02:TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E
chr1 228635 228645 M0082_1.02:TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E
chr1 331068 331078 M0082_1.02:TFAP2A;TFAP2B;TFAP2C;TFAP2D;TFAP2E
```

### Sample Code
```
# prepare_elements
prepare_elements_regions_tf <- prepare_elements("regions.bed", 
                                                "tf_bindingsites.bed", 
                                                 window_size = 50000)

# find_drivers
results <- find_drivers(prepare_elements_regions_tf, 
                        "all_mutations_WGS.txt")

# write results to a table
write.table(results_, 
            “results.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)
```

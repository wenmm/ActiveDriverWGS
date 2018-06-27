# ActiveDriverWGS

ActiveDriverWGS is a driver discovery tool for analysis of whole genome sequencing data. It analyzes the mutational burden of SNVs and short INDELs in functionally defined regions of interest and retrieves regions which are significantly mutated compared to a background window. 

For more information, please refer to the ActiveDriverWGS publication. https://www.biorxiv.org/content/early/2017/12/19/236802

### Data Input
##### NOTE: Currently ActiveDriverWGS only works on regions/mutations from GRCh37. Stay tuned for an update.
ActiveDriverWGS requires three input files:

#### 1. Mutations (all_mutations_WGS.txt)
6 Column Tab Separated File (Chromosome, Start, End, Reference, Alternate, Patient ID)

```
chr17	7577058	7577058	C	A	patient_1
chr17	7577094	7577094	G	A	patient_2
chr17	7578239	7578239	C	A	patient_3
chr17	7579312	7579312	C	G	patient_4
chr17	7577127	7577127	C	G	patient_5
chr17	7578190	7578190	T	C	patient_6
chr17	7578534	7578534	C	A	patient_7
chr17	7577106	7577106	G	T	patient_8
chr17	7577022	7577022	G	A	patient_9
chr17	7578265	7578265	A	G	patient_10
chr17	7577114	7577114	C	T	patient_11
chr17	7577120	7577120	C	T	patient_12
chr17	7578457	7578457	C	G	patient_13
chr17	7577022	7577022	G	A	patient_14
chr17	7577508	7577508	T	C	patient_15
chr17	7578403	7578403	C	A	patient_16
chr17	7577106	7577106	G	A	patient_17
chr17	7578507	7578507	G	C	patient_18
chr17	7578404	7578404	A	T	patient_19
chr17	7579414	7579414	C	T	patient_20
chr17	7577018	7577018	C	T	patient_21
chr17	7577538	7577538	C	T	patient_22
chr17	7578461	7578461	C	A	patient_23
chr17	7577538	7577538	C	T	patient_24
chr17	7577545	7577545	T	C	patient_25
chr17	7578208	7578208	T	C	patient_26
chr17	7579312	7579312	C	G	patient_27
chr17	7577538	7577538	C	T	patient_28
chr17	7577551	7577551	C	T	patient_29
chr17	7578271	7578271	T	A	patient_30
chr17	7577124	7577124	C	T	patient_31
chr17	7577538	7577538	C	T	patient_32
chr17	7577548	7577548	C	T	patient_33
chr17	7579311	7579311	C	A	patient_34
chr17	7577536	7577536	T	C	patient_35
chr17	7578366	7578366	C	G	patient_36
chr17	7578177	7578177	C	A	patient_37
chr17	7578479	7578479	G	A	patient_38
chr17	7578555	7578555	C	T	patient_39
chr17	7578406	7578406	C	T	patient_40
chr17	7567709	7567709	C	G	patient_41
chr17	7578389	7578389	G	A	patient_42
chr17	7578212	7578212	G	A	patient_43
chr17	7577097	7577097	C	T	patient_44
chr17	7577570	7577570	C	T	patient_45
chr17	7578449	7578449	C	T	patient_46
chr17	7579377	7579377	G	A	patient_47
chr17	7574003	7574003	G	A	patient_48
chr17	7578461	7578461	C	A	patient_49
chr17	7577530	7577530	T	G	patient_50
chr17	7578467	7578467	T	G	patient_51
chr17	7578535	7578535	T	C	patient_52
chr17	7577094	7577094	G	A	patient_53
chr17	7577120	7577120	C	T	patient_54
chr17	7574035	7574035	T	A	patient_55
chr17	7577610	7577610	T	C	patient_56
chr17	7578416	7578416	C	A	patient_57
chr17	7577094	7577094	G	A	patient_58
chr17	7577094	7577094	G	A	patient_59
chr17	7578551	7578551	A	G	patient_60
chr17	7577120	7577120	C	T	patient_61
chr17	7578394	7578394	T	C	patient_62
chr17	7578236	7578236	A	G	patient_63
chr17	7578205	7578205	C	T	patient_64
chr17	7577123	7577123	A	C	patient_65
chr17	7578407	7578407	G	C	patient_66
chr17	7578190	7578190	T	C	patient_67
chr17	7577559	7577559	G	T	patient_68
chr17	7577105	7577105	G	C	patient_69
chr17	7577124	7577124	C	A	patient_70
chr17	7577609	7577609	C	A	patient_71
chr17	7578536	7578536	T	C	patient_72
chr17	7577556	7577556	C	T	patient_73
chr17	7577108	7577108	C	A	patient_74
chr17	7577082	7577082	C	T	patient_75
chr17	7578524	7578524	G	C	patient_76
chr17	7569712	7569712	G	T	patient_77
chr17	7577117	7577117	A	G	patient_78
chr17	7577548	7577548	C	T	patient_79
chr17	7579485	7579485	C	A	patient_80
chr17	7577120	7577120	C	T	patient_81
chr17	7574003	7574003	G	A	patient_82
chr17	7577539	7577539	G	A	patient_83
chr17	7578442	7578442	T	C	patient_84
chr17	7578395	7578395	G	A	patient_85
chr17	7578457	7578457	C	A	patient_86
chr17	7578535	7578535	T	C	patient_87
chr17	7574003	7574003	G	A	patient_88
chr17	7578271	7578271	T	C	patient_89
chr17	7578275	7578275	G	A	patient_90
chr17	7577120	7577120	C	T	patient_91
chr17	7577539	7577539	G	A	patient_92
chr17	7576855	7576855	G	A	patient_93
chr17	7576855	7576855	G	A	patient_94
chr17	7579369	7579369	G	C	patient_95
chr17	7577559	7577559	G	C	patient_96
chr17	7578527	7578527	A	G	patient_97
chr17	7577559	7577559	G	C	patient_98
chr17	7574035	7574035	T	G	patient_99
chr17	7578190	7578190	T	C	patient_100
chr17	7577085	7577085	C	T	patient_101
chr17	7576852	7576852	C	T	patient_102
chr17	7578406	7578406	C	T	patient_103
chr17	7577570	7577570	C	T	patient_104
chr17	7578479	7578479	G	A	patient_105
chr17	7574018	7574018	G	A	patient_106
chr17	7577114	7577114	C	T	patient_107
chr17	7578449	7578449	C	T	patient_108
chr17	7578275	7578275	G	A	patient_109
chr17	7578263	7578263	G	A	patient_110
chr17	7578406	7578406	C	T	patient_111
chr17	7577094	7577094	G	A	patient_112
chr17	7578534	7578534	C	G	patient_113
chr17	7578442	7578442	T	C	patient_114
chr17	7577124	7577124	C	T	patient_115
chr17	7578406	7578406	C	T	patient_116
chr17	7577117	7577117	A	C	patient_117
chr17	7578514	7578514	T	-	patient_118
chr17	7577558	7577558	G	-	patient_119
chr17	7577087	7577087	G	-	patient_120
chr17	7577153	7577153	C	-	patient_121
chr17	7578397	7578398	-	G	patient_122
chr17	7578169	7578183	AACCAGACCTCAGGC	-	patient_123
chr17	7577017	7577018	AC	TA	patient_124
chr17	7577594	7577595	AC	-	patient_125
chr17	7579476	7579477	GA	-	patient_126
chr17	7579315	7579316	-	C	patient_127
chr17	7579420	7579420	G	-	patient_128
chr17	7578409	7578437	CTCACAACCTCCGTCATGTGCTGTGACTG	-	patient_129
chr17	7574004	7574005	-	A	patient_130
chr17	7573981	7573982	-	CCAAGG	patient_131
chr17	7578384	7578401	GCAGCGCTCATGGTGGGG	-	patient_132
chr17	7572933	7572934	-	G	patient_133
chr17	7577036	7577036	G	-	patient_134
chr17	7579401	7579401	A	-	patient_135
chr17	7578232	7578235	AAAT	-	patient_136
chr17	7578298	7578308	AATCAGTGAGG	-	patient_137
chr17	7577145	7577146	-	T	patient_138
chr17	7578384	7578401	GCAGCGCTCATGGTGGGG	-	patient_139
chr17	7576872	7576872	C	-	patient_140
chr17	7576904	7576904	G	-	patient_141
chr17	7573220	7573223	TTTA	-	patient_142
chr17	7577605	7577606	-	AACCT	patient_143
chr17	7579470	7579470	C	-	patient_144
chr17	7621497	7621497	T	C	patient_145
chr17	7540272	7540272	G	T	patient_146
chr17	7558658	7558658	G	A	patient_147
chr17	7529893	7529893	C	T	patient_148
chr17	7542958	7542958	C	A	patient_149
chr17	7615333	7615333	C	A	patient_150
chr17	7522960	7522960	G	T	patient_151
chr17	7542032	7542032	C	T	patient_152
chr17	7549794	7549794	T	C	patient_153
chr17	7519200	7519200	A	T	patient_154
chr17	7553740	7553740	G	C	patient_155
chr17	7596010	7596010	C	T	patient_156
chr17	7524594	7524594	C	T	patient_157
chr17	7552067	7552067	T	C	patient_158
chr17	7556256	7556256	C	T	patient_159
chr17	7545359	7545359	G	T	patient_160
chr17	7593025	7593025	G	A	patient_161
chr17	7545803	7545803	C	T	patient_162
chr17	7592703	7592703	G	A	patient_163
chr17	7543742	7543742	C	G	patient_164
chr17	7608922	7608922	G	C	patient_165
chr17	7527414	7527414	C	T	patient_166
chr17	7547095	7547095	G	A	patient_167
chr17	7600411	7600411	G	T	patient_168
chr17	7587152	7587152	C	G	patient_169
chr17	7615385	7615385	C	G	patient_170
chr17	7556843	7556843	C	A	patient_171
chr17	7542899	7542900	-	T	patient_172
chr17	7611545	7611546	CG	AT	patient_173
chr17	7517223	7517223	A	-	patient_174
chr17	7590689	7590692	GCTT	-	patient_175
```

#### 2. Regions of Interest (regions.bed)
Regions of interest can be coding or noncoding should be in a BED12 format

```
chr17	7565256	7579912	gc19_pc.cds::gencode::TP53::ENSG00000141510.11	0	-	7565256	7579912	0	14	76,39,82,107,48,33,74,137,110,113,184,279,22,74,	0,4267,7670,8670,11280,11368,11596,11762,12242,12920,13114,14055,14443,14582,

```

#### 3. Sites of Interest (tf_bindingsites.bed)
Sites of interest (e.g. Transcription factor binding sites) which may be active sites in regions of interest should be in a BED4 format

```
chr17	7565257	7565267	site1
chr17	7566256	7566266	site2
chr17	7577087	7577097	site3
chr17	7578244	7578254	site4
chr17	7579901	7579911	site5
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
write.table(results, 
            “results.txt", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE)
```

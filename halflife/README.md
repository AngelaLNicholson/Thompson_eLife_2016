# halflife

1) Plot measured mRNA half-lives vs. translation efficiency:
python halflife_v_TE.py <Miller_halflives> <Presnyak_halflives> <Neymotin_halflives> <mergedAttributes> <outname>

file descriptions:

<Miller_halflives> Miller et al., 2011, text file S1

<Neymotin_halflives> Neymotin et al., Table S5, in csv format

<Presnyak_halflives> Presnyak et al., 2015, Table S1, in csv format

<mergedAttributes> csv file output of collectAtts.py, e.g. ORF_length_comparisons/1411014_filtered_final.csv

<outname> prefix for output file or graph

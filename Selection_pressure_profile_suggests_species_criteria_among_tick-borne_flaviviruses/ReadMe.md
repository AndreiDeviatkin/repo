## Selection pressure profile suggests species criteria among tick-borne flaviviruses

Data and code for "Selection pressure profile suggests species criteria among tick-borne flaviviruses" paper by Andrei Deviatkin, Yulia Aleshina, Galina Karganova, and Alexander Lukashev.

Data:
- `data/codon_based_alignments` - codon-based alignments of ORF, NS3, and E for mosquito-borne flaviviruses (MBFV) and tick-borne flaviviruses (TBFV).
- `data/codon_based_alignments` - amino acid alignments of ORF, NS3, and E for MBFV and TBFV.
- `data/pairwise_dnds_tables`  - tables with pairwise Nei-Gojobori and Yang-Nielsen dN/dS values

Scripts:
- `scripts/distance_visualization_intra_inter.r ` - visualization of uncorrected pairwise distances for nucleotide and amino acid sequences.
- `scripts/computeNeiGojobori.py` - calculation of Nei-Gojobori dN/dS.
- `scripts/dnds_hist.r` - visualization of Nei-Gojobori dN/dS.
- `scripts/parse_yn00.py` - parsing of yn00.exe output to produce a table with pairwise Yang-Nielsen dN/dS values.
- `scripts/dnds_hist_sup.r` - visualization of pairwise Yang-Nielsen dN/dS.
- `scripts/dnds_tree.r` - cladogram based on pairwise Nei-Gojobori dN/dS ratio using neighbour-joining algorithm.


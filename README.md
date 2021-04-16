# BioPredictions
A package for simple estimations and predictions for biomolecules. Included are:
- Isoelectric point
- Molecular weight
- GRAVY index
- Protease cleavage sites

# Usage
```julia
julia> using BioSequences, ProteinStats
julia> seq = aa"MTPSIRQPLRLRRLPATVSPGGETDAMEYRELPQPQPIPSDGLAEAASPNRLLGYLLLHWPMVLILGSMLGAGMAYLAYTLIPAKYTTYAMIRVALVPPSVSGFQNEEAARNDFLTCLKTQTQLIKSHFVLNAAIRDPAIAELPMIRSQVDPVAFLQDEVRVEYTDNS";
julia> weight(seq)
18567.54638
julia> pI(seq)
5.496
julia> cleave(seq)
(positions = [6, 10, 12, 13, 30, 51, 85, 93, 111, 119, 126, 136, 147, 161], peptides = LongAminoAcidSeq[MTPSIR, QPLR, LR, R, LPATVSPGGETDAMEYR, ELPQPQPIPSDGLAEAASPNR, LLGYLLLHWPMVLILGSMLGAGMAYLAYTLIPAK, YTTYAMIR, VALVPPSVSGFQNEEAAR, NDFLTCLK, TQTQLIK, SHFVLNAAIR, DPAIAELPMIR, SQVDPVAFLQDEVR, VEYTDNS])
julia> cleave(seq, BioPredictions.proteinaseK)
(positions = [2, 5, 9, 11, 14, 16, 17, 18, 23, 24  …  153, 154, 155, 156, 159, 160, 162, 163, 164, 165], peptides = LongAminoAcidSeq[MT, PSI, RQPL, RL, RRL, PA, T, V, SPGGE, T  …  A, F, L, QDE, V, RV, E, Y, T, DNS])
julia> gravy(seq)
0.015476190476190442
```

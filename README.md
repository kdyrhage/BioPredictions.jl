# ProteinStats.jl
A package for calculating and estimating simple statistics for biomolecules. Included are:
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
(positions = [6, 10, 12, 13, 30, 51, 85, 93, 111, 119, 126, 136, 147, 161], peptides = LongAminoAcidSeq[MTPSIR, QPLR, LR, R, LPATVSPGGETDAMEYR, ELPQPQPIPSDGLAEAASPNR, LLGYLLLHWPMVLILGSMLGAGMAYLAYTLIPAK, YTTYAMIR, VALVPPSVSGFQNEEAAR, NDFLTCLK, TQTQLIK, SHFVLNAAIR, DPAIAELPMIR, SQVDPVAFLQDEVR, RVEYTDNS])
julia> gravy(seq)
0.015476190476190442
```

# ProteinStats.jl
A package for calculating simple protein statistics. Included are:
- Isoelectric point
- Molecular weight
- GRAVY index
- Protease cleavage sites

# Usage
```julia
using BioSequences, ProteinStats
seq = aa"MTPSIRQPLRLRRLPATVSPGGETDAMEYRELPQPQPIPSDGLAEAASPNRLLGYLLLHWPMVLILGSMLGAGMAYLAYTLIPAKYTTYAMIRVALVPPSVSGFQNEEAARNDFLTCLKTQTQLIKSHFVLNAAIRDPAIAELPMIRSQVDPVAFLQDEVRVEYTDNS"
weight(seq)
pI(seq)
cleave(seq)
gravy(seq)
```

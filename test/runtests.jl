using Test
using ProteinStats
using BioSequences

seq = LongAminoAcidSeq([aa for aa in alphabet(AminoAcid)])

@testset "Gravy index" begin
    @test gravy(seq) == -0.35
end

@testset "Molecular weight" begin
    @test isapprox(weight(seq),
                   weight(seq; ignoreunknown = true) + sum(map(aa -> ProteinStats.aminoacidweights[aa], [AA_X, AA_Z, AA_B, AA_J])))
end

@testset "Isoelectric point" begin
    @test pI(seq) == 6.6
    @test pI(seq, "Solomon") == 7.14
    @test pI(seq, collect(1:9)) == 4.5
end

@testset "Cleavage sites" begin
    @test cleave(seq[1:20]) == (positions = [2, 12], peptides = LongAminoAcidSeq[aa"AR", aa"NDCQEGHILK", aa"MFPSTWYV"])
    @test cleave(seq[1:20], ProteinStats.proteinaseK) == (positions = [1, 7, 10, 11, 14, 17, 18, 19], peptides = LongAminoAcidSeq[aa"A", aa"RNDCQE", aa"GHI", aa"L", aa"KMF", aa"PST", aa"W", aa"Y", aa"V"])
    @test cleave(seq[1:20], ProteinStats.lysC) == (positions = [12], peptides = LongAminoAcidSeq[aa"ARNDCQEGHILK", aa"MFPSTWYV"])
    @test cleave(seq[1:20], ProteinStats.proline_endopeptidase) == (positions = Int[], peptides = LongAminoAcidSeq[seq[1:20]])
end

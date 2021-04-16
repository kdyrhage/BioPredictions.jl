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
    @test cleave(seq) == (positions = [2, 12], peptides = LongAminoAcidSeq[aa"AR", aa"NDCQEGHILK", aa"KMFPSTWYVOUBJZX"])
end

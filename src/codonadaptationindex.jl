function codon_frequencies(chrs)
    s = dna""
    for gene in @genes(chrs, CDS, iscomplete(gene))
        s *= sequence(gene)
    end
    cf = countmap(each_codon(s))
end

function optimal_codons(cf)
    maxf(v) = maximum(map(c -> get(cf, DNACodon(c), 0.0), v))
    oc = Dict(DNACodon("TTT") => maxf(["TTT", "TTC"]), # Phenylalanine
              DNACodon("TTC") => maxf(["TTT", "TTC"]),
              DNACodon("TTA") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]), # Leucine
              DNACodon("TTG") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
              DNACodon("CTT") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
              DNACodon("CTC") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
              DNACodon("CTA") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
              DNACodon("CTG") => maxf(["TTA", "TTG", "CTT", "CTC", "CTA", "CTG"]),
              DNACodon("ATT") => maxf(["ATT", "ATC", "ATA", "ATG"]), # Isoleucine
              DNACodon("ATC") => maxf(["ATT", "ATC", "ATA", "ATG"]),
              DNACodon("ATA") => maxf(["ATT", "ATC", "ATA", "ATG"]),
              DNACodon("ATG") => maxf(["ATT", "ATC", "ATA", "ATG"]),
              DNACodon("GTT") => maxf(["GTT", "GTC", "GTA", "GTG"]), # Valine
              DNACodon("GTC") => maxf(["GTT", "GTC", "GTA", "GTG"]),
              DNACodon("GTA") => maxf(["GTT", "GTC", "GTA", "GTG"]),
              DNACodon("GTG") => maxf(["GTT", "GTC", "GTA", "GTG"]),
              DNACodon("TCT") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]), # Serine
              DNACodon("TCC") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
              DNACodon("TCA") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
              DNACodon("TCG") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
              DNACodon("AGT") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
              DNACodon("AGC") => maxf(["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]),
              DNACodon("CCT") => maxf(["CCT", "CCC", "CCA", "CCG"]), # Proline
              DNACodon("CCC") => maxf(["CCT", "CCC", "CCA", "CCG"]),
              DNACodon("CCA") => maxf(["CCT", "CCC", "CCA", "CCG"]),
              DNACodon("CCG") => maxf(["CCT", "CCC", "CCA", "CCG"]),
              DNACodon("ACT") => maxf(["ACT", "ACC", "ACA", "ACG"]), # Threonine
              DNACodon("ACC") => maxf(["ACT", "ACC", "ACA", "ACG"]),
              DNACodon("ACA") => maxf(["ACT", "ACC", "ACA", "ACG"]),
              DNACodon("ACG") => maxf(["ACT", "ACC", "ACA", "ACG"]),
              DNACodon("GCT") => maxf(["GCT", "GCC", "GCA", "GCG"]), # Alanine
              DNACodon("GCC") => maxf(["GCT", "GCC", "GCA", "GCG"]),
              DNACodon("GCA") => maxf(["GCT", "GCC", "GCA", "GCG"]),
              DNACodon("GCG") => maxf(["GCT", "GCC", "GCA", "GCG"]),
              DNACodon("TAT") => maxf(["TAT", "TAC"]), # Tyrosine
              DNACodon("TAC") => maxf(["TAT", "TAC"]),
              DNACodon("TAA") => maxf(["TAA", "TAG", "TGA"]), # Stop
              DNACodon("TAG") => maxf(["TAA", "TAG", "TGA"]),
              DNACodon("TGA") => maxf(["TAA", "TAG", "TGA"]),
              DNACodon("CAT") => maxf(["CAT", "CAC"]), # Histidine
              DNACodon("CAC") => maxf(["CAT", "CAC"]),
              DNACodon("CAA") => maxf(["CAA", "CAG"]), # Glutamine
              DNACodon("CAG") => maxf(["CAA", "CAG"]),
              DNACodon("AAT") => maxf(["AAT", "AAC"]), # Asparagine
              DNACodon("AAC") => maxf(["AAT", "AAC"]),
              DNACodon("AAA") => maxf(["AAA", "AAG"]), # Lysine
              DNACodon("AAG") => maxf(["AAA", "AAG"]),
              DNACodon("GAT") => maxf(["GAT", "GAC"]), # Aspartic acid
              DNACodon("GAC") => maxf(["GAT", "GAC"]),
              DNACodon("GAA") => maxf(["GAA", "GAG"]), # Glutamic acid
              DNACodon("GAG") => maxf(["GAA", "GAG"]),
              DNACodon("TGT") => maxf(["TGT", "TGC"]), # Cysteine
              DNACodon("TGC") => maxf(["TGT", "TGC"]),
              DNACodon("TGG") => maxf(["TGG"]), # Tryptophan
              DNACodon("CGT") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]), # Arginine
              DNACodon("CGC") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
              DNACodon("CGA") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
              DNACodon("CGG") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
              DNACodon("AGA") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
              DNACodon("AGG") => maxf(["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"]),
              DNACodon("GGT") => maxf(["GGT", "GGC", "GGA", "GGG"]), # Glycine
              DNACodon("GGC") => maxf(["GGT", "GGC", "GGA", "GGG"]),
              DNACodon("GGA") => maxf(["GGT", "GGC", "GGA", "GGG"]),
              DNACodon("GGG") => maxf(["GGT", "GGC", "GGA", "GGG"]))
end

function relative_adaptiveness(codon, cf, oc)
    get(cf, codon, 0.0) / get(oc, codon, 0.0)
end

function cai(seq::LongDNA{N}, cf, oc) where N
    codons = each_codon(seq)
    w = fill(0.0, length(codons))
    for (i, codon) in enumerate(codons)
        w[i] = relative_adaptiveness(codon, cf, oc)
    end
    isempty(w) ? 0.0 : geomean(w)
end
cai(gene, cf, oc) = cai(sequence(gene), cf, oc)

function cai(chrs)
    cf = codon_frequencies(chrs)
    oc = optimal_codons(cf)
    genes = @genes(chrs, CDS, iscomplete(gene))
    result = fill(0.0, length(genes))
    Threads.@threads for (i, gene) in collect(enumerate(genes))
        result[i] = cai(gene, cf, oc)
    end
    result
end


function weightfactors(chrs)
    dict = Dict()
    ngenes = 0
    for gene in @genes(chrs, CDS, iscomplete(gene))
        codons = each_codon(sequence(gene))
        for codon in codons
            get!(dict, codon, 0.0)
            dict[codon] += 1.0
        end
        ngenes += 1
    end
    for (codon, count) in dict
        dict[codon] = count / ngenes
    end
    return dict
end

function gcai(seq::LongDNA{N}, cf, oc, wf) where N
    codons = each_codon(seq)
    w = fill(0.0, length(codons))
    for (i, codon) in enumerate(codons)
        w[i] = relative_adaptiveness(codon, cf, oc)
    end
    isempty(w) ? 0.0 : geomean(w)
end
gcai(gene, cf, oc, wf) = cai(sequence(gene), cf, oc, wf)

function gcai(chrs)
    cf = codon_frequencies(chrs)
    oc = optimal_codons(cf)
    wf = weightfactors(chrs)
    genes = @genes(chrs, CDS, iscomplete(gene))
    result = fill(0.0, length(genes))
    Threads.@threads for (i, gene) in collect(enumerate(genes))
        result[i] = cai(gene, cf, oc, wf)
    end
    result
end
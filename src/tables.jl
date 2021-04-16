"""
`Dict` containing pK values from different sources. Values are ordered [NTerm, CTerm, C, D, E, H, K, R, Y].
Available options are "EMBOSS", "Solomon", "Rodwell", "IPC protein", and "IPC peptide".
Values were collected from http://isoelectric.org/theory.html
"""
const pKtable = Dict{String, Vector{Float64}}(
    "EMBOSS" => [8.6, 3.6, 8.5, 3.9, 4.1, 6.5, 10.8, 12.5, 10.1],
    "Solomon" => [9.6, 2.4, 8.3, 3.9, 4.3, 6.0, 10.5, 12.5, 10.1],
    "Rodwell" => [8.0, 3.1, 8.33, 3.68, 4.25, 6.0, 11.5, 11.5, 10.07],
    "IPC protein" => [9.094, 2.869, 7.555, 3.872, 4.412, 5.637, 9.052, 11.84, 10.85],
    "IPC peptide" => [9.564, 2.383, 8.297, 3.887, 4.317, 6.018, 10.517, 12.503, 10.071]
)


"""
Amino acid weights (source: http://education.expasy.org/student_projects/isotopident/htdocs/aa-list.html)
"""
const aminoacidweights = Dict(
    AA_A => 71.0788,
    AA_R => 156.1875,
    AA_N => 114.1038,
    AA_D => 115.0886,
    AA_C => 103.1388,
    AA_Q => 128.1307,
    AA_E => 129.1155,
    AA_G => 57.0519,
    AA_H => 137.1411,
    AA_I => 113.1594,
    AA_L => 113.1594,
    AA_K => 128.1741,
    AA_M => 131.1926,
    AA_F => 147.1766,
    AA_P => 97.1167,
    AA_S => 87.0782,
    AA_T => 101.1051,
    AA_W => 186.2132,
    AA_Y => 163.1760,
    AA_V => 99.1326,
    AA_Z => 146.6385,
    AA_J => 131.175,
    AA_B => 132.6115,
    AA_U => 168.064,
    AA_O => 255.313,
    AA_X => 118.88603)


"""
Hydropathy indices (source: https://www.ncbi.nlm.nih.gov/pubmed/7108955)
"""
const hydropathy = Dict(
    AA_A => 1.8,
    AA_R => -4.5,
    AA_N => -3.5,
    AA_D => -3.5,
    AA_C => 2.5,
    AA_Q => -3.5,
    AA_E => -3.5,
    AA_G => -0.4,
    AA_H => -3.2,
    AA_I => 4.5,
    AA_L => 3.8,
    AA_K => -3.9,
    AA_M => 1.9,
    AA_F => 2.8,
    AA_P => -1.6,
    AA_S => -0.8,
    AA_T => -0.7,
    AA_W => -0.9,
    AA_Y => -1.3,
    AA_V => 4.2)

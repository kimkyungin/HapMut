local aux = require"auxlib"
local chr2num = aux.chr2num
local make_ref = aux.make_ref
local make_seq = aux.make_seq
local extract_site = aux.extract_site -- zero based genomic position
local extract_read = aux.extract_read -- zero based genomic position
local posterior = aux.posterior
 
local blib = require"blocklib"
local make_hapblock = blib.mkhapblock
local make_blockread = blib.mkblockread
local assign_prob = blib.assignpr

-- local genome = "human_g1k_v37.fasta"
-- local chr = "21"
-- local nbam = "pe_30x_e1p_002h_m00.sorted.bam"
-- local tbam = "pe_30x_e1p_002h_m06.sorted.bam"

genome, chr, nbam, tbam = arg[1], arg[2], arg[3], arg[4]
local chr = chr2num(chr)

local ref = make_ref(genome, chr)
local normal = make_seq(nbam, chr, ref)
local tumor = make_seq(tbam, chr, ref)
local snp, snp_pos = extract_site(tumor.bam, tumor.chr, normal.info, "biallelic")
local hapblock = make_hapblock(snp, snp_pos)
local blockread = make_blockread(snp, hapblock)

local list = {}
local t = normal.somatic
for k, v in pairs(tumor.somatic) do
   if not t[k] then
      list[#list + 1] = k
   end
end
table.sort(list)

local out, out_pos = extract_read(tumor.bam, tumor.chr, list, ref)
assign_prob(out, out_pos, blockread)      

local dna = {A = {}, C = {}, G = {}, T = {}}
dna.A = {"C", "G", "T"}
dna.C = {"G", "T", "A"}
dna.G = {"T", "A", "C"}
dna.T = {"A", "C", "G"}

for i = 1, #out_pos do
   local pos = out_pos[i] 
   local somatic = posterior(out[pos])
   if somatic then
      local ref_allele = string.char(ref[pos])
      io.write(pos .. "\t") -- 0-based genomic position
      io.write(ref_allele .. "\t")
      io.write(dna[ref_allele][somatic.imax] .. "\t")
      io.write(table.concat(somatic.prob, "\t"))
      io.write("\n")
   end
end



-- path = "/Users/arrti/Documents/HapBlock2/C/ucharray.so"
-- uc = package.loadlib(path, "luaopen_ucharray")()

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

-- local genome = "/Users/arrti/Public/human_g1k_v37.fasta"
-- genome = "/Users/arrti/Documents/LongRead/db/human_g1k_v37.fasta"
-- local chr = "21"
-- local nbam = "/Users/arrti/Public/pe_30x_e1p_002h_m00.sorted.bam"
-- local tbam = "/Users/arrti/Public/pe_30x_e1p_002h_m06.sorted.bam"

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

for i = 1, #out_pos do
   local pos = out_pos[i] 
   local t = posterior(out[pos])
   if t then
      io.write(pos .. "\t")
      for _, v in ipairs(t) do
	 io.write(v .. "\t")
      end
      io.write("\n")
   end
end



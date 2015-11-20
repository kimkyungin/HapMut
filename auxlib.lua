local aux = {}

-- package.cpath = package.cpath .. ";./?.so"
local make_array = require"ucharray".new
local scan = require"readscan".scan
local decode = require"readscan".decode

function aux.chr2num (x)
   if x == "X" then chr = 23
   elseif x == "Y" then chr = 24
   else chr = tonumber(x) end
   return chr
end

function aux.num2chr (x)
   if x == 23 then chr = "X"
   elseif x == 24 then chr = "Y"
   else chr = tostring(x) end
   return chr
end

local function file_exists (name)
   local f = io.open(name, "r")
   if f ~= nil then 
      io.close(f)
      return true
   else
      return false
   end
end

function aux.make_ref (genome, chr)
   local genomeinfo = genome:gsub(".*/", "./"):gsub(".fasta", "_chr_size.lua")
   if not file_exists(genomeinfo) then
      print("No " .. genomeinfo .. " found")
      print("Should generate " .. genomeinfo .. " first")
      print("Please run the following command on your terminal:")
      print("   $ lua get-chr-size.lua " .. genome)
      assert(false)
   end
   dofile(genomeinfo)
   local tag = ">" .. aux.num2chr(chr)
   return make_array(chr_size[chr], genome, tag)
end

function aux.make_seq (bam, chr, ref)
   if not file_exists(bam) then
      error("No BAM file")
   end
   local info = make_array(#ref)
   scan(bam, chr, ref, info)
   local s, b, m = decode(info)
   return {bam=bam, chr=chr, ref=ref, info=info, biallelic=b, monoallelic=m, somatic=s}
end

aux.extract_site = require"snpextract".extract
aux.extract_read = require"readextract".extract
aux.posterior = require"mutlib".posterior

return aux

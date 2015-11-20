local genome = arg[1]
local finfo = genome:gsub(".*/", "./"):gsub(".fasta", "_chr_size.lua")
local f = io.open(finfo, "r")
if f ~= nil then
   f:close()
   print(finfo .. " is alread generated")
else
   local fp = io.open(genome, "r")
   local chr = {}
   local line = fp:read("*line")
   local tag = ">"
   while line do
      if line:find(tag) then
	 local count = 0
	 local name = line:match(tag .. "(%S+)%s+")
	 while true do
	    line = fp:read("*line")
	    if not line or line:find(tag) then break end
	    count = count + line:len()
	 end
	 if not name:find("%.") then
	    if name == "X" then chr[23] = count
	    elseif name == "Y" then chr[24] = count
	    elseif name:find("%d+") then chr[tonumber(name)] = count
	    end
	 end
      end      
   end
   fp:close()
   
   local fp= io.open(finfo, "w")
   fp:write("chr_size = {}\n")
   for i, v in ipairs(chr) do
      fp:write("chr_size[" .. tostring(i) .. "] = " .. v .. "\n")
   end
   fp:close()
end

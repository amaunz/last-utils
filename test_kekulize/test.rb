require 'openbabel'

def analysis (m,p)
	p.match(m)
	hits = p.get_umap_list
	print "Found #{hits.size} instances."
	if hits.size>0
	   puts " Here are the atom indices:"
	else
	   print "\n"
	end

	hits.each_with_index do |hit, index|
	   print "  Hit #{index}: [ "
	   hit.each do |atom_index|
	       print "#{atom_index} "
	   end
	   puts "]"
	end
end

c=OpenBabel::OBConversion.new
c.set_in_format 'smi'
m=OpenBabel::OBMol.new
p=OpenBabel::OBSmartsPattern.new

# all bonds aromatic by default
#c.read_string m, "c1ccccc1"
smi="n1(c2ccccc2)n(c(C)cc1=O)C"
c.read_string m, smi
puts "#{smi}: #{m.num_atoms} atoms."

### NEED TO COMMENT OUT FOR LOWER BLOCK TO WORK ###


###### SWITCH TO KEKULE ######
puts "kekulize"
m.kekulize
##############################


puts "[#6]-[#6](=[#6]-[#6]-[#7]-[#6])(-[#7])"
p.init("[#6]-[#6](=[#6]-[#6]-[#7]-[#6])(-[#7])")
analysis(m,p)

puts "[#6]-,:[#6](=,:[#6]-,:[#6]-,:[#7]-,:[#6])(-,:[#7])"
p.init("[#6]-,:[#6](=,:[#6]-,:[#6]-,:[#7]-,:[#6])(-,:[#7])")
analysis(m,p)

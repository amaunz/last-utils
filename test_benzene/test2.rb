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
c.read_string m, "c1ccccc1"

### NEED TO COMMENT OUT FOR LOWER BLOCK TO WORK ###

# matches all bonds
#p "[#6][#6]"
#p.init("[#6][#6]")
#analysis(m,p)

# matches no bonds
#p "[#6]-[#6]"
#p.init("[#6]-[#6]")
#analysis(m,p)

# matches no bonds
#p "[#6]=[#6]"
#p.init("[#6]=[#6]")
#analysis(m,p)


###### SWITCH TO KEKULE ######
p "kekulize and set aromatic"
m.kekulize
m.set_aromatic_perceived
##############################

# matches only single bonds
p "[#6][#6]"
p.init("[#6]-[#6]")
analysis(m,p)


# matches only single bonds
p "[#6]-[#6]"
p.init("[#6]-[#6]")
analysis(m,p)

# matches only double bonds
p "[#6]=[#6]"
p.init("[#6]=[#6]")
analysis(m,p)

p "[#6]-[#6]=[#6]"
p.init("[#6]-[#6]=[#6]")
analysis(m,p)

p "[#6]=[#6]-[#6]"
p.init("[#6]=[#6]-[#6]")
analysis(m,p)

p "[#6]=[#6]-[#6]=[#6]"
p.init("[#6]=[#6]-[#6]=[#6]")
analysis(m,p)

p "[#6]=[#6]-[#6]=[#6]"
p.init("[#6]=[#6]-[#6]=[#6]")
analysis(m,p)


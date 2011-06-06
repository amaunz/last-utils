# Host application for LAST-UTILS, to be called from the command line
# See http://last-pm.maunz.de
#
# Author::    Andreas Maunz (mailto:andreas@maunz.de)

begin 
  require File.dirname(__FILE__) +  "/lu.rb"
rescue
  puts "last-utils.rb: lu.rb not found!"
  exit false
end

# Main
STDOUT.sync = true

status=false
if $*.size<1 || $*.size>4
    status=true
end

lu = LU.new

case $*[0]
when '1' 
    if !status
	no_aromatic=false
	aromatic_wc=false
	if $*.size>=3
		no_aromatic=true unless $*[2]!="nna"
		if $*.size>=4
			aromatic_wc=true unless $*[3]!="wcb"
		end
	end
        dom = lu.read(nil,no_aromatic,aromatic_wc)
        lu.smarts(dom, $*[1])
    end
when '2'
    if !status
        ENV["FMINER_NR_HITS"].nil? ? nr_hits=false : nr_hits=true 
        lu.match_file($*[1], nr_hits)
    end
when '3'
    exit lu.ob_test
else
    status=true
end

if status
    puts "Usage: #{$0} 1 [ msa | nls | nop ] [nna] [wcb] < /path/to/graphmlfile.graphml > /path/to/smartsfile.smarts" 
    puts "       #{$0} 2 /path/to/smifile.smi < /path/to/smartsfile.smarts > /path/to/lastpmfile.lastpm" 
    puts "       #{$0} 3" 
    puts "       cmd=1 : convert GraphML to SMARTS."
    puts "           msa : All level max opt."
    puts "           nls : Next level opt."
    puts "           nop : No opt."
    puts 
    puts "           nna: Disable explicit node annotation for aromatic/aliphatic."
    puts "              Useful when you previously used the -a option in fminer (if not set all atoms will be marked as aromatic (if applicable) or aliphatic)"
    puts "           wcb: Enable aromatic wildcarding on bonds."
    puts 
    puts "       cmd=2 : match SMARTS to SMILES file and create LASTPM file."
    puts "       cmd=3 : test mode."
    exit
end

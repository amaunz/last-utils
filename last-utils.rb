require File.dirname(__FILE__) +  "/lu.rb"

# Main
STDOUT.sync = true


status=false
if $*.size<1 || $*.size>3
    status=true
end

lu = LU.new

case $*[0]
when '1' 
    if !status
	aromatic=true
	if $*.size==3
		aromatic=false unless $*[2]!="a"
	end
        dom = lu.read(nil,aromatic)
        lu.smarts(dom, $*[1])
    end
when '2'
    if !status
        lu.match_file($*[1])
    end
when '3'
    lu.demo
else
    status=true
end

if status
    puts "Usage: #{$0} 1 [ msa | nls | nop ] [a] < /path/to/graphmlfile.graphml > /path/to/smartsfile.smarts" 
    puts "       #{$0} 2 /path/to/smifile.smi < /path/to/smartsfile.smarts > /path/to/lastpmfile.lastpm" 
    puts "       #{$0} 3" 
    puts "       cmd=1 : convert GraphML to SMARTS."
    puts "           msa : All level max opt."
    puts "           nls : Next level opt."
    puts "           nop : No opt."
    puts "           a: Disable explicit node annotation for aromatic/aliphatic."
    puts "              Useful when you previously used the -a option in fminer (or all atoms will be marked by aliphatic '&A')."
    puts "       cmd=2 : match SMARTS to SMILES file and create LASTPM file."
    puts "       cmd=3 : demo mode."
    exit
end

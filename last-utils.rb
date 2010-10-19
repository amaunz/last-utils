#!/usr/bin/env ruby
require 'lu.rb'

# Main
STDOUT.sync = true


status=false
if $*.size==0 || $*.size>2
    status=true
end

lu = LU.new

case $*[0]
when '1' 
    !if !status
        dom = lu.read
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
    puts "Usage: #{$0} 1 [ msa | nls | nop ] < /path/to/graphmlfile.graphml > /path/to/smartsfile.smarts" 
    puts "       #{$0} 2 /path/to/smifile.smi < /path/to/smartsfile.smarts > /path/to/lastpmfile.lastpm" 
    puts "       cmd=1 : convert GraphML to SMARTS."
    puts "           msa : All level max opt."
    puts "           nls : Next level opt."
    puts "           nop : No opt."
    puts "       cmd=2 : match SMARTS to SMILES file and create LASTPM file."
    exit
end

#!/usr/bin/env ruby

require('openbabel')
#include OpenBabel   <= may be used for direct access to OB namespace, i.e. w/o "OpenBabel::". Below, I want namespaces for clarity.

%w[pp rubygems hpricot].each { |x| require x }

class LUGraph
  attr_accessor :nodes, :edges
  def initialize(nodes, edges)
    @nodes = nodes
    @edges = edges
  end

  def to_smarts(backw,f,opt)

      mand_branches=0
      opt_branches = @edges[f].size
      mand_t = []
      mand_e = []
      opt_t = []
      opt_e = []
      @edges[f].each do |t,e|
          if e.opt_e==0 || (e.opt_e==opt)
              mand_branches+=1
              opt_branches-=1
              mand_t << t
              mand_e << e
          else
              opt_t << t
              opt_e << e
          end
      end

      print "["
      @nodes[f].to_smarts
      do_branch=(opt_branches>1)
      if do_branch
          # 1) Uncomment next line to not force at least one branch (c.f. not. 17.11.09)
          # print "$(" ; @nodes[f].to_smarts ; print ")," 
          # 1) end
          print ";"
          for i in 0..opt_t.size-1
              t = opt_t[i]
              e = opt_e[i]
              print "$("
                  print "["                   # self
                  @nodes[f].to_smarts         # 
                  print "]"                   # 

                  print "("                   # branch
                  e.to_smarts                 # 
                  to_smarts(f,t,e.opt_e)      # 
                  print ")"                   # 
                  if backw!=f
                      print "(" unless mand_branches==0   # backw
                      print "["                           # backw
                      @nodes[backw].to_smarts             # 
                      print "]"                           # 
                      print ")" unless mand_branches==0   # 
                  end
                  if mand_branches>0
                      print "["               # forw
                      @nodes[t].to_smarts     # 
                      print "]"               # 
                  end
              print ")"
              print "," unless i==opt_t.size-1
          end
      else
          # 2) Uncomment next block to make single optional edges mandatory (c.f. not. 17.11.09)
          #if opt_branches==1
          #    @edges[f].each do |t,e|
          #        if !mand_t.include?(t)
          #            print "("
          #            e.to_smarts
          #            to_smarts(t,e.opt_e)
          #            print ")"
          #        end
          #    end
          #end
          # 2) end
      end
      print "]"
      print "(~*)" unless !do_branch

      for i in 0..mand_t.size-1
          t = mand_t[i]
          e = mand_e[i]
          print "(" unless i==mand_t.size-1
          e.to_smarts
          to_smarts(f,t,e.opt_e)
          print ")" unless i==mand_t.size-1
      end

  end


#        @edges[f].each do |t,e|
#
#                # handle deleted
#                if (e.del_e == 0)
#
#                    # handle opt open 
#                    opt_flag=false
#                    (opt_flag=true; print "["; opt=e.opt_e) unless e.opt_e<=opt 
#                    print "(" unless @edges[f].size==1
#
#                    # recurse
#                    e.to_smarts
#                    to_smarts(t,opt)
#                    print ")" unless @edges[f].size==1
#
#                    # handle opt close
#                    print "]" unless !opt_flag
#
#                end # end handle deleted
#            end
#        end
end

class LUNode
  attr_accessor :id, :lab_n, :weight_n, :del_n, :opt_n
  def initialize(id)
    @id = id
  end
  def to_smarts
    l_s = @lab_n.split
    if l_s.size > 1 
        for i in 0..l_s.size-1 do
            print "##{l_s[i]}"
            print "," unless i==l_s.size-1
        end
    else 
        print "##{@lab_n}"
    end
  end
end

class LUEdge
  attr_accessor :source, :target, :lab_e, :weight_e, :del_e, :opt_e
  def initialize(source, target)
    @source = source
    @target = target
  end
  def to_s
    puts "'#{@source} #{@target} #{@lab_e} #{@weight_e} #{@del_e} #{@opt_e}'"
  end
  def to_smarts
     l_e = @lab_e.split
    if l_e.size > 1 
        for i in 0..l_e.size-1
            case l_e[i]
                when '1' then print "-"
                when '2' then print "="
                when '3' then print "#"
                else 
            end
            print "," unless i==l_e.size-1
        end
    else 
        case @lab_e
            # Uncomment the next line for explicit aliphatic bindings in notation
            #when '1' then print "-"
            when '2' then print "="
            when '3' then print "#"
            else 
        end
    end
  end
end


def read
    # Store activites and hops seperately
    activities = {}
    hops = {}

    xml = STDIN.read
    doc, lu_graphs = Hpricot::XML(xml), []

    # For each graph
    graphs = {}
    (doc/:graph).each do |g|
        id=g.attributes['id'].to_i

        # For each data tag
        (g/:data).each do |d|
            key = d.attributes['key']
            assoc = case key
                when 'act' then activities
                when 'hops' then hops
                else nil
            end
            assoc[id] = d.inner_html.to_i unless assoc.nil?
        end

        # For each node tag
        graph_nodes = {}
        (g/:node).each do |n|
            node = LUNode.new(n.attributes['id'].to_i)
            (n/:data).each do |data|
                slot = data.inner_html
                case data.attributes['key']
                    when 'lab_n' then node.lab_n = slot
                    when 'weight_n' then node.weight_n = slot.to_i
                    when 'del_n' then node.del_n = slot.to_i
                    when 'opt_n' then node.opt_n = slot.to_i
                    else nil
                end
            end
            graph_nodes[n.attributes['id'].to_i]=node
        end

        # For each edge tag
        graph_edges = Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }
        (g/:edge).each do |e|
            edge = LUEdge.new(e.attributes['source'].to_i, e.attributes['target'].to_i)
            (e/:data).each do |data|
                slot = data.inner_html
                case data.attributes['key']
                    when 'lab_e' then edge.lab_e = slot
                    when 'weight_e' then edge.weight_e = slot.to_i
                    when 'del_e' then edge.del_e = slot.to_i
                    when 'opt_e' then edge.opt_e = slot.to_i
                    else nil
                end
            end
            graph_edges[e.attributes['source'].to_i][e.attributes['target'].to_i] = edge
        end
        
        graphs[id] = LUGraph.new(graph_nodes, graph_edges)
        
        #begin
        #    cnt=0
        #    graphs[id].edges.each do |f,_cmp_|
        #        cnt = cnt + graphs[id].edges[f].size
        #    end
        #    puts "Graph '#{id}' has '#{cnt}' edges and '#{graphs[id].nodes.size}' nodes."
        #end

    end

    return {:grps => graphs, :acts => activities, :hops => hops}
end

def smarts(dom)
    dom[:grps].sort{|a,b| a[0]<=>b[0]}.each do |id, g| 
        print "#{id}\t"
        g.to_smarts(0,0,0)
        print "\n"
    end
end


def match (smiles, smarts, verbose=true)
    c=OpenBabel::OBConversion.new
    c.set_in_format 'smi'
    m=OpenBabel::OBMol.new
    c.read_string m, smiles
    m.set_aromatic_perceived

    p=OpenBabel::OBSmartsPattern.new
    if !p.init(smarts)
        puts "Error! Smarts pattern invalid."
        exit
    end
   
    if verbose  
        p.match(m)
        hits = p.get_umap_list
        print "Found #{hits.size} instances of the SMARTS pattern '#{smarts}' in the SMILES string '#{smiles}'."
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

    p.match(m,true)
end

def match_file (file)
    smarts = STDIN.read
    File.open(file, "r") do |infile|
        while (line = infile.gets)
            result=""
            result << "#{line.split[0]} "
            smarts.each do |s|
                result << "#{s.split[0]}," unless !match(line.split[1],s.split[1],false)
            end
            result.chop!
            puts "#{result}\n"
        end
    end
end

# Main
STDOUT.sync = true

status=false
if $*.size==0 || $*.size>1
    status=true
end

case $*[0]
when '1' 
    if !status
        dom = read
        smarts(dom)
    end
when '2'
    if status && $*.size==2
        status=false
        match_file($*[1])
    end
when '3'
    demo
else
    status=true
end

if status
    puts "Usage: #{$0} cmd [/path/to/smifile.smi] < file" 
    puts "  cmd=1 : convert GraphML to SMARTS"
    puts "  cmd=2 : match SMARTS to SMILES file 'smifile.smi'"
    puts "  Output goes to $stdout."
    exit
end




# Demonstrates different SMARTS patterns
def demo
    # This pattern is equivalent to  [#7][#7][#6] :(
    puts
    match("NNC",                "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    match("N(=O)NC",            "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    match("N(=O)NCN",           "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    match("N(=O)NC(N)N",        "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    match("N(=O)NC=O",          "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    match("N(=O)NC(N)=O",       "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (yes)

    # This pattern is better (seems to demand a branch), but actually demands no branch since no explict degree is forced, which allows "folding" :(
    puts
    match("NNC",            "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)
    match("NN(C)C",         "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)
    match("NN(C)C(N)",      "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    match("NN(C)C(=O)",     "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    match("NN(CC)(C)C(=O)", "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    match("NNC(=O)",        "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)

    # This enhanced pattern enforces 1-step environments around the current node (f) and requires at least one branch :)
    # Recursive definition: 
    # p:
    # ------------------
    #          t: #int:n
    #      local: [t:self](p:branch)([t:back])[t:forw]
    # branchnode: [t:self;$(local:env),...](~*)
    #       node: [t]
    # ------------------
    puts
    match("NNC",                     "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # no (no)
    match("NN(C)C",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # no (no)
    match("NNC(N)",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # no (no)
    match("NN(C)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # yes (yes)
    match("NN(CC)(C)C(N)",           "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # yes (yes)
    match("NN(CC)(C)C(N)(=O)",       "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # yes (yes)
    match("NN(CC)C(=C)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # no (no)
    match("NN(CN)C(=N)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # yes (yes) <- two embeddings is correct
    match("NN(N)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])](~*)[#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])](~*)") # no (no)
end

#!/usr/bin/ruby1.8

require('openbabel')

#include OpenBabel   <= may be used for direct access to OB namespace, i.e. w/o "OpenBabel::". Below, I use namespaces for clarity.
require('set')

%w[pp rubygems hpricot].each { |x| require x }

class LUGraph
  attr_accessor :nodes, :edges
  def initialize(nodes, edges)
    @nodes = nodes
    @edges = edges
  end

  # LAST-SMARTS
  #   to_smarts(    nil,      0,0,  0,   1,<pa>)  <== DEFAULT
  def to_smarts(backw_e,backw_n,f,opt,init,mode)
      opt = @edges[0][1].weight unless !init # initialize opt level to weight of first edge
      mand_branches=0
      opt_branches = @edges[f].size
      different_opts = {}
      mand_t = []
      mand_e = []
      opt_t = []
      opt_e = []
      @edges[f].each do |t,e|
          if ((mode == 'ade' || mode == 'nde' || mode == 'ode') || e.del == 0)
              if (mode == 'nop' || mode == 'ode' || e.weight >= opt)
                  mand_branches+=1
                  opt_branches-=1
                  mand_t << t
                  mand_e << e
              else
                  opt_t << t
                  opt_e << e
              end
          else
              opt_branches-=1
          end
      end
      opt_e.each do |e|
          if different_opts.has_key?(e.weight)
            different_opts[e.weight]+=1
          else
            different_opts[e.weight]=1
          end
      end
      nr_b = 0
      if (mode == 'nla' || mode == 'nde')
        nr_b = different_opts.max[1] unless different_opts.size == 0 # EITHER nla variant (DEFAULT) OR
      elsif (mode == 'ala' || mode == 'ade')
        nr_b = different_opts.values.max unless different_opts.size == 0 # ala variant
      end

      if opt_branches != opt_t.size
        puts "Error! O: #{opt_branches} #{opt_t.size}"
        exit
      end
      if mand_branches != mand_t.size
        puts "Error! M: #{mand_branches} #{mand_t.size}"
        exit
      end
      if (((nr_b == 0) && (opt_branches > 0)) || (nr_b > opt_branches))
        puts "Error!"
        exit
      end

      print "["
      @nodes[f].to_smarts
      do_branch=(opt_branches>1)
      if do_branch
          # 1) Uncomment next line to not force at least one branch (c.f. not. 17.11.09)
          # print "$(" ; @nodes[f].to_smarts ; print ")," 
          # 1) end
          print ";"
          for i in 0..opt_branches-1
              t = opt_t[i]
              e = opt_e[i]
              print "$("

                  print "["                   # self
                  @nodes[f].to_smarts         # 
                  print "]"                   # 

                  print "("                   # branch
                  e.to_smarts                 # 
                  to_smarts(e,f,t,e.weight,0,mode)    # recursive: LAST-SMARTS
                  print ")"                   # 

                  if backw_n!=f && !backw_e.nil?          # backw
                      print "(" unless mand_branches==0   #
                      backw_e.to_smarts                   #
                      print "["                           #
                      @nodes[backw_n].to_smarts           #
                      print "]"                           #
                      print ")" unless mand_branches==0   #
                  end

                  for j in 0..mand_branches-1               # forw
                      print "(" unless j==mand_branches-1   #
                      mand_e[j].to_smarts                   #
                      print "["                             #
                      @nodes[mand_t[j]].to_smarts           #
                      print "]"                             #
                      print ")" unless j==mand_branches-1   #
                  end

              print ")"
              print "," unless i==opt_branches-1
          end
      end
      print "]"
      for i in 1..nr_b
        print "(~*)" unless !do_branch
      end

      # 2) Uncomment next block to make single optional edges mandatory (c.f. not. 17.11.09)
      if !do_branch && opt_branches==1
          t = opt_t[0]
          e = opt_e[0]
          print "(" unless mand_branches==0
          e.to_smarts
          to_smarts(e,f,t,e.weight,0,mode)
          print ")" unless mand_branches==0
      end
      # 2) end

      for i in 0..mand_branches-1
          t = mand_t[i]
          e = mand_e[i]
          print "(" unless i==mand_branches-1
          e.to_smarts
          to_smarts(e,f,t,e.weight,0,mode)
          print ")" unless i==mand_branches-1
      end
  end
end

class LUNode
  attr_accessor :id, :lab_n
  def initialize(id)
    @id = id
  end
  def to_smarts
    l_s = @lab_n.split
    if l_s.size > 1 
        for i in 0..l_s.size-1 do
            if l_s[i].to_i > 150 # aromatic carbon support
                print "##{l_s[i].to_i-150}&a"  
            else
                print "##{l_s[i]}"  
            end
            print "," unless i==l_s.size-1
        end
    else 
        if @lab_n.to_i > 150 # aromatic carbon support
            print "##{@lab_n.to_i-150}&a"  
        else
            print "##{@lab_n}"
        end
    end
  end
end

class LUEdge
  attr_accessor :source, :target, :lab_e, :weight, :del, :opt
  def initialize(source, target)
    @source = source
    @target = target
  end
  def to_s
    puts "'#{@source} #{@target} #{@lab_e} #{@weight} #{@del} #{@opt}'"
  end
  def to_smarts
     l_e = @lab_e.split
    if l_e.size > 1 
        for i in 0..l_e.size-1
            case l_e[i]
                when '1' then print "-"
                when '2' then print "="
                when '3' then print "#"
                when '4' then print ":"
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
            when '4' then print ":"
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
                    when 'weight' then edge.weight = slot.to_i
                    when 'del' then edge.del = slot.to_i
                    when 'opt' then edge.opt = slot.to_i
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

def smarts(dom, mode)
    dom[:grps].sort{|a,b| a[0]<=>b[0]}.each do |id, g| 
        #if g.edges[0][1].del == 0
            print "#{id}\t"
            print "#{dom[:acts][id]}\t"
            g.to_smarts(nil,0,0,0,1,mode)
            print "\n"
        #end
    end
end


def match (smiles, smarts, verbose=true)
    c=OpenBabel::OBConversion.new
    c.set_in_format 'smi'
    m=OpenBabel::OBMol.new
    c.read_string m, smiles
    m.kekulize

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
    smarts = STDIN.readlines

    # build act hash
    act_hash={}
    $stderr.puts "Building activity database..."
    File.open(file.sub(/.smi/, '.class')) do |cfile| # class
        while (line = cfile.gets)
            line_arr = line.split
            act_hash[line_arr.first]=line_arr.last
        end
    end   

    smi_arr=[]
    $stderr.puts "Reading instances..."
    File.open(file, "r") do |infile| # smi
        while (line = infile.gets)
            smi_arr << line
        end
    end

    $stderr.print "Processing smarts... "
    smarts.each do |s|
        result_hash={}
        string_result=""
        
        smi_arr.each do |smi|
            result_hash[smi.split.first] = act_hash[smi.split.first] unless !match(smi.split.last,s.split.last,false)
            #$stderr.print "#{result_hash.size} "
        end

        string_actives = " "
        string_inactives = " "

        result_hash.each do |id, act|
            if act == "1"
                string_actives << id << " "
            elsif act == "0"
                string_inactives << id << " "
            end
        end

        string_result << "\"#{s.split.last}\"\t["
        string_result << string_actives << "] ["
        string_result << string_inactives << "]"

        puts "#{string_result}"
        nom = string_actives.split.size
        den = (string_actives.split.size + string_inactives.split.size)
        $stderr.print "#{nom}/#{den}(#{s.split[1]}) "
    end
    $stderr.puts
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
    # See README for the recursive definition: 
    puts
    match("NNC",                     "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    match("NN(C)C",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    match("NNC(N)",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    match("NN(C)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    match("NN(C)C-N",                "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") 
    match("NN(CC)(C)C(N)",           "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    match("NN(CC)(C)C(N)(=O)",       "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    match("NN(CC)C(=C)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    match("NN(CN)C(=N)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes) <- two embeddings is correct
    match("NN(N)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
end



# Main
STDOUT.sync = true


status=false
if $*.size==0 || $*.size>2
    status=true
end

case $*[0]
when '1' 
    if !status
        dom = read
        smarts(dom, $*[1])
    end
when '2'
    if !status
        match_file($*[1])
    end
when '3'
    demo
else
    status=true
end

if status
    puts "Usage: #{$0} 1 [ ala | ade | nla | nde | nop | ode ] < /path/to/graphmlfile.graphml > /path/to/smartsfile.smarts" 
    puts "       #{$0} 2 /path/to/smifile.smi < /path/to/smartsfile.smarts > /path/to/lastpmfile.lastpm" 
    puts "       cmd=1 : convert GraphML to SMARTS."
    puts "           ala (ade) : All level max opt (no deletions)."
    puts "           nla (nde) : Next level opt (no deletions)."
    puts "           nop (ode) : No opt (no deletions)."
    puts "       cmd=2 : match SMARTS to SMILES file and create LASTPM file."
    exit
end

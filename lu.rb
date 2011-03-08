require('openbabel')

#include OpenBabel   <= may be used for direct access to OB namespace, i.e. w/o "OpenBabel::". Below, I use namespaces for clarity.
require('set')

%w[pp rubygems nokogiri].each { |x| require x }

class LUGraph
  attr_accessor :nodes, :edges
  def initialize(nodes, edges)
    @nodes = nodes
    @edges = edges
  end

  # LAST-SMARTS
  #   to_smarts(    nil,      0,0,  0,   1,<pa>)  <== DEFAULT
  def to_smarts(backw_e,backw_n,f,opt,init,mode,s=$stdout)
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
    if (mode == 'nls' || mode == 'nde')
      nr_b = different_opts.max[1] unless different_opts.size == 0 # EITHER nls variant (DEFAULT) OR
    elsif (mode == 'msa' || mode == 'ade')
      nr_b = different_opts.values.max unless different_opts.size == 0 # msa variant
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

    s.print "["
    @nodes[f].to_smarts(s)
    do_branch=(opt_branches>1)
    if do_branch
      # 1) Uncomment next line to not force at least one branch (c.f. not. 17.11.09)
      # s.print "$(" ; @nodes[f].to_smarts ; s.print ")," 
      # 1) end
      s.print ";"
      for i in 0..opt_branches-1
        t = opt_t[i]
        e = opt_e[i]
        s.print "$("

        s.print "["                   # self
        @nodes[f].to_smarts(s)         # 
        s.print "]"                   # 

        s.print "("                   # branch
        e.to_smarts(s)                 # 
        to_smarts(e,f,t,e.weight,0,mode,s)    # recursive: LAST-SMARTS
        s.print ")"                   # 

        if backw_n!=f && !backw_e.nil?          # backw
          s.print "(" unless mand_branches==0   #
          backw_e.to_smarts(s)                   #
          s.print "["                           #
          @nodes[backw_n].to_smarts(s)           #
          s.print "]"                           #
          s.print ")" unless mand_branches==0   #
        end

        for j in 0..mand_branches-1               # forw
          s.print "(" unless j==mand_branches-1   #
          mand_e[j].to_smarts(s)                #
          s.print "["                             #
          @nodes[mand_t[j]].to_smarts(s)           #
          s.print "]"                             #
          s.print ")" unless j==mand_branches-1   #
        end

        s.print ")"
        s.print "," unless i==opt_branches-1
      end
    end
    s.print "]"
    for i in 1..nr_b
      s.print "(~*)" unless !do_branch
    end

    # 2) Uncomment next block to make single optional edges mandatory (c.f. not. 17.11.09)
    if !do_branch && opt_branches==1
      t = opt_t[0]
      e = opt_e[0]
      s.print "(" unless mand_branches==0
      e.to_smarts(s)
      to_smarts(e,f,t,e.weight,0,mode,s)
      s.print ")" unless mand_branches==0
    end
    # 2) end

    for i in 0..mand_branches-1
      t = mand_t[i]
      e = mand_e[i]
      s.print "(" unless i==mand_branches-1
      e.to_smarts(s)
      to_smarts(e,f,t,e.weight,0,mode,s)
      s.print ")" unless i==mand_branches-1
    end
  end
end

class LUNode
  attr_accessor :id, :lab_n
  def initialize(id)
    @id = id
  end
  def to_smarts(s=$stdout)
    l_s = @lab_n.split
    if l_s.size > 1 
      for i in 0..l_s.size-1 do
        if l_s[i].to_i > 150 # aromatic carbon support
          s.print "##{l_s[i].to_i-150}"  
          s.print "&a"
        else
          s.print "##{l_s[i]}"  
          s.print "&A"
        end
        s.print "," unless i==l_s.size-1
      end
    else 
      if @lab_n.to_i > 150 # aromatic carbon support
        s.print "##{@lab_n.to_i-150}"  
        s.print "&a"
      else
        s.print "##{@lab_n}"
        s.print "&A"
      end
    end
  end
end

class LUEdge
  attr_accessor :source, :target, :lab_e, :weight, :del, :opt, :no_aromatic
  def initialize(source, target, no_aromatic=true)
    @source = source
    @target = target
    @no_aromatic = no_aromatic
  end
  def to_s
    puts "'#{@source} #{@target} #{@lab_e} #{@weight} #{@del} #{@opt}'"
  end
  def to_smarts(s=$stdout)
    l_e = @lab_e.split
    if l_e.size > 1 
      aromatized = false
      for i in 0..l_e.size-1
        case l_e[i]
        when '1' then 
          s.print "-"
        when '2' then 
          s.print "="
        when '3' then 
          s.print "#"
        when '4' then 
          s.print ":" 
          aromatized = true
        else 
        end
        s.print "," unless i==l_e.size-1
      end
      s.print ",:" unless no_aromatic || aromatized # always allow aromatic contexts
    else 
      case @lab_e
        # Uncomment the next line for explicit aliphatic bindings in notation
      when '1' then 
        s.print "-" 
      when '2' then 
        s.print "=" 
      when '3' then
        s.print "#" 
      when '4' then 
        s.print ":" 
        aromatized = true
      else
      end
      s.print ",:" unless no_aromatic || aromatized
    end
  end
end

class LU
  def read(xml=nil)
    # Store activites and hops seperately
    activities = {}
    hops = {}

    xml = STDIN.read unless !xml.nil?
    doc, lu_graphs = Nokogiri::XML(xml), []

    # For each graph
    graphs = {}


    ## (doc/:graph).each do |g|
    ##    id=g.attributes['id'].to_i
    doc.xpath('//graphml:graph', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |g|
      id=g['id'].to_i

      # For each data tag
      #(g/:data).each do |d|
      #    key = d.attributes['key']
      #    assoc = case key
      #        when 'act' then activities
      #        when 'hops' then hops
      #        else nil
      #    end
      #    assoc[id] = d.inner_html.to_i unless assoc.nil?
      #end
      #

      g.xpath('graphml:data', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |d|
        key = d['key']
        assoc = case key 
          when 'act' then activities
          when 'hops' then hops
          else nil
        end
        assoc[id] = d.text.to_i unless assoc.nil?
      }

      # For each node tag
      #graph_nodes = {}
      #(g/:node).each do |n|
      #  node = LUNode.new(n.attributes['id'].to_i)
      #  (n/:data).each do |data|
      #    slot = data.inner_html
      #    case data.attributes['key']
      #    when 'lab_n' then node.lab_n = slot
      #    else nil
      #    end
      #  end
      #  graph_nodes[n.attributes['id'].to_i]=node
      #end

      graph_nodes = {}
      g.xpath('graphml:node', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |n|
        nid = n['id'].to_i
        node = LUNode.new(nid)
        n.xpath('graphml:data', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |d|
          case d['key']
            when 'lab_n' then node.lab_n = d.text
            else nil
          end
        }
        graph_nodes[nid]=node
      }


      # For each edge tag
      #graph_edges = Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }
      #(g/:edge).each do |e|
      #  edge = LUEdge.new(e.attributes['source'].to_i, e.attributes['target'].to_i) # AM LAST: pass 'false' as 3rd arg to enable aromatic wildcarding
      #  (e/:data).each do |data|
      #    slot = data.inner_html
      #    case data.attributes['key']
      #    when 'lab_e' then edge.lab_e = slot
      #    when 'weight' then edge.weight = slot.to_i
      #    when 'del' then edge.del = slot.to_i
      #    when 'opt' then edge.opt = slot.to_i
      #    else nil
      #    end
      #  end
      #  graph_edges[e.attributes['source'].to_i][e.attributes['target'].to_i] = edge
      #end
      
      graph_edges = Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }
      g.xpath('graphml:edge', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |e|
        id1 = e['source'].to_i
        id2 = e['target'].to_i
        edge = LUEdge.new(id1, id2)
        e.xpath('graphml:data', {"graphml"=>"http://graphml.graphdrawing.org/xmlns"}).each { |d|
            case d['key']
               when 'lab_e' then edge.lab_e = d.text
               when 'weight' then edge.weight = d.text.to_i
               when 'del' then edge.del = d.text.to_i
               when 'opt' then edge.opt = d.text.to_i
               else nil
            end
        }
        graph_edges[id1][id2] = edge
      }

      graphs[id] = LUGraph.new(graph_nodes, graph_edges)

      #begin
      #    cnt=0
      #    graphs[id].edges.each do |f,_cmp_|
      #        cnt = cnt + graphs[id].edges[f].size
      #    end
      #    puts "Graph '#{id}' has '#{cnt}' edges and '#{graphs[id].nodes.size}' nodes."
      #end

    }

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

def smarts_rb(dom, mode)
  s = StringIO.new
  dom[:grps].sort{|a,b| a[0]<=>b[0]}.each do |id, g| 
    g.to_smarts(nil,0,0,0,1,mode,s)
    s.print(" ")
  end
  s.string.split # return array
end

def match (smiles, smarts, verbose=true)
  c=OpenBabel::OBConversion.new
  c.set_in_format 'smi'
  m=OpenBabel::OBMol.new
  c.read_string m, smiles
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

    fminer_output=false
    if ENV['FMINER_LAZAR'].size
            fminer_output=true
    end
    string_result << "\"" unless fminer_output
    string_result << "#{s.split.last}"
    string_result << "\"" unless fminer_output
    string_result << "\t["
    string_actives.chomp!(" ") if fminer_output
    string_result << string_actives
    string_result << "] [" unless fminer_output
    string_result << string_inactives << "]"

    puts "#{string_result}"
    nom = string_actives.split.size
    den = (string_actives.split.size + string_inactives.split.size)
    $stderr.print "#{nom}/#{den}(#{s.split[1]}) "
  end
  $stderr.puts
end

def match_rb (smiles,smarts) # AM LAST-PM: smiles= array id->smi
  result={}
  smarts.each do |s|
    ids=[]
    smiles.each_index do |id|
      if (id>1) 
        ids << id unless !match(smiles[id],s,false)
      end
    end
    result[s] = ids
  end
  result
end

def match_rb_hash (smiles,smarts) # AM LAST-PM: smiles= hash id->smi
  result={}
  smarts.each do |s|
    ids=[]
    smiles.each do |id,v|
      if (id>1) 
        ids << id unless !match(v,s,false)
      end
    end
    result[s] = ids
  end
  result
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
end


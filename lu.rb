# The library provides conversion from LAST-PM output to LAST-SMARTS 
# and matching of LAST-SMARTS to compounds.
# See http://last-pm.maunz.de
#
# Author::    Andreas Maunz (mailto:andreas@maunz.de)

begin 
  for g in [ "openbabel", "set", "yaml", "pp", "rubygems", "nokogiri"]
    require(g)
  end
rescue
  puts "\nlu.rb: Could not find gem '#{g}!'"
  exit false
end



# Represents a Graph
# Contains nodes and edges, where the labels can be multimaps

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
      puts "\nError! O: #{opt_branches} #{opt_t.size}"
      exit
    end
    if mand_branches != mand_t.size
      puts "\nError! M: #{mand_branches} #{mand_t.size}"
      exit
    end
    if (((nr_b == 0) && (opt_branches > 0)) || (nr_b > opt_branches))
      puts "\nError!"
      exit
    end

    s.print "["
    if @nodes[f].nil?
      puts "\nlu.rb: Node '#{f}' not found!"
      exit false
    else
      @nodes[f].to_smarts(s)
    end
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
        if @nodes[f].nil?
          puts "\nlu.rb: Node '#{f}' not found!"
          exit false
        else
          @nodes[f].to_smarts(s)
        end
        s.print "]"                   # 

        s.print "("                   # branch
        if e.nil?
          puts "\nlu.rb: Edge not found!"
          exit false
        else
          e.to_smarts(s)                 # 
        end
        to_smarts(e,f,t,e.weight,0,mode,s)    # recursive: LAST-SMARTS
        s.print ")"                   # 

        if backw_n!=f && !backw_e.nil?          # backw
          s.print "(" unless mand_branches==0   #

          if backw_e.nil?
            puts "\nlu.rb: Edge not found!"
            exit false
          else
            backw_e.to_smarts(s)                 # 
          end

          s.print "["                           #

          if @nodes[backw_n].nil?
            puts "\nlu.rb: Node '#{f}' not found!"
            exit false
          else
            @nodes[backw_n].to_smarts(s)
          end

          s.print "]"                           #
          s.print ")" unless mand_branches==0   #
        end

        for j in 0..mand_branches-1               # forw
          s.print "(" unless j==mand_branches-1   #
          mand_e[j].to_smarts(s)                #
          s.print "["                             #
          if @nodes[mand_t[j]].nil?
            puts "\nlu.rb: Node '#{j}' not found!"
            exit false
          else
            @nodes[mand_t[j]].to_smarts(s)
          end
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
      if e.nil?
        puts "\nlu.rb: Edge not found!"
        exit false
      else
        e.to_smarts(s)
      end
      to_smarts(e,f,t,e.weight,0,mode,s)
      s.print ")" unless mand_branches==0
    end
    # 2) end

    for i in 0..mand_branches-1
      t = mand_t[i]
      e = mand_e[i]
      s.print "(" unless i==mand_branches-1
      if e.nil?
        puts "\nlu.rb: Edge not found!"
        exit false
      else
        e.to_smarts(s)
      end
      to_smarts(e,f,t,e.weight,0,mode,s)
      s.print ")" unless i==mand_branches-1
    end
  end
end



# A Node
# Can be flagged as 'not aromatic' which is useful for
#  - Deliberate ambiguous notation, even if mining was done with arom. perc.
#  - Kekule mining (LAST-PM was used with '-a'

class LUNode
  attr_accessor :id, :lab_n, :no_aromatic
  def initialize(id, no_aromatic=false)
    @id = id
    @no_aromatic = no_aromatic
  end
  def to_smarts(s=$stdout)
    l_s = @lab_n.split
    if l_s.size > 1 
      for i in 0..l_s.size-1 do
        if l_s[i].to_i > 150 # aromatic carbon support
          s.print "##{l_s[i].to_i-150}"  
          s.print "&a" unless @no_aromatic
        else
          s.print "##{l_s[i]}"  
          s.print "&A" unless @no_aromatic
        end
        s.print "," unless i==l_s.size-1
      end
    else 
      if @lab_n.to_i > 150 # aromatic carbon support
        s.print "##{@lab_n.to_i-150}"  
        s.print "&a" unless @no_aromatic
      else
        s.print "##{@lab_n}"
        s.print "&A" unless @no_aromatic
      end
    end
  end
end



# An Edge
# Can be flagged as 'aromatic wildcarding' which is useful for
#  - Kekule mining (LAST-PM was used with '-a'

class LUEdge
  attr_accessor :source, :target, :lab_e, :weight, :del, :opt, :no_aromatic
  def initialize(source, target, aromatic_wc=false)
    @source = source
    @target = target
    @aromatic_wc = aromatic_wc
  end
  def to_s
    puts "\n'#{@source} #{@target} #{@lab_e} #{@weight} #{@del} #{@opt}'"
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
          puts "\nlu.rb: Edge type unknown!"
          exit false
        end
        s.print "," unless i==l_e.size-1
      end
      s.print ",:" if @aromatic_wc && !aromatized # always allow aromatic contexts
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
        puts "\nlu.rb: Edge type unknown!"
        exit false
      end
      s.print ",:" if @aromatic_wc && !aromatized
    end
  end
end


# SAX parser for reading GraphML
# Avoids memory overload

class XMLHandler < Nokogiri::XML::SAX::Document
  attr_accessor :graphs, :activities, :hops
  def initialize(no_aromatic, aromatic_wc)
    @no_aromatic = no_aromatic
    @aromatic_wc = aromatic_wc 
    @g_id = 0
    @n_id = 0
    @e_f = 0
    @e_t = 0
    @node = nil # intermediate holder
    @edge = nil #

    @graphs = {}
    @activities = {}
    @hops = {}

    @in_graph=false
    @graph_act=false
    @graph_hops=false
    @node_lab=false
    @edge_lab=false
    @edge_weight=false
    @edge_del=false

    @graph_nodes={}
    @graph_edges=Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }

    @in_node=false
    @in_edge=false
  end

  def start_element(name, attrs = [])
    case name
    when 'graph' then
      if @in_graph || @in_node || @in_edge
        puts "\nlu.rb: invalid structure (graph in graph or node or edge)!" 
        exit false
      end
      @in_graph=true
      @g_id=attrs[1].to_i  # get graph id
    when 'node'  then 
      if @in_graph then 
        if @in_node 
          puts "\nlu.rb: invalid structure (node in node)!" 
          exit false
        end
        if @in_edge
          puts "\nlu.rb: invalid structure (node in edge)!" 
          exit false
        end
        @in_node=true
        @n_id=attrs[1].to_i
      else
        puts "\nlu.rb: invalid structure (node outside graph)!" 
        exit false
      end # get node id
    when 'edge'  then
      if @in_graph then
        if @in_edge
          puts "\nlu.rb: invalid structure (edge in edge)!" 
          exit false
        end
        if @in_node
          puts "\nlu.rb: invalid structure (edge in node)!" 
          exit false
        end
        @in_edge=true;  
        @e_f=attrs[1].to_i
        @e_t=attrs[3].to_i   
      else
        puts "\nlu.rb: invalid structure (edge outside graph)!" 
        exit false
      end # get edge nodes
    when 'data' then 
      if @in_graph && attrs[1]=='act'   then  
        @graph_act=true   
      end # get graph act
      if @in_graph && attrs[1]=='hops'  then  
        @graph_hops=true  
      end  # get graph hops
      if @in_node && attrs[1]=='lab_n'  then  
        @node_lab=true   
      end # get graph act
      if @in_edge && attrs[1]=='lab_e'  then  
        @edge_lab=true    
      end  # get edge elements
      if @in_edge && attrs[1]=='weight' then  
        @edge_weight=true 
      end # 
      if @in_edge && attrs[1]=='del'    then  
        @edge_del=true    
      end  # 
    end
  end

  def characters(str)
    if @in_graph
      if @graph_act then 
        @activities[@g_id]=str.to_i
        @graph_act=false
      end  # OK, non-hierarchical
      if @graph_hops then 
        @hops[@g_id]=str.to_i;
        @graph_hops=false
      end  #
      if @in_node then 
        if @node.nil? then 
          @node=LUNode.new(@n_id,@no_aromatic) 
        end
        if @node_lab then 
          @node.lab_n=str
          @node_lab=false
        end
      end      # hierarchical.
      if @in_edge then   # hierarchical, non-trivial case.
        if @edge.nil? then 
          @edge=LUEdge.new(@e_f,@e_t,@aromatic_wc) 
        end
        if @edge_lab then 
          @edge.lab_e=str
          @edge_lab=false
        end
        if @edge_weight then 
          @edge.weight=str.to_i
          @edge_weight=false
        end
        if @edge_del then 
          @edge.del=str.to_i
          @edge_del=false
        end
      end
    end
  end

  def end_element(name)
    case name
    when 'graph' then 
      @graphs[@g_id] = LUGraph.new(@graph_nodes, @graph_edges)
      @in_graph=false;# create new graph here 
      @graph_nodes={}
      @graph_edges=Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }
    when 'node' then  
      if @graph_nodes[@n_id].nil? 
        @graph_nodes[@n_id]=@node 
      else
        puts "\nlu.rb: Can not insert node '#{@n_id}', since it exists."
        exit false
      end
      @in_node=false
      @node=nil  # store away
    when 'edge' then  
      @graph_edges[@e_f][@e_t]=@edge
      @in_edge=false
      @edge=nil  # 
    end
  end

end



class LU
  def read(xml=nil, no_aromatic=false, aromatic_wc=false)
    handler = XMLHandler.new(no_aromatic, aromatic_wc)
    parser = Nokogiri::XML::SAX::Parser.new(handler)
    parser.parse(STDIN)
    return {:grps => handler.graphs, :acts => handler.activities, :hops => handler.hops}
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
      puts "\nError! Smarts pattern invalid."
      exit
    end

    if verbose  
      p.match(m)
      hits = p.get_umap_list
      print "Found #{hits.size} instances of the SMARTS pattern '#{smarts}' in the SMILES string '#{smiles}'."
      if hits.size>0
        puts "\n Here are the atom indices:"
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
    File.open(file.sub(/.smi/, '.class').sub(/.nob/,'')) do |cfile| # class
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
  def ob_test 
    status = true
    # This pattern is equivalent to  [#7][#7][#6] :(
    puts
    status = false unless match("NNC",                "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    status = false unless match("N(=O)NC",            "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    status = false unless match("N(=O)NCN",           "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    status = false unless match("N(=O)NC(N)N",        "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    status = false unless match("N(=O)NC=O",          "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (no)
    status = false unless match("N(=O)NC(N)=O",       "[$([#7]),$([#7]=[#8])][$([#7]),$([#7][#6][#6]),$([#7][#6])][$([#6]),$([#6][#7]),$([#6]~[#7,#8])]")    # yes (yes)

    # This pattern is better (seems to demand a branch), but actually demands no branch since no explict degree is forced, which allows "folding" :(
    puts
    status = false unless match("NNC",            "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)
    status = false unless match("NN(C)C",         "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)
    status = false unless match("NN(C)C(N)",      "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    status = false unless match("NN(C)C(=O)",     "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    status = false unless match("NN(CC)(C)C(=O)", "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (yes)
    status = false unless match("NNC(=O)",        "[#7][$([#7][#6][#6]),$([#7][#6])][$([#6][#7]),$([#6]~[#7,#8])]")                         # yes (no)

    # This enhanced pattern enforces 1-step environments around the current node (f) and requires at least one branch :)
    # See README for the recursive definition: 
    puts
    status = false if match("NNC",                     "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    status = false if match("NN(C)C",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    status = false if match("NNC(N)",                  "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    status = false unless match("NN(C)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    status = false unless match("NN(C)C-N",                "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes 
    status = false unless match("NN(CC)(C)C(N)",           "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    status = false unless match("NN(CC)(C)C(N)(=O)",       "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes)
    status = false if match("NN(CC)C(=C)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    status = false unless match("NN(CN)C(=N)",             "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # yes (yes) <- two embeddings is correct
    status = false if match("NN(N)C(=N)",              "[#7][#7;$([#7]([#6][#6])([#7])[#6]),$([#7]([#6])([#7])[#6])][#6;$([#6]([#7])[#7]),$([#6](-,=[#7,#8])[#7])]") # no (no)
    puts "\nLAST-UTILS: Tests failed!" unless status
    status
  end
end


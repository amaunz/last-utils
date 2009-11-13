#!/usr/bin/env ruby

%w[pp rubygems hpricot].each { |x| require x }

class LUGraph
  attr_accessor :nodes, :edges
  def initialize(nodes, edges)
    @nodes = nodes
    @edges = edges

    def to_smarts(f,opt)
        if !@nodes[f].nil? 
            @nodes[f].to_smarts 
        else 
            puts "No Node #{f}."
            pp @nodes
        end
        @edges[f].each do |t,e|

                # handle deleted
                if (e.del_e == 0)

                    # handle opt open 
                    opt_flag=false
                    (opt_flag=true; print "["; opt=e.opt_e) unless e.opt_e<=opt 
                    print "(" unless @edges[f].size==1

                    # recurse
                    e.to_smarts
                    to_smarts(t,opt)
                    print ")" unless @edges[f].size==1

                    # handle opt close
                    print "]" unless !opt_flag

                end # end handle deleted
            end
        end
    end
end

class LUNode
  attr_accessor :id, :lab_n, :weight_n, :del_n, :opt_n
  def initialize(id)
    @id = id
  end
  def to_smarts
    l_s = @lab_n.split
    if l_s.size > 1 
        print "["
        for i in 0..l_s.size-1 do
            print "##{l_s[i]}"
            if i<l_s.size-1 
                print "," 
            end
        end
        print "]"
    else 
        print "[##{@lab_n}]"
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
        print "~" 
    else 
        case @lab_e
            when '1' then print "-"
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
        g.to_smarts(0,0)
        print "\n"
    end
end


# Main
STDOUT.sync = true
dom = read
smarts(dom)

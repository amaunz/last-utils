#!/usr/bin/env ruby

%w[pp rubygems hpricot].each { |x| require x }

class LUGraph
  attr_accessor :nodes, :edges
  def initialize(nodes, edges)
    @nodes = nodes
    @edges = edges
  end
end

class LUNode
  attr_accessor :id, :lab_n, :weight_n, :del_n, :opt_n
  def initialize(id)
    @id = id
  end
end

class LUEdge
  attr_accessor :source, :target, :lab_e, :weight_e, :del_e, :opt_e
  def initialize(source, target)
    @source = source
    @target = target
  end
end


# Main
if ARGV.size==0 || ARGV.size > 1
    puts "Usage: " << $0 << " <xmlfile>"
    exit
end

def read
    # Store activites and hops seperately
    activities = {}
    hops = {}

    xml = File.read(ARGV[0])
    doc, lu_graphs = Hpricot::XML(xml), []

    # For each graph
    graphs = {}
    (doc/:graph).each do |g|
        id=g.attributes['id']

        # For each data tag
        (g/:data).each do |d|
            key = d.attributes['key']
            assoc = case key
                when 'act' then activities
                when 'hops' then hops
                else nil
            end
            assoc[id] = d.inner_html unless assoc.nil?
        end

        # For each node tag
        graph_nodes = {}
        (g/:node).each do |n|
            node = LUNode.new(n.attributes['id'])
            (n/:data).each do |data|
                slot = case data.attributes['key']
                    when 'lab_n' then node.lab_n
                    when 'weight_n' then node.weight_n
                    when 'del_n' then node.del_n
                    when 'opt_n' then node.opt_n
                    else nil
                end
                slot = data.inner_html unless slot.nil? 
            end
            graph_nodes[n.attributes['id']]=node
        end

        # For each edge tag
        graph_edges = Hash.new{ |h,k| h[k]=Hash.new(&h.default_proc) }
        (g/:edge).each do |e|
            edge = LUEdge.new(e.attributes['source'], e.attributes['target'])
            (e/:data).each do |data|
                slot = case data.attributes['key']
                    when 'lab_e' then edge.lab_e
                    when 'weight_e' then edge.weight_e
                    when 'del_e' then edge.del_e
                    when 'opt_e' then edge.opt_e
                    else nil
                end
                slot = data.inner_html unless slot.nil?
            end
            graph_edges[e.attributes['source']][e.attributes['target']] = edge
        end
        
        graphs[id] = LUGraph.new(graph_nodes, graph_edges)
    end

    return {:grps => graphs, :acts => activities, :hops => hops}
end

dom = read
#pp dom[:acts]
#pp dom[:hops]
#pp dom[:grps]

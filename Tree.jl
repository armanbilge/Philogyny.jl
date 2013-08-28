# Tree.jl
# 
# Philogyny: Phylogenetics for the love of Julia
#

abstract Tree
abstract Node

type SimpleTree <: Tree
    id::String
    root::SimpleNode
    nodes::Vector{SimpleNode}
    node_count::Uint
    external_node_count::Uint
    internal_node_count::Uint
    in_edit::Bool
    traits::Dict{String,Any}
end

function SimpleTree(id::String, root::SimpleNode)
    node_count = get_node_count(node)
    nodes = Array(SimpleNode, count)
    for node in @task preorder_traversal(node)
        nodes[node.number] = node
    end
    SimpleTree(id, root, nodes, node_count, (node_count + 1) / 2,
               (node_count - 1) / 2, false, Dict{String,Any}())
end

type SimpleNode <: Node
    number::Uint
    parent::SimpleNode
    children::Vector{SimpleNode}
    height::Float64
    rate::Float64
    taxon::Taxon
    traits::Dict{String,Any}
end

type Taxon
    id::String
end

function is_root(tree::Tree, node::Node)
    tree.root == node
end

function is_external(node::Node)
    size(node.children)[1] == 0
end

function is_internal(node::Node)
    size(node.children)[1] > 0
end

function get_child_count(node::Node)
    size(node.children)[1]
end

function get_node_count(node::Node)
    x * get_leaf_count(node) - 1
end

function get_leaf_count(node::Node)
    child_count = get_child_count(node)
    if child_count == 0
        return 1
    end
    leaf_count = 0
    for child in node.children
        leaf_count += get_leaf_count(child)
    end
    return leaf_count
end

function preorder_traverse(tree::Tree)
    for node in @task preorder_traverse(tree.root)
        produce(node)
    end
end

function preorder_traverse(node::Node)
    produce(node)
    if not is_external(node)
        for child in node.children
            for child_node in @task preorder_traverse(node::Node)
                produce(child_node)
            end
        end
    end
end

function postorder_traverse(tree::Tree)
    for node in @task postorder_traverse(tree.root)
        produce(node)
    end
end

function postorder_traverse(node::Node)
    if is_internal(node)
        for child in node.children
            for child_node in @task preorder_traverse(node::Node)
                produce(child_node)
            end
        end
    end
    produce(node)
end


/*
 * FRICO - Friendly Radiation Integral COde
 *
 * Copyright (c) 2026, Matteo Cicuttin - IV3IWE
 * Politecnico di Torino
 * Dipartimento di Scienze Matematiche "G. L. Lagrange"
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <print>
#include <iostream>
#include <fstream>
#include <optional>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>
#include <queue>
#include <stack>

/* This is a trivial graph implementation. It is actually the code that I show
 * to my students, so don't expect too much. Assumptions:
 *  - no labels on edges or vertices, only connectivity
 *  - no isolated nodes, there's no add_node() precisely for this reason
 *  - DFS/BFS are classical and return trees, not forests. Remember this if
 *    your graph has more than one CC.
 *  - Performance? What is this? Adjacency list is stored as an unordered_map
 *    where neighbours are placed in sets. Easy, no issues with noncontinuous
 *    numbering or duplicated nodes, not optimal.
 */

namespace frico {

template<typename NodeIdT>
class undirected_graph {

public:
    class undirected_edge {
        NodeIdT   from_;
        NodeIdT   to_;

    public:
        undirected_edge(NodeIdT from, NodeIdT to)
            : from_(std::min(from, to)), to_(std::max(from, to))
        {}

        auto operator<=>(const undirected_edge& other) const = default;

        NodeIdT from() const {
            return from_;
        }

        NodeIdT to() const {
            return to_;
        }
    };

    using edge_type = undirected_edge;

private:
    std::unordered_map<NodeIdT, std::set<NodeIdT>>  adj_;

    mutable bool                            edgemap_dirty_;
    mutable std::vector<undirected_edge>    edgemap_;

    void recompute_edgemap() const {
        edgemap_.clear();
        for (const auto& [from, neighs] : adj_) {
            for (const auto& to : neighs) {
                edgemap_.push_back({from, to});
            }
        }
                    
        std::sort(edgemap_.begin(), edgemap_.end());
        edgemap_.erase(
            std::unique(edgemap_.begin(), edgemap_.end()),
            edgemap_.end()
        );
        edgemap_dirty_ = false;
    }

public:

    undirected_graph()
        : edgemap_dirty_(false)
    {}

    undirected_graph(const undirected_graph&) = default;
    undirected_graph(undirected_graph&&) = default;

    /* Return all neighbours of 'node', throws if 'node' does not exist */
    const std::set<NodeIdT>& neighbours(const NodeIdT& node) const {
        return adj_.at(node);
    }

    /* Return all neighbours of 'node', throws if 'node' does not exist */
    std::set<NodeIdT>& neighbours(const NodeIdT& node) {
        return adj_.at(node);
    }

    /* Return all nodes in the graph */
    std::vector<NodeIdT> all_nodes() const {
        std::vector<NodeIdT> ret;
        ret.reserve(adj_.size());
        for (auto& [from, neighs] : adj_) {
            ret.push_back(from);
        }
        return ret;
    }

    /* Return all edges. This perhaps should return by copy,
     * as the edgemap can change. Beware of this. */
    const std::vector<undirected_edge>& all_edges() const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }

        return edgemap_;
    }

    /* Add an edge between 'u' and 'v' */
    bool add_edge(const NodeIdT& u, const NodeIdT& v) {
        auto itor = adj_.find(u);
        if (itor != adj_.end()) {
            const auto& neighs_u = adj_[u];
            if ( neighs_u.find(v) != neighs_u.end() ) {
                std::cout << u << " " << v << std::endl;
                return false;
            } 
        }

        adj_[u].insert(v);
        adj_[v].insert(u);
        edgemap_dirty_ = true;
        return true;
    }

    /* Add an edge */
    bool add_edge(const undirected_edge& e) {
        return add_edge(e.from(), e.to());
    }

    /* Number of nodes in the graph */
    size_t num_nodes() const {
        return adj_.size();
    }

    /* Number of edges in the graph */
    size_t num_edges() const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }

        return edgemap_.size();
    }

    using EdgeIdT = size_t;

    /* Given an edge, return its number. Returns an empty optional
     * if the edge does not exist. */
    std::optional<EdgeIdT>
    edge_number(const undirected_edge& edge) const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }

        auto itor = std::lower_bound(edgemap_.begin(), edgemap_.end(), edge);
        if (itor != edgemap_.end() and *itor == edge)
            return std::distance(edgemap_.begin(), itor);
        
        return {};
    }

    /* Given an edge between 'u' and 'v', return its number. Returns an
     * empty optional if the edge does not exist. */
    std::optional<EdgeIdT>
    edge_number(const NodeIdT& u, const NodeIdT& v) const {
        return edge_number({u, v});
    }

    /* Return the edge in position 'eid' */
    undirected_edge edge_at(const EdgeIdT& eid) const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }

        assert(eid < edgemap_.size());
        return edgemap_[eid];
    }

    bool save_to_graphviz(const std::string& filename) const {
        std::ofstream ofs(filename);
        if (not ofs.is_open()) {
            std::println(stderr, "Can't open '{}'", filename);
            return false;
        }

        std::println(ofs, "strict graph G {{");
        for (const auto& edge : all_edges()) {
            std::println(ofs, "\t {} -- {} ;", edge.from(), edge.to());
        }
        std::println(ofs, "}}");
        return true;
    }

    auto begin() const { return adj_.begin(); }
    auto end() const { return adj_.end(); }

    /* Those two are pretty dangerous and should
     * not exist in their current form */
    auto edges_begin() const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }
        return edgemap_.begin();
    }
    auto edges_end() const {
        if (edgemap_dirty_) {
            recompute_edgemap();
        }
        return edgemap_.end();
    }

    /* Given two graphs G1 and G2, compute G = G1 \ G2. */
    undirected_graph operator-(const undirected_graph& other) const {
        undirected_graph ret;

        for (auto& [from, neighs] : adj_) {
            auto othadj_itor = other.adj_.find(from);
            if (othadj_itor == other.adj_.end()) {
                continue;
            }
            auto& retadj = ret.adj_[from];
            const auto& myadj = neighs;
            const auto& othadj = (*othadj_itor).second;
            std::set_difference(myadj.begin(), myadj.end(),
                othadj.begin(), othadj.end(), 
                    std::inserter(retadj, retadj.end()));
        }

        ret.recompute_edgemap();
        return ret;
    }
};

template<typename NodeIdT>
struct graph_edges {
    const undirected_graph<NodeIdT>& G_;
    graph_edges(const undirected_graph<NodeIdT>& G) : G_(G) {}
    auto begin() const { return G_.edges_begin(); }
    auto end() const { return G_.edges_end(); }
};

namespace detail {

enum node_color {
    white = 0,
    gray = 1,
    black = 2
};

/* Recursive DFS */
template<typename NodeIdT>
static void
dfs_aux(const undirected_graph<NodeIdT>& G, undirected_graph<NodeIdT>& dfstree,
    std::unordered_map<NodeIdT, node_color>& nodecolors, const NodeIdT& node)
{
    nodecolors[node] = node_color::gray;

    const auto& neighs = G.neighbours(node);
    for (const auto& n : neighs) {
        auto it = nodecolors.find(n);
        if (it != nodecolors.end() && it->second == node_color::white) {
            continue;
        }
        
        dfstree.add_edge(node, n);
        dfs_aux(G, dfstree, nodecolors, n);
    }

    nodecolors[node] = node_color::black;
}

/* Push and pop adapters */
template<typename NodeIdT>
void push(std::stack<NodeIdT>& stk, const NodeIdT& val)
{
    stk.push(val);
}

template<typename NodeIdT>
void push(std::queue<NodeIdT>& queue, const NodeIdT& val)
{
    queue.push(val);
}

template<typename NodeIdT>
NodeIdT pop(std::stack<NodeIdT>& stk)
{
    NodeIdT ret = stk.top();
    stk.pop();
    return ret;
}

template<typename NodeIdT>
NodeIdT pop(std::queue<NodeIdT>& queue)
{
    NodeIdT ret = queue.front();
    queue.pop();
    return ret;
}

/* Graph visit: does either BFS or DFS depending on VisitStorage */
template<typename NodeIdT, typename VisitStorage>
static undirected_graph<NodeIdT>
graph_visit(const undirected_graph<NodeIdT>& G,
    const NodeIdT& source)
{
    VisitStorage vs;
    undirected_graph<NodeIdT> visit_tree;
    std::unordered_map<NodeIdT, node_color> nodecolors;
    for (const auto& [from, neighs] : G) {
        nodecolors[from] = node_color::white;
    }
    
    nodecolors[source] = node_color::gray;

    push(vs, source);
    while ( not vs.empty() ) {
        NodeIdT u = pop(vs);

        for (const auto& v : G.neighbours(u)) {
            auto it = nodecolors.find(v);
            if (it != nodecolors.end() && it->second == node_color::white) {
                nodecolors[v] = node_color::gray;
                push(vs, v);
                visit_tree.add_edge(u,v);
            }
        }
        nodecolors[u] = node_color::black;
    }

    return visit_tree;
}

} // namespace detail

/* BFS wrapper */
template<typename NodeIdT>
static undirected_graph<NodeIdT>
bfs(const undirected_graph<NodeIdT>& G,
    const NodeIdT& source)
{
    return detail::graph_visit<NodeIdT, std::queue<NodeIdT>>(G, source);
}

/* DFS wrapper */
template<typename NodeIdT>
static undirected_graph<NodeIdT>
dfs(const undirected_graph<NodeIdT>& G,
    const NodeIdT& source)
{
    return detail::graph_visit<NodeIdT, std::stack<NodeIdT>>(G, source);
}

/* Recursive DFS wrapper */
template<typename NodeIdT>
undirected_graph<NodeIdT>
dfs_recursive(const undirected_graph<NodeIdT>& G, const NodeIdT& startnode)
{
    undirected_graph<NodeIdT> dfstree;
    std::unordered_map<NodeIdT, detail::node_color> nodecolors;
    for (const auto& [from, neighs] : G) {
        nodecolors[from] = detail::node_color::white;
    }
    detail::dfs_aux(G, dfstree, nodecolors, startnode);
    return dfstree;
}

/* Return a BFS tree */
template<typename NodeIdT>
undirected_graph<NodeIdT>
bfs_tree(const undirected_graph<NodeIdT>& G)
{
    if ( std::distance(G.begin(), G.end()) > 0 ) {
        NodeIdT firstnode = (*G.begin()).first;
        return bfs(G, firstnode);
    }
    return {};
}

/* Return a DFS tree */
template<typename NodeIdT>
undirected_graph<NodeIdT>
dfs_tree(const undirected_graph<NodeIdT>& G)
{
    if ( std::distance(G.begin(), G.end()) > 0 ) {
        NodeIdT firstnode = (*G.begin()).first;
        return dfs(G, firstnode);
    }
    return {};
}

namespace detail {
/* Given a DFS tree, find the path between 'curr' and 'to'.
 * 'dfsT' must be a tree.
 */
template<typename NodeIdT>
bool
find_path(const undirected_graph<NodeIdT>& dfsT, const NodeIdT& curr,
    const NodeIdT& to, std::set<NodeIdT>& visited,
    std::vector<NodeIdT>& path)
{
    visited.insert(curr);
    path.push_back(curr);

    if (curr == to) {
        return true;
    }

    for (const auto& n : dfsT.neighbours(curr)) {
        if ( !visited.count(n) ) {
            if (find_path(dfsT, n, to, visited, path)) {
                return true;
            }
        }
    }

    path.pop_back();
    return false;
}
}

/* Given a DFS tree dfsT and an edge (u,v) from the cotree, return the
 * list of nodes of the fundamental cycle associated to (u,v) */
template<typename NodeIdT>
std::vector<NodeIdT>
cycle_from_cotree_edge(const undirected_graph<NodeIdT>& dfsT,
                       const NodeIdT& u, const NodeIdT& v)
{
    std::set<NodeIdT> visited;
    std::vector<NodeIdT> path;
    detail::find_path(dfsT, u, v, visited, path);
    path.push_back(u);
    return path;
}

/* Get all the fundamental cycles of a graph using precalculated
 * DFS tree and cotree. Returns a list of edges.
 */
template<typename NodeIdT>
std::vector<std::vector<typename undirected_graph<NodeIdT>::edge_type>>
fundamental_cycles(const undirected_graph<NodeIdT>& dfs_tree,
    const undirected_graph<NodeIdT>& dfs_cotree)
{
    using edge_type = typename undirected_graph<NodeIdT>::edge_type;

    std::vector<std::vector<edge_type>> cycles;
    for (auto& edge : graph_edges(dfs_cotree) ) {
        auto nodes = cycle_from_cotree_edge(dfs_tree, edge.from(), edge.to());
        std::vector<edge_type> cycle;
        for (size_t i = 1; i < nodes.size(); i++) {
            cycle.push_back( {nodes[i-1], nodes[i]} );
        }
        cycles.push_back( std::move(cycle) );
    }

    return cycles;
}

/* Compute all the fundamental cycles of a graph */
template<typename NodeIdT>
std::vector<std::vector<typename undirected_graph<NodeIdT>::edge_type>>
fundamental_cycles(const undirected_graph<NodeIdT>& G)
{
    auto dfsT = dfs_tree(G);
    auto dfsC = G - dfsT;
    return fundamental_cycles(dfsT, dfsC);
}

template<typename NodeIdT, typename T>
using nmap = std::unordered_map<NodeIdT, T>;

template<typename NodeIdT>
void
dijkstra(const undirected_graph<NodeIdT>& G,
    nmap<NodeIdT, size_t>& dist, nmap<NodeIdT, NodeIdT>& prev, 
    const NodeIdT& src)
{   
    std::priority_queue<NodeIdT, std::vector<NodeIdT>, std::greater<NodeIdT>> Q;

    for (auto& [vtx, neighs] : G) {
        dist[vtx] = std::numeric_limits<size_t>::max();
        prev[vtx] = std::numeric_limits<NodeIdT>::max();
    }

    dist[src] = 0;
    Q.push(src);

    while (not Q.empty()) {
        NodeIdT u = Q.top();
        Q.pop();

        const auto& neighs_u = G.neighbours(u);
        for (auto& v : neighs_u) {
            size_t alt = dist[u] + 1;
            if (alt < dist[v]) {
                dist[v] = alt;
                prev[v] = u;
                Q.push(v);
            }
        }
    }
}

} // namespace frico

template<typename NodeIdT>
std::ostream&
operator<<(std::ostream& os,
    const typename frico::undirected_graph<NodeIdT>::edge_type& e) {
    os << "(" << e.from() << ", " << e.to() << ")";
    return os; 
}
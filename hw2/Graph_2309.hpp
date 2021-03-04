#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edge_map. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V, typename E>
class Graph {
private:

    /** Struct that holds all data of a node. */
    struct NodeInfo;
    /** Struct that holds all data of an edge. */
    struct EdgeInfo;

public:

    //
    // PUBLIC TYPE DEFINITIONS
    //

    /** Type of this graph. */
    using graph_type = Graph;

    /** Predeclaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;

    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;

    /** Type of node iterators, which iterate over all graph nodes. */
    class NodeIterator;
    /** Synonym for NodeIterator */
    using node_iterator = NodeIterator;

    /** Type of edge iterators, which iterate over all graph edge_map. */
    class EdgeIterator;
    /** Synonym for EdgeIterator */
    using edge_iterator = EdgeIterator;

    /** Type of incident iterators, which iterate incident edge_map to a node. */
    class IncidentIterator;
    /** Synonym for IncidentIterator */
    using incident_iterator = IncidentIterator;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    /** Type of node value.
        Return type of Graph::Node::value() */
    typedef V node_value_type;
    /** Type of edge value.
        Return type of Graph::Edge::value() */
    typedef E edge_value_type;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        nodes.clear();
        edges.clear();
        node_uids.clear();
        edge_uids.clear();
        edge_map.clear();
        adj.clear();
        nnodes = 0;
        nedges = 0;
    }

    /** Default destructor */
    ~Graph() = default;

    //
    // NODES
    //

    /** @class Graph::Node
     * @brief Class representing the graph's nodes.
     *
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private totally_ordered<Node> {
    public:
        /** Construct an invalid node.
         *
         * Valid nodes are obtained from the Graph class, but it
         * is occasionally useful to declare an @i invalid node, and assign a
         * valid node to it later. For example:
         *
         * @code
         * Graph::node_type x;
         * if (...should pick the first node...)
         *   x = graph.node(0);
         * else
         *   x = some other node using a complicated calculation
         * do_something(x);
         * @endcode
         */
        Node() {
        }

        /** Return this node's position. */
        Point& position() {
            return graph->nodes.at(graph->node_uids.at(ind)).position;
        }

        /** Return this node's position. */
        const Point& position() const {
            return graph->nodes.at(graph->node_uids.at(ind)).position;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type& index() {
            return ind;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        const size_type& index() const {
            return ind;
        }

        // HW1: YOUR CODE HERE
        /** Return this node's value. */
        node_value_type& value() {
            return graph->nodes.at(graph->node_uids.at(ind)).value;
        }

        /** Return this node's value. */
        const node_value_type& value() const {
            return graph->nodes.at(graph->node_uids.at(ind)).value;
        }

        /** Update this node's value in place.
         * @param[in] val The value to be assigned to this node
         * @post value() == val
         */
        void update_value(node_value_type val) {
            graph->nodes.at(graph->node_uids.at(ind)).value = val;
        }

        /** Return this node's degree. */
        size_type degree() const {
            return graph->nodes.at(graph->node_uids.at(ind)).degree;
        }

        /** Return this node's adjacency list stored as a unordered set. */
        std::unordered_set<size_type>& adj() const {
            return graph->adj.at(graph->node_uids.at(ind));
        }

        /** Return an incident iterator whose dereference is the first edge
         * incident to this node.
         */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph, graph->node_uids.at(ind), adj().begin());
        }

        /** Return an incident iterator whose dereference is the last edge
         * incident to this node.
         */
        incident_iterator edge_end() const {
            return IncidentIterator(graph, graph->node_uids.at(ind), adj().end());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return (graph == n.graph) && (ind == n.ind);
        }

        /** Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node& n) const {
            if (graph == n.graph) {
                return (ind < n.ind);
            }
            return (graph < n.graph);
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        graph_type* graph;
        size_type ind;

        /** Construct a (valid) node given graph and index. */
        Node(graph_type* g, const size_type i) : graph(g), ind(i) {};
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return this->nnodes;
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] p The new node's position
     * @param[in] v The new node's value, default to default value of
     *              node_value_type
     * @param[out] A Node with position p, value v or default value, index old
     *             num_nodes(), degree zero, and empty adjacency list
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     * @post result_node.degree() == 0
     * @post result_node.adj() == std::unordered_set<size_type>()
     * @post A new NodeInfo is stored in graph to hold all data of result_node
     *
     * Complexity: O(1) amortized operations.
     */
//    Node add_node(Point& p, const node_value_type& v = node_value_type()){
//        node_type n = Node(this, this->nnodes);
//        node_value_type v_copy = v;
//        NodeInfo nd = NodeInfo(&p, n, v_copy);
//        this->nodes.push_back(nd);
//        this->adj.push_back(std::unordered_set<size_type>());
//        ++this->nnodes;
//        return n;
//    }

    Node add_node(Point p, const node_value_type& v = node_value_type()){
        node_type n = Node(this, num_nodes());
        node_value_type v_copy = v;
        NodeInfo nd = NodeInfo(&p, n, v_copy);
        nodes.push_back(nd);
        node_uids.push_back(nodes.size() - 1);
        adj.push_back(std::unordered_set<size_type>());
        ++nnodes;
        return n;
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        return (this == n.graph) && (nnodes > n.index());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(const size_type i) const {
        assert(i < this->num_nodes());
        return this->nodes.at(node_uids.at(i)).node;
    }

    /** Remove a node from this graph
     *
     * @param n Node to be removed
     * @post all edges incident to the node are removed
     * @post num_nodes() = num_nodes() - 1 if removal succeeded;
     *       num_nodes() = num_nodes() otherwise
     * @post node indices rearranged
     * @return 1 if a node is removed; 0 if no node is removed
     */
    size_type remove_node(const Node& n){
        if (has_node(n)) {
            size_type i = n.index();
            for (auto it = n.edge_begin(); it != n.edge_end();) {
                it = remove_edge(it);
            }
            std::swap(node_uids.at(i), node_uids.back());
            nodes.at(node_uids.at(i)).node.index() = i;
            node_uids.pop_back();
            --nnodes;
            return 1;
        }
        return 0;
    }

    /** Remove a node from this graph
     *
     * @param n_it Node iterator pointing to the edge to be removed
     * @post num_nodes() = num_nodes() - 1 if removal succeeded;
     *       num_nodes() = num_nodes() otherwise
     * @post node indices rearranged
     * @return Node iterator e_it pointing to the next valid node or equal to
     *         node_end() of this graph
     */
    node_iterator remove_node(node_iterator n_it) {
        remove_node(*n_it);
        return n_it;
    }

    //
    // EDGES
    //

    /** @class Graph::Edge
     * @brief Class representing the graph's edge_map.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge() {};

        /** Return a node of this Edge */
        Node node1() const {
            return graph->nodes.at(node1_uid).node;
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return graph->nodes.at(node2_uid).node;
        }

        /** Return this edge's index, a number in the range [0, num_edges()). */
        size_type& index() {
            return this->ind;
        }

        const size_type& index() const {
            return this->ind;
        }

        /** Return this edge's value. */
        edge_value_type& value() {
            return graph->edges.at(graph->edge_uids.at(ind)).value;
        }

        /** Return this edge's value. */
        const edge_value_type& value() const {
            return graph->edges.at(graph->edge_uids.at(ind)).value;
        }

        /** Update this edge's value in place.
         *
         * @param val New value of this edge.
         */
        void update_value(edge_value_type val) {
            graph->edges.at(graph->edge_uids.at(ind)).value = val;
        }

        /** Return this edge's (current) length. */
        double length() const {
            return norm(node1().position() - node2().position());
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edge_map represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            if ((graph == e.graph) && (ind == e.ind)) {
                if (((node1_uid == e.node1_uid) && (node2_uid == e.node2_uid)) ||
                    ((node1_uid == e.node2_uid) && (node2_uid == e.node1_uid))) {
                    return true;
                }
            }
            return false;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            if (graph == e.graph) {
                return (ind < e.ind);
            }
            return (graph < e.graph);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        graph_type* graph;
        size_type node1_uid;    //node_uid
        size_type node2_uid;    //node_uid
        size_type ind;    // edge.index()

        /** Construct an edge given graph, indices of two nodes connected by
         * this edge, and edge index. */
        Edge(graph_type* g, size_type i1, size_type i2, size_type i) :
                graph(g), node1_uid(i1), node2_uid(i2), ind(i) {};
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     * Complexity: O(1).
     */
    size_type num_edges() const {
        return this->nedges;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     * Complexity: O(1).
     */
    Edge edge(size_type i) const {
        assert(i < this->nedges);
        return edges.at(edge_uids.at(i)).edge;
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        assert(has_node(a) && has_node(b));
//        assert(has_node(a) && has_node(b) && a.index() != b.index());
        return adj.at(node_uids.at(a.index())).count(node_uids.at(b.index()));
    }

    /** Add an edge to the graph, or return the current edge if it already exists.
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges() + 1.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b, const edge_value_type& v =
                  edge_value_type()) {
        size_type i1 = a.index();
        size_type i2 = b.index();
        size_type ui1 = node_uids.at(i1);
        size_type ui2 = node_uids.at(i2);
        if (has_edge(a, b)) {
            size_type i = edge_map[std::minmax(ui1, ui2)];
            return Edge(this, ui1, ui2, edges.at(i).edge.index());
        }

        assert(has_node(a) && has_node(b) && i1 != i2);
        edge_type e = Edge(this, ui1, ui2, num_edges());
        edge_value_type v_copy = v;
        EdgeInfo einfo = EdgeInfo(e, v_copy);
        edge_uids.push_back(edges.size());
        edge_map[std::minmax(ui1, ui2)] = edges.size();
        edges.push_back(einfo);
        adj.at(ui1).insert(ui2);
        adj.at(ui2).insert(ui1);
        ++nodes.at(ui1).degree;
        ++nodes.at(ui2).degree;
        ++nedges;
        return e;
    }

    /** Remove an edge from this graph
     *
     * @param a, b Endpoints of the edge to be removed
     * @post num_edges() = num_edges() - 1 if removal succeeded;
     *       num_edges() = num_edges() otherwise
     * @post adjacency lists of both endpoints updated
     * @post edge map updated
     * @post edge indices rearranged
     * @return 1 if an edge is removed; 0 if no edge is removed
     */
    size_type remove_edge(const Node& a, const Node& b) {
        if (has_edge(a, b)) {
            size_type ui1 = node_uids.at(a.index());
            size_type ui2 = node_uids.at(b.index());
            size_type eui = edge_map[std::minmax(ui1, ui2)];
            size_type ei = edges.at(eui).edge.index();
            std::swap(edge_uids.at(ei), edge_uids.back());
            edges.at(edge_uids.at(ei)).edge.index() = ei;
            edge_uids.pop_back();
            adj.at(ui1).erase(ui2);
            adj.at(ui2).erase(ui1);
            --nodes.at(ui1).degree;
            --nodes.at(ui2).degree;
            --nedges;
            return 1;
        }
        return 0;
    }

    /** Remove an edge from this graph
     *
     * @param e Edge to be removed
     * @post num_edges() = num_edges() - 1 if removal succeeded;
     *       num_edges() = num_edges() otherwise
     * @post adjacency lists of both endpoints updated
     * @post edge map updated
     * @post edge indices rearranged
     * @return 1 if an edge is removed; 0 if no edge is removed
     */
    size_type remove_edge(const Edge& e) {
        return remove_edge(e.node1(), e.node2());
    }

    /** Remove an edge from this graph
     *
     * @param e_it Edge iterator pointing to the edge to be removed
     * @post num_edges() = num_edges() - 1 if removal succeeded;
     *       num_edges() = num_edges() otherwise
     * @post adjacency lists of both endpoints updated
     * @post edge map updated
     * @post edge indices rearranged
     * @return Edge iterator e_it pointing to the next valid edge or equal to
     *         edge_end() of this graph
     */
    edge_iterator remove_edge(edge_iterator e_it){
        Edge e = *e_it;
        remove_edge(e.node1(), e.node2());
        return e_it;
    }

    /** Remove an edge from this graph
     *
     * @param i_it Incident iterator pointing to the edge to be removed
     * @post num_edges() = num_edges() - 1 if removal succeeded;
     *       num_edges() = num_edges() otherwise
     * @post adjacency lists of both endpoints updated
     * @post edge map updated
     * @post edge indices rearranged
     * @return Incident iterator i_it pointing to the next valid edge or equal to
     *         edge_end() of the corresponding node
     */
    incident_iterator remove_edge(incident_iterator i_it){
        Edge e = *i_it;
        if (has_edge(e.node1(), e.node2())) {
            size_type ui1 = node_uids.at(e.node1().index());
            size_type ui2 = node_uids.at(e.node2().index());
            size_type eui = edge_map[std::minmax(ui1, ui2)];
            size_type ei = edges.at(eui).edge.index();
            std::swap(edge_uids.at(ei), edge_uids.back());
            edges.at(edge_uids.at(ei)).edge.index() = ei;
            edge_uids.pop_back();
            if (i_it.node_ind == ui1){
                i_it.adj_it = adj.at(ui1).erase(i_it.adj_it);
                adj.at(ui2).erase(ui1);
            } else {
                adj.at(ui1).erase(ui2);
                i_it.adj_it = adj.at(ui2).erase(i_it.adj_it);
            }
            --nodes.at(ui1).degree;
            --nodes.at(ui2).degree;
            --nedges;
        }
        return i_it;
    }

    /** Remove all nodes and edge_map from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     * @post nodes, edges, edge_map, adj are empty
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        nodes.clear();
        edges.clear();
        node_uids.clear();
        edge_uids.clear();
        edge_map.clear();
        adj.clear();
        nnodes = 0;
        nedges = 0;
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : public std::iterator<std::forward_iterator_tag, Node>,
                         private equality_comparable<NodeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Node;                     // Element type
        using pointer = Node*;                    // Pointers to elements
        using reference = Node&;                    // Reference to elements
        using difference_type = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid NodeIterator. */
        NodeIterator() {
        }

        /** Returns the current node of this iterator.
         *
         * @return A node this iterator currently points to
         */
        Node operator*() const {
            return graph->nodes.at(ind).node;
        }

        /** Increment this iterator to point to the next node (or right after
         * the last node if this iterator is on the last node).
         *
         * @return The incremented-by-one NodeIterator
         */
        NodeIterator& operator++() {
            ++ind;
            return *this;
        }

        /** Test if a given node iterator is equal to this node iterator.
         *
         * Two NodeIterators are equal if they have the same graph pointer and
         * the same (node) index.
         *
         * @return true if the given iterator is equal to this iterator
         */
        bool operator==(const NodeIterator& iter) const {
            return (graph == iter.graph) && (ind == iter.ind);
        }

    private:
        friend class Graph;

        const graph_type* graph;
        size_type ind;

        /** Constructs a node iterator given a graph pointer and node index.
         *
         * @param g A pointer to the graph this iterator works for
         * @param i Index of the current node of this iterator
         */
        NodeIterator(const graph_type* g, const size_type i) :
                     graph(g), ind(i) {};
    };

    /** Returns a node iterator whose current node is the first node.
     *
     * @return A node iterator for this graph and with node index 0
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /** Returns a node iterator that points to right after the last node.
    *
    * @return A node iterator for this graph and with node index num_nodes(),
    *         i.e., index of last node + 1
    */
    node_iterator node_end() const {
        return NodeIterator(this, nnodes);
    }

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edge_map incident to a node. A forward iterator. */
    class IncidentIterator : private equality_comparable<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Edge;                     // Element type
        using pointer = Edge*;                    // Pointers to elements
        using reference = Edge&;                    // Reference to elements
        using difference_type = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() {
        }

        /** Returns the current edge of this iterator.
         *
         * @return An edge this iterator currently points to
         */
        Edge operator*() const {
            auto endpoints = std::minmax(node_ind, *adj_it);
            size_type ui = graph->edge_map.find(endpoints)->second;
            return graph->edges.at(ui).edge;
        }

        /** Increment this iterator to point to the next edge (or right after
         * the last edge if this iterator is on the last edge of its node),
         * where an edge is indexed by the corresponding neighbor node.
         *
         * @return The incremented-by-one IncidentIterator
         */
        IncidentIterator& operator++() {
            ++adj_it;
            return *this;
        }

        /** Test if a given incident iterator is equal to this node iterator.
         *
         * Two IncidentIterators are equal if they have the same graph pointer,
         * the same node index, and the same edge (neighbor) index.
         *
         * @return true if the given iterator is equal to this iterator
         */
        bool operator==(const IncidentIterator& iter) const {
            return (graph == iter.graph) && (node_ind == iter.node_ind)
                   && (adj_it == iter.adj_it);
        }

    private:
        friend class Graph;

        const graph_type* graph;
        size_type node_ind;
        std::unordered_set<size_type>::iterator adj_it;

        /** Constructs an incident iterator given a graph pointer, node index,
         * and edge (neighbor) index.
         *
         * @param g A pointer to the graph that this iterator works for
         * @param ni The index of the node in *g that this iterator works for
         * @param ai Index of the current edge (neighbor node) of this iterator
         */
        IncidentIterator(const graph_type* g, size_type ni,
                         std::unordered_set<size_type>::iterator ai)
                         : graph(g), node_ind(ni), adj_it(ai) {};
    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edge_map. A forward iterator. */
    class EdgeIterator : private equality_comparable<EdgeIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Edge;                     // Element type
        using pointer = Edge*;                    // Pointers to elements
        using reference = Edge&;                    // Reference to elements
        using difference_type = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {
        }

        /** Returns the current edge of this iterator.
         *
         * @return A edge this iterator currently points to
         */
        Edge operator*() const {
            return graph->edges.at(graph->edge_uids.at(ind)).edge;
        }

        /** Increment this iterator to point to the next edge (or right after
         * the last edge if this iterator is on the last edge).
         *
         * @return The incremented-by-one EdgeIterator
         */
        EdgeIterator& operator++() {
            ++ind;
            return *this;
        }

        /** Test if a given edge iterator is equal to this edge iterator.
         *
         * Two EdgeIterators are equal if they have the same graph pointer and
         * the same (edge) index.
         *
         * @return true if the given iterator is equal to this iterator
         */
        bool operator==(const EdgeIterator& iter) const {
            return (graph == iter.graph) && (ind == iter.ind);
        }

    private:
        friend class Graph;

        const graph_type* graph;
        size_type ind;

        /** Constructs an edge iterator given a graph pointer and node index.
         *
         * @param g A pointer to the graph this iterator works for
         * @param i Index of the current edge of this iterator
         */
        EdgeIterator(const graph_type* g, const size_type i) : graph(g),
                                                               ind(i) {};
    };

    /** Returns an edge iterator whose current edge is the first edge.
     *
     * @return An edge iterator for this graph and with edge index 0
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** Returns an edge iterator that points to right after the last edge.
    *
    * @return An edge iterator for this graph and with edge index num_edges(),
    *         i.e., index of last edge + 1
    */
    edge_iterator edge_end() const {
        return EdgeIterator(this, nedges);
    }


private:

    /** The struct stores all the data of a node. */
    struct NodeInfo {
        Point position;
        node_type node;
        node_value_type value;
        size_type degree;

        /** Construct a NodeInfo object given position, node (the corresponding
         * Node object), node value, and degree. Degree defaults to zero.  */
        NodeInfo(Point* p, node_type& n, node_value_type& v,
                 size_type d = 0) :
                position(*p), node(n), value(v), degree(d) {};
    };

    /** The struct stores all the data of an edge. */
    struct EdgeInfo {
        edge_type edge;
        edge_value_type value;

        /** Construct an EdgeInfo object given edge (the corresponding
         * Edge object), and edge value.  */
        EdgeInfo(edge_type& e, edge_value_type& v) : edge(e), value(v) {};
    };

    std::vector<NodeInfo> nodes;    // indexed by node_uid
    std::vector<EdgeInfo> edges;    // indexed by edge_uid
    std::vector<size_type> node_uids;    // indexed by node.index()
    std::vector<size_type> edge_uids;    // indexed by edge.index()
    std::map<std::pair<size_type, size_type>, size_type> edge_map;    // maps a
    // pair of node_uid to an edge_uid
    std::vector<std::unordered_set<size_type>> adj;    // indexed by node_uid
    size_type nnodes;
    size_type nedges;

    /** Test whether this graph and a given graph are equal. */
    bool operator==(const graph_type* g) const {
        return (this == g);
    }

};

#endif // CME212_GRAPH_HPP

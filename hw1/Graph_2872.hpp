#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edge_map. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph {
private:

    /** Struct that holds all data of a node, including position, node (the
     * corresponding Node object), value, degree. */
    struct NodeData;

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
    using node_value_type = V;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() {
        this->nodes.clear();
        this->edges.clear();
        this->edge_map.clear();
        this->nnodes = 0;
        this->nedges = 0;
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
        const Point& position() const {
            return this->graph->nodes.at(ind).position;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return this->ind;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;
        /** Return this node's value. */
        const node_value_type& value() const {
            return this->graph->nodes.at(ind).value;
        }

        /** Return this node's value. */
        node_value_type& value() {
//            return const_cast<node_value_type&>( static_cast<const Node&> (*this).value());
            return this->graph->nodes.at(ind).value;
        }

        /** Update this node's value in place.
         * @param[in] val The value to be assigned to this node
         * @post value() == val
         */
        void update_value(node_value_type val) {
            this->graph->nodes.at(ind).value = val;
        }

        /** Return this node's degree. */
        size_type degree() const {
            return this->graph->nodes.at(ind).degree;
        }

        /** Return this node's adjacency list stored as a set. */
        std::set<size_type> adj() const {
            return this->graph->adj[ind];
        }

        /** Return an incident iterator whose dereference is the first edge
         * incident to this node.
         */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph, ind, 0);
        }

        /** Return an incident iterator whose dereference is the last edge
         * incident to this node.
         */
        incident_iterator edge_end() const {
            return IncidentIterator(graph, ind, degree());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return (this->graph == n.graph) && (this->ind == n.ind);
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
            if (this->graph == n.graph) {
                return (this->ind < n.ind);
            }
            return (this->graph < n.graph);
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        graph_type* graph;
        size_type ind;

        /** Construct a (valid) node given graph and index. */
        Node(graph_type* g, const size_type i) {
            this->graph = g;
            this->ind = i;
        }
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
     * @post result_node.adj() == std::set<size_type>()
     * @post A new NodeData is stored in graph to hold all data of result_node
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& p, const node_value_type& v = node_value_type()){
        node_type n = Node(this, this->nnodes);
        node_value_type v_copy = v;
        NodeData nd = NodeData(&p, n, v_copy);
        this->nodes.push_back(nd);
        this->adj.push_back(std::set<size_type>());
        ++this->nnodes;
        return n;
    }


    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        return (this == n.graph) && (this->nnodes > n.index());
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(const size_type i) const {
        assert(i < this->num_nodes());
        return this->nodes.at(i).node;
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {};

        /** Return a node of this Edge */
        Node node1() const {
            return this->graph->node(this->ind1);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return this->graph->node(this->ind2);
        }

        /** Return this edge's index, a number in the range [0, num_edges()). */
        size_type index() const {
            return this->ind;
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edge_map represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            if ((this->graph == e.graph) && (this->ind == e.ind)) {
                if (((this->ind1 == e.ind1) && (this->ind2 == e.ind2)) ||
                    ((this->ind1 == e.ind2) && (this->ind2 == e.ind1))) {
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
            if (this->graph == e.graph) {
                return (this->ind < e.ind);
            }
            return (this->graph < e.graph);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        graph_type* graph;
        size_type ind1;
        size_type ind2;
        size_type ind;

        /** Construct an edge given graph, indices of two nodes connected by
         * this edge, and edge index. */
        Edge(graph_type* g, size_type i1, size_type i2, size_type i) {
            this->graph = g;
            this->ind1 = i1;
            this->ind2 = i2;
            this->ind = i;
        }
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
        return this->edges.at(i);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        assert(this->has_node(a) && this->has_node(b) && a.ind != b.ind);
        return this->adj[a.ind].count(b.ind);
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
    Edge add_edge(const Node& a, const Node& b) {
        size_type i1 = a.ind;
        size_type i2 = b.ind;
        if (has_edge(a, b)) {
            size_type i = edge_map[std::minmax(i1, i2)];
            return Edge(this, i1, i2, i);
        }

        assert(this->has_node(a) && this->has_node(b) && i1 != i2);
        edge_type e = Edge(this, i1, i2, this->nedges);
        this->edge_map[std::minmax(i1, i2)] = this->nedges;
        this->edges.push_back(e);
        this->adj[i1].insert(i2);
        this->adj[i2].insert(i1);
        ++this->nedges;
        ++this->nodes.at(i1).degree;
        ++this->nodes.at(i2).degree;
        return e;
    }

    /** Remove all nodes and edge_map from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     * @post nodes, edges, edge_map, adj are empty
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        this->nodes.clear();
        this->edges.clear();
        this->edge_map.clear();
        this->adj.clear();
        this->nnodes = 0;
        this->nedges = 0;
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
         * @return A edge this iterator currently points to
         */
        Edge operator*() const {
            auto endpoints = std::minmax(node_ind, adj_ind);
            size_type i = graph->edge_map.find(endpoints)->second;
            return graph->edge(i);
        }

        /** Increment this iterator to point to the next edge (or right after
         * the last edge if this iterator is on the last edge of its node),
         * where an edge is indexed by the corresponding neighbor node.
         *
         * @return The incremented-by-one IncidentIterator
         */
        IncidentIterator& operator++() {
            ++adj_ind;
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
                   && (adj_ind == iter.adj_ind);
        }

    private:
        friend class Graph;

        const graph_type* graph;
        const size_type node_ind;
        size_type adj_ind;  // index of a neighbor node in the adjacency list

        /** Constructs an incident iterator given a graph pointer, node index,
         * and edge (neighbor) index.
         *
         * @param g A pointer to the graph that this iterator works for
         * @param ni The index of the node in *g that this iterator works for
         * @param ai Index of the current edge (neighbor node) of this iterator
         */
        IncidentIterator(const graph_type* g, const size_type ni, size_type ai)
                : graph(g), node_ind(ni), adj_ind(ai) {}
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
            return graph->edges.at(ind);
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
    struct NodeData {
        const Point position;
        node_type node;
        node_value_type value;
        size_type degree;

        /** Construct a NodeData object given position, node (the corresponding
         * Node object), value, degree. Degree defaults to zero.  */
        NodeData(const Point* p, node_type n, node_value_type& v,
                 size_type d = 0) :
                position(*p), node(n), value(v), degree(d) {};
    };

    std::vector<NodeData> nodes;
    std::vector<edge_type> edges;
    std::map<std::pair<size_type, size_type>, size_type> edge_map;   // maps a
    // pair of node indices to an edge index
    std::vector<std::set<size_type>> adj;  // adjacency lists
    size_type nnodes;
    size_type nedges;

    /** Test whether this graph and a given graph are equal. */
    bool operator==(const graph_type* g) const {
        return (this == g);
    }

};

#endif // CME212_GRAPH_HPP

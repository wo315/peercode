#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <set>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


 /** @class Graph
  * @brief A template for 3D undirected graphs.
  *
  * Users can add and retrieve nodes and edges. Edges are unique (there is at
  * most one edge between any pair of distinct nodes).
  */
class Graph {
private:

    // HW0: YOUR CODE HERE
    // Use this space for declarations of important internal types you need
    // later in the Graph's definition.
    // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
    // code here. Just use the space if you need it.)
    
    // Predeclare internal node struct
    struct internal_node;

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

    /** Type of edge iterators, which iterate over all graph edges. */
    class EdgeIterator;
    /** Synonym for EdgeIterator */
    using edge_iterator = EdgeIterator;

    /** Type of incident iterators, which iterate incident edges to a node. */
    class IncidentIterator;
    /** Synonym for IncidentIterator */
    using incident_iterator = IncidentIterator;

    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty graph. */
    Graph() : nodes_(), size_(0), next_uid_(0) {
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
    class Node {
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
            uid_ = 0;
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
        }

        /** Return this node's position. */
        const Point& position() const {
            return graph_->nodes_[uid_].point;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return uid_;
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // node_value_type& value();
        // const node_value_type& value() const;
        // size_type degree() const;
        // incident_iterator edge_begin() const;
        // incident_iterator edge_end() const;

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            bool equal_nodes = false;
            if (uid_ == n.uid_ && graph_ == n.graph_) {
                equal_nodes = true;
            }
            return equal_nodes;
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
            bool less_than = false;
            if (uid_ < n.uid_) {
                less_than = true;
            }
            return less_than;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        Graph* graph_;
        size_type uid_;

        // Private constructor
        Node(const Graph* graph, size_type uid)
            : graph_(const_cast<Graph*>(graph)), uid_(uid) {
        }
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        return size_;
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position) {
        // Set the point information and UID for the new node
        internal_node new_node;
        new_node.point = position;
        new_node.uid = next_uid_;
        nodes_.emplace(next_uid_, new_node);

        // Update node UIDs vector and Node/Edge dictionary
        nodes_uids_.push_back(next_uid_);
        node_edge_dictionary_.emplace(next_uid_, std::set<size_type> {});

        // Update graph information
        size_++;
        next_uid_++;

        return Node(this, next_uid_ - 1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        bool node_exists = false;
        if ((n.graph_ == this) && (n.index() < size())) {
            node_exists = true;
        }
        return node_exists;
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        return Node(this, nodes_uids_[i]);
    }

    //
    // EDGES
    //

    /** @class Graph::Edge
     * @brief Class representing the graph's edges.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
            uid_a_ = 0;
            uid_b_ = 1;
        }

        /** Return a node of this Edge */
        Node node1() const {
            return Node(graph_, uid_a_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            return Node(graph_, uid_b_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            bool equal_nodes = false;
            if ((uid_a_ == e.uid_a_ && uid_b_ == e.uid_b_) ||
                (uid_a_ == e.uid_b_ && uid_b_ == e.uid_a_)) {
                equal_nodes = true;
            }
            return equal_nodes;
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            bool less_than = false;
            if (
                (uid_a_ + uid_b_ < e.uid_a_ + e.uid_b_) ||
                ((uid_a_ + uid_b_ == e.uid_a_ + e.uid_b_) &&
                    (std::min(uid_a_, uid_b_) < std::min(e.uid_a_, e.uid_b_))
                    )) {
                less_than = true;
            }
            return less_than;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        Graph* graph_;
        size_type uid_a_;
        size_type uid_b_;

        // Private constructor
        Edge(const Graph* graph, size_type uid_a, size_type uid_b)
            : graph_(const_cast<Graph*>(graph)), uid_a_(uid_a), uid_b_(uid_b) {
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return edges_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        return Edge(this, edges_[i].first, edges_[i].second);
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        bool edge_exists = false;

        // Obtain Node a's edges and check for Node b
        std::set<size_type> a_edges = node_edge_dictionary_.at(a.uid_);
        auto edge_a_b = a_edges.find(b.uid_);
        if (edge_a_b != a_edges.end()) {
            edge_exists = true;
        }

        return edge_exists;
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
        // Check if edge already exists
        if (has_edge(a, b) == false) {
            edges_.push_back(std::make_pair(a.uid_, b.uid_));

            // Update edge dictionary
            node_edge_dictionary_[a.uid_].emplace(b.uid_);
            node_edge_dictionary_[b.uid_].emplace(a.uid_);
        }
        return Edge(this, a.uid_, b.uid_);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // Node data
        nodes_.clear();
        nodes_uids_.clear();
        size_ = 0;
        next_uid_ = 0;

        // Edge data
        edges_.clear();
        node_edge_dictionary_.clear();
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator {
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

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Node operator*() const
        // NodeIterator& operator++()
        // bool operator==(const NodeIterator&) const

    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator {
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

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator {
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

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // EdgeIterator& operator++()
        // bool operator==(const EdgeIterator&) const

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const

private:
    // Nodes
    struct internal_node {
        Point point;
        size_type uid;
    };

    std::map<size_type, internal_node> nodes_;
    std::vector<size_type> nodes_uids_;
    size_type size_;
    size_type next_uid_;

    // Note: the nodes_uids_ std::vector serves to accommodate for
    // Node re-indexing in future implementations.

    // Edges
    std::vector<std::pair<size_type, size_type>> edges_;
    std::map<size_type, std::set<size_type>> node_edge_dictionary_;

};

#endif // CME212_GRAPH_HPP


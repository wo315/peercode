#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

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
    Graph() {
        // HW0: YOUR CODE HERE
    }

    /** Default destructor */
    ~Graph() = default;

    //
    // NODES
    //
    // -------------------------------------------------------------------------
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
        Node()
            : graph_(), uid_(), point_() {
            // HW0: YOUR CODE HERE
        }

        Node(Graph* graph, size_type uid, const Point* point)
            : graph_(graph),
              uid_(uid),
              point_(point) {}

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            return *point_;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
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
            // HW0: YOUR CODE HERE
            return (n.graph_ == graph_) && (n.uid_ == uid_);
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
            // HW0: YOUR CODE HERE
            return uid_ < n.uid_;
        }

       private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects

        // Pointer back to the graph
        graph_type* graph_;
        // This node's unique identification number
        size_type uid_;
        const Point* point_;
    };  // class node
    // -------------------------------------------------------------------------
    /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
    size_type size() const {
        // HW0: YOUR CODE HERE
        return nodes_.size();
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
        // HW0: YOUR CODE HERE
        auto new_node = new Node(this, num_nodes(), new Point(position));
        nodes_.push_back(new_node);
        return *new_node;
    }

    /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        assert(n.index() < num_nodes());
        return true;
    }

    /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(i < num_nodes());
        return *nodes_[i];
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
        Edge()
            : node1_(), node2_(), nid_() {
            // HW0: YOUR CODE HERE
        }

        Edge(const Node* node1, const Node* node2, size_type nid)
            : node1_(node1), node2_(node2), nid_(nid) {}

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return *node1_;  // Invalid Node
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return *node2_;  // Invalid Node
        }

        /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
        bool operator==(const Edge& e) const {
            return ((*node1_ == e.node1()) && (*node2_ == e.node2())) || ((*node2_ == e.node1()) && (*node1_ == e.node2()));
        }

        /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
        bool operator<(const Edge& e) const {
            return nid_ < e.nid_;
        }

       private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects
        const Node* node1_;
        const Node* node2_;
        size_type nid_;
    };  //clase Edge
    // -------------------------------------------------------------------------

    /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        return edges_.size();
    }

    /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(i < num_edges());
        return *edges_[i];  // Invalid Edge
    }

    /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        Edge in_edge(&a, &b, 0);
        for (size_t i = 0; i < edges_.size(); i++) {
            if (in_edge == *edges_[i]) {
                return true;
            }
        }
        return false;
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
        // HW0: YOUR CODE HERE
        auto in_edge = new Edge(&a, &b, edges_.size());
        for (size_t i = 0; i < edges_.size(); i++) {
            if (*in_edge == *edges_[i]) {
                return *edges_[i];
            }
        }
        edges_.push_back(in_edge);
        return *in_edge;
    }

    /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
    void clear() {
        // HW0: YOUR CODE HERE
        auto edgs_num = edges_.size();
        for (size_t i = 0; i < edgs_num; i++) {
            edges_.pop_back();  //only clear the vector of edge,don't care about the node's memory
        }
        auto nodes_num = nodes_.size();
        for (size_t i = 0; i < nodes_num; i++) {
            delete nodes_[nodes_num - 1 - i]->point_;  // releas the point's memory first
            delete nodes_[nodes_num - 1 - i];          // releas the node's memory first
            nodes_.pop_back();
        }
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator {
       public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Node;                            // Element type
        using pointer = Node*;                              // Pointers to elements
        using reference = Node&;                            // Reference to elements
        using difference_type = std::ptrdiff_t;             // Signed difference
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
        using value_type = Edge;                            // Element type
        using pointer = Edge*;                              // Pointers to elements
        using reference = Edge&;                            // Reference to elements
        using difference_type = std::ptrdiff_t;             // Signed difference
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
        using value_type = Edge;                            // Element type
        using pointer = Edge*;                              // Pointers to elements
        using reference = Edge&;                            // Reference to elements
        using difference_type = std::ptrdiff_t;             // Signed difference
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
    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.
    std::vector<Node*> nodes_;
    std::vector<Edge*> edges_;
};

#endif  // CME212_GRAPH_HPP

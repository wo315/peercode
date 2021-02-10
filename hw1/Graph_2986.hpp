#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <set>
#include <unordered_map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


 /** @class Graph
  * @brief A template for 3D undirected graphs.
  *
  * Users can add and retrieve nodes and edges. Edges are unique (there is at
  * most one edge between any pair of distinct nodes).
  */
template <typename V>
class Graph {
private:
   
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

    /** User-specified node values */
    using node_value_type = V;

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
    Graph() : nodes_({}), size_(0), next_uid_(0), edges_({}), node_edge_dictionary_() {
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
    class Node : private totally_ordered<Node>{
    public:
        /** 
         * @brief Construct an invalid node.
         * @return A Node object associated to an empty graph
         */
        Node() {
            uid_ = 0;
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
        }

        /** 
         * @brief Return this node's position
         * @pre The node's uid_ represents a valid position in the graph's node
         * dictionary
         * @return a reference to a Point containing the node's position information
         */
        const Point& position() const {
            return graph_->nodes_[uid_].point;
        }

        /** 
         * @brief Return this node's index. 
         * @return The node's index of type size_type, a number in the range [0, graph_size)
         * 
         * Complexity: O(1)
         */
        size_type index() const {
            return uid_;
        }

        /**
         * @brief Get the Node's user-specified value 
         * @return a reference to the value in the graph's node dictionary
         * 
         * Complexity: O(1)
         */
        node_value_type& value() {
            return graph_->nodes_[uid_].value_;
        }
        
        /**
         * @brief Get the Node's user-specified value 
         * @return a reference to the (const) value in the graph's node dictionary
         * 
         * Complexity: O(1)
         */
        const node_value_type& value() const {
            return graph_->nodes_[uid_].value_;
        }
        
        /**
         * @brief Return the number of incident edges of the Node
         * @return a value of type size_type representing the number of adjacent edges
         * 
         * Complexity: O(n.degree())
         */
        size_type degree() const {
            // Handle invalid nodes
            if (graph_->nodes_.size() == 0) {
                assert(0);
            }
            else {
				size_type num_edges = 0;
				for (auto it = this->edge_begin(); it != this->edge_end(); ++it) {
					num_edges++;
				}
				return num_edges;
            }
        }
        
        /**
         * @brief Return an iterator to the beginning of the Node's adjacent edges
         * @return an IncidentIterator pointing to the start of the edges
         */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph_, uid_, 0);
        }

        /**
         * @brief Return an iterator to the end of the Node's adjacent edges
         * @return an IncidentIterator pointing to one past the end of the edges
         */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_, uid_,
                graph_->node_edge_dictionary_[uid_].size());
        }

        /** 
         * @brief Test whether this node and @a n are equal.
         * @return a boolean value, true if the nodes share the same index and graph
         * 
         * Complexity: O(1)
         */
        bool operator==(const Node& n) const {
            bool equal_nodes = false;
            if (uid_ == n.uid_ && graph_ == n.graph_) {
                equal_nodes = true;
            }
            return equal_nodes;
        }

        /** 
         * @brief Test whether this node is less than @a n in a global order.
         * @return a boolean value, true if the node's index is less @a n's index
         * 
         * Complexity: O(1)
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

        /* Private Node constructor */
        Node(const Graph* graph, size_type uid)
            : graph_(const_cast<Graph*>(graph)), uid_(uid) { }
    };

    /** 
     * @brief Return the number of nodes in the graph.
     * @return a value of type size_type representing the number of nodes
     * Complexity: O(1).
     */
    size_type size() const {
        return size_;
    }

    /** Synonym for size(). */
    size_type num_nodes() const {
        return size();
    }

    /** 
     * @brief Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @param[in] node_value_type the type of the nodes' user-specified value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& = node_value_type()) {
        // Set the point information and UID for the new node
        internal_node new_node;
        new_node.point = position;
        new_node.uid = next_uid_;
        new_node.value_ = node_value_type();
        nodes_.push_back(new_node);

        // Update Node-Edge dictionary
        node_edge_dictionary_.emplace(next_uid_, std::vector<size_type> {});

        // Update graph information
        size_++;
        next_uid_++;

        return Node(this, next_uid_ - 1);
    }

    /** 
     * @brief Determine if a Node belongs to this Graph
     * @param[in] n the node to be checked
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

    /** 
     * @brief Return the node with index @a i.
     * @param[in] i the index of the desired node
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        return Node(this, i);
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
    class Edge : private totally_ordered<Edge>{
    public:
        /** 
         * @brief Construct an invalid Edge.
         * @return an Edge object associated to an empty graph and two nonexistent nodes
         */
        Edge() {
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
            uid_a_ = 0;
            uid_b_ = 1;
        }

        /** 
         * @brief Return a node of this Edge
         * @return a Node object associated to this edge
         * 
         * Complexity: O(1)
         */
        Node node1() const {
            return Node(graph_, uid_a_);
        }

        /** 
         * @brief Return the other node of this Edge 
         * @return a Node object associated to this edge
         * 
         * Complexity: O(1)
         */
        Node node2() const {
            return Node(graph_, uid_b_);
        }

        /** 
         * @brief Test whether this edge and @a e are equal.
         * @param[in] e, the other edge to be verified for equality
         * @return a bool value, true if the edges represent the same 
         *      undirected edge between two nodes.
         * 
         * Complexity: O(1)
         */
        bool operator==(const Edge& e) const {
            bool equal_nodes = false;
            if (((uid_a_ == e.uid_a_ && uid_b_ == e.uid_b_) ||
                (uid_a_ == e.uid_b_ && uid_b_ == e.uid_a_))
                && (graph_ == e.graph_)) {
                equal_nodes = true;
            }
            return equal_nodes;
        }

        /** 
         * @brief Test whether this edge is less than @a e in a global order.
         * @param[in] e, the other edge to be verified for inequality
         * @return a bool value (of no geometric significance) comparing the edges
         * 
         * Complexity: O(1)
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

        /* Private Edge constructor */
        Edge(const Graph* graph, size_type uid_a, size_type uid_b)
            : graph_(const_cast<Graph*>(graph)), uid_a_(uid_a), uid_b_(uid_b) {
        }
    };

    /** 
     * @brief Return the total number of edges in the graph.
     * @return a value of type size_type representing the total number of edges
     * 
     * Complexity: O(1)
     */
    size_type num_edges() const {
        return edges_.size();
    }

    /** 
     * @brief Return the edge with index @a i. 
     * @param[in] i, the index of the desired edge
     * @return an Edge object
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: O(1)
     */
    Edge edge(size_type i) const {
        return Edge(this, edges_[i].first, edges_[i].second);
    }

    /** 
     * @brief Test whether two nodes are connected by an edge.
     * @param[in] a, an incident node
     * @param[in] b, another incident node
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: O(a.degree())
     */
    bool has_edge(const Node& a, const Node& b) const {
        bool edge_exists = false;

        // Obtain Node a's edges and check for Node b
        std::vector<size_type> a_edges = node_edge_dictionary_.at(a.uid_);
        auto it = std::find(a_edges.begin(), a_edges.end(), b.uid_);
        if (it != a_edges.end()) {
            edge_exists = true;
        }

        return edge_exists;
    }

    /** 
     * @brief Add an edge to the graph, or return the current edge if it already exists.
	 * @param[in] a, an incident node
	 * @param[in] b, another incident node
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
     *       Else,                        new num_edges() == old num_edges() + 1.
     *
     * Complexity: O(a.degree())
     */
    Edge add_edge(const Node& a, const Node& b) {
        // Check if edge already exists
        if (has_edge(a, b) == false) {
            edges_.push_back(std::make_pair(a.uid_, b.uid_));

            // Update edge dictionary
            node_edge_dictionary_[a.uid_].push_back(b.uid_);
            node_edge_dictionary_[b.uid_].push_back(a.uid_);
        }
        return Edge(this, a.uid_, b.uid_);
    }

    /** 
     * @brief Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Complexity: O(1)
     */
    void clear() {
        // Node data
        nodes_.clear();
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

        /** 
         * @brief Construct an invalid NodeIterator.
         * @return a NodeIterator object associated to an empty graph
         */
        NodeIterator() {
            node_i_ = 0;
			Graph empty_graph = Graph();
			graph_ = &empty_graph;
        }

        /**
         * @brief Define increment operator on NodeIterator
         * @return a NodeIterator pointing to the following node in the container
         */
        NodeIterator& operator++() {
            node_i_++;
            return *this;
        }

        /**
         * @brief Define dereference operator on NodeIterator
         * @return the Node object of the graph and index pointed to by the NodeIterator
         */
		Node operator*() const {
            return Node(graph_, node_i_);
		}

        /**
         * @brief Define equality operation on NodeIterators
         * @param[in] node_iterator_b, another NodeIterator object to be compared
         * @return a bool value, true if the NodeIterators share the same graph
         *      and point to the same node index in the graph
         */
        bool operator==(const NodeIterator& node_iterator_b) const {
            bool equal_node = false;
            if ((graph_ == node_iterator_b.graph_) && 
                (node_i_ == node_iterator_b.node_i_)) {
                equal_node = true;
            }
            return equal_node;
        }

        /**
         * @brief Define inequality operation on NodeIterators
         * @param[in] node_iterator_b, another NodeIterator object to be compared
		 * @return a bool value, true if the NodeIterators do not share the same graph
		 *      or do not point to the same node index in the graph
         */
        bool operator!=(const NodeIterator& node_iterator_b) const {
            return !(*this == node_iterator_b);
        }

    private:
        friend class Graph;
        // Attributes for node location
        Graph* graph_;
        size_type node_i_;

        /* Private NodeIterator constructor that can be accessed by the Graph class */
        NodeIterator(const Graph* graph, size_type node_i) :
            graph_(const_cast<Graph*>(graph)), node_i_(node_i) {};
    };

    /**
     * @brief NodeIterator pointing to the beginning of the graph's node container
     * @return a NodeIterator object
     * 
     * Complexity: O(1)
     */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }

    /**
     * @brief NodeIterator pointing to one past the end of the graph's node container
     * @return a NodeIterator object
     * 
     * Complexity: O(1)
     */
    node_iterator node_end() const {
        return NodeIterator(this, size_);
    }

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

        /** 
         * @brief Construct an invalid IncidentIterator. 
         * @return an IncidentIterator object associated to an empty graph
         */
        IncidentIterator() {
			Graph empty_graph = Graph();
			graph_ = &empty_graph;
            node_i_ = 0;
            adj_node_j_ = 0;
        }

        /**
         * @brief Define dereference operator on IncidentIterator
         * @return the Edge object of the graph and node indexes pointed to by the
         *      iterator
         * 
         * Complexity: O(1)
         */
        Edge operator*() const {
            return Edge(graph_, node_i_, node_edges_[adj_node_j_]);
        }

        /**
         * @brief Define increment operator on IncidentIterator
         * @return an IncidentIterator pointing to the following node in the
         *      incident node container
         */
        IncidentIterator& operator++() {
			adj_node_j_++;
			return *this;
        }

        /**
         * @brief Define equality operation on IncidentIterator objects
         * @param[in] incident_it_b, an IncidentIterator to be compared
         * @return a bool value, true if the IncidentIterators share the same graph
         *      and the same node indexes
         */
        bool operator==(const IncidentIterator& incident_it_b) const {
            bool iterator_equals = false;
            if ((graph_ == incident_it_b.graph_) && (node_i_ == incident_it_b.node_i_) && 
                (adj_node_j_ == incident_it_b.adj_node_j_)) {
                iterator_equals = true;
            }
            return iterator_equals;
        }

        /**
         * @brief Define inequality operation on IncidentIterator objects
		 * @param[in] incident_it_b, an IncidentIterator to be compared
		 * @return a bool value, true if the IncidentIterators do not share the same graph
		 *      or have differing same node indexes
         */
		bool operator!=(const IncidentIterator& incident_it_b) const {
            return !(*this == incident_it_b);
		}

    private:
        friend class Graph;

		// Attributes for node to which edges are incident
		Graph* graph_;
		size_type node_i_;
        size_type adj_node_j_;
        std::vector<size_type> node_edges_ = graph_->node_edge_dictionary_[node_i_];

		/* Private IncidentIterator constructor that can be accessed by the Graph class */
        IncidentIterator(const Graph* graph, size_type node_i, size_type adj_node_j) :
			graph_(const_cast<Graph*>(graph)), node_i_(node_i), adj_node_j_(adj_node_j) {};

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

        /**
         * @brief Construct an invalid EdgeIterator.
         * @return an EdgeIterator object associated to an empty graph
         */
        EdgeIterator() {
			Graph empty_graph = Graph();
			graph_ = &empty_graph;
			edge_i_ = 0;
        }

        /**
         * @brief Define dereference operator on EdgeIterator objects
         * @return the Edge object of the graph and edge index being pointed to
         */
        Edge operator*() const {
            return Edge(graph_, graph_->edges_[edge_i_].first, 
                graph_->edges_[edge_i_].second);
        }

        /**
         * @brief Define increment operator on EdgeIterator objects
         * @return an EdgeIterator object pointing to the next edge in the graph's
         *      edge container
         */
        EdgeIterator& operator++() {
            edge_i_++;
            return *this;
        }

        /**
         * @brief Define equality operator on EdgeIterators
         * @param[in] edge_iter, an EdgeIterator object to be compared
         * @return a bool value, true if the iterators share the same graph and
         *      edge index
         */
        bool operator==(const EdgeIterator& edge_iter) const {
            bool equal_iter = false;
            if ((graph_ == edge_iter.graph_) && (edge_i_ == edge_iter.edge_i_)) {
                equal_iter = true;
            }
            return equal_iter;
        }

        /**
         * @brief Define inequality operator on EdgeIterators
         * @param[in] edge_iter, an EdgeIterator object to be compared
         * @return a bool value, true if the iterators do not share the same graph
         *      or have differing edge indexes
         */
		bool operator!=(const EdgeIterator& edge_iter) const {
            return !(*this == edge_iter);
		}

    private:
        friend class Graph;

		// Attributes for edge location
		Graph* graph_;
		size_type edge_i_;

		/* Private EdgeIterator constructor that can be accessed by the Graph class */
		EdgeIterator(const Graph* graph, size_type edge_i) :
			graph_(const_cast<Graph*>(graph)), edge_i_(edge_i) {};

    };

    /**
     * @brief EdgeIterator pointing to the beginning of a graph's edge container
     * @return an EdgeIterator object
     * 
     * Complexity: O(1)
     */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /**
     * @brief EdgeIterator pointing to one past the end of a graph's edge container
     * @return an EdgeIterator object
     * 
     * Complexity: O(1)
     */
    edge_iterator edge_end() const {
        return EdgeIterator(this, edges_.size());
    }

private:
    /* Define the internal_node struct which contains information for each node */
    struct internal_node {
        Point point;
        size_type uid;
        V value_;
    };

    /* Node container and graph size information */
    std::vector<internal_node> nodes_;
    size_type size_;
    size_type next_uid_;

    /* Edge container and Node-Edge dictionary */
    std::vector<std::pair<size_type, size_type>> edges_;
    std::unordered_map<size_type, std::vector<size_type>> node_edge_dictionary_;

};

#endif // CME212_GRAPH_HPP


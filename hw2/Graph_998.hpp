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
template <typename V, typename E>
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

    /** Edge values */
    typedef E edge_value_type;

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
    Graph() : nodes_({}), i2u_({}), next_uid_(0), edges_({}),
        edge_values_({}), node_edge_dictionary_() {
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
    private:
        bool valid() const {
            return uid_ >= 0 && uid_ < graph_->nodes_.size()
                && graph_->nodes_[uid_].idx_ < graph_->i2u_.size()
                && graph_->i2u_[graph_->nodes_[uid_].idx_] == uid_;
        }

    public:
        /**
         * @brief Construct an invalid node.
         * @return A Node object associated to an empty graph
         */
        Node() {
            idx_ = 0;
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
        }

        /**
         * @brief Return this node's position
         * @pre The node's idx_ represents a valid position in the graph's i2u_ vector
         * @return a reference to a Point containing the node's position information
         */
        const Point& position() const {
            assert(this->valid());
            return graph_->nodes_.at(uid_).point;
        }

        /**
         * @brief Return this node's modifiable position
         * @pre The node's idx_ represents a valid position in the graph's i2u_ vector
         * @return a modifiable reference to a Point containing the node's position information
         */
        Point& position() {
            assert(this->valid());
            return graph_->nodes_.at(uid_).point;
        }

        /**
         * @brief Return this node's index.
         * @return The node's index of type size_type, a number in the range [0, graph_size)
         *
         * Complexity: O(1)
         */
        size_type index() const {
            return idx_;
        }

        /**
         * @brief Get the Node's user-specified value
         * @return a reference to the value in the graph's node dictionary
         *
         * Complexity: O(1)
         */
        node_value_type& value() {
            assert(this->valid());
            return graph_->nodes_.at(uid_).value_;
        }

        /**
         * @brief Get the Node's user-specified value
         * @return a reference to the (const) value in the graph's node dictionary
         *
         * Complexity: O(1)
         */
        const node_value_type& value() const {
            assert(this->valid());
            return graph_->nodes_.at(uid_).value_;
        }

        /**
         * @brief Return the number of incident edges of the Node
         * @return a value of type size_type representing the number of adjacent edges
         *
         * Complexity: O(n.degree())
         */
        size_type degree() const {
            assert(this->valid());
            size_type num_edges = 0;
            for (auto it = this->edge_begin(); it != this->edge_end(); ++it) {
                num_edges++;
            }
            return num_edges;
        }

        /**
         * @brief Return an iterator to the beginning of the Node's adjacent edges
         * @return an IncidentIterator pointing to the start of the edges
         */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph_, idx_, 0);
        }

        /**
         * @brief Return an iterator to the end of the Node's adjacent edges
         * @return an IncidentIterator pointing to one past the end of the edges
         */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_, idx_,
                graph_->node_edge_dictionary_.at(uid_).size());
        }

        /**
         * @brief Test whether this node and @a n are equal.
         * @return a boolean value, true if the nodes share the same index, uid and graph
         *
         * Complexity: O(1)
         */
        bool operator==(const Node& n) const {
            bool equal_nodes = false;
            if (uid_ == n.uid_ && graph_ == n.graph_ && idx_ == n.idx_) {
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
        size_type idx_;
        size_type uid_;

        /* Private Node constructor */
        Node(const Graph* graph, size_type idx)
            : graph_(const_cast<Graph*>(graph)), idx_(idx) {
            uid_ = graph_->i2u_.at(idx_);
        }
    };

    /**
     * @brief Return the number of nodes in the graph.
     * @return a value of type size_type representing the number of nodes
     * Complexity: O(1).
     */
    size_type size() const {
        return this->i2u_.size();
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
    Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
        // Set the point information and index for the new node
        internal_node new_node;
        new_node.point = position;
        new_node.value_ = val;
        new_node.idx_ = this->i2u_.size();

        // Update the nodes_ vector with the internal_node
        nodes_.push_back(new_node);

        // Update the i2u_ vector with the node's UID
        i2u_.push_back(next_uid_);

        // Update Node-Edge dictionary
        node_edge_dictionary_.emplace(next_uid_, std::vector<size_type> {});

        // Update edge values dictionary
        edge_values_.emplace(next_uid_, std::unordered_map<size_type, edge_value_type>{});

        // Update graph information
        next_uid_++;

        return Node(this, new_node.idx_);
    }

    /**
     * @brief Determine if a Node belongs to this Graph
     * @param[in] n the node to be checked
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        return n.valid() && n.graph_ == this;
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
    class Edge : private totally_ordered<Edge> {
    public:
        /**
         * @brief Construct an invalid Edge.
         * @return an Edge object associated to an empty graph and two nonexistent nodes
         */
        Edge() {
            Graph empty_graph = Graph();
            graph_ = &empty_graph;
            idx_a_ = 0;
            idx_b_ = 1;
        }

        /**
         * @brief Return a node of this Edge
         * @return a Node object associated to this edge
         *
         * Complexity: O(1)
         */
        Node node1() const {
            return Node(graph_, idx_a_);
        }

        /**
         * @brief Return the other node of this Edge
         * @return a Node object associated to this edge
         *
         * Complexity: O(1)
         */
        Node node2() const {
            return Node(graph_, idx_b_);
        }

        /**
         * @brief Return a uniform Edge length
         * @return double
         */
        double length() const {
            return 0.025;
        }

        /**
         * @brief Get the Edge's user-specified value
         * @return a reference to the value in the graph's edge value dictionary
         *
         * Complexity: O(1)
         */
        edge_value_type& value() {
            if (uid_a_ < uid_b_) {
                return graph_->edge_values_.at(uid_a_).at(uid_b_);
            }
            else {
                return graph_->edge_values_.at(uid_b_).at(uid_a_);
            }

        }

        /**
         * @brief Get the Edge's user-specified value
         * @return a reference to the (const) value in the graph's edge value dictionary
         *
         * Complexity: O(1)
         */
        const edge_value_type& value() const {
            if (uid_a_ < uid_b_) {
                return graph_->edge_values_.at(uid_a_).at(uid_b_);
            }
            else {
                return graph_->edge_values_.at(uid_b_).at(uid_a_);
            }
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
            if (this->graph_ == e.graph_) {
                if (((uid_a_ + uid_b_ < e.uid_a_ + e.uid_b_) || ((uid_a_ + uid_b_ == e.uid_a_ + e.uid_b_) &&
                        (std::min(uid_a_, uid_b_) < std::min(e.uid_a_, e.uid_b_))))) {
                    less_than = true;
                }
            } else {
                less_than = true;
            }
            return less_than;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        Graph* graph_;
        size_type idx_a_;
        size_type idx_b_;
        size_type uid_a_;
        size_type uid_b_;

        /* Private Edge constructor */
        Edge(const Graph* graph, size_type idx_a, size_type idx_b)
            : graph_(const_cast<Graph*>(graph)), idx_a_(idx_a), idx_b_(idx_b) {
            uid_a_ = graph_->i2u_[idx_a_];
            uid_b_ = graph_->i2u_[idx_b_];
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
        size_type uid_a = edges_[i].first;
        size_type uid_b = edges_[i].second;
        return Edge(this, this->nodes_[uid_a].idx_, this->nodes_[uid_b].idx_);
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

            // Update edge value dictionary
            if (a.uid_ < b.uid_) {
                edge_values_.at(a.uid_).emplace(b.uid_, edge_value_type{});
            }
            else {
                edge_values_.at(b.uid_).emplace(a.uid_, edge_value_type{});
            }
        }
        return Edge(this, a.idx_, b.idx_);
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
        i2u_.clear();
        next_uid_ = 0;

        // Edge data
        edges_.clear();
        node_edge_dictionary_.clear();
        edge_values_.clear();
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
        // Attributes for node location (graph and node index as we only iterate
        // through active nodes)
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
        return NodeIterator(this, i2u_.size());
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
            node_edges_ = {};
        }

        /**
         * @brief Define dereference operator on IncidentIterator
         * @return the Edge object of the graph and node indexes pointed to by the
         *      iterator
         *
         * Complexity: O(1)
         */
        Edge operator*() const {
            // Grab the uid and idx of the adjacent node
            size_type uid_j = node_edges_.at(adj_node_j_);
            size_type idx_j = graph_->nodes_.at(uid_j).idx_;
            return Edge(graph_, node_i_, idx_j);
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
        size_type node_i_; // Idx of the main node
        size_type adj_node_j_; // Local index of the adjacent node's position in the adjacent node vector
        std::vector<size_type> node_edges_;

        /* Private IncidentIterator constructor that can be accessed by the Graph class */
        IncidentIterator(const Graph* graph, size_type node_i, size_type adj_node_j) :
            graph_(const_cast<Graph*>(graph)), node_i_(node_i), adj_node_j_(adj_node_j) {
            // Grab the uid of Node with index node_i
            size_type uid_i = graph_->i2u_.at(node_i_);
            node_edges_ = graph_->node_edge_dictionary_.at(uid_i);
        };

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
            // Grab the uids of the adjacent edges
            size_type uid_a = graph_->edges_.at(edge_i_).first;
            size_type uid_b = graph_->edges_.at(edge_i_).second;
            return Edge(graph_, graph_->nodes_.at(uid_a).idx_,
                graph_->nodes_.at(uid_b).idx_);
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

    /* -------------------- Node and Edge removal ------------------------- */

    /**
     * @brief Remove a node from the graph as specified by a Node object. This
     *      removes the node's connecting edges as well.
     * @param[in] @n a Node object indicating the node to be removed from the graph
     * @pre Node @n is a valid node of the graph
     * @return a size_type indicating if the node was removed (1) or not (0)
     * @post If 1 is returned, new i2u_.size() == old i2u_.size() - 1
     *       Else if 0 is returned, new i2u_.size() == old i2u_.size()
     * @post new nodes_.size() == old nodes_.size()
     * @post if 1 is returned, new num_edges() == old num_edges() - old n.degree()
     *       Else if 0 is returned, new num_edges() == old num_edges()
     * @post for 0 <= i < graph.size(), nodes_[i2u_[i]].idx < graph.size()
     * @post If 1 is returned, new node_edge_dictionary_.at(n.uid_).size() = 0
     * @post If 1 is returned, new edge_values_.at((*n_it).uid_).size() = 0
     * 
     * Complexity: O(n.degree())
     * 
     */
    size_type remove_node(const Node& n) {
        if (n.idx_ < i2u_.size()) {

            size_type idx = n.idx_;

            // Identify node edges
            std::vector<Edge> edges_to_remove = {};
            for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
                edges_to_remove.push_back(*it);
            }

            // Delete edges
            for (unsigned int i = 0; i < edges_to_remove.size(); i++) {
                this->remove_edge(edges_to_remove.at(i));
            }

            // Update i2u_ and nodes_
            if (idx == i2u_.size() - 1) {
                i2u_.pop_back();
            }
            else {
                std::swap(i2u_.at(idx), i2u_.at(i2u_.size() - 1));
                i2u_.pop_back();

                // Modify swapped node
                size_type j_uid = i2u_.at(idx);
                nodes_.at(j_uid).idx_ = idx;
            }

            return 1;
        } else {
            return 0;
        }
    }

    /**
     * @brief Remove a node from the graph as specified by a node iterator. This
     *      removes the node's connecting edges as well.
     * @param[in] @n_it a NodeIterator object pointing to the Node to be deleted
     * @pre @n_it points to a valid node of the graph
     * @return a NodeIterator object pointing to the next Node in the graph
     * @post new i2u_.size() == old i2u_.size() - 1
     * @post new nodes_.size() == old nodes_.size()
     * @post num_edges() == old num_edges() - old (*n_it).degree()
     * @post for 0 <= i < graph.size(), nodes_[i2u_[i]].idx < graph.size()
     * @post new node_edge_dictionary_.at((*n_it).uid_).size() = 0
     * @post new edge_values_.at((*n_it).uid_).size() = 0
     * 
     * Complexity: O(n.degree())
     */
    node_iterator remove_node(node_iterator n_it) {
        Node n = *n_it;
        size_type idx = n.idx_;
        size_type node_removal = this->remove_node(n);
        return NodeIterator(this, idx);
    }

    /**
     * @brief Remove an edge defined by two adjacent nodes.
     * @param[in] @a, @b two Node objects whose connecting Edge we want to remove.
     * @pre @a and @b represent to valid nodes of the graph.
     * @return a size_type indicating if an Edge connecting @a and @b was removed (1)
     *      or not (0)
     * @post If 1 is returned, new num_edges() == old num_edges() - 1
     *       Else new num_edges() == old num_edges()
     * @post new edge_values_.at(a.uid_).at(b.uid_) = new edge_values_.at(b.uid_).at(a.uid_) = edge_value_type{}
     * @post b.uid_ not in new node_edge_dictionary_.at(a.uid_) and
     *          a.uid_ not in new node_edge_dictionary_.at(b.uid_)
     *
     * Complexity: O(num_edges() + a.degree() + b.degree()) < O(num_edges() + num_nodes())
     */
    size_type remove_edge(const Node& a, const Node& b) {
        // Verify edge existence
        int edge_loc = -1;
        for (unsigned int i = 0; i < edges_.size(); i++) {
            if ((edges_[i].first == a.uid_ && edges_[i].second == b.uid_) ||
                (edges_[i].first == b.uid_ && edges_[i].second == a.uid_)) {
                edge_loc = i;
            }
        }
           
        // Remove edge from containers if it is valid
        if (edge_loc == -1) { 
            return 0;
        } else {
            // Update edges_
            std::swap(edges_.at(edge_loc), edges_.at(edges_.size() - 1));
            edges_.pop_back();

            // Update edge_values_
            if (a.uid_ < b.uid_) {
                edge_values_.at(a.uid_).at(b.uid_) = edge_value_type{};
            }
            else {
                edge_values_.at(b.uid_).at(a.uid_) = edge_value_type{};
            }
            
            // Update node_edge_dictionary
            std::vector<size_type>& a_edges = node_edge_dictionary_.at(a.uid_);
            for (auto it = a_edges.begin(); it != a_edges.end(); ) {
                if (*it == b.uid_) {
                    it = a_edges.erase(it);
                }
                else {
                    ++it;
                }
            }

            std::vector<size_type>& b_edges = node_edge_dictionary_.at(b.uid_);
            for (auto it = b_edges.begin(); it != b_edges.end(); ) {
                if (*it == a.uid_) {
                    it = b_edges.erase(it);
                }
                else {
                    ++it;
                }
            }
            
            return 1;
        }
    }

    /**
     * @brief Remove an edge from the graph as specified by an Edge object
     * @param[in] @ed The Edge object representing the edge we want to remove from the graph.
     * @pre @ed is a valid edge of the graph
     * @return a size_type indicating if the edge represented by @ed was removed (1)
     *      or not (0)
     * @post If 1 is returned, new num_edges() == old num_edges() - 1
     *       Else new num_edges() == old num_edges()
     * @post new edge_values_.at(ed.node1().uid_).at(ed.node2().uid_) = 
                new edge_values_.at(ed.node2().uid_).at(ed.node1().uid_) = edge_value_type{}
     * @post ed.node2().uid_ not in new node_edge_dictionary_.at(ed.node1().uid_) and
     *          ed.node1().uid_ not in new node_edge_dictionary_.at(ed.node2().uid_)
     * 
     * Complexity: O(num_edges() + a.degree() + b.degree()) < O(num_edges() + num_nodes())
     */
    size_type remove_edge(const Edge& ed) {
        return remove_edge(ed.node1(), ed.node2());
    }

    /**
     * @brief Remove an edge from the graph as specified by an edge iterator
     * @param[in] @e_it an EdgeIterator pointing to the Edge representing the edge we want to remove
     * @pre @e_it is a valid EdgeIterator of the graph
     * @return an EdgeIterator pointing to the next Edge of the graph
     * @post new num_edges() == old num_edges() - 1
     * @post new edge_values_.at((*e_it).node1().uid_).at((*e_it).node2().uid_) =
                new edge_values_.at((*e_it).node2().uid_).at((*e_it).node1().uid_) = edge_value_type{}
     * @post (*e_it).node2().uid_ not in new node_edge_dictionary_.at((*e_it).node1().uid_) and
     *          (*e_it).node1().uid_ not in new node_edge_dictionary_.at((*e_it).node2().uid_)
     *
     * Complexity: O(num_edges() + a.degree() + b.degree()) < O(num_edges() + num_nodes())
     */
    edge_iterator remove_edge(edge_iterator e_it) {
        Edge ed = *e_it;
        size_type ed_idx = e_it.edge_i_;
        this->remove_edge(ed);
        return EdgeIterator(this, ed_idx);
    }

private:
    /* Define the internal_node struct which contains information for each node */
    struct internal_node {
        Point point;
        V value_;
        size_type idx_;
    };

    /* Node containers */
    // Used to store internal_node for any node added, even if later removed.
    std::vector<internal_node> nodes_;

    // Store the currently 'active' set of nodes
    std::vector<size_type> i2u_;

    /* Graph size information */
    size_type next_uid_;

    /* Edge containers */
    std::vector<std::pair<size_type, size_type>> edges_;
    std::unordered_map<size_type, std::unordered_map<size_type, edge_value_type>> edge_values_;
    std::unordered_map<size_type, std::vector<size_type>> node_edge_dictionary_;

};

#endif // CME212_GRAPH_HPP

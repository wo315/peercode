#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/**
 * @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/**
 * @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
   private:
    // Private node representation
    struct internal_node;
    struct internal_edge;

   public:
    /** Type of this graph. */
    using graph_type = Graph;

    /** Predeclaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;
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

    /**
     * Type of indexes and sizes.
     * Return type of Graph::Node::index(), Graph::num_nodes(),
     * Graph::num_edges(), and argument type of Graph::node(size_type)
     */
    using size_type = unsigned;
    using edge_id_type = size_type;
    using node_id_type = size_type;

    /** Construct an empty graph. */
    Graph()
        : next_uid_(0),
          node_vec_(),
          node_uid_map_(),
          edge_vec_(),
          edge_uid_map_(),
          node_to_edges_() {}

    /** Default destructor */
    ~Graph() = default;

    /**
     * @class Graph::Node
     * @brief Class representing the graph's nodes.
     *
     * Node objects are used to access information about the Graph's nodes.
     */
    class Node : private totally_ordered<Node> {
       public:
        /** Useful for testing/hashing purposes. */
        node_id_type get_uid() const { return this->uid_; }

        /**
         * Construct an invalid node.
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
        Node() : graph_(nullptr), uid_(-1) {}

        /** Return this node's position. */
        const Point& position() const {
            return graph_->node_uid_map_.at(uid_).point;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            for (size_type i = 0; i < graph_->size(); i++) {
                if (graph_->node_vec_.at(i).uid == uid_) {
                    return i;
                }
            }
            assert(false);
        }

        /** Read only access to a Node's @a value */
        const node_value_type& value() const {
            const internal_node& node = graph_->node_uid_map_.at(uid_);
            return node.value;
        }

        /** Write (and read) access to a Node's @a value */
        node_value_type& value() {
            internal_node& node = graph_->node_uid_map_.at(uid_);
            return node.value;
        }

        /**
         * Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            return graph_ == n.graph_ && index() == n.index();
        }

        /**
         * Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node& n) const {
            if (graph_ == n.graph_) {
                return index() < n.index();
            }
            return graph_ < n.graph_;
        }

        /** Number of edges incident to a Node */
        size_type degree() const {
            return graph_->node_to_edges_.at(uid_).size();
        }

        /** Start an iterator of incident edges */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph_, uid_, 0);
        }

        /** Last iterator in the sequence. Useful for equality checks */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_, uid_, degree());
        }

       private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;

        Graph* graph_;
        node_id_type uid_;

        Node(const Graph* graph, const node_id_type uid)
            : graph_(const_cast<Graph* const>(graph)), uid_(uid) {}
    };

    /**
     * Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const { return this->node_vec_.size(); }

    /** Synonym for size(). */
    size_type num_nodes() const { return this->size(); }

    /**
     * Add a node to the graph, returning the added node.
     * @param[in] position The new node's position
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point position,
                  const node_value_type& value = node_value_type()) {
        internal_node node = {next_uid_, position, value};
        node_vec_.push_back(node);
        node_uid_map_.insert({next_uid_, node});
        node_to_edges_.insert({next_uid_, std::vector<internal_edge>()});
        next_uid_++;
        return Node(this, next_uid_ - 1);
    }

    /**
     * Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        auto search = node_uid_map_.find(n.uid_);
        return search != node_uid_map_.end();
    }

    /**
     * Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(const size_type i) const {
        internal_node node = node_vec_.at(i);
        return Node(this, node.uid);
    }

    /**
     * @class Graph::Edge
     * @brief Class representing the graph's edges.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge> {
       public:
        /** Construct an invalid Edge. */
        Edge() : graph_(nullptr), edge_uid_(-1), m_one_(-1), m_two_(-1) {}

        /** Return a node of this Edge */
        Node node1() const { return get_node_by_uid(m_one_); }

        /** Return the other node of this Edge */
        Node node2() const { return get_node_by_uid(m_two_); }

        /**
         * Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            return (this->m_one_ == e.m_one_ && this->m_two_ == e.m_two_) ||
                   (this->m_one_ == e.m_two_ && this->m_two_ == e.m_one_);
        }

        /**
         * Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            return this->edge_uid_ < e.edge_uid_;
        }

       private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;

        Graph* graph_;
        edge_id_type edge_uid_;
        node_id_type m_one_;
        node_id_type m_two_;

        /**
         * Privately construct an Edge by specifying an edge id, and two node
         * ids
         */
        Edge(const Graph* graph, const edge_id_type uid, const node_id_type one,
             const node_id_type two)
            : graph_(const_cast<Graph* const>(graph)),
              edge_uid_(uid),
              m_one_(one),
              m_two_(two) {}

        /** Construct a node with the corresponding uid */
        Node get_node_by_uid(const node_id_type uid) const {
            internal_node node = graph_->node_uid_map_.at(uid);
            return Node(graph_, node.uid);
        }
    };

    /**
     * Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const { return edge_vec_.size(); }

    /**
     * Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(const size_type i) const {
        edge_id_type uid = this->edge_vec_.at(i);
        auto it = this->edge_uid_map_.find(uid);
        internal_edge edge = it->second;
        return Edge(this, uid, edge.one, edge.two);
    }

    /**
     * Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        std::vector<internal_edge> edges = node_to_edges_.at(a.uid_);
        for (auto const& edge : edges) {
            if (b.uid_ == edge.one) {
                assert(a.uid_ == edge.two);
                return true;
            }
            if (b.uid_ == edge.two) {
                assert(a.uid_ == edge.one);
                return true;
            }
        }
        return false;
    }

    /**
     * Add an edge to the graph, or return the current edge if it already
     * exists.
     * @pre @a a and @a b are distinct valid nodes of this graph
     * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
     * @post has_edge(@a a, @a b) == true
     * @post If old has_edge(@a a, @a b):
     *          new num_edges() == old num_edges().
     *       else:
     *          new num_edges() == old num_edges() + 1.
     *
     * Can invalidate edge indexes -- in other words, old edge(@a i) might not
     * equal new edge(@a i). Must not invalidate outstanding Edge objects.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge add_edge(const Node& a, const Node& b) {
        // this essentially reimplements 'has_edge' except also grabs edge uid
        std::vector<internal_edge> edges = node_to_edges_.at(a.uid_);
        for (auto const& e : edges) {
            if (b.uid_ == e.one || b.uid_ == e.two) {
                assert(a.uid_ == e.two || a.uid_ == e.one);
                return Edge(this, e.uid, a.uid_, b.uid_);
            }
        }
        // If we reach here, we are adding a new edge
        edge_id_type edge_uid = next_uid_;
        next_uid_++;

        edge_vec_.push_back(edge_uid);
        internal_edge edge = {edge_uid, a.uid_, b.uid_};
        edge_uid_map_.insert({edge_uid, edge});

        // add the edge to BOTH node lists
        assert(node_to_edges_.find(a.uid_) != node_to_edges_.end());
        node_to_edges_.at(a.uid_).push_back(edge);
        assert(node_to_edges_.find(b.uid_) != node_to_edges_.end());
        node_to_edges_.at(b.uid_).push_back(edge);

        return Edge(this, edge_uid, a.uid_, b.uid_);
    }

    /**
     * Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        node_vec_.clear();
        node_uid_map_.clear();
        edge_vec_.clear();
        edge_uid_map_.clear();
        node_to_edges_.clear();
    }

    /**
     * @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator.
     */
    class NodeIterator : private totally_ordered<NodeIterator> {
       public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Node;                 // Element type
        using pointer = Node*;                   // Pointers to elements
        using reference = Node&;                 // Reference to elements
        using difference_type = std::ptrdiff_t;  // Signed difference
        using iterator_category =
            std::forward_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid NodeIterator. */
        NodeIterator() : m_graph(nullptr), m_i(-1) {}

        Node operator*() const { return m_graph->node(m_i); }

        /** Advance the iterator to the next Node in the graph */
        NodeIterator& operator++() {
            m_i++;
            return *this;
        }

        /**
         * @brief Determine if two iterators are equal
         * 
         * Equality is determined by graph equality AND
         * node index equality.
         */
        bool operator==(const NodeIterator& other) const {
            return m_graph == other.m_graph && m_i == other.m_i;
        }

       private:
        friend class Graph;
        friend class EdgeIterator;

        Graph* m_graph;
        size_type m_i;

        NodeIterator(const Graph* graph, size_type i)
            : m_graph(const_cast<Graph*>(graph)), m_i(i) {}
    };

    node_iterator node_begin() const { return NodeIterator(this, 0); }
    node_iterator node_end() const { return NodeIterator(this, size()); }

    /**
     * @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator.
     */
    class IncidentIterator : private totally_ordered<IncidentIterator> {
       public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Edge;                 // Element type
        using pointer = Edge*;                   // Pointers to elements
        using reference = Edge&;                 // Reference to elements
        using difference_type = std::ptrdiff_t;  // Signed difference
        using iterator_category =
            std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator() : m_graph(nullptr), m_node(-1), m_j(-1) {}

        /**
         * @brief Dereference the iterator to retrieve an Edge.
         * 
         * IncidentIterator's are associated with a defining Node,
         * so when we dereference the associated node will always be
         * the result of edge.node1() of the returned edge.
         */
        Edge operator*() const {
            internal_edge edge = m_graph->node_to_edges_.at(m_node).at(m_j);
            node_id_type other = m_node == edge.one ? edge.two : edge.one;
            return Edge(m_graph, edge.uid, m_node, other);
        }

        IncidentIterator& operator++() {
            m_j++;
            return *this;
        }

        bool operator==(const IncidentIterator& other) const {
            return m_graph == other.m_graph && m_node == other.m_node &&
                   m_j == other.m_j;
        }

       private:
        friend class Graph;

        Graph* m_graph;
        // Making this const would break EdgeIterator since it erases assignment
        node_id_type m_node;
        size_type m_j;

        IncidentIterator(const Graph* graph, const node_id_type node,
                         size_type index)
            : m_graph(const_cast<Graph*>(graph)), m_node(node), m_j(index) {}
    };

    /**
     * @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator.
     * 
     * An @a EdgeIterator provides a wrapper around a @a Graph::NodeIterator
     * and a @a Graph::IncidentIterator.
     */
    class EdgeIterator : private totally_ordered<EdgeIterator> {
       public:
        // These type definitions let us use STL's iterator_traits.
        using value_type = Edge;                 // Element type
        using pointer = Edge*;                   // Pointers to elements
        using reference = Edge&;                 // Reference to elements
        using difference_type = std::ptrdiff_t;  // Signed difference
        using iterator_category =
            std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid EdgeIterator. */
        EdgeIterator() {}

        Edge operator*() const { return *m_incident_it; }

        /**
         * @brief Move the iterator to the next Edge.
         * 
         * We have to be careful that we don't advance the 
         * underlying @a NodeIterator to an invalid index since
         * we dereference it. We also have to make sure that we don't
         * repeat edges so there are two conditions that must be met.
         */
        EdgeIterator& operator++() {
            // always start by advancing the incident iterator
            ++m_incident_it;

            // we keep moving onward until both conditions are false
            while (m_incident_it == (*m_node_it).edge_end() ||
                   m_seen.find((*m_incident_it).edge_uid_) != m_seen.end()) {
                // first we make sure our incident_iterator is indexable
                while (m_incident_it == (*m_node_it).edge_end()) {
                    // by advancing the node
                    ++m_node_it;
                    if (m_node_it == m_graph->node_end()) {
                        // and returning a sentinel if we're at the end
                        m_incident_it = IncidentIterator();
                        return *this;
                    }
                    // otherwise we reset the IncidentIterator to zero
                    node_id_type node_id = (*m_node_it).uid_;
                    m_incident_it = IncidentIterator(m_graph, node_id, 0);
                }
                // but we also can't return seen edges
                if (m_seen.find((*m_incident_it).edge_uid_) != m_seen.end()) {
                    ++m_incident_it;
                }
            }
            // mark the newest edge id as seen
            m_seen.insert((*m_incident_it).edge_uid_);
            return *this;
        }

        bool operator==(const EdgeIterator& other) const {
            return m_graph == other.m_graph && m_node_it == other.m_node_it &&
                   m_incident_it == other.m_incident_it;
        }

       private:
        friend class Graph;

        Graph* const m_graph;
        node_iterator m_node_it;
        incident_iterator m_incident_it;
        std::unordered_set<edge_id_type> m_seen;

        EdgeIterator(const Graph* graph, size_type node_index,
                     size_type edge_index)
            : m_graph(const_cast<Graph* const>(graph)),
              m_node_it(graph, node_index),
              m_incident_it(graph, (*m_node_it).uid_, edge_index),
              m_seen{(*m_incident_it).edge_uid_} {}

        EdgeIterator(const Graph* graph, size_type node_index)
            : m_graph(const_cast<Graph* const>(graph)),
              m_node_it(graph, node_index),
              m_incident_it(),
              m_seen() {}
    };

    edge_iterator edge_begin() const { return EdgeIterator(this, 0, 0); }
    edge_iterator edge_end() const { return EdgeIterator(this, num_nodes()); }

   private:
    struct internal_node {
        const node_id_type uid;
        const Point point;
        node_value_type value;
    };

    struct internal_edge {
        const edge_id_type uid;
        const node_id_type one;
        const node_id_type two;
    };

    // is used for both node and edge ids, so they are
    // indeed unique
    size_type next_uid_;
    std::vector<internal_node> node_vec_;
    std::unordered_map<node_id_type, internal_node> node_uid_map_;

    std::vector<edge_id_type> edge_vec_;
    std::unordered_map<edge_id_type, internal_edge> edge_uid_map_;

    // maps node ids to a set of incident edge ids. both nodes
    // get mapped to the edge, so naively iterating through this
    // would yield duplicate edges
    std::unordered_map<node_id_type, std::vector<internal_edge>> node_to_edges_;
};

#endif  // CME212_GRAPH_HPP

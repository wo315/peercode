#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

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

    struct node_p;
    struct edge_p;

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
    Graph()
            :vec_nodes_(), size_nodes_(0), next_uidNode_(0), vec_edges_(),
             size_edges_(0),next_uidEdge_(0){
        // HW0: YOUR CODE HERE
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
            // HW0: YOUR CODE HERE
        }

        /** Return this node's position. */
        const Point& position() const {
            // HW0: YOUR CODE HERE
            return fetch().point_;
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return fetch().uid_;
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
            return ((this->graph_ == n.graph_) && (this->uidNode_ == n.uidNode_) ) ? true : false ;
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
            return ((this->graph_ == n.graph_) &&  this->uidNode_ < n.uidNode_) ? true : false;

        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects

        // Pointer back to the Graph container
        graph_type * graph_;
        // This node's unique identification number
        size_type uidNode_;
        /** Private Constructor */
        Node(const graph_type* graph, size_type uidNode)
                : graph_(const_cast<graph_type*>(graph)), uidNode_(uidNode) {
        };

        node_p& fetch() const {
            return graph_->vec_nodes_[uidNode_];
        }

    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        // HW0: YOUR CODE HERE
        return size_nodes_;
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
        vec_nodes_.emplace_back(next_uidNode_,position);
        ++size_nodes_;
        ++next_uidNode_;
        // Returns a Node that points to the new element
        return Node(this, next_uidNode_-1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE

        // added this ckeck just in case, but shouldn't be necessary
        bool b_larger_graph = (n.uidNode_ < this->size());

        return ((this == n.graph_) && b_larger_graph) ? true : false;
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(i < size());
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
    class Edge {
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HW0: YOUR CODE HERE
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(fetch().node1_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(fetch().node2_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            //HW0: YOUR CODE HERE
            bool same_graph = (this->graph_ == e.graph_);
            bool same_start = (node1() == e.node1()) ;
            bool same_end = (node2() == e.node2()) ;

            return ( same_graph && same_start && same_end ) ? true : false ;

        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            //HW0: YOUR CODE HERE
            return (this->uidEdge_ < e.uidEdge_) ? true : false;
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Edge objects

        // Pointer back to the Graph container

        graph_type * graph_;
        // This node's unique identification number
        size_type uidEdge_;
        /** Private Constructor */
        Edge(const graph_type* graph, size_type uidEdge)
                : graph_(const_cast<graph_type*>(graph)), uidEdge_(uidEdge) {
        };

        edge_p& fetch() const {
            return graph_->vec_edges_[uidEdge_];
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        return size_edges_;
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        assert(i < num_edges());
        return Edge(this, i);
    }


    /** sub function to test whether to nodes are connected by an edge (called by has_edge)
    * @pre @a a and @a b are valid nodes of this graph, with a.id < b.id
    * @return True if for some @a i, edge(@a i) connects @a a and @a b.
    *
    */
    int index_has_edge(const Node& a, const Node& b) const {
        for(unsigned i =0; i< this->vec_edges_.size(); ++i){
            Edge cur_edge = Edge(this,i);
            if((cur_edge.node1() == a) && (cur_edge.node2() == b)){
                return (int) i;
            }
        }
        return -1;
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        if(b<a){
            return (index_has_edge(b,a)>-1);
        }else{
            return (index_has_edge(a,b)>-1);
        }
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

        /*If we can modify the type returned by has_edge we could make this faster */

        if(has_edge(a,b)){
            if(b<a){
                size_type idx_edge = (unsigned int) index_has_edge(b,a);
                return Edge(this,idx_edge);
            }else{
                size_type idx_edge = (unsigned int) index_has_edge(a,b);
                return Edge(this,idx_edge);
            }
        }else{
            //(a < b) ? vec_edges_.emplace_back(next_uidEdge_,a,b) : vec_edges_.emplace_back(next_uidEdge_,b,a);
            if(a<b){
                vec_edges_.emplace_back(next_uidEdge_,a.index(),b.index());
            }else{
                vec_edges_.emplace_back(next_uidEdge_,b.index(),a.index());
            }
            ++size_edges_;
            ++next_uidEdge_;
            // Returns a Node that points to the new element
            return Edge(this, next_uidEdge_-1);
        }

        // this last line should never be called
        return Edge(this,num_edges()-1);
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        vec_nodes_   = std::vector<node_p>();
        size_nodes_ = 0 ;
        next_uidNode_ = 0;
        vec_edges_  = std::vector<edge_p>();
        size_edges_ = 0 ;
        next_uidEdge_ = 0;
    }

    //
    // Node Iterator
    //

    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Node;                     // Element type
        using pointer           = Node*;                    // Pointers to elements
        using reference         = Node&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
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
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
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
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
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

    // Internal type for graph nodes
    struct node_p {
        size_type uid_;
        Point point_;

        node_p(size_type uid , const Point& point)
                :uid_(uid),point_(point) {}
    };
    // Internal type for graph edges
    struct edge_p {
        size_type uid_;
        const size_type node1_;
        const size_type node2_;

        edge_p(size_type uid , const size_type node1,const size_type node2)
                :uid_(uid),node1_(node1),node2_(node2) {}
    };

    std::vector<node_p> vec_nodes_;
    size_type size_nodes_;
    size_type next_uidNode_;
    std::vector<edge_p> vec_edges_;
    size_type size_edges_;
    size_type next_uidEdge_;

    // Disable copy and assignment of a SimpleSet
    Graph(const Graph&) = delete;
    Graph& operator=(const Graph&) = delete;

};

#endif // CME212_GRAPH_HPP

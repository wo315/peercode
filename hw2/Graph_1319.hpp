#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
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
    /** Predeclaration of Node value type. */
    using node_value_type = V;

    /** Predeclaration of Edge type. */
    class Edge;
    /** Synonym for Edge (following STL conventions). */
    using edge_type = Edge;
    /** Predeclaration of Edge value type. */
    using edge_value_type = E;

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
            // HW0: YOUR CODE HERE
            this->graph_ = nullptr;
            this->idx_ = -1;
        }

        /** Return this node's position. */
        Point& position() {
            // HW0: YOUR CODE HERE
            return this->graph_->Points_[this->idx_];
        }

        /** Return this node's position. */
        const Point& position() const{
            return this->graph_->Points_[this->idx_];
        }

        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            // HW0: YOUR CODE HERE
            return this->graph_->nu2i_[this->idx_];
        }

        // HW1: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Return this node's value. */
        node_value_type& value(){
            return this->graph_->Values_[this->idx_];
        }
        const node_value_type& value() const{
            return this->graph_->Values_[this->idx_];
        }

        /** return the number of indicent edges. */
        size_type degree() const {
            if(!(this->graph_->AdjLst_.size() > this->idx_))
                return 0;
            return this->graph_->AdjLst_[idx_].size();
        }

        /** Start of the incident iterator. */
        incident_iterator edge_begin() const {
            return IncidentIterator(graph_, idx_, 0);
        }

        /** End of incident iterator. */
        incident_iterator edge_end() const {
            return IncidentIterator(graph_, idx_, degree());
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            return (this->graph_ == n.graph_) and (this->idx_ == n.idx_);
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
            return ((graph_ == n.graph_) and idx_ < n.idx_)
                or (long long)graph_ < (long long)n.graph_;
        }

    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Node
        // that will not be visible to users, but may be useful within Graph.
        // i.e. Graph needs a way to construct valid Node objects

        // Pointer back to the Graph container
        Graph* graph_;
        // This element's unique identification number
        size_type idx_;
        // Constructor given position and index
        Node(const Graph* graph, size_type index)
            :graph_(const_cast<Graph*>(graph)), idx_(index){
        }
        bool valid() const {
            return idx_ < graph_->nu2i_.size() and graph_->nu2i_[idx_] < graph_->ni2u_.size()
                and graph_->ni2u_[graph_->nu2i_[idx_]] == idx_;
        }
    };

    /** Return the number of nodes in the graph.
     *
     * Complexity: O(1).
     */
    size_type size() const {
        // HW0: YOUR CODE HERE
        return ni2u_.size();
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
    Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
        // HW0: YOUR CODE HERE
        Points_.push_back(position);
        Values_.push_back(val);
        nu2i_.push_back(ni2u_.size());
        ni2u_.push_back(nu2i_.size() - 1);
        return Node(this, nu2i_.size() - 1);
    }

    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        // HW0: YOUR CODE HERE
        return n.graph_ == this and n.valid();
    }

    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        // HW0: YOUR CODE HERE
        return Node(this, ni2u_[i]);        // Invalid node
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
    class Edge : totally_ordered<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HW0: YOUR CODE HERE
            graph_ = nullptr;
            nidx1_ = nidx2_ = -1;
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(this->graph_->nu2i_[nidx1_]);      // Invalid Node
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(this->graph_->nu2i_[nidx2_]);      // Invalid Node
        }

        /** Return the length of this Edge */
        double length() const {
            return norm(node1().position() - node2().position());
        }

        /** Return modifiable edge value reference */
        edge_value_type& value(){
            return graph_->EdgeValues_[graph_->Edgemap_[*this]];
        }

        /** Return constant edge value reference */
        const edge_value_type& value() const{
            return graph_->EdgeValues_[graph_->Edgemap_[*this]];
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            //HW0: YOUR CODE HERE
            return (graph_ == e.graph_) and ((nidx1_ == e.nidx1_ and nidx2_ == e.nidx2_)
                                             or (nidx1_ == e.nidx2_ and nidx2_ == e.nidx1_));
        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            //HW0: YOUR CODE HERE
            size_type a, b, c, d;
            if((long long)graph_ < (long long)e.graph_)
                return true;
            if(nidx1_ < nidx2_)     { a = nidx1_; b = nidx2_; }
            else                    { a = nidx2_; b = nidx1_; }
            if(e.nidx1_ < e.nidx2_) { c = e.nidx1_; d = e.nidx2_; }
            else                    { c = e.nidx2_; d = e.nidx1_; }
            return a < c or (a == c and b < d);
        }

    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        // HW0: YOUR CODE HERE
        // Use this space to declare private data members and methods for Edge
        // that will not be visible to users, but may be useful within Graph.p
        // i.e. Graph needs a way to construct valid Edge objects


        // sizeof(Graph::Edge) = 16
        // two indices
        size_type nidx1_, nidx2_;
        // pointer to Graph
        Graph* graph_;
        // Constructor
        Edge(const Node& a, const Node& b){
            graph_ = a.graph_;
            nidx1_ = a.idx_;
            nidx2_ = b.idx_;
        }
    };

    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        // HW0: YOUR CODE HERE
        // O(1)
        return ei2u_.size();
    }

    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        // HW0: YOUR CODE HERE
        // O(1)
        return Edgevec_[ei2u_[i]];
    }

    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        // HW0: YOUR CODE HERE
        // O(log(num_edges()))
        Edge e(a, b);
        return Edgemap_.find(e) != Edgemap_.end();
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
    Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
        // HW0: YOUR CODE HERE
        // O(log(num_edges()))
        Edge e(a, b);
        size_type maxidx = std::max(a.idx_, b.idx_);
        while(AdjLst_.size() < maxidx + 1)
            AdjLst_.push_back(std::vector<size_type>());
        if(!has_edge(a, b)){
            AdjLst_[a.idx_].push_back(b.idx_);
            AdjLst_[b.idx_].push_back(a.idx_);
            Edgemap_[e] = eu2i_.size();
            Edgevec_.push_back(e);
            EdgeValues_.push_back(val);
            eu2i_.push_back(ei2u_.size());
            ei2u_.push_back(eu2i_.size() - 1);
        }
        return e;
    }

    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        // HW0: YOUR CODE HERE
        Points_.clear();
        Values_.clear();
        Edgevec_.clear();
        Edgemap_.clear();
        AdjLst_.clear();
        nu2i_.clear();
        ni2u_.clear();
        eu2i_.clear();
        ei2u_.clear();
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

        /** Dereference operator */
        Node operator*() const {
            return Node(graph_, graph_->ni2u_[idx_]);
        }

        /** Increment opreator */
        NodeIterator& operator++() {
            idx_++;
            return *(this);
        }

        /** Check equality */
        bool operator==(const NodeIterator& iter) const {
            return graph_ == iter.graph_ and
                idx_ == iter.idx_;
        }

        /** Check inequality */
        bool operator!=(const NodeIterator& iter) const {
            return graph_ != iter.graph_ or
                idx_ != iter.idx_;
        }

    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        Graph* graph_;
        size_type idx_;
        NodeIterator(const Graph* graph, size_type index):
            graph_(const_cast<Graph*>(graph)), idx_(index){}
    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Return a NodeIterator pointing to the first Node. */
    node_iterator node_begin() const {
        return NodeIterator(this, 0);
    }
    /** Return a NodeIterator pointing to one past the last Node. */
    node_iterator node_end() const {
        return NodeIterator(this, num_nodes());
    }

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
        IncidentIterator() : graph_(nullptr), idx1_(-1), idx2_(-1) {
        }

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:

        /** Dereference operator. */
        Edge operator*() const {
            return Edge(Node(graph_, idx1_), Node(graph_, graph_->AdjLst_[idx1_][idx2_]));
        }
        IncidentIterator& operator++(){
            idx2_++;
            return *(this);
        }
        bool operator==(const IncidentIterator& iter) const {
            return graph_ == iter.graph_ and idx1_ == iter.idx1_ and
                idx2_ == iter.idx2_;
        }
        bool operator!=(const IncidentIterator& iter) const {
            return graph_ != iter.graph_ or idx1_ != iter.idx1_ or
                idx2_ != iter.idx2_;
        }

    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
        Graph* graph_;
        const size_type idx1_;
        size_type idx2_;
        IncidentIterator(const Graph* graph, size_type idx1, size_type idx2) :
            graph_(const_cast<Graph*>(graph)), idx1_(idx1), idx2_(idx2){}
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

        /** Dereference. */
        Edge operator*() const {
            return graph_->Edgevec_[graph_->ei2u_[idx_]];
        }

        /** Increment operator. */
        EdgeIterator& operator++() {
            idx_++;
            return *(this);
        }

        /** Check equality. */
        bool operator==(const EdgeIterator& iter) const {
            return graph_ == iter.graph_ and idx_ == iter.idx_;
        }

        /** Check inequality. */
        bool operator!=(const EdgeIterator& iter) const {
            return graph_ != iter.graph_ or idx_ != iter.idx_;
        }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        Graph* graph_;
        size_type idx_;
        EdgeIterator(const Graph* graph, size_type idx):
            graph_(const_cast<Graph*>(graph)), idx_(idx) {}
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Start of edges in Graph. */
    edge_iterator edge_begin() const {
        return EdgeIterator(this, 0);
    }

    /** End of edges in Graph. */
    edge_iterator edge_end() const {
        return EdgeIterator(this, num_edges());
    }

    /** Remove node and its incident edges.
     *  @param[in] n    Node
     *  @return         0 (Node is not valid) or 1 (Node is valid and successfully removed
     *  @post           graph.has_node(n) == false
     *  @post           no edge in graph is incident to n
     *  @post           graph.num_nodes is reduced by at most 1
     *  Complexity      O(max_degree * (max_degree + log(num_edges)))
     */
    size_type remove_node(const Node& n){
        // check if Node is in the graph
        if(!has_node(n))
            return 0;
        // remove one edge at a time
        while(n.degree() != 0){
            auto iter = n.edge_begin();
            Edge e = *iter;
            remove_edge(e);
        }
        // update ni2u_ and nu2i_
        size_type uid = ni2u_.back();
        ni2u_[nu2i_[n.idx_]] = uid;
        ni2u_.pop_back();
        nu2i_[uid] = nu2i_[n.idx_];
        return 1;
    }

    /** Override remove_node
     *  @param[in]    n_it    node_iterator pointing to a Node
     *  @return       node_iterator pointing to the beginning of the graph
     *  @post         the same as remove_node(const Node&) above
     *  Complexity    the same as remove_node(const Node&) above
     */
    //--functionality_1
    //--node_begin() can be a node that is already removed
    //--START
    node_iterator remove_node(node_iterator n_it){
        Node n = *n_it;
        remove_node(n);
        return node_begin();
    }
    //--END

    /** Remove an edge
     *  @param[in]    e    Edge
     *  @return       0 (Edge is not in the graph) or 1 (Edge is in the graph and successfully removed)
     *  @post         graph.has_edge(e) == false
     *  @post         graph.num_edges() is reduced by at most 1
     *  Complexity    O(log(num_edges) + max_degree)
     */
    size_type remove_edge(const Edge& e){
        // check if the edge is in the graph
        if(!has_edge(e.node1(), e.node2()))
            return 0;
        // get uid of edge
        size_type uid1 = Edgemap_[e];
        // erase edge information from Edgemap_
        Edgemap_.erase(e);
        // erase edge information from AdjLst_
        for(auto iter = AdjLst_[e.nidx1_].begin(); iter != AdjLst_[e.nidx1_].end(); ++iter){
            if((*iter) == e.nidx2_){
                AdjLst_[e.nidx1_].erase(iter);
                break;
            }
        }
        for(auto iter = AdjLst_[e.nidx2_].begin(); iter != AdjLst_[e.nidx2_].end(); ++iter){
            if((*iter) == e.nidx1_){
                AdjLst_[e.nidx2_].erase(iter);
                break;
            }
        }
        // update eu2i_ and ei2u_
        size_type uid2 = ei2u_.back();
        ei2u_[eu2i_[uid1]] = uid2;
        ei2u_.pop_back();
        eu2i_[uid2] = eu2i_[uid1];
        return 1;
    }

    /** Override remove_edge
     *  @param[in]    n1    Node
     *  @param[in]    n2    Node
     *  @return       the same as remove_edge(const Edge&)
     *  @post         the same as remove_edge(const Edge&)
     *  Complexity    the same as remove_edge(const Edge&)
     */
    size_type remove_edge(const Node& n1, const Node n2){
        Edge e(n1, n2);
        return remove_edge(e);
    }

    /** Override remove_edge
     *  @param[in]    e_it    edge_iterator pointing to an Edge
     *  @return       edge_iterator pointing to the beginning of the graph
     *  @post         the same as remove_edge(const Edge&)
     *  Complexity    the same as remove_edge(const Edge&)
     */
    edge_iterator remove_edge(edge_iterator e_it){
        Edge e = *e_it;
        remove_edge(e);
        return edge_begin();
    }

private:

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.

    // Point information
    std::vector<Point> Points_;
    // Value information
    std::vector<node_value_type> Values_;
    // map node uid to idx
    std::vector<size_type> nu2i_;
    // map node idx to uid
    std::vector<size_type> ni2u_;

    // Two containers to store edges information.
    // map for log(n) find (map Edge to index)
    std::map<Edge, size_type> Edgemap_;
    // vector for constant indexing (map index to Edge)
    std::vector<Edge> Edgevec_;
    // Edge information
    std::vector<edge_value_type> EdgeValues_;
    // adjacency list
    std::vector<std::vector<size_type>> AdjLst_;
    // map edge uid to idx
    std::vector<size_type> eu2i_;
    // map edge idx to uid
    std::vector<size_type> ei2u_;
};

#endif // CME212_GRAPH_HPP

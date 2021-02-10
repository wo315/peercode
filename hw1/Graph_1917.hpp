#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include<set>

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

    /** Type of nodes value
     */
    using node_value_type = V;

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
    class Node: private totally_ordered<Node> {
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

        /** Return this node's value, a real number (I don't know whether or not the
         * values can be smaller than 0). */
        node_value_type& value(){
            return fetch().value_;
        }

        /** Same method as above except it can be called on const object and you can only
        * read the value.
        * Return this node's value, a real number (I don't know whether or not the
        * values can be smaller than 0). */
        const node_value_type& value() const{
            return fetch().value_;
        }

        /** Return this node's degree, a number in the range [0, graph_size).
         * This represent the number of incident edges for this node*/
        size_type degree() const{
            return fetch().incident_edges_.size();
        }

        /** Iterator pointing at the start of the incident edges container of the current node
        * @pre Should be called if graph initialized
        * @return an Iterator over incident_edge of the cur node  (std) starting at the beginning of it
        */
        incident_iterator edge_begin() const
        {
            typename std::set<size_type>::const_iterator begin = fetch().incident_edges_.begin();
            return IncidentIterator(this->graph_,this->index(),begin);
        }

        /** Iterator pointing at the end of the incident edges container of the current node.
         * @pre Should be called if graph initialized
         * @post Used for stl function over iterators to end call
         * @return an Iterator over incident_edge of the cur node (std) pointing at the end
        */
        incident_iterator edge_end() const
        {
            typename std::set<size_type>::const_iterator end = fetch().incident_edges_.end();
            return IncidentIterator(this->graph_,this->index(),end);
        }

        /** Test whether this node and @a n are equal.
         *
         * Equal nodes have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            // HW0: YOUR CODE HERE
            return ((this->graph_ == n.graph_) && (this->uidNode_ == n.uidNode_));
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
            assert(this->graph_ == n.graph_);
            return (this->uidNode_ < n.uidNode_);

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
     * @param[in] value The new node's value
     * @post new num_nodes() == old num_nodes() + 1
     * @post result_node.index() == old num_nodes()
     *
     * Complexity: O(1) amortized operations.
     */
    Node add_node(const Point& position, const node_value_type& node_val= node_value_type()){
        // HW0: YOUR CODE HERE
        std::set<size_type> empty_incident =  std::set<size_type>();
        vec_nodes_.emplace_back(next_uidNode_,position,node_val,empty_incident);
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

        return ((this == n.graph_) && b_larger_graph);
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
    class Edge: private totally_ordered<Edge> {
    public:
        /** Construct an invalid Edge. */
        Edge() {
            // HW0: YOUR CODE HERE
        }

        /** Return a node of this Edge */
        Node node1() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(uidNode1_);
        }

        /** Return the other node of this Edge */
        Node node2() const {
            // HW0: YOUR CODE HERE
            return this->graph_->node(uidNode2_);
        }

        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            //HW0: YOUR CODE HERE
            bool same_graph = (this->graph_ == e.graph_);
            bool direction_1 = (node1() == e.node1() && node2() == e.node2()) ;
            bool direction_2 = (node1() == e.node2() && node2() == e.node1()) ;

            return (same_graph && (direction_1 || direction_2));

        }

        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            //HW0: YOUR CODE HERE
            assert(this->graph_ == e.graph_);
            return (this->uidEdge_ < e.uidEdge_);
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
        // The edge's first node
        size_type uidNode1_;
        // The edge's second node
        size_type uidNode2_;
        /** Private Constructor */
        Edge(const graph_type* graph, size_type uidEdge, size_type uidNode1, size_type uidNode2)
                : graph_(const_cast<graph_type*>(graph)), uidEdge_(uidEdge),
                  uidNode1_(uidNode1), uidNode2_(uidNode2){
        };

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
        size_type idN1 = vec_edges_[i].node1_;
        size_type idN2 = vec_edges_[i].node2_;

        return Edge(this, i, idN1, idN2);
    }


    /** sub function to test whether two nodes are connected by an edge (called by has_edge/add_edge)
    * @pre @a a and @a b are valid nodes of this graph
    * @return integer i if for some @a i, edge(@a i) connects @a a and @a b.
    *
    */
    int index_has_edge(const Node& a, const Node& b) const {
        std::set<size_type> inc_ed_a = vec_nodes_[a.index()].incident_edges_;
        std::set<size_type>::iterator it;
        for (it = inc_ed_a.begin(); it != inc_ed_a.end(); ++it) {
            size_type id_edge = *it;
            size_type f_node = this->vec_edges_[id_edge].node1_;
            size_type s_node = this->vec_edges_[id_edge].node2_;

            // I can also directly test the equality of the nodes
            // but I just want to stay with this for the moment
            // and in case of next assignements changing the structure
            Edge cur_edge = Edge(this,id_edge,f_node,s_node);

            if(((cur_edge.node1() == a) && (cur_edge.node2() == b)) ||
               ((cur_edge.node1() == b) && (cur_edge.node2() == a))){
                return (int) id_edge;
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
        size_type idx_a = a.index();
        size_type idx_b = b.index();

        int find_edge = index_has_edge(b,a);
        if(find_edge>-1){
            size_type idx_edge = (unsigned int) find_edge;
            return Edge(this,idx_edge,idx_a,idx_b);
        }else{
            vec_edges_.emplace_back(next_uidEdge_,idx_a,idx_b);

            // we add the edge to the incident vector of each end node
            vec_nodes_[idx_a].add_inci_edge(next_uidEdge_);
            vec_nodes_[idx_b].add_inci_edge(next_uidEdge_);

            ++size_edges_;
            ++next_uidEdge_;
            // Returns a Node that points to the new element
            return Edge(this, next_uidEdge_-1,idx_a,idx_b);
        }
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
    class NodeIterator : private equality_comparable<NodeIterator>{
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


        /** Increment operator for node iterator.
       * @pre Should be called on a valid node iterator
       * @return a node iterator which incremented by one the previous
       */
        NodeIterator& operator++()
        {
            it_node_++;
            return *this;
        }

        /** Equality of iterator.
        * @pre Should be called on two valids node iterator
        * @return whether ot not the 2 iterator are equals
        */
        //Defines equality between two iterators
        bool operator==(const NodeIterator& node_iterator) const
        {
            return(it_node_ == node_iterator.it_node_ && graph_==node_iterator.graph_);
        }

        /** Dereference iterator.
        * @pre Should be called on a valid node iterator
        * @return The node associated to the current step in the NodeIterator
        */
        value_type operator*() const
        {
            return graph_->node((*it_node_).uid_);
        }


    private:
        friend class Graph;
        // HW1 #2: YOUR CODE HERE
        graph_type * graph_;
        typename std::vector<node_p>::const_iterator it_node_;

        /** Private Constructor */
        NodeIterator(const graph_type* graph, typename std::vector<node_p>::const_iterator it_node)
                : graph_(const_cast<graph_type*>(graph)), it_node_(it_node) {
        };

    };

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_iterator node_begin() const
    // node_iterator node_end() const

    /** Iterator pointing at the start of the vec_nodes (which contains all nodes of the graph).
    * @pre Should be called if graph initialized
    * @return an Iterator over vec_nodes (std) starting at the beginning of it
    */
    NodeIterator node_begin() const
    {
        typename std::vector<node_p>::const_iterator begin = this->vec_nodes_.begin();
        return NodeIterator(this,begin);
    }

    /** Iterator pointing at the end of the vec_nodes (which contains all nodes of the graph).
    * @pre Should be called if graph initialized
    * @post Used for stl function over iterators to end call
    * @return an Iterator over vec_nodes (std) pointing at the end
    */
    NodeIterator node_end() const
    {
        typename std::vector<node_p>::const_iterator end = this->vec_nodes_.end();
        return NodeIterator(this,end);
    }

    //
    // Incident Iterator
    //

    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator: private equality_comparable<IncidentIterator> {
    public:
        // These type definitions let us use STL's iterator_traits.
        using value_type        = Edge;                     // Element type
        using pointer           = Edge*;                    // Pointers to elements
        using reference         = Edge&;                    // Reference to elements
        using difference_type   = std::ptrdiff_t;           // Signed difference
        using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

        /** Construct an invalid IncidentIterator. */
        IncidentIterator(){
        }

        // HW1 #3: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        // Edge operator*() const
        // IncidentIterator& operator++()
        // bool operator==(const IncidentIterator&) const

        /** Increment operator for incident iterator.
        * @pre Should be called on a valid incident iterator
        * @return a incident iterator which incremented by one the previous
        */
        IncidentIterator& operator++()
        {
            it_incident_++;
            return *this;
        }

        /** Equality of iterator.
         * @pre Should be called on two valids incident iterator
         * @return whether ot not the 2 iterator are equals
         */
        //Defines equality between two iterators
        bool operator==(const IncidentIterator& incident_iterator) const
        {
            return(it_incident_ == incident_iterator.it_incident_
                   && id_node_ == incident_iterator.id_node_
                   && graph_==incident_iterator.graph_);
        }

        /** Dereference iterator.
         * @pre Should be called on a valid incident iterator
         * @pre Be sure that graph.node(i).index() = i otherwise there could be conflict
         *      between nodes lookup
         * @return The node associated to the current step in the IncidentIterator
         */
        value_type operator*() const
        {
            edge_p edge = graph_->vec_edges_[*it_incident_];
            size_type id_edge = edge.uid_;
            size_type id_node1 =edge.node1_;
            size_type id_node2 =edge.node2_;
            return (id_node_==id_node1) ? Edge(graph_,id_edge,id_node_,id_node2) : Edge(graph_,id_edge,id_node_,id_node1);
        }


    private:
        friend class Graph;
        // HW1 #3: YOUR CODE HERE
        graph_type* graph_;
        size_type id_node_;
        typename std::set<size_type>::const_iterator it_incident_;

        /** Private Constructor */
        IncidentIterator(const graph_type* graph, const size_type id_node,
                         typename std::set<size_type >::const_iterator it_incident)
                : graph_(const_cast<graph_type*>(graph)), id_node_(id_node), it_incident_(it_incident) {
        };

    };

    //
    // Edge Iterator
    //

    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator: private equality_comparable<EdgeIterator> {
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
        /** Increment operator for edge iterator.
        * @pre Should be called on valid edge iterator
        * @return a edge iterator which incremented by one the previous
        */
        EdgeIterator& operator++()
        {
            it_edge_++;
            return *this;
        }

        /** Equality of iterator.
        * @pre Should be called on two valids edge iterator
        * @return whether ot not the 2 iterator are equals
        */
        //Defines equality between two iterators
        bool operator==(const EdgeIterator& edge_iterator) const
        {
            return(it_edge_ == edge_iterator.it_edge_ && graph_==edge_iterator.graph_);
        }

        /** Dereference iterator.
        * @pre Should be called on valid edge iterator
        * @return The edge associated to the current step in the NodeIterator
        */
        value_type operator*() const
        {
            return graph_->edge((*it_edge_).uid_);
        }

    private:
        friend class Graph;
        // HW1 #5: YOUR CODE HERE
        graph_type * graph_;
        typename std::vector<edge_p>::const_iterator it_edge_;

        /** Private Constructor */
        EdgeIterator(const graph_type* graph, typename std::vector<edge_p>::const_iterator it_edge)
                : graph_(const_cast<graph_type*>(graph)), it_edge_(it_edge) {
        };
    };

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // edge_iterator edge_begin() const
    // edge_iterator edge_end() const

    /** Iterator pointing at the start of the vec_edges (which contains all edges of the graph).
  * @pre Should be called if graph initialized
  * @return an Iterator over vec_edges (std) starting at the beginning of it
  */
    EdgeIterator edge_begin() const
    {
        typename std::vector<edge_p>::const_iterator begin = this->vec_edges_.begin();
        return EdgeIterator(this,begin);
    }

    /** Iterator pointing at the end of the vec_edges (which contains all edges of the graph).
    * @pre Should be called if graph initialized
    * @post Used for stl function over iterators to end call
    * @return an Iterator over vec_edges (std) pointing at the end
    */
    EdgeIterator edge_end() const
    {
        typename std::vector<edge_p>::const_iterator end = this->vec_edges_.end();
        return EdgeIterator(this,end);
    }

private:

    // HW0: YOUR CODE HERE
    // Use this space for your Graph class's internals:
    //   helper functions, data members, and so forth.

    // Internal type for graph nodes
    struct node_p {
        size_type uid_;
        Point point_;
        node_value_type value_;
        std::set<size_type> incident_edges_;

        node_p(size_type uid , const Point& point, const node_value_type& value,
               const std::set<size_type> incident_edges)
                :uid_(uid),point_(point),value_(value),incident_edges_(incident_edges) {};

        void add_inci_edge(size_type uid_edge){
            incident_edges_.insert(uid_edge);
        }

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

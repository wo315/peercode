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
template <typename V,typename E>
class Graph {
private:
    
    // Declare the internal structs:
    struct internal_node;
    struct internal_edge;
    struct linked_e_n; //linked node and edge to the node
    
    
public:
    
    //
    // PUBLIC TYPE DEFINITIONS
    //
    typedef V node_value_type;
    typedef E edge_value_type;
    
    /** Type of this graph. */
    using graph_type = Graph;
    
    
    /** Predeclaration of Node type. */
    class Node;
    /** Synonym for Node (following STL conventions). */
    using node_type = Node;
    //    /** Synonym for Node value. */
    //    using node_value_type = V;
    
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
    
    // CONSTRUCTORS AND DESTRUCTOR
    
    /** Construct an empty graph. */
    Graph() {
        // Initialize number of nodes and edges in the graph to 0.
        num_nodes_ = 0;
        num_edges_ = 0;
    }
    
    /** Default destructor */
    ~Graph() = default;
    
    /**************************** NODES ***************************************/
    
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
        
        /** HW2: Return this node's position. */
        Point& position() const {
            //pointer -> vec<internal_node> --> internal_node --> Point
            return pointer_graph_-> all_nodes_[node_idx_].node_point_;
        }
        
        
        /** Return this node's index, a number in the range [0, graph_size). */
        size_type index() const {
            return pointer_graph_-> all_nodes_[node_idx_].node_idx_;
        }
        
        
        /** @brief Return this node's value
         *
         * @return value for the current node
         *
         * This function accesses all internal_node via pointer_graph_,
         * then use the node_idx_ to find the current node and get its value.
         */
        node_value_type& value(){
            //pointer -> vec<internal_node> --> internal_node --> node_value_
            return pointer_graph_->all_nodes_[node_idx_].node_value_;
        }
        
        
        /** @brief Return this node's value, which cannot be changed
         *
         * @return value for the current node
         *
         * This function accesses all internal_node via pointer_graph_,
         * then use the node_idx_ to find the current node and get its value.
         */
        const node_value_type& value() const{
            return pointer_graph_->all_nodes_[node_idx_].node_value_;
        }
        
        
        /** @brief Return this node's degree, which is the number of incident edges
         *
         * @return degree of the current node
         *
         * This function use the vector<linked_e_n> all_linked_e_n_, which contains
         * all linked nodes and edges for a given node,
         * through this we can get the degree of each node.
         */
        size_type degree() const{
            //pointer -> vec<internal_node> --> internal_node --> all_linked_e_n_
            return pointer_graph_->all_nodes_[node_idx_].all_linked_e_n_.size();
        }
        
        
        /** @brief Return the start of the incident iterator pointing to edge that linked to the node
         *
         * @return Incident iterator pointing to the start of edges that linked to the node
         * The inicident iterator here helps us iterate over all edges that linked to the node
         */
        incident_iterator edge_begin() const{
            return IncidentIterator(pointer_graph_, this, 0);
        }
        
        
        /** @brief Return the end of the incident iterator pointing to edge that linked to the node
         *
         * @return Incident iterator pointing to the end of edges that linked to the node
         * The inicident iterator here helps us iterate over all edges that linked to the node
         */
        incident_iterator edge_end() const{
            return IncidentIterator(pointer_graph_, this, this->degree());
        }
        
        
        /** @brief Test whether this node and @a n are equal.
         *
         *  @return boolean value
         *  Equal nodes should have the same graph and the same index.
         */
        bool operator==(const Node& n) const {
            (void) n;          // Quiet compiler warning
            // check whether it's in same graph and has valid node's index
            return ((this->pointer_graph_ == n.pointer_graph_) && (this->node_idx_ == n.node_idx_));
        }
        
        /** @brief Test whether this node is less than @a n in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any geometric meaning.
         *
         * The node ordering relation must obey trichotomy: For any two nodes x
         * and y, exactly one of x == y, x < y, and y < x is true.
         */
        bool operator<(const Node& n) const {
            (void) n;           // Quiet compiler warning
            // first check graph
            if(this->pointer_graph_ < n.pointer_graph_){
                return true;
            }
            // then check if in the same graph and compare node's index
            else{
                return ((this->pointer_graph_ == n.pointer_graph_) && (this->node_idx_ < n.node_idx_));
            }
        }
        
    private:
        // Allow Graph to access Node's private member data and functions.
        friend class Graph;
        // Node's pointer back to the Graph
        graph_type* pointer_graph_;
        // Node's index
        size_type node_idx_;
        // Private constructor of Node
        Node(const graph_type* g, size_type nidx){
            pointer_graph_ = const_cast<graph_type*>(g);
            node_idx_ = nidx;
        }
    };
    
    /** @brief Return the number of nodes in the graph.
     *
     *  Complexity: O(1).
     */
    size_type size() const {
        return num_nodes_;
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
    
    Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
        /* step1: initialize a vector<linked_e_n> to store all nodes that linked to this new node later
         step2: create an internal_node and add it to the vector that store all nodes
         step3: update the number of nodes in the graph*/
        (void) position;      // Quiet compiler warning
        size_type idx = num_nodes_;
        std::vector<linked_e_n> linked{};
        internal_node new_node = internal_node(position, idx, linked, node_value);
        all_nodes_.push_back(new_node);
        node_uid_.push_back(idx);
        num_nodes_ += 1;
        return Node(this, idx);
    }
    
    /** Determine if a Node belongs to this Graph
     * @return True if @a n is currently a Node of this Graph
     *
     * Complexity: O(1).
     */
    bool has_node(const Node& n) const {
        (void) n;            // Quiet compiler warning
        return node(n.index()) == n;
    }
    
    /** Return the node with index @a i.
     * @pre 0 <= @a i < num_nodes()
     * @post result_node.index() == i
     *
     * Complexity: O(1).
     */
    Node node(size_type i) const {
        (void) i;             // Quiet compiler warning
        assert(i>=0);
        return Node(this, node_uid_[i]);
    }
    
    
    /** Remove node from the graph.
     * @return 1 if graph contains node n and it is successfully removed; else false
     * @param[in] Node n
     * @post has_node(n) is false
     * @post If original graph has node n, new num_nodes() == old num_nodes() - 1.
     *       else,  new num_nodes() == old num_nodes().
     * Complexity: O(num_nodes())
     */
    size_type remove_node(const Node& n){
        if (has_node(n)){
            //remove all connected edges
            for (auto it = n.edge_begin(); it != n.edge_end(); it = n.edge_begin()){
                remove_edge(*n.edge_begin());
            }
            //remove node from vector all_nodes_ and from vector node_uid_
            size_type n_idx = all_nodes_[n.node_idx_].node_idx_;
            size_type temp_uid = node_uid_.back();
            all_nodes_[temp_uid].node_idx_ = n_idx;
            node_uid_[n_idx] = temp_uid;
            node_uid_.pop_back();
            //update num_nodes
            num_nodes_--;
            return 1;
        }
        return 0;
    }
    
    
    /** Remove node from the graph.
     * @return true if graph contains node *n_it and it is successfully removed.
     * @param[in] node_iterator n_it
     * @post has_node(*n_it) is false
     * @post If node removed, new num_nodes() == old num_nodes() - 1.
     *       else,  new num_nodes() == old num_nodes().
     * Complexity: O(num_nodes())
     */
    node_iterator remove_node(node_iterator n_it){
        remove_node(*n_it);
        return n_it;
    }
    
    
    /************************************ EDGES ***************************************/
    
    /** @class Graph::Edge
     * @brief Class representing the graph's edges.
     *
     * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
     * are considered equal if they connect the same nodes, in either order.
     */
    class Edge : private totally_ordered<Edge>{
    public:
        /** Construct an invalid Edge. */
        Edge() {
        }
        
        /** Return a node of this Edge */
        Node node1() const {
            if (pointer_graph_ == nullptr){
                return Node(); // Invalid Node
            }
            return Node(pointer_graph_, node1_idx_);
        }
        
        /** Return the other node of this Edge */
        Node node2() const {
            if (pointer_graph_ == nullptr){
                return Node(); // Invalid Node
            }
            return Node(pointer_graph_, node2_idx_);
        }
        
        /** Test whether this edge and @a e are equal.
         *
         * Equal edges represent the same undirected edge between two nodes.
         */
        bool operator==(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            // check if the edges have same nodes
            bool check_nodes = (node1() == e.node1() && node2() == e.node2());
            bool check_nodes_rev = (node1() == e.node2() && node2() == e.node1());
            // check if the edges are in the same graph
            bool check_graph = (pointer_graph_ == e.pointer_graph_);
            // check if the edges have same index
            bool check_edge = (edge_idx_ == e.edge_idx_);
            return (check_graph && (check_nodes || check_nodes_rev) && check_edge);
        }
        
        /** Test whether this edge is less than @a e in a global order.
         *
         * This ordering function is useful for STL containers such as
         * std::map<>. It need not have any interpretive meaning.
         */
        bool operator<(const Edge& e) const {
            (void) e;           // Quiet compiler warning
            if (pointer_graph_ < e.pointer_graph_){
                return true;
            }else{
                return ((pointer_graph_ == e.pointer_graph_) && (this->edge_idx_ < e.edge_idx_));
            }
        }
        
        /** Returns edge's value. */
        edge_value_type& value() {
            return pointer_graph_->all_edges_[edge_idx_].edge_value_;
        }
        
        const edge_value_type& value() const {
            return pointer_graph_->all_edges_[edge_idx_].edge_value_;
        }
        
    private:
        // Allow Graph to access Edge's private member data and functions.
        friend class Graph;
        graph_type* pointer_graph_;
        size_type edge_idx_;
        size_type node1_idx_;
        size_type node2_idx_;
        // Private constructor of Edge
        Edge(const graph_type* g, size_type eidx, size_type n1, size_type n2){
            pointer_graph_ = const_cast<graph_type*>(g);
            edge_idx_ = eidx;
            node1_idx_ = n1;
            node2_idx_ = n2;
        }
    };
    
    /** Return the total number of edges in the graph.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    size_type num_edges() const {
        return num_edges_;
    }
    
    /** Return the edge with index @a i.
     * @pre 0 <= @a i < num_edges()
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    Edge edge(size_type i) const {
        (void) i;             // Quiet compiler warning
        assert(i>=0);
        assert(i< num_edges_);
        size_type j = edge_uid_[i];
        size_type nid1 = all_edges_[j].node1_idx_;
        size_type nid2 = all_edges_[j].node2_idx_;
        return Edge(this,j, nid1, nid2);
    }
    
    /** Test whether two nodes are connected by an edge.
     * @pre @a a and @a b are valid nodes of this graph
     * @return True if for some @a i, edge(@a i) connects @a a and @a b.
     *
     * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
     */
    bool has_edge(const Node& a, const Node& b) const {
        (void) a; (void) b;   // Quiet compiler warning
        /* step1: find the vector that contains all linked nodes and edges to Node a
         step2: loop over the vector, check if it contains node b's index */
        std::vector<linked_e_n> all_linked = all_nodes_[a.node_idx_].all_linked_e_n_;
        for (size_type i=0; i< all_linked.size(); i++){
            if (all_linked[i].node_idx_ == b.node_idx_){
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
    Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value = edge_value_type()) {
        (void) a, (void) b;   // Quiet compiler warning
        assert(has_node(a));
        assert(has_node(b));
        assert(!(a==b));
        
        size_type eidx;
        bool cond = false; //check if the edge already exists
        std::vector<linked_e_n> all_linked = all_nodes_[a.node_idx_].all_linked_e_n_;
        for (size_type i=0; i< all_linked.size(); i++){
            if (all_linked[i].node_idx_ == b.node_idx_){
                cond = true;
                eidx = all_linked[i].edge_idx_;
                break;
            }
        }
        
        //if the graph not contain edge(a,b)
        if(cond == false){
            //add the edge to the vector all_edges_
            size_type eidx = num_edges_;
            internal_edge new_edge = internal_edge(a.node_idx_, b.node_idx_,eidx, edge_value);
            all_edges_.push_back(new_edge);
            
            //update linked struct linked to node a and node b
            all_nodes_[a.node_idx_].all_linked_e_n_.push_back(linked_e_n(b.node_idx_,eidx));
            all_nodes_[b.node_idx_].all_linked_e_n_.push_back(linked_e_n(a.node_idx_,eidx));
            
            edge_uid_.push_back(eidx);
            //update number of edges
            num_edges_ += 1;
            return Edge(this, eidx, a.node_idx_, b.node_idx_);
        }
        else{
            //if the edge already exists
            //return Edge();
            return Edge(this, eidx, a.node_idx_, b.node_idx_);
        }
        
    }
    
    /** Remove all nodes and edges from this graph.
     * @post num_nodes() == 0 && num_edges() == 0
     *
     * Invalidates all outstanding Node and Edge objects.
     */
    void clear() {
        num_nodes_ = 0;
        num_edges_ = 0;
        all_edges_.clear();
        all_nodes_.clear();
        node_uid_.clear();
        edge_uid_.clear();
    }
    
    
    /** Remove edges from the graph.
     * @return true if graph contains edge between node n1, n2 and it is successfully removed.
     * @param[in] Node n1, Node2
     * @pre has_node(n1) is true, has_node(n2) is true
     * @post has_node(n1) is false, has_node(n2) is false
     * @post If original graph has node n, new num_edges() == old num_edges() - 1.
     *       else,  new num_edges() == old num_edges().
     * Complexity: O(num_nodes() + num_edges())
     */
    size_type remove_edge(const Node& n1, const Node& n2){
        if(has_edge(n1,n2)){
            unsigned idx1 = n1.node_idx_; //node1 index
            unsigned idx2 = n2.node_idx_; //node2 index
            std::vector<linked_e_n> linked1 = all_nodes_[idx1].all_linked_e_n_; //all nodes and edges linked to n1
            
            for (size_type i=0; i< linked1.size(); i++){
                if (linked1[i].node_idx_ == idx2){
                    linked_e_n e = linked1[i];
                    //erase edge and node2 from the n1's linked vector
                    all_nodes_[idx1].all_linked_e_n_[i] = all_nodes_[idx1].all_linked_e_n_.back();
                    all_nodes_[idx1].all_linked_e_n_.pop_back();
                    //erase edge and node1 from the n2's linked vector
                    for (size_type j = 0; j < n2.degree(); j++) {
                        if (all_nodes_[idx2].all_linked_e_n_[j].node_idx_ == idx1){
                            all_nodes_[idx2].all_linked_e_n_[j] = all_nodes_[idx2].all_linked_e_n_.back();
                            all_nodes_[idx2].all_linked_e_n_.pop_back();
                            break;
                        }
                    }
                    //remove edge from vector all_edges_ and from vector edge_uid_
                    size_type e_idx = all_edges_[e.edge_idx_].edge_idx_;
                    size_type temp_uid = edge_uid_.back();
                    all_edges_[temp_uid].edge_idx_ = e_idx;
                    edge_uid_[e_idx] = temp_uid;
                    edge_uid_.pop_back();
                    //update num_edges
                    num_edges_--;
                    return 1;
                }
            }
        }
        return 0;
    }
    
    
    /** Remove edges from the graph.
     * @return true if graph contains edge e and it is successfully removed.
     * @param[in] Edge e
     * @pre has_edge(e) is true
     * @post has_edge(e) is false
     * @post If original graph has node n, new num_edges() == old num_edges() - 1.
     *       else,  new num_edges() == old num_edges().
     * Complexity: O(num_nodes() + num_edges())
     */
    size_type remove_edge(const Edge& e){
        return remove_edge(e.node1(), e.node2());
    }
    
    /** Remove edges from the graph.
     * @return true if the edge that edge_iterator at is successfully removed.
     * @param[in] edge_iterator e_it
     * @pre has_edge(*e_it) is true
     * @post has_edge(*e_it) is false
     * @post If original graph has node n, new num_edges() == old num_edges() - 1.
     *       else,  new num_edges() == old num_edges().
     * Complexity: O(num_nodes() + num_edges())
     */
    edge_iterator remove_edge(edge_iterator e_it){
        return remove_edge(*e_it);
    }
    
    /**************************** Node Iterator ********************************************/
    
    /** @class Graph::NodeIterator
     * @brief Iterator class for nodes. A forward iterator. */
    class NodeIterator : private totally_ordered<NodeIterator> {
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
        
        /** Return the node that the iterator is pointing to */
        Node operator*() const{
            return pointer_graph_->node(node_iter_idx_);
        }
        
        /** Return the iterator that is moving to next node.*/
        NodeIterator& operator++(){
            ++node_iter_idx_;
            return *this;
        }
        
        /** Test whether this node iterator and node_iterator2 are equal.
         * Equal node iterators should have the same graph and the same node_iter index.*/
        bool operator==(const NodeIterator& iterator2) const{
            //same graph
            bool cond1 = (pointer_graph_ == iterator2.pointer_graph_);
            //same iterator
            bool cond2 = (node_iter_idx_ == iterator2.node_iter_idx_);
            return cond1 && cond2;
        }
        
    private:
        friend class Graph;
        graph_type* pointer_graph_;
        size_type node_iter_idx_;
        NodeIterator(const graph_type* g, size_type idx){
            pointer_graph_ = const_cast<graph_type*>(g);
            node_iter_idx_ = idx;
        }
        
    };
    
    /** @brief Return node iterator pointing to the start of nodes
     *
     * @return Node iterator pointing to the start of nodes.
     * This can be used to iterate over all the nodes.
     */
    node_iterator node_begin() const{
        return NodeIterator(this, 0);
    }
    
    /** @brief Return the node iterator pointing to the end of nodes
     *
     * @return Node iterator pointing to the end of nodes
     * This can be used to iterate over all the nodes
     */
    node_iterator node_end() const{
        return NodeIterator(this, this->num_nodes());
    }
    
    /**************************** Incident Iterator  ********************************/
    
    /** @class Graph::IncidentIterator
     * @brief Iterator class for edges incident to a node. A forward iterator. */
    class IncidentIterator : private totally_ordered<IncidentIterator>{
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
        
        /** Returns an edge that connects to the node
         *
         * Based on the node_idx of linked node and edge_idx_,
         * find the edge that connects to the current node.
         * @post Edge that connect this incident linked to the nodd
         */
        Edge operator*() const{
            size_type n1 = node_->node_idx_;
            //pointer_graph_ --> vector<internalnode> --> internalnode --> vector<linked_e_n> --> linked_e_n --> node_idx_
            size_type n2 = pointer_graph_-> all_nodes_[n1].all_linked_e_n_[incident_iter_idx_].node_idx_;
            size_type e = pointer_graph_-> all_nodes_[n1].all_linked_e_n_[incident_iter_idx_].edge_idx_;
            return Edge(pointer_graph_, e, n1, n2);
        }
        
        /** Move the incident iterator to the next one */
        IncidentIterator& operator++(){
            ++incident_iter_idx_;
            return *this;
        }
        
        /** Check if two incident iterators are equal
         * Equal node iterators should have the same graph and the same node_iter index. */
        bool operator==(const IncidentIterator& iter2) const{
            //same graph
            bool cond1 = (pointer_graph_ == iter2.pointer_graph_);
            //same node
            bool cond2 = (node_ == iter2.node_);
            //same iterator
            bool cond3 = (incident_iter_idx_ == iter2.incident_iter_idx_);
            return (cond1 && cond2 && cond3);
        }
        
    private:
        friend class Graph;
        graph_type* pointer_graph_;
        node_type* node_;
        size_type incident_iter_idx_;
        
        //constructor
        IncidentIterator(const graph_type* g,const node_type* n,const size_type idx){
            pointer_graph_ = const_cast<graph_type*>(g);
            node_ = const_cast<node_type*>(n);
            incident_iter_idx_ = idx;
        }
    };
    
    /********************************** Edge Iterator  ************************************/
    
    /** @class Graph::EdgeIterator
     * @brief Iterator class for edges. A forward iterator. */
    class EdgeIterator : private totally_ordered<EdgeIterator>{
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
        
        /** Return the edge with its edge_iter_idx_.*/
        Edge operator*() const{
            return pointer_graph_->edge(edge_iter_idx_);
        }
        
        /** Move the edge iterator to the next one.*/
        EdgeIterator& operator++(){
            ++edge_iter_idx_;
            return *this;
        }
        
        /** Test whether the two edge iterators are equal.
         *
         * Equal edge iterators should have the same graph and the same index.
         */
        bool operator==(const EdgeIterator& iter2) const{
            //same graph
            bool cond1 = (pointer_graph_ == iter2.pointer_graph_);
            //same edge iterator
            bool cond2 = (edge_iter_idx_ == iter2.edge_iter_idx_);
            return cond1 && cond2;
        }
        
    private:
        friend class Graph;
        graph_type* pointer_graph_;
        size_type edge_iter_idx_;
        EdgeIterator(const graph_type* g, size_type idx){
            pointer_graph_ = const_cast<graph_type*>(g);
            edge_iter_idx_ = idx;
        }
    };
    
    /** Return an iterator pointing at the start of the edges */
    edge_iterator edge_begin() const{
        return EdgeIterator(this, 0);
    }
    
    /** Return an iterator pointing at the end of the edges */
    edge_iterator edge_end() const{
        return EdgeIterator(this, this->num_edges());
    }
    
private:
    /*Private members for class Graph*/
    size_type num_nodes_; //number of nodes in the graph
    size_type num_edges_; //number of edges in the graph
    
    struct linked_e_n{
        size_type node_idx_; //node that linked to
        size_type edge_idx_; //edge that used to link
        
        //constructor for linked_e_n
        linked_e_n(const size_type nidx, const size_type eidx){
            node_idx_ = nidx;
            edge_idx_ = eidx;
        }
    };
    
    struct internal_node{
        Point node_point_; //node's point
        size_type node_idx_; //node's index
        std::vector<linked_e_n> all_linked_e_n_; //all lineked edges and nodes
        node_value_type node_value_; //value of the node
        
        //hw1: constructor for internal_node
        internal_node(Point p, size_type idx, std::vector<linked_e_n> linked, node_value_type node_value){
            node_point_ = p;
            node_idx_ = idx;
            all_linked_e_n_ = linked;
            node_value_ = node_value;
        }
    };
    
    struct internal_edge{
        size_type node1_idx_; //node1's index
        size_type node2_idx_; //node2's index
        size_type edge_idx_; //edge's index
        edge_value_type edge_value_; //value of the node
        
        //constructor for internal_edge
        internal_edge(size_type n1, size_type n2, size_type eidx, edge_value_type v){
            node1_idx_ = n1;
            node2_idx_ = n2;
            edge_idx_ = eidx;
            edge_value_ = v;
        }
    };
    
    // use vector<internal_nodes> & vector<internal_edges> to store all nodes and edges
    std::vector<internal_node> all_nodes_;
    std::vector<internal_edge> all_edges_;
    
    // map nodes/edges' index to unique id
    std::vector<unsigned int> node_uid_;
    std::vector<unsigned int> edge_uid_;
};

#endif // CME212_GRAPH_HPP


#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>            
#include <set>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/* util func for swapping arbitrary elements (taken from lec notes)*/
template < class T>
void swap (T &a, T &b) {
    T tmp { static_cast <T&& >(a )};
    a = static_cast <T&& >(b);
    b = static_cast <T&& >(tmp );
}

/* util func for swapping arbitrary elements (taken from lec notes)*/
template < class T>
void swap_graph_hpp (T &a, T &b) {
    T tmp { static_cast <T&& >(a )};
    a = static_cast <T&& >(b);
    b = static_cast <T&& >(tmp );
}


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

    public:

        //
        // PUBLIC TYPE DEFINITIONS
        //

        /** Type of templated variable */
        using node_value_type = V;

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
            num_nodes_ = 0;
            num_edges_ = 0;
            std::vector<Node> nodes_;
            std::vector<Edge> edges_;
            std::vector<Point*> point_pointers_;
            std::map<size_type, std::set<size_type>> node_to_edges_;
            std::map<size_type, std::set<size_type>> node_to_nodes_;
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
                    uid_=-1;
                    graph_=nullptr;    
                }

                /** Return this node's position. */
                const Point& position() const {
                    // HW0: YOUR CODE HERE
                    size_type indx = graph_->node_uid_to_indx_[uid_];
                    return *(graph_->point_pointers_)[indx];
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    // HW0: YOUR CODE HERE
                    return graph_->node_uid_to_indx_[uid_];
                }

                // HW1: YOUR CODE HERE
                //Supply definitions AND SPECIFICATIONS for:
                node_value_type& value() {
                    size_type node_indx = this->index();
                    return (graph_->node_values_[node_indx]);
                }
                const node_value_type& value() const {
                    size_type node_indx = this->index();
                    return (graph_->node_values_[node_indx]);
                }
                size_type degree() const {
                    return graph_->node_to_nodes_[uid_].size();
                }
                incident_iterator edge_begin() const {
                    IncidentIterator ie {graph_, uid_, 0};
                    return ie;
                }

                 incident_iterator edge_end() const {
                    size_type incident_edges = this->graph_->node_to_edges_[this->uid_].size();
                    IncidentIterator ie {graph_, uid_, incident_edges};
                    return ie;
                 }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    // HW0: YOUR CODE HERE
                    if (this->uid_==n.uid_ && this->graph_==n.graph_)
                        return true; 
                    return false;
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
                    if (uid_<n.uid_)
                        return true; 
                    return false;
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                // HW0: YOUR CODE HERE
                // Use this space to declare private data members and methods for Node
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Node objects
                graph_type* graph_;
                size_type uid_;

                /* private constructor for when the associated graph and point is known 
                 * This handles uid increment 
                 * But the caller is responsible for adding the Node to graph and 
                 * keeping track of the node count 
                 */

                Node(graph_type* graph, size_type uid) {
                    graph_ = graph;
                    uid_ = uid;
                }
        };

        /** Return the number of nodes in the graph.
         *
         * Complexity: O(1).
         */
        size_type size() const {
            // HW0: YOUR CODE HERE
            return this->num_nodes_;
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
        Node add_node(const Point& position,  const node_value_type& val= node_value_type()) {
            // HW0: YOUR CODE HERE
            // Create new point object to hold data
            Point* new_point = new Point(position.x, position.y, position.z);
            point_pointers_.push_back(new_point);
            // Create new Node
            size_type new_uid = next_node_uid_;
            next_node_uid_++;
            node_uid_to_indx_[new_uid] = nodes_.size(); 
            Node n(this, new_uid);
            // Add extra node_value_type
            node_values_.push_back(val);
            // Add node 
            nodes_.push_back(n);
            num_nodes_++;

            return n; 
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         */
        bool has_node(const Node& n) const {
            // HW0: YOUR CODE HERE
            if (n.graph_ == this) {
                return true;
            }
            return false;
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
            // HW0: YOUR CODE HERE
            if (i>= num_nodes())
                throw std::invalid_argument("Index too large."); 
            return nodes_[i];        
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
                    node1_ = new Node();
                    node2_ = new Node();
                }

                /**
                 * Acts in place
                 */
                void swap_edge_nodes_inplace_() {
                    swap_graph_hpp(this->node1_, this->node2_);
                }

                /** Return a node of this Edge */
                Node node1() const {
                    // HW0: YOUR CODE HERE
                    return *node1_; 
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    // HW0: YOUR CODE HERE
                    return *node2_;
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    //HW0: YOUR CODE HERE
                    if ( 
                            // e.node1().graph_ == e.node2().graph_
                            //
                            (
                             (
                              e.node1() == this->node1() && e.node2() == this->node2()
                             ) ||
                             (
                              e.node1() == this->node2() && e.node2() == this->node1()
                             )
                            )
                       ) return true;
                    return false;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    //HW0: YOUR CODE HERE
                    if (this->node1() < e.node1())
                        return true;
                    if (this->node2() < e.node2())
                        return true;
                    return false;
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // HW0: YOUR CODE HERE
                // Use this space to declare private data members and methods for Edge
                // that will not be visible to users, but may be useful within Graph.
                // i.e. Graph needs a way to construct valid Edge objects

                // Note that its a non-const pointer to const node
                const Graph::Node* node1_;
                const Graph::Node* node2_;
                size_type uid_=0;
                graph_type* graph_;

                /* Private constructor
                 * Like the private Node contructor it does create a new uid and
                 * increments the graph.next_edgeuid_ counter. But it does not add 
                 * the edges to the graph container, or update the size, and this
                 * must be handled by the caller.
                 */
                Edge(const Node* a, const Node* b, graph_type* graph, size_type uid)  {
                    node1_ = a;
                    node2_ = b;
                    uid_ = uid;
                    graph_ = graph;
                }
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const {
            // HW0: YOUR CODE HERE
            return num_edges_;
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const {
            // HW0: YOUR CODE HERE
            if (i>= num_edges_)
                throw std::invalid_argument("Index too large.");

            return edges_[i];
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const {
            // HW0: YOUR CODE HERE
            if (a.graph_!=this || b.graph_!=this)
                return false;
            // check the node map
            if (node_to_nodes_.count(a.uid_) && node_to_nodes_.at(a.uid_).count(b.uid_)>0) {
                return true;
            }
            return false;
        }
        
        size_type find_edge_index_(const Node& a, const Node& b) { 
            std::set<size_type> a_edges = node_to_edges_.at(a.uid_);
            std::set<size_type> b_edges = node_to_edges_.at(b.uid_);
            for (auto it=a_edges.begin(); it!=a_edges.end(); it++) {
                if (b_edges.count(*it)>0)
                    return edge_uid_to_indx_[*it];
            }
            return -1;
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
            // Check this isn't a duplicate edge
            if (this->has_edge(a,b)) {
                size_type indx = find_edge_index_(a,b);
                Edge e = edges_[indx];
                // Ensure that a is node1() in the returned edge
                if (e.node1() != a) 
                    e.swap_edge_nodes_inplace_();
                return e;
            }

            if (a.graph_ != b.graph_)
                throw "Nodes are not from the same graph";
                
            // create new edge id and do other accounting
            size_type new_uid = next_edge_uid_; 
            next_edge_uid_++;
            edge_uid_to_indx_[new_uid] = edges_.size(); 
            Edge e {&a, &b, this, new_uid};
            edges_.push_back(e);
            num_edges_++;
            // record the edge that these nodes are in 
            node_to_edges_[a.uid_].insert(e.uid_);
            node_to_edges_[b.uid_].insert(e.uid_);
            node_to_nodes_[a.uid_].insert(b.uid_);
            node_to_nodes_[b.uid_].insert(a.uid_);

            return e; 
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
            // HW0: YOUR CODE HERE
            for (size_type i=0; i<num_nodes_; i++) {
                nodes_[i].graph_ = nullptr;
                nodes_[i].uid_= -1;
            }
            nodes_.clear();
            num_nodes_=0;
            for (size_type i=0; i<num_edges_; i++) {
                edges_[i].graph_ = nullptr;
                edges_[i].uid_= -1;
            }
            edges_.clear();
            num_edges_=0;

        }

        //
        // Node Iterator
        //

        /** @class Graph::NodeIterator
         * @brief Iterator class for nodes. A forward iterator. */
        class NodeIterator: private totally_ordered<NodeIterator>  {
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
                Node operator*() const{
                    Node n = graph_->nodes_[node_indx_];
                    return n;
                }
                NodeIterator& operator++() {
                    node_indx_++;
                    return *this;
                }
                bool operator==(const NodeIterator& n) const {
                    return (n.node_indx_ == this->node_indx_ 
                            && n.graph_==this->graph_);
                }

            private:
                friend class Graph;
                // HW1 #2: YOUR CODE HERE
                const graph_type* graph_;
                size_type node_indx_;
                NodeIterator(const graph_type* graph, size_type node_indx) 
                    //: graph_(graph), node_indx_(node_indx) {};
                  { graph_ = graph; node_indx_=node_indx;
                  }
        };

        // HW1 #2: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        node_iterator node_begin() const {
            NodeIterator ni(this, 0);
            return ni;
        }
        node_iterator node_end() const {
            size_type num_nodes = nodes_.size();
            NodeIterator ni {this, num_nodes};
            return ni;
        }
        
        //
        // Incident Iterator
        //

        /** @class Graph::IncidentIterator
         * @brief Iterator class for edges incident to a node. A forward iterator. */
        class IncidentIterator : private totally_ordered<IncidentIterator> {
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
                Edge operator*() const {
                    graph_type* graph_=this->graph_;
                    // get iterator to the start of the set of neighbouring node uids
                    auto it = graph_->node_to_edges_[central_node_uid_].begin();
                    // iterate the number of spaces based on index
                    for (size_type j=0; j<indx_; j++)
                        it++;
                    // Get the edge indx to return by looking up the uid
                    size_type edge_indx = graph_->edge_uid_to_indx_[*it];
                    Edge e = graph_->edges_[edge_indx];
                    // check that the nodes are in the right ordering, if not then swap
                    if (e.node1().uid_ != central_node_uid_) 
                        e.swap_edge_nodes_inplace_();
                    return e;    
                }
                IncidentIterator& operator++() {
                    indx_++;
                    return *this;
                }
                bool operator==(const IncidentIterator& ie) const {
                    return  (ie.graph_==this->graph_
                        && ie.central_node_uid_==this->central_node_uid_
                        && ie.indx_==this->indx_);
                }

            private:
                friend class Graph;
                // HW1 #3: YOUR CODE HERE
                graph_type* graph_;
                size_type central_node_uid_; //
                size_type indx_;
                
                // Constructor 
                IncidentIterator(graph_type* graph, size_type central_node_uid, size_type indx) 
                    : graph_(graph), central_node_uid_(central_node_uid), indx_(indx) {};
        };

        //
        // Edge Iterator
        //

        /** @class Graph::EdgeIterator
         * @brief Iterator class for edges. A forward iterator. */
        class EdgeIterator : totally_ordered<EdgeIterator> {
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
                Edge operator*() const {
                    return graph_->edges_[indx_];
                }
                EdgeIterator& operator++() {
                    indx_++;
                    return *this;
                }
                bool operator==(const EdgeIterator& ei) const {
                    return (ei.graph_==this->graph_
                            && ei.indx_==this->indx_);
                }

            private:
                friend class Graph;
                // HW1 #5: YOUR CODE HERE
                const graph_type* graph_;
                size_type indx_;

                // Constructor
                EdgeIterator(const graph_type* graph, size_type indx) 
                    : graph_(graph), indx_(indx) {};
        };

        // HW1 #5: YOUR CODE HERE
        // Supply definitions AND SPECIFICATIONS for:
        
        edge_iterator edge_begin() const {
            EdgeIterator ei {this, 0};
            return ei;
        }
        edge_iterator edge_end() const {
            size_type num_edges = edges_.size(); 
            EdgeIterator ei {this, num_edges};
            return ei;
        }

    private:
        // HW0: YOUR CODE HERE
        // Use this space for your Graph class's internals:
        //   helper functions, data members, and so forth.
        std::vector<Node> nodes_;
        std::vector<Edge> edges_;
        std::vector<const Point*> point_pointers_;
        // lookup node and edge position in vector
        std::map<size_type, size_type> node_uid_to_indx_;
        std::map<size_type, size_type> edge_uid_to_indx_;
        // map of nodes and the edge number they're in. And to nodes they have an edge with
        std::map<size_type, std::set<size_type>> node_to_edges_;
        std::map<size_type, std::set<size_type>> node_to_nodes_;

        // container to store node_values. Index will match nodes and pointer
        std::vector<node_value_type> node_values_;

        // track num elements
        size_type num_nodes_;
        size_type num_edges_;
        
        // track next_uids that will only be incremented
        size_type next_node_uid_;
        size_type next_edge_uid_;
};

#endif // CME212_GRAPH_HPP

#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>


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

  struct node_element;
  
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
  
  /** Type of value of node **/
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
  Graph() : nodes_(), adjacency_list_(),edge_ptrs_(){
    size_ = size_type();
    num_edges_ = size_type();
    next_uid_ = size_type();
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
    }
    
    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return uid_;
    }

    /** Returns a reference to the node's value on the graph, or an empty reference
      * if the node is not found.
      */
    node_value_type& value() {
      return graph_->nodes_[uid_].value;        
    }
    
    /** Returns a const reference to the node's value on the graph, or an empty reference
      * if the node is not found.
      */
    const node_value_type& value() const {
      return graph_->nodes_[uid_].value;  
    }
    
    /** Returns the number of nodes connected to this node.
      * @return degree of node, or -1 if node cannot be found on graph
      */
    size_type degree() const {
      auto reachable = graph_->adjacency_list_.find(*this);
      if (reachable != graph_->adjacency_list_.end()) {
        return reachable->second.size();
      }
      return -1;
    }
    
    /** Returns Incident Iterator over edges connected to this node.
      * @return Incident Iterator pointing to first connected node
      */
    incident_iterator edge_begin() const{
      auto nodes = graph_->adjacency_list_.find(*this);
      return incident_iterator(this, nodes->second.cbegin());
    }
    
    /** Returns Iterator at end of edges connected to this node.
      * @return Incident Iterator pointing to end of connected nodes
      */
    incident_iterator edge_end() const {
      auto nodes = graph_->adjacency_list_.find(*this);
      return incident_iterator(this, nodes->second.cend());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return graph_==n.graph_ && index()==n.index();
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
      if (graph_ !=n.graph_) return (uintptr_t) graph_ + (uintptr_t) uid_ < (uintptr_t) n.graph_ + (uintptr_t) n.uid_;
      return uid_ < n.uid_;
    
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;    
    /** Pointer to container graph **/
    Graph* graph_;
    /** Unique Identification **/
    size_type uid_;
    /** Private Constructor */
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    node_element new_element = {position, value, next_uid_};
    nodes_.push_back(new_element);
    std::set<Node> empty;
    adjacency_list_[Node(this,next_uid_)] = empty;
    next_uid_++;
    size_++;
    return Node(this, next_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.graph_==this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i<size_) {
      return Node(this, nodes_[i].uid);
    }
    return Node();        // Invalid node
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
    /** Construct an invalid Edge. */
    Edge() {
    }

    /** Return a node of this Edge */
    Node node1() const {
      return (*node_1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return (*node_2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return ((*node_1_==*e.node_1_ && *node_2_ == *e.node_2_) || (*node_2_==*e.node_1_ && *node_1_ == *e.node_2_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Using the Cantor Pairing Function 
      return (!(*this==e)) &&((node_1_->index()+node_2_->index()+1)*(node_1_->index()+node_2_->index())/2+node_2_->index()) <
          ((e.node_1_->index()+e.node_2_->index()+1)*(e.node_1_->index()+e.node_2_->index())/2+e.node_2_->index());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    /** Nodes **/
    Node* node_1_;
    Node* node_2_;
    /** Private Constructor **/
    Edge(const Node* node_1, const Node* node_2):  node_1_(const_cast<Node*>(node_1)), node_2_(const_cast<Node*>(node_2)){
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
    return (edge_ptrs_[i]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    auto node_a = adjacency_list_.find(a);
    if (node_a != adjacency_list_.end()) {
      return node_a->second.find(b) != node_a->second.end();
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
    Edge new_element = Edge(&a,&b);
    if(has_edge(a,b)) return new_element;
    adjacency_list_[a].insert(b);
    adjacency_list_[b].insert(a);
    edge_ptrs_.push_back(new_element);
    num_edges_++;
    return new_element;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    adjacency_list_.clear();
    edge_ptrs_.clear();
    size_ = 0;
    num_edges_=0;
    next_uid_=0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private equality_comparable<NodeIterator> {
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

    /** Returns the node pointed to by the iterator.
      * @pre NodeIterator does not already point to the end of the list
      */
    Node operator*() const {
      return Node(graph_, (*iter_).uid);
    }
    
    /** Increments the iterator over list of nodes.
      * @pre NodeIterator does not already point to the end of the list
      * @post NodeIterator points to next node
      * @return Incremented NodeIterator
      */
    NodeIterator& operator++() {
      ++iter_;
      return *this;
    }
    
    /** Checks if two NodeIterators are equal.
      * @pre other node is not a null reference
      */
    bool operator==(const NodeIterator& other) const {
      return iter_ == other.iter_ && graph_ == other.graph_;
    }

   private:
    friend class Graph;
    /** Underlying iterating element **/
    typename std::vector<node_element>::const_iterator iter_;
    /** Container graph**/
    Graph* graph_;
    /** Constructor from Graph **/
    NodeIterator(const Graph* graph, typename std::vector<node_element>::const_iterator it) : iter_(it), graph_(const_cast<Graph*>(graph)){}
  };

  /** Returns iterator pointing to first node on graph */
  node_iterator node_begin() const {
    return node_iterator(this, nodes_.begin());
  }
  /** Returns iterator pointing to last node on graph */
  node_iterator node_end() const {
    return node_iterator(this, nodes_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
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

    /** Returns the edge pointed to by the iterator.
      * @pre iterator does not point to end of the list
      */
    Edge operator*() const {
      return Edge(origin_, &(*iter_));
    }
    
    /** Increments the iterator over list of connected edges.
      * @pre IncidentIterator does not already point to the end of the list
      * @post IncidentIterator points to next node-connected edge
      * @return Incremented IncidentIterator
      */
    IncidentIterator& operator++() {
      ++iter_;
      return *this;
    }
    
    /** Checks if two IncidentIterators are equal.
      * @pre other iterator is not a null reference
      */
    bool operator==(const IncidentIterator& other) const {
      return origin_==other.origin_ && iter_==other.iter_;
    }

   private:
    friend class Graph;
    /** Original Node which all iterated edges connected to */
    Node* origin_;
    /** Iterator over nodes edges */
    typename std::set<Node>::const_iterator iter_;
    /** Private Constructor */
    IncidentIterator(const Node * node,typename std::set<Node>::const_iterator itr) :  origin_(const_cast<Node*>(node)), iter_(itr) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
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

    /** Returns the edge pointed to by the iterator.
      * @pre iterator does already point to the end of the list
      */
    Edge operator*() const {
      return *itr_;
    }
    
    /** Increments the iterator over list of connected edges.
      * @pre EdgeIterator does not already point to the end of the list
      * @post EdgeIterator points to next edge on graph
      * @return Incremented EdgeIterator
      */
    EdgeIterator& operator++() {
      ++itr_;
      return *this;
    }
    /** Checks if two EdgeIterators are equal.
      * @pre other iterator is not a null reference
      */
    bool operator==(const EdgeIterator& other) const {
      return graph_ == other.graph_ && itr_ == other.itr_;
    }

   private:
    friend class Graph;
    /** Iterator over all edges */
    typename std::vector<Edge>::const_iterator itr_;
    /** Owner Graph of edges */
    Graph* graph_;
    /** Private Constructor */
    EdgeIterator(const Graph* graph, typename std::vector<Edge>::const_iterator itr) : itr_(itr), graph_(const_cast<Graph*>(graph)){}
  };

  /** Returns iterator pointing to the first edge on the graph */
  edge_iterator edge_begin() const {
    return edge_iterator(this, edge_ptrs_.begin());
  }
  /** Returns iterator pointing to the last edge on the graph */
  edge_iterator edge_end() const {
    return edge_iterator(this, edge_ptrs_.end());
  }

 private:
  struct node_element {
    Point point;
    node_value_type value;
    size_type uid;
  };
  /** Vector of Nodes **/
  std::vector<node_element> nodes_;
  /** Adjacency List **/
  std::map<Node,std::set<Node>> adjacency_list_;
  /** Vector of Unique Edges **/
  std::vector<Edge> edge_ptrs_;
  /** Number of Nodes **/
  size_type size_;
  /** Number of Edges **/
  size_type num_edges_;
  /** Next Unique Id, For now, maintaining that @a next_uid_ is the node index**/
  size_type next_uid_;

};

#endif // CME212_GRAPH_HPP

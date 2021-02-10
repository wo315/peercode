#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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
  //internal node that contains all info on node
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  /** type of template node value */
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
  Graph()
            : size_(0), num_distinct_edges_(0), internal_nodes_(), internal_edges_()  // init internal_nodes and size_
  {
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
  class Node  {
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
      const Point &position() const {
          return fetch().position;
      }

      /** Return this node's index, a number in the range [0, graph_size). */
      size_type index() const {
          return fetch().index;
      }

      /** Return this node's value, in template variable */
      node_value_type& value() {
          return fetch().node_value;
      }

      /** Const return this node's value, in template variable */
      const node_value_type& value() const {
          const node_value_type node_value = fetch().node_value;
          return node_value;
      }

      /** Return the number of incident edges */
      size_type degree() const {
          return fetch().incident_edges_.size();
      }

      /** Start of the incident iterator. */
      incident_iterator edge_begin() const{
          return IncidentIterator(0, this);
      }

      /** End of the incident iterator. */
      incident_iterator edge_end() const{
          return IncidentIterator(this->degree(), this);
      }

      /** Test whether this node and @a n are equal.
       *
       * Equal nodes have the same graph and the same index.
       */
      bool operator==(const Node &n) const {
          if ((fetch().index == n.index()) and (graph_ == n.graph_)) { // n.graph_ i'm a little unsure about
              return true;
          }
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
      bool operator<(const Node &n) const {
          if (fetch().index < n.index()) {
              return true;
          }
          return false;
      }

  private:
      // Allow Graph to access Node's private member data and functions.
      friend class Graph;

      //pointer to graph this node belongs to
      graph_type* graph_;

      //node's unique id
      size_type uid_;

      /** Private Constructor */
      Node(const graph_type *graph, size_type uid)
              : graph_(const_cast<graph_type *>(graph)), uid_(uid) {
      }

      /** Helper method to return the appropriate element.
       * This uses vector indexing
       */
      internal_node& fetch() const {
          return *(graph_->internal_nodes_.at(uid_));
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
       * @param[in] node value, the new node's value
       * @post new num_nodes() == old num_nodes() + 1
       * @post result_node.index() == old num_nodes()
       *
       * Complexity: O(1) amortized operations.
       */
      Node add_node(const Point& position, const node_value_type& node_value = node_value_type()) {
          //add new node to internal nodes
          internal_nodes_.push_back(new internal_node(size_, position, node_value));
          size_++;
          return {this, size_ - 1};        // Invalid node
      }

      /** Determine if a Node belongs to this Graph
       * @return True if @a n is currently a Node of this Graph
       *
       * Complexity: O(1).
       */
      bool has_node(const Node &n) const {
          if (n.index() < size_) {
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
          return {this, i};        // proxy to node
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
    /** Construct an invalid Edge. */
    Edge() {}

    /** Return this edge's index, a number in the range [0, num_edges). */
    size_type index() const {
      return fetch_index();
    }

    /** Return a node of this Edge */
    Node node1() const {
            return *node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
            return *node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->graph_ == e.graph_)) {
          if ((this->node1_ == e.node1_) and (this->node2_ == e.node2_)) {
              return true;
          }

          if ((this->node2_ == e.node1_) and (this->node1_ == e.node2_)) {
              return true;
          }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if ((this->node1_ == e.node2_) and (this->node2_ == e.node1_)){
          return false;
      }

      if (this->node1_ == e.node1_){
          return (this->node2_ < e.node2_);
      }
      return (this->node1_ < e.node1_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    //pointer to graph this node belongs to
    graph_type* graph_;

    //edges's node1 and node2
    const Node* node1_;
    const Node* node2_;

    /** Private Constructor */
    Edge(const graph_type* graph, const Node* node1, const Node* node2)
          : graph_(const_cast<graph_type *>(graph)), node1_(node1), node2_(node2) {
    }

    /** Helper method to return the appropriate edge.
    * This uses map key
    */
    size_type fetch_index() const {
      return this->graph_->internal_edges_.at(Edge(graph_, node1_, node2_));
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_distinct_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edges_vec_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return internal_edges_.count(Edge(this, &a, &b))>=1 | internal_edges_.count(Edge(this, &b, &a))>=1;
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
    //find the existing index
    if (this->has_edge(a,b)){
        return {this, &a, &b};
    }
    else{
        //add the edge to our map of internal edges and vector of edges
        internal_edges_.insert({Edge(this, &a, &b), edges_vec_.size()});
        edges_vec_.push_back(Edge(this, &a, &b));
        num_distinct_edges_++;

        //add the edge to the incident list of both nodes
        internal_nodes_[a.index()]->incident_edges_.push_back(Edge(this, &a, &b));
        internal_nodes_[b.index()]->incident_edges_.push_back(Edge(this, &b, &a));
        return {this, &a, &b};        // Proxy to edge
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      //clear nodes
      this->internal_nodes_.clear();
      this->size_ = 0;

      //clear edges
      this->num_distinct_edges_ = 0;
      this->internal_edges_.clear();
  }

  //
  // Node Iterator
  //

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

    /** dereference operator*/
    Node operator* () const
    {
        return(Node(this->graph, iterator_index));
    }

    /** increment (++) operator*/
    NodeIterator operator++(){
        iterator_index++;
        return *this;
    }

    bool operator==(const NodeIterator& other_iterator) const{
        return (iterator_index == other_iterator.iterator_index & graph == other_iterator.graph);
    }

   private:
    friend class Graph;

    int iterator_index;
    const graph_type* graph;

    /** private constructor that can be accessed by graph class */
    NodeIterator(int iterator_index_, const graph_type* graph_)
        :iterator_index(iterator_index_), graph(graph_) {}
  };

  /** get's the first iterator object */
  node_iterator node_begin() const{
      return NodeIterator(0, this);
  }

  /** get's the end (e.g. last+1) iterator object */
  node_iterator node_end() const{
      return NodeIterator(this->internal_nodes_.size(), this);
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

    //dereference operator
    Edge operator*() const{
        return (node->graph_->internal_nodes_[node->index()]->incident_edges_[current_edge]);
    }

    //increment operator that goes 1 beyond incident edges
    IncidentIterator& operator++(){
        if (current_edge < node->graph_->internal_nodes_[node->index()]->incident_edges_.size()){
            ++current_edge;
        }
        return (*this);
    }

    //Equal if current edge is the same and iterator over the same node
    bool operator==(const IncidentIterator& other) const{
        return ((current_edge == other.current_edge) & (node == other.node));
    }

   private:
    friend class Graph;
    //pointer to edge
    size_type current_edge; //pointer to edge element within internal_node.incident_edges

    //the node of this iterator
    const Node* node; //contains index (so we can use internal_nodes_[node.index()].incident_edges

    //create a private constructor
    IncidentIterator(int current_edge_, const Node* node_)
        : current_edge(current_edge_), node(node_) {};
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
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

    /** dereference operator*/
    Edge operator* () const
    {
      return(this->graph->edges_vec_[edge_index]);
    }

    /** increment (++) operator*/
    EdgeIterator& operator++(){
      edge_index++;
      return *this;
    }

    bool operator==(const EdgeIterator& other_iterator) const{
      return (edge_index == other_iterator.edge_index & graph == other_iterator.graph);
    }


   private:
    friend class Graph;

    // keeps the index of the edge
    int edge_index;

    // keeps the graph
    const graph_type* graph;

    /** private constructor that can be accessed by graph class */
    EdgeIterator(int edge_index_, const graph_type* graph_)
          :edge_index(edge_index_), graph(graph_) {}
  };

  /** get's the first iterator object */
  edge_iterator edge_begin() const{
      return EdgeIterator(0, this);
  }

  /** get's the end (e.g. last+1) iterator object */
  edge_iterator edge_end() const{
      return EdgeIterator(this->edges_vec_.size(), this);
  }

 private:
  //private object that contains all info on node with constructor
  struct internal_node {
      size_type index;   // index of the node
      Point position;    // position of the node
      node_value_type node_value; // node value

      //vector of edges that incident this node
      std::vector<Edge> incident_edges_;

      internal_node(size_type index_, Point position_, const node_value_type& node_value_ = node_value_type())
                :index(index_), position(position_), node_value(const_cast<node_value_type &>(node_value_)) {
      };
  };

  //number of nodes
  size_type size_;

  //number of distinct edges
  size_type num_distinct_edges_;

  //vector that contains pointers to internal node objects
  std::vector<internal_node*> internal_nodes_;

  //map that contains edges as keys and index as values
  std::map<Edge, size_type> internal_edges_;

  //vector of edges for faster edge(i)
  std::vector<Edge> edges_vec_;
};

#endif // CME212_GRAPH_HPP

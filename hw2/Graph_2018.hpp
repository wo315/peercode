
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
template <typename V, typename E = int>
class Graph {
 private:

  // HW0: YOUR CODE HERE


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  typedef V node_value_type;
  typedef E edge_value_type;
  //using node_value_type = V;
  //using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
    : nodes_(), size_(0), edge_node1_(), edge_node2_(), size_edge_(0), 
      adj_list_(), adj_list_index_(), nodes_values_(), edges_values_(), 
      node_removal_(), num_nodes_rem_(), edge_removal_(), num_edges_rem_(),
      node_pop_(), edge_pop_(), edge_node1_pop_(), edge_node2_pop_()
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return g_->nodes_[uid_];
    }

    Point& position() {
      // HW0: YOUR CODE HERE
      return g_->nodes_[uid_];
    }

    /** Return this node's user facing index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_user_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    //* Store value of node by indexing. */
    node_value_type& value() {
       return g_->nodes_values_[uid_];
    }

    //* Store value of node as const by indexing. */ 
    const node_value_type& value() const {
      return g_->nodes_values_[uid_];
    }

    //* Store the degree of a node, i.e., number of incident edges. */
    size_type degree() const {
      std::vector<int> edge_indices =  g_->adj_list_index_[uid_];
      size_type deg = 0;
      for(auto ii = edge_indices.begin(); ii != edge_indices.end(); ++ii) {
        if(g_->edge_removal_[*ii] != -1) {
          deg += 1;
        }
      }
      return deg;
    }

    //* Iterator to start at edge with index 0 for vertices incident to node. */
    IncidentIterator edge_begin() const {
      std::vector<size_type> edge_indices =  g_->adj_list_index_[uid_];
      size_type i=0;
      for(auto ii = edge_indices.begin(); ii != edge_indices.end(); ++ii) {
        if(g_->edge_removal_[*ii] != -1) {
          return IncidentIterator(g_, uid_, g_->adj_list_index_[uid_], g_->adj_list_[uid_], g_->edge_removal_, g_->node_removal_, i);
        }
        i += 1;
      }
      return edge_end();
    }

    //* Iterator to end at edge with index one beyond the last edge for vertices incident to node. */
    IncidentIterator edge_end() const {
      std::vector<size_type> edge_indices =  g_->adj_list_index_[uid_];
      return IncidentIterator(this->g_, uid_, edge_indices, g_->adj_list_[uid_], g_->edge_removal_, g_->node_removal_, edge_indices.size());
    }    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (this->g_ == n.g_ && this->uid_ == n.uid_){
        return true;
      }
      else{
        return false;
      }
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
      if (this->uid_ < n.uid_) {
        return true;
      }
      if (this->g_ != n.g_ && n.uid_ == this->uid_){
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    //Pointer back to graph container
    Graph* g_;
    //This Node's unique id
    size_type uid_;
    //This Node's user facing indx
    size_type uid_user_;
    //Private Constructor
    Node(const Graph* g, size_type uid, size_type uid_user)
      : g_(const_cast<Graph*>(g)), uid_(uid), uid_user_(uid_user) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& val_= node_value_type()) {
    // HW0: YOUR CODE HERE
    //Append to set of nodes in Graph
    nodes_.emplace_back(position);
    node_removal_.emplace_back(size_);
    node_pop_.emplace_back(num_nodes_rem_ + size_);
    //Update size and index
    size_ += 1;
    //Update adjacency list with node andn index of node
    adj_list_.push_back(std::vector<size_type> ());
    adj_list_index_.push_back(std::vector<size_type> ());
    //Update nodes values
    nodes_values_.push_back(val_);
    return Node(this, num_nodes_rem_ + size_ - 1, size_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.g_ && node_removal_[n.uid_] != -1) {
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
    assert(i<num_nodes());
    return Node(this, node_pop_[i], node_removal_[node_pop_[i]]);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(g__, uid_1, uid_user1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(g__, uid_2, uid_user2_);
    }

    //* Store value of edge by indexing. */
    edge_value_type& value() {
       return g__->edges_values_[uid_edge_];
    }

    //* Store value of edge as const by indexing. */
    const edge_value_type& value() const {
      return g__->edges_values_[uid_edge_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (this->g__ == e.g__ && this->uid_edge_ == e.uid_edge_) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (this->uid_edge_ < e.uid_edge_) {
        return true;
      }
      if (this->g__ != e.g__ && e.uid_edge_ == this->uid_edge_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* g__;
    //Pointer to nodes 1 and 2
    size_type uid_1;
    size_type uid_2;
    //User Pointers to nodes 1 and 2
    size_type uid_user1_;
    size_type uid_user2_;
    //Pointer to edge
    size_type uid_edge_;
    //Private Constructor
    Edge(const Graph* g, size_type uid1, size_type uid2, size_type uid_edge, size_type uid_user1, size_type uid_user2)
      : g__(const_cast<Graph*>(g)), uid_1(uid1), uid_2(uid2), uid_edge_(uid_edge), uid_user1_(uid_user1), uid_user2_(uid_user2) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return size_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    return Edge(this, edge_node1_pop_[i], edge_node2_pop_[i], edge_pop_[i], node_removal_[edge_node1_pop_[i]], node_removal_[edge_node2_pop_[i]]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b)) {
      std::vector<size_type> small_vec = adj_list_[a.uid_];
      for (size_type j = 0; j<small_vec.size(); ++j) {
        if (small_vec[j] == b.uid_ && edge_removal_[adj_list_index_[a.uid_][j]] != -1) {
          return true;
        }
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val_ = edge_value_type()) {
    // HW0: YOUR CODE HERE
    if (has_node(a) && has_node(b) && !(a==b)) {
      if (!(has_edge(a, b))) {
        edge_node1_.push_back(a.uid_);
        edge_node2_.push_back(b.uid_);
        adj_list_[a.uid_].push_back(b.uid_);
        adj_list_[b.uid_].push_back(a.uid_);
        adj_list_index_[a.uid_].push_back(size_edge_ + num_edges_rem_);
        adj_list_index_[b.uid_].push_back(size_edge_ + num_edges_rem_);
        edges_values_.push_back(val_);
        edge_removal_.emplace_back(size_edge_);
        edge_pop_.emplace_back(size_edge_+num_edges_rem_);
        edge_node1_pop_.emplace_back(a.uid_);
        edge_node2_pop_.emplace_back(b.uid_);
        size_edge_ += 1;
        return Edge(this, a.uid_, b.uid_, num_edges_rem_ + size_edge_ - 1, a.uid_user_, b.uid_user_);
      }
      else {
        std::vector<size_type> small_vec = adj_list_[a.uid_];
        for (size_type k = 0; k<small_vec.size(); ++k) {
          if (small_vec[k] == b.uid_) {
            size_type tem = adj_list_index_[a.uid_][k];
            return Edge(this, a.uid_, b.uid_, tem, a.uid_user_, b.uid_user_);
          }
        }
      }
      
    }
    return Edge();        // Invalid Edge
  }

  /** Remove edge from graph if it exists otherwise return 0.
  * @pre @a a and @a b are distinct valid nodes of this graph
  * @pre has_edge(@a a, @a b) == true
  * @post has_edge(@a a, @a b) == false
  * @post the total number of edges decreases by 1
  * Complexity: at most O(num_nodes())
  */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_node(a) && has_node(b) && !(a==b)) {
      if (has_edge(a, b)) {
        for (int i = 0; i<adj_list_[a.uid_].size(); ++i) {
          if (adj_list_[a.uid_][i] == b.uid_) {
            size_type curr_edge_idx = adj_list_index_[a.uid_][i];
            size_type temp = edge_pop_[edge_pop_.size()-1];
            std::swap(edge_pop_[edge_removal_[curr_edge_idx]], edge_pop_[edge_pop_.size()-1]);
            std::swap(edge_node1_pop_[edge_removal_[curr_edge_idx]], edge_node1_pop_[edge_node1_pop_.size()-1]);
            std::swap(edge_node2_pop_[edge_removal_[curr_edge_idx]], edge_node2_pop_[edge_node2_pop_.size()-1]);
            edge_removal_[temp] = edge_removal_[curr_edge_idx];
            edge_pop_.pop_back();
            edge_node1_pop_.pop_back();
            edge_node2_pop_.pop_back();
            edge_removal_[curr_edge_idx] = -1;
          }
        }
      size_edge_ -= 1;
      num_edges_rem_ += 1;
      return 1;
      }
    }
    return 0;
  }

  /** Remove edge from graph if it exists otherwise return 0.
  * @pre @a a is a valid edge of this graph
  * @pre has_edge(@a e.node1(), @a e.node2()) == true
  * @post has_edge(@a e.node1(), @a e.node2()) == false
  * @post the total number of edges decreases by 1
  * Complexity: at most O(num_nodes())
  */
  size_type remove_edge(const Edge& e) {
    if (this == e.g__) {
      remove_edge(e.node1(), e.node2());
      return 1;
    }
    return 0;
  }

  /** Remove node and all incident edges from graph, if node in graph otherwise return 0.
  * @pre @a a is a valid node of this graph
  * @pre has_node(@a a) == true
  * @post has_node(@a a) == false
  * @post all edges incident to @a are removed
  * @post the total number of nodes decreases by 1
  * Complexity: at most O(num_nodes())
  */
  size_type remove_node(const Node& n) {
    //std::cout << n.index() << std::endl;
    if (has_node(n)) {
      for (size_type i = 0; i<adj_list_[n.uid_].size(); ++i) {
        //std::cout << edge_removal_[adj_list_index_[n.uid_][i]] << num_edges() << std::endl;
        if (edge_removal_[adj_list_index_[n.uid_][i]] != -1){
          remove_edge(EdgeIterator(this, edge_removal_[adj_list_index_[n.uid_][i]]));
        }
      }
      size_type temp = node_pop_[node_pop_.size()-1];
      std::swap(node_pop_[node_removal_[n.uid_]], node_pop_[node_pop_.size()-1]);
      node_removal_[temp] = node_removal_[n.uid_];
      node_pop_.pop_back();
      node_removal_[n.uid_] = -1;
      size_ -= 1;
      num_nodes_rem_ += 1;
      return 1;
    }
    return 0;
  }

  /** Change the node_iterator to reflect a removed node, and output the new iterator.
  * Swap the node to be removed with the last node and return the corresponding iterator
  * In the call to size_type remove_node(Node) we will pop off this last node, which will then
  * change the end element of the iterator
  */
  node_iterator remove_node(node_iterator n_it) {
    node_type n_ = *n_it;
    size_type node_removed = remove_node(n_);
    if (node_removed == 0){return n_it;}
    if (num_nodes() == 0){return node_end();}

    return n_it;
  }

  /** Change the edge_iterator to reflect a removed edge, and output the new iterator.
  * Swap the edge to be removed with the last edge and return the corresponding iterator
  * In the call to size_type remove_edge(edge) we will pop off this last edge, which will then
  * change the end element of the iterator
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    node_type node1 = (*e_it).node1();
    node_type node2 = (*e_it).node2();
    size_type edge_removed = remove_edge(node1, node2);
    if (edge_removed == 0){return e_it;}
    if (num_edges() == 0){return edge_end();}

    return e_it;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    size_ = 0;
    edge_node1_.clear();
    edge_node2_.clear(); 
    size_edge_ = 0;
    adj_list_.clear();
    adj_list_index_.clear();
    nodes_values_.clear();
    edges_values_.clear();
    node_removal_.clear();
    edge_removal_.clear();
    num_nodes_rem_ = 0;
    num_edges_rem_ = 0;
    node_pop_.clear();
    edge_pop_.clear();
    edge_node1_pop_.clear();
    edge_node2_pop_.clear();
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

    // HW1 #2: YOUR CODE HERE

   /** Increment a node by one. */
   NodeIterator& operator++() {
     ++node_counter;
     return *this;
   }

   /** Check if two nodes are equal by checking their index and graph. */
   bool operator==(const NodeIterator& node_iter) const {
     return g_node_iter == node_iter.g_node_iter && node_iter.node_counter == node_counter;
   }

   /** Dereference the pointer to a node. */
   Node operator*() const {
     return g_node_iter->node(node_counter);
   }

   private:
    friend class Graph;
    Graph* g_node_iter;
    size_type node_counter;
    NodeIterator(const Graph* g_node_iter_, size_type node_counter_) {
      node_counter = node_counter_;
      g_node_iter = const_cast<Graph*>(g_node_iter_);
    }
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE

  /** Iterator to start at the node with index 0. */
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /** Iterator to end at node with index one beyond the last node. */
  NodeIterator node_end() const {
    return NodeIterator(this, size());
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

   /** Dereference the pointer to an edge. */
   Edge operator*() const {
     return Edge(g_inc_iter, node_real_id, connect_nodes[edge_node_ptr], edge_indices[edge_node_ptr], 
                 node_removal[node_real_id], node_removal[connect_nodes[edge_node_ptr]]);
   }

   /** Increment an incident edge by one. */
   IncidentIterator& operator++() {
      ++edge_node_ptr;
      if (edge_node_ptr < connect_nodes.size()){
        while (edge_removal[edge_indices[edge_node_ptr]] == -1) {
          ++edge_node_ptr;
          if (edge_node_ptr == connect_nodes.size()){break;}
        }
      }
      return *this;
   }

   /** Check if two edges are equal by checking if the pointer, graphs and nodes are the same. */
   bool operator==(const IncidentIterator& inc_iter) const {
      return inc_iter.edge_node_ptr == edge_node_ptr && g_inc_iter == inc_iter.g_inc_iter && node_real_id == inc_iter.node_real_id;
   }

   private:
    friend class Graph;
    Graph* g_inc_iter;
    size_type edge_node_ptr;
    std::vector<size_type> edge_indices;
    std::vector<size_type> connect_nodes;
    std::vector<float> edge_removal;
    size_type node_real_id;
    std::vector<float> node_removal;
    IncidentIterator(const Graph* g_inc_iter_, size_type node_real_id_, std::vector<size_type> edge_indices_,
                     std::vector<size_type> connect_nodes_, std::vector<float> edge_removal__, std::vector<float> node_removal__, size_type edge_node_ptr_) {
      edge_node_ptr = edge_node_ptr_;
      g_inc_iter = const_cast<Graph*>(g_inc_iter_);
      edge_indices = edge_indices_;
      connect_nodes = connect_nodes_;
      edge_removal = edge_removal__;
      node_removal = node_removal__;
      node_real_id = node_real_id_;
    }
    // HW1 #3: YOUR CODE HERE
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

   /** Dereference the pointer to an edge. */
   Edge operator*() const {
     return g_edge_iter->edge(edge_counter);
   }

   /** Increment an edge by one. */
   EdgeIterator& operator++() {
     ++edge_counter;
     return *this;
   }

   /** Check if two edges are equal by comparing the graph and edge index. */
   bool operator==(const EdgeIterator& edge_val) const {
     return edge_counter == edge_val.edge_counter && g_edge_iter == edge_val.g_edge_iter;
   }

   private:
    friend class Graph;
    Graph* g_edge_iter;
    int edge_counter;
    EdgeIterator(const Graph* g_edge_iter_, size_type edge_counter_) {
      edge_counter = edge_counter_;
      g_edge_iter = const_cast<Graph*>(g_edge_iter_);
    }
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE

  /** Iterator to start at the edge with index 0. */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /** Iterator to end at edge with index one beyond the last edge. */
  EdgeIterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  std::vector<Point> nodes_; //Store points in graph
  size_type size_; //Store number of nodes
  std::vector<size_type> edge_node1_; //Store first node in edge
  std::vector<size_type> edge_node2_; //Store second node in edge
  size_type size_edge_; //Store number of edges
  std::vector<std::vector<size_type>> adj_list_; //Store adjacent nodes
  std::vector<std::vector<size_type>> adj_list_index_; //Store adjacent edges
  std::vector<node_value_type> nodes_values_; //Store values of nodes
  std::vector<edge_value_type> edges_values_; //Store values of edges
  std::vector<float> node_removal_; //Store indices of user facing nodes
  size_type num_nodes_rem_; //Store number of removed nodes
  std::vector<float> edge_removal_; //Store "indices" of user facing edges
  size_type num_edges_rem_; //Store number of removed edges
  std::vector<size_type> node_pop_; //Store popped array of nodes
  std::vector<size_type> edge_pop_; //Store popped array of edges
  std::vector<size_type> edge_node1_pop_; //Store popped array of adjacent nodes1
  std::vector<size_type> edge_node2_pop_; //Store popped array of adjacent nodes2
};

#endif // CME212_GRAPH_HPP

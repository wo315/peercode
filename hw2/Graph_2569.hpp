#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @struct pair_hash
 * @brief provided by https://stackoverflow.com/questions/32685540/why-cant-i-compile-an-unordered-map-with-a-pair-as-key
 * useful for unordered_map hash function.
 */
struct pair_hash
{
  template <class T1, class T2>
  std::size_t operator()(const std::pair<T1, T2> &p) const
  {
    auto h1 = std::hash<T1>{}(p.first);
    auto h2 = std::hash<T2>{}(p.second);

    // Mainly for demonstration purposes, i.e. works but is overly simple
    // In the real world, use sth. like boost.hash_combine
    return h1 ^ h2;
  }
};

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph
{
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
  /** Type of user defined value carried by the nodes. */
  using node_value_type = V;

  /** Type of user defined value carried by the edges. */
  using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<node_value_type, edge_value_type>;

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

  /** Type of the iterator on the nodes_ attribute of the graph.
   * Support the implementation of the node_iterator */
  using node_map_iterator = typename std::unordered_map<size_type, node_type>::const_iterator;
  using index_iterator = typename std::unordered_map<size_type, size_type>::const_iterator;

  /** Type of the iterator on the edges_ attribute of the graph.
   * Support the implementation of the edge_iterator */
  using edge_map_iterator = typename std::unordered_map<size_type, edge_type>::const_iterator;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : nodes_{},
        points_{}, edges_{}, nodes_to_edges_{}, size_(0),
        num_edges_(0), next_node_index_(0), next_edge_index_(0)
  {
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
  class Node : private totally_ordered<Node>
  {
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
    Node()
    {
      // HW0: YOUR CODE HERE
    }

    /**
     * @brief Provide a reference to the position of the node
     * @return Node's position (reference)
     */
    const Point &position() const
    {
      // HW0: YOUR CODE HERE
      node_type &node = fetch();
      return node.graph_->points_.at(node.index_);
      // return Point();
    }

    /**
     * @brief Provide a reference to the position of the node.
     * @return Node's position (reference).
     */
    Point &position()
    {
      node_type &node = fetch();
      return node.graph_->points_.at(node.index_);
    }

    /**
     * @brief Provide the index to the node.
     * @post Result is in [0, graph_size).
     * @return Node's index.
     */
    size_type index() const
    {
      // HW0: YOUR CODE HERE
      return this->graph_->hidden_index_to_node_index_.at(fetch().index_);
    }

    /**
     * @brief Provide the hidden index to the node.
     * @post Result is not necessarily in [0, graph_size). In fact, as
     * soon as a node gets removed, the hidden index and the index are desynchronized.
     * @return Node's hidden index.
     */
    size_type hidden_index() const
    {
      return fetch().index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /**
     * @brief Computes the degree of a node.
     * 
     * @return Number of adjacent nodes to the current node
     * 
     * @post node.degree() == size({other_node for other_node in Nodes(graph) if (node, other_node) in Edges(graph)})
     */
    size_type degree() const
    {
      return (*(fetch().graph_)).nodes_to_incidents_.at(index_).size();
    }

    /**
     * @brief Provides an iterator over the edges of a node.
     * 
     * @return incident_iterator, which iterates over all the adjacent nodes of a node.
     * 
     * @post number_adjacent_nodes = 0; 
     *       for(auto ei = node.edge_begin(); ei != node.edge_end(); ++ei){
     *          number_adjacent_nodes += 1;
     *          (*ei).node1() == node;
     *          (node, (*ei).node2()) is in Edges(graph);
     *        }
     *        number_adjacent_nodes == node.degree();
     */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(graph_, graph_->nodes_to_incidents_.at(index_).begin(), &fetch());
    }

    /**
     * @brief Provides an end iterator over the edges of a node.
     * 
     * @return incident_iterator, which is equal to the end iterator.
     * 
     * @post number_adjacent_nodes = 0; 
     *       for(auto ei = node.edge_begin(); ei != node.edge_end(); ++ei){
     *          number_adjacent_nodes += 1;
     *          (*ei).node1() == node;
     *          (node, (*ei).node2()) is in Edges(graph);
     *        }
     *        number_adjacent_nodes == node.degree();
     */
    incident_iterator edge_end() const
    {
      return IncidentIterator(graph_, graph_->nodes_to_incidents_.at(index_).end(), &fetch());
    }

    /** @brief Outputs the node's value.
     * 
     * @return node_value_type&
     * 
     * @post Result == node.value_;
     */
    node_value_type &value()
    {
      node_type &node = fetch();
      return node.graph_->node_values_.at(node.index_);
    }

    /** @brief Outputs the node's value.
     * 
     * @return const node_value_type&
     * 
     * @post Result == node.value_;
     */
    const node_value_type &value() const
    {
      node_type &node = fetch();
      return node.graph_->node_values_.at(node.index_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      // HW0: YOUR CODE HERE
      node_type &node = fetch();
      return node.graph_ == n.graph_ && node.index_ == n.index_;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      // HW0: YOUR CODE HERE
      return fetch().index_ < n.index_;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type *graph_;
    size_type index_;

    /** @brief Initializes a node object.
     * 
     * @param[in] graph, pointer to the graph of the node.
     * @param[in] index, index of the node in its graph.
     * @param[in] value, value of the node, of type node_value_type.
     */
    Node(graph_type *graph, size_type index) : graph_(graph), index_(index) {}

    /** @brief fetch the node object content.
     * 
     * @return A value of type node_type.
     * 
     * @post node == node.fetch()
     */
    node_type &fetch() const { return graph_->nodes_.at(index_); }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    // HW0: YOUR CODE HERE
    return this->size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type &value = node_value_type())
  {
    // HW0: YOUR CODE HERE
    node_type node = Node(this, this->next_node_index_);
    std::unordered_set<size_type> incidents;
    this->nodes_.insert({this->next_node_index_, node});
    this->points_.insert({this->next_node_index_, position});
    this->node_values_.insert({this->next_node_index_, value});
    this->nodes_to_incidents_.insert({this->next_node_index_, incidents});
    this->hidden_index_to_node_index_.insert({this->next_node_index_, this->hidden_index_to_node_index_.size()});
    this->node_index_to_hidden_index_.insert({this->node_index_to_hidden_index_.size(), this->next_node_index_});

    this->next_node_index_++;
    this->size_++;

    return node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    // HW0: YOUR CODE HERE
    index_iterator node_found = this->hidden_index_to_node_index_.find(n.index_);
    if (node_found == this->hidden_index_to_node_index_.end())
      return false;
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    // HW0: YOUR CODE HERE
    return this->nodes_.at(node_index_to_hidden_index_.at(i));
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
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      // HW0: YOUR CODE HERE
      return *(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      // HW0: YOUR CODE HERE
      return *(node2_);
    }

    /**
     * @brief Provide the index to the edge.
     * @post Result is in [0, graph_size).
     * @return Edge's index.
     */
    size_type index() const
    {
      // HW0: YOUR CODE HERE
      return this->graph_->hidden_index_to_edge_index_.at(fetch().index_);
    }

    /**
     * @brief Provide the value to the edge.
     * @return Edge's value, reference.
     */
    edge_value_type &value()
    {
      edge_type &edge = fetch();
      return edge.graph_->edge_values_.at(edge.index_);
    }

    /** @brief Outputs the node's value.
     * 
     * @return const edge_value_type&
     * 
     * @post Result == node.value_;
     */
    const edge_value_type &value() const
    {
      edge_type &edge = fetch();
      return edge.graph_->edge_values_.at(edge.index_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      // HW0: YOUR CODE HERE
      edge_type edge = fetch();
      return ((*(edge.node1_) == *(e.node1_)) && (*(edge.node2_) == *(e.node2_))) ||
             ((*(edge.node1_) == *(e.node2_)) && (*(edge.node2_) == *(e.node1_)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      // HW0: YOUR CODE HERE
      return index_ < e.index_ || graph_ < e.graph_;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type *graph_;
    size_type index_;
    node_type *node1_;
    node_type *node2_;

    /**
     * @brief Edge constructor.
     */
    Edge(graph_type *graph, size_type index, const node_type *node1,
         const node_type *node2)
        : graph_(graph), index_(index), node1_(const_cast<node_type *>(node1)),
          node2_(const_cast<node_type *>(node2)) {}

    edge_type &fetch() const { return graph_->edges_[index_]; }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    // HW0: YOUR CODE HERE
    return this->num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    // HW0: YOUR CODE HERE

    return this->edges_.at(edge_index_to_hidden_index_.at(i));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    // HW0: YOUR CODE HERE
    auto edge_found =
        this->nodes_to_edges_.find(std::minmax(a.index_, b.index_));
    if (edge_found == this->nodes_to_edges_.end())
      return false;
    return true;
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
  Edge add_edge(const Node &a, const Node &b)
  {
    if (has_edge(a, b))
    {
      auto edge_found =
          this->nodes_to_edges_.find(std::minmax(a.index_, b.index_));
      return Edge(this, edge_found->second, &a, &b);
    }
    edge_type edge = Edge(this, this->next_edge_index_, &(this->nodes_.at(a.index_)), &(this->nodes_.at(b.index_)));
    this->edges_.insert({this->next_edge_index_, edge});
    this->nodes_to_edges_.insert(
        {std::minmax(a.index_, b.index_), this->next_edge_index_});
    this->edge_values_.insert({this->next_edge_index_, edge_value_type()});
    this->nodes_to_incidents_[a.index_].insert(this->next_edge_index_);
    this->nodes_to_incidents_[b.index_].insert(this->next_edge_index_);
    this->hidden_index_to_edge_index_.insert({this->next_edge_index_, this->hidden_index_to_edge_index_.size()});
    this->edge_index_to_hidden_index_.insert({this->edge_index_to_hidden_index_.size(), this->next_edge_index_});

    this->next_edge_index_++;
    this->num_edges_++;
    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    // HW0: YOUR CODE HERE
    this->nodes_.clear();
    this->edges_.clear();
    this->points_.clear();
    this->nodes_to_edges_.clear();
    this->node_values_.clear();
    this->hidden_index_to_node_index_.clear();
    this->node_index_to_hidden_index_.clear();
    this->edge_values_.clear();
    this->hidden_index_to_edge_index_.clear();
    this->edge_index_to_hidden_index_.clear();
    this->nodes_to_incidents_.clear();
    this->size_ = 0;
    this->next_edge_index_ = 0;
    this->next_node_index_ = 0;
    this->num_edges_ = 0;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** @brief dereference operator.
     * 
     * @return node at the current iterator position.
     * 
     * @post Result is in Nodes(graph)
     */
    Node operator*() const
    {
      return graph_->nodes_.at(iterator_->second);
    }

    /** @brief incrementation operator.
     * 
     * @return node at the current iterator position.
     * 
     * @post Result is in Nodes(graph).
     * @post *Result is changes each time we use ++.
     */
    NodeIterator &operator++()
    {
      iterator_++;
      return *this;
    }

    /** @brief equality operator.
     * 
     * @return true if the two iterators are equal, otherwise false.
     * 
     * @post node_iterator == node_iterator returns true.
     */
    bool operator==(const NodeIterator &node_iterator) const
    {
      return iterator_ == node_iterator.iterator_;
    }

    /** @brief inequality operator.
     * 
     * @return true if the two iterators are different, otherwise true.
     * 
     * @post node_iterator != node_iterator returns false.
     */
    bool operator!=(const NodeIterator &node_iterator) const
    {
      return iterator_ != node_iterator.iterator_;
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    index_iterator iterator_;
    graph_type *graph_;

    /** @brief Constructor of Node iterator.
     * 
     * @param[in] @a index_iterator An iterator over graph.nodes_.
     * @param[in] @a graph a pointer to the graph.
     */
    NodeIterator(const index_iterator &iterator, graph_type *graph) : iterator_(iterator), graph_(graph) {}
    NodeIterator(index_iterator &iterator, graph_type *graph) : iterator_(iterator), graph_(graph) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** @brief Provides an iterator over the nodes of the graph.
   * 
   * @return node_iterator.
   * 
   * @post number_nodes = 0
   *       for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
   *          number_nodes++;
   *          *ni is in Nodes(graph);
   *       }
   *       number_nodes == graph.num_nodes();
   */
  node_iterator node_begin() const
  {
    return NodeIterator(this->node_index_to_hidden_index_.begin(), const_cast<graph_type *>(this));
  }

  /** @brief Provides an iterator over the nodes of the graph.
   * 
   * @return node_iterator.
   * 
   * @post number_nodes = 0
   *       for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
   *          number_nodes++;
   *          *ni is in Nodes(graph);
   *       }
   *       number_nodes == graph.num_nodes();
   */
  node_iterator node_begin()
  {
    return NodeIterator(this->node_index_to_hidden_index_.begin(), this);
  }

  /** @brief Provides an end iterator for the nodes of the graph.
   * 
   * @return end node_iterator.
   * 
   * @post number_nodes = 0
   *       for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
   *          number_nodes++;
   *          *ni is in Nodes(graph);
   *       }
   *       number_nodes == graph.num_nodes();
   */
  node_iterator node_end() const
  {
    return NodeIterator(this->node_index_to_hidden_index_.end(), const_cast<graph_type *>(this));
  }

  /** @brief Provides an end iterator for the nodes of the graph.
   * 
   * @return end node_iterator.
   * 
   * @post number_nodes = 0
   *       for(auto ni = graph.node_begin(); ni != graph.node_end(); ++ni){
   *          number_nodes++;
   *          *ni is in Nodes(graph);
   *       }
   *       number_nodes == graph.num_nodes();
   */
  node_iterator node_end()
  {
    return NodeIterator(this->node_index_to_hidden_index_.end(), this);
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** @brief dereference operator.
     * 
     * @return edge at the current iterator position.
     * 
     * @post Result is in Edges(graph)
     * @post Result.node1() == node
     */
    Edge operator*() const
    {
      Edge edge = (*graph_).edge(graph_->hidden_index_to_edge_index_.at(*iterator_));
      if (edge.node1() == *node_)
      {
        return Edge(graph_, *iterator_, node_, &(graph_->nodes_.at(edge.node2().hidden_index())));
      }
      else
      {
        return Edge(graph_, *iterator_, node_, &(graph_->nodes_.at(edge.node1().hidden_index())));
      }
    }

    /**
     * @brief increment operator.
     * 
     * @return current iterator.
     * 
     * @post
     * for (auto ei = node.edge_begin(); ei != node.edge_end(); ++ei) {
     *  Edge e = *ei;
     *  Node n1 = e.node1();
     *  node == n1;
     * }
     */
    IncidentIterator &operator++()
    {
      iterator_++;
      return *this;
    }

    /** @brief equality operator.
     * 
     * @return true if the two iterators are equal, otherwise false.
     * 
     * @post incident_iterator == incident_iterator returns true.
     */
    bool operator==(const IncidentIterator &incident_iterator) const
    {
      return iterator_ == incident_iterator.iterator_;
    }

    /** @brief inequality operator.
     * 
     * @return true if the two iterators are different, otherwise false.
     * 
     * @post incident_iterator == incident_iterator returns false.
     */
    bool operator!=(const IncidentIterator &incident_iterator) const
    {
      return iterator_ != incident_iterator.iterator_;
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type *graph_;
    std::unordered_set<size_type>::const_iterator iterator_;
    node_type *node_;

    IncidentIterator(graph_type *graph, const std::unordered_set<size_type>::const_iterator &iterator, node_type *node) : graph_(graph), iterator_(iterator), node_(node) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** @brief dereference operator.
     * 
     * @return edge at the current iterator position.
     * 
     * @post Result is in Edges(graph)
     */
    Edge operator*() const
    {
      return graph_->edges_.at(iterator_->second);
    }

    /** @brief incrementation operator.
     * 
     * @return edge at the current iterator position.
     * 
     * @post Result is in Edges(graph).
     * @post *Result is changes each time we use ++.
     */
    EdgeIterator &operator++()
    {
      iterator_++;
      return *this;
    }

    /** @brief equality operator.
     * 
     * @return true if the two iterators are equal, otherwise false.
     * 
     * @post edge_iterator == edge_iterator returns true.
     */
    bool operator==(const EdgeIterator &edge_iterator) const
    {
      return iterator_ == edge_iterator.iterator_;
    }

    /** @brief inequality operator.
     * 
     * @return true if the two iterators are different, otherwise false.
     * 
     * @post edge_iterator == edge_iterator returns false.
     */
    bool operator!=(const EdgeIterator &edge_iterator) const
    {
      return iterator_ != edge_iterator.iterator_;
    }

  private:
    friend class Graph;
    index_iterator iterator_;
    graph_type *graph_;

    /** @brief Constructor.
     * 
     * @param[in] iterator over graph.edges_.
     * 
     */
    EdgeIterator(const index_iterator &iterator, graph_type *graph) : iterator_(iterator), graph_(graph) {}
  };

  /** @brief Provides an iterator over the edges of the graph.
   * 
   * @return edge_iterator.
   * 
   * @post number_edges = 0
   *       for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
   *          number_edges++;
   *          *ei is in Edges(graph);
   *       }
   *       number_edges == graph.num_edges();
   */
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this->edge_index_to_hidden_index_.begin(), const_cast<graph_type *>(this));
  }

  /** @brief Provides an iterator over the edges of the graph.
   * 
   * @return edge_iterator.
   * 
   * @post number_edges = 0
   *       for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
   *          number_edges++;
   *          *ei is in Edges(graph);
   *       }
   *       number_edges == graph.num_edges();
   */
  edge_iterator edge_begin()
  {
    return EdgeIterator(this->edge_index_to_hidden_index_.begin(), this);
  }

  /** @brief Provides an iterator over the edges of the graph.
   * 
   * @return edge_iterator.
   * 
   * @post number_edges = 0
   *       for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
   *          number_edges++;
   *          *ei is in Edges(graph);
   *       }
   *       number_edges == graph.num_edges();
   */
  edge_iterator edge_end() const
  {
    return EdgeIterator(this->edge_index_to_hidden_index_.end(), const_cast<graph_type *>(this));
  }

  /** @brief Provides an iterator over the edges of the graph.
   * 
   * @return edge_iterator.
   * 
   * @post number_edges = 0
   *       for(auto ei = graph.edge_begin(); ei != graph.edge_end(); ++ei){
   *          number_edges++;
   *          *ei is in Edges(graph);
   *       }
   *       number_edges == graph.num_edges();
   */
  edge_iterator edge_end()
  {
    return EdgeIterator(this->edge_index_to_hidden_index_.end(), this);
  }

  /** @brief Removes a node from the graph. O(num nodes())
   * 
   * @param[in,out] @a node, the node to remove.
   * 
   * @return 1 if the node was removed, and zero otherwise.
   * 
   * @post invalidates iterators over the nodes of the graph, 
   * over the edges (as we also remove the adjacent edges to the removed node),
   * over the incident edges of a node.
   * Also, invalidates the argument @a node, and its adjacent edges.
   * 
   * @post graph.num_nodes() has -1 node.
   * @post graph.num_edges() has -(number of incidents edges to @a node).
   */
  size_type remove_node(const Node &node)
  {
    try
    {
      size_type hidden_index = node.index_;
      size_type node_index = node.index();
      auto edge_it = node.edge_begin();
      while (edge_it != node.edge_end())
      {

        remove_edge(*edge_it);
        edge_it = node.edge_begin();
      }

      this->hidden_index_to_node_index_.erase(hidden_index);
      this->node_index_to_hidden_index_.erase(node_index);

      this->nodes_to_incidents_.erase(hidden_index);

      if (node_index != this->node_index_to_hidden_index_.size())
      {
        size_type last_hidden_index = this->node_index_to_hidden_index_.at(this->node_index_to_hidden_index_.size());
        size_type last_node_index = this->node_index_to_hidden_index_.size();

        this->node_index_to_hidden_index_.erase(last_node_index);
        this->hidden_index_to_node_index_.erase(last_hidden_index);

        this->node_index_to_hidden_index_.insert({node_index, last_hidden_index});
        this->hidden_index_to_node_index_.insert({last_hidden_index, node_index});
      }

      this->size_--;
      return 1;
    }
    catch (const std::exception &e)
    {
      return 0;
    }
  }

  /** @brief Removes a node from the graph. O(num nodes())
   * 
   * @param[in,out] @a n_it, an iterator to the node to remove.
   * 
   * @return A valid iterator over the nodes of the graph.
   * 
   * @post invalidates iterators over the nodes of the graph, 
   * over the edges (as we also remove the adjacent edges to the removed node),
   * over the incident edges of a node.
   * Also, invalidates the argument @a node, and its adjacent edges.
   * 
   * @post graph.num_nodes() has -1 node.
   * @post graph.num_edges() has -(number of incidents edges to @a node).
   */
  node_iterator remove_node(node_iterator n_it)
  {
    *n_it;
    remove_node(*n_it);
    return this->node_begin();
  }

  /** @brief Removes an edge from the graph. O(num nodes() + num edges())
   * 
   * @param[in] @a node1 and @a node2, two nodes that define the edge to remove.
   * 
   * @return 1 if the node was removed, and zero otherwise.
   * 
   * @post invalidates iterators over the edges of the graph, 
   * over the incident edges of a node.
   * Also, invalidates the edge that was removed.
   * 
   * @post graph.num_nodes() has same number of nodes than before.
   * @post graph.num_edges() has -1 edge.
   */
  size_type remove_edge(const Node &node1, const Node &node2)
  {
    if (has_edge(node1, node2))
    {
      size_type hidden_index = this->nodes_to_edges_.at(std::minmax(node1.index_, node2.index_));
      return remove_edge(this->edges_.at(hidden_index));
    }
    else
    {
      return 0;
    }
  }

  /** @brief Removes an edge from the graph. O(num nodes() + num edges())
   * 
   * @param[in,out] @a edge, the edge to remove.
   * 
   * @return 1 if the node was removed, and zero otherwise.
   * 
   * @post invalidates iterators over the edges of the graph, 
   * over the incident edges of a node.
   * Also, invalidates the edge that was removed.
   * 
   * @post graph.num_nodes() has same number of nodes than before.
   * @post graph.num_edges() has -1 edge.
   */
  size_type remove_edge(const Edge &edge)
  {
    try
    {
      size_type hidden_index = edge.index_;
      size_type edge_index = edge.index();

      this->hidden_index_to_edge_index_.erase(hidden_index);
      this->edge_index_to_hidden_index_.erase(edge_index);
      this->nodes_to_edges_.erase(std::minmax(edge.node1().index_, edge.node2().index_));

      auto &node_incidents1 = this->nodes_to_incidents_.at(edge.node1().index_);
      node_incidents1.erase(hidden_index);
      auto &node_incidents2 = this->nodes_to_incidents_.at(edge.node2().index_);
      node_incidents2.erase(hidden_index);
      if (edge_index != this->edge_index_to_hidden_index_.size())
      {
        size_type last_hidden_index = this->edge_index_to_hidden_index_.at(this->edge_index_to_hidden_index_.size());
        size_type last_edge_index = this->edge_index_to_hidden_index_.size();
        this->edge_index_to_hidden_index_.erase(last_edge_index);
        this->hidden_index_to_edge_index_.erase(last_hidden_index);

        this->edge_index_to_hidden_index_.insert({edge_index, last_hidden_index});
        this->hidden_index_to_edge_index_.insert({last_hidden_index, edge_index});
      }

      this->num_edges_--;
      return 1;
    }
    catch (...)
    {
      return 0;
    }
  }

  /** @brief Removes an edge from the graph. O(num nodes() + num edges())
   * 
   * @param[in,out] @a e_it, an iterator to the edge to remove.
   * 
   * @return A valid iterator over the edges of the graph.
   * 
   * @post invalidates iterators over the edges of the graph, 
   * over the incident edges of a node.
   * Also, invalidates the edge that was removed.
   * 
   * @post graph.num_nodes() has same number of nodes than before.
   * @post graph.num_edges() has -1 edge.
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    remove_edge(*e_it);
    return this->edge_begin();
  }

private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::unordered_map<size_type, node_type> nodes_;
  std::unordered_map<size_type, Point> points_;
  std::unordered_map<size_type, node_value_type> node_values_;
  std::unordered_map<size_type, size_type> hidden_index_to_node_index_;
  std::unordered_map<size_type, size_type> node_index_to_hidden_index_;
  std::unordered_map<size_type, edge_type> edges_;
  std::unordered_map<size_type, edge_value_type> edge_values_;
  std::unordered_map<size_type, size_type> hidden_index_to_edge_index_;
  std::unordered_map<size_type, size_type> edge_index_to_hidden_index_;
  std::unordered_map<std::pair<size_type, size_type>, size_type, pair_hash> nodes_to_edges_;
  std::unordered_map<size_type, std::unordered_set<size_type>> nodes_to_incidents_;
  size_type size_;
  size_type num_edges_;
  size_type next_node_index_;
  size_type next_edge_index_;
};

#endif // CME212_GRAPH_HPP

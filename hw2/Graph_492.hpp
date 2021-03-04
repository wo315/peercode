#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph
{
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

private:
  /** Custom structure of data to store with Nodes */
  struct EdgeInfo
  {
    unsigned int n1id;
    unsigned int n2id;
    E value;
  };

  /** Custom structure of data to store with Edges */
  struct NodeInfo
  {
    Point p_;
    node_value_type value_;
    size_type idx_;
  };

  std::vector<NodeInfo> nodes_;                                                          //vector of points/nodes, indexed by node uid
  std::vector<size_type> i2u_;                                                           //indexed by node idx
  std::unordered_map<unsigned int, EdgeInfo> edges_;                                     //maps of edges with indexes as keys
  std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>> n2e_; //maps of nodes to connected nodes and edge index

public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
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

    //Construct an invalid node
    Node() : graph_(nullptr), nid_(0)
    {
    }

    /** Return this node's position. */
    Point &position()
    {
      return this->graph_->nodes_[nid_].p_;
    }

    const Point &position() const
    {
      return this->graph_->nodes_[nid_].p_;
    }

    /** Return this node's user facing index, a number in the range [0, graph_size). */
    size_type index() const
    {
      assert(this->graph_ != nullptr);
      size_type idx = this->graph_->nodes_[nid_].idx_;
      assert(this->graph_->i2u_[idx] == nid_);
      return idx;
    }

    /** Retuen the value of this node */
    node_value_type &value()
    {
      return this->graph_->nodes_[nid_].value_;
    };
    const node_value_type &value() const
    {
      return this->graph_->nodes_[nid_].value_;
    };

    /* Return the number incidents/edges of a node*/
    size_type degree() const
    {
      return this->graph_->n2e_[this->nid_].size();
    };

    /* Begin and end for inciedent iterator */
    incident_iterator edge_begin() const
    {
      std::unordered_map<size_type, size_type>::const_iterator cit = this->graph_->n2e_[this->nid_].begin();
      return IncidentIterator(this->graph_, *this, cit);
    };

    incident_iterator edge_end() const
    {
      std::unordered_map<size_type, size_type>::const_iterator cit = this->graph_->n2e_[this->nid_].end();
      return IncidentIterator(this->graph_, *this, cit);
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      assert(n.graph_ != nullptr and this->graph_ != nullptr);
      return (this->graph_ == n.graph_ and this->nid_ == n.nid_);
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
      assert(n.graph_ != nullptr and this->graph_ != nullptr);
      if (graph_ < n.graph_)
      {
        return true;
      }
      else if (graph_ == n.graph_ and nid_ < n.nid_)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type *graph_; //pointer to graph the node belongs to
    size_type nid_;     //node index in the nodes_ vector
    Node(const graph_type *graph, size_type nid)
        : graph_(const_cast<graph_type *>(graph)), nid_(nid)
    {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return i2u_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return this->size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type &value = node_value_type())
  {
    size_type idx = this->size();
    i2u_.push_back(nodes_.size());
    nodes_.push_back(NodeInfo{position, value, idx});
    return Node(this, i2u_[idx]);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    if (this == n.graph_ and n.index() < size())
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    if (i < i2u_.size())
    {
      return Node(this, i2u_[i]);
    }
    else
    {
      return Node(); // Invalid node
    }
  }

  /** Remove a node from the graph, returning 1 if removed
    * and 0 node doesn't exist.
   * @param[in] n node to be removed
   * @post new num_nodes() == old num_nodes() - 1
   * @post g.node(i).index == i for all i with 0<=i<g.num_nodes()
   * @post g.node(n.idex())==n
   *
   * Complexity: O(num_edges())amortized operations.
   */
  size_type remove_node(const Node &n)
  {
    if (has_node(n))
    {
      unsigned int nid = n.nid_;
      size_type last_idb=0;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
      {
        Edge e = (*it);
        unsigned int idb = e.n2id_;
        if (it != n.edge_begin())
        {
          remove_edge_with_ids(nid, last_idb);
        }
        last_idb = idb;
      }

      if (n.degree() != 0)
      {
        remove_edge_with_ids(nid, last_idb);
      }

      n2e_.erase(nid);

      unsigned int idx = n.index();
      i2u_[idx] = std::move(i2u_.back());
      this->nodes_[i2u_[idx]].idx_ = idx;
      i2u_.pop_back();
      return 1;
    }
    else
    {
      return 0;
    }
  }

   /** Remove a node from the graph, returning 1 if removed
    * and 0 node doesn't exist.
   * @param[in] n_it node iterator pointing to the memory address of the node
   * to be removed
   * @post new num_nodes() == old num_nodes() - 1
   * @post g.node(i).index == i for all i with 0<=i<g.num_nodes()
   * @post g.node(n.idex())==n for all nodes
   *
   * Complexity: O(num_edges())amortized operations.
   */
  node_iterator remove_node(node_iterator n_it)
  {
    Node n = *n_it;
    n_it++;
    remove_node(n);
    return n_it;
  }

    /** Remove an edge from the graph, returning 1 after removal.
   * @param[in] ida index of the one node connecting the edge
   * @param[in] idb index of the other node connecting the edge
   * @post new num_edges() == old num_edges() - 1
   *
   * Complexity: O(1) amortized operations.
   */
  size_type remove_edge_with_ids(size_type ida, size_type idb)
  {
    unsigned int eidx = n2e_.at(ida).at(idb);
    edges_.erase(eidx);
    n2e_.at(ida).erase(idb);
    n2e_.at(idb).erase(ida);
    return 1;
  }

  /** Remove an edge from the graph, returning 1 after removal.
   * @param[in] a one node connecting the edge
   * @param[in] b the other node connecting the edge
   * @post new num_edges() == old num_edges() - 1
   *
   * Complexity: O(1) amortized operations.
   */
  size_type remove_edge(const Node &a, const Node &b)
  {
    unsigned int ida = a.nid_;
    unsigned int idb = b.nid_;
    if (has_edge(a, b))
    {
      return remove_edge_with_ids(ida, idb);
    }
    else
    {
      return 0;
    }
  }

  /** Remove an edge from the graph, returning iterator
   * pointing to the next element after removal.
   * @param[in] e_it edge_iterator pointing to the memory location 
   * of the edge to be removed
   * @post new num_edges() == old num_edges() - 1
   *
   * Complexity: O(1) amortized operations.
   */
  edge_iterator remove_edge(edge_iterator e_it)
  {
    Edge e = *e_it;
    e_it++;
    unsigned int ida = e.n1id_;
    unsigned int idb = e.n2id_;
    remove_edge_with_ids(ida, idb);
    return e_it;
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
    Edge() : graph_(nullptr), n1id_(0), n2id_(0), eid_(0)
    {
    }

    /** Return the value of Edge. */
    edge_value_type &value()
    {
      return this->graph_->edges_[eid_].value;
    }

    const edge_value_type &value() const
    {
      return this->graph_->edges_[eid_].value;
    }

    /** Computes the length of an edge.*/
    double length() const
    {
      Point p1 = this->graph_->nodes_[n1id_].p_;
      Point p2 = this->graph_->nodes_[n2id_].p_;
      double len = norm(p1 - p2);
      return len;
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return Node(this->graph_, n1id_);
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return Node(this->graph_, n2id_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      assert(e.graph_ != nullptr and this->graph_ != nullptr);
      if (graph_ == e.graph_ and n1id_ == e.n1id_ and n2id_ == e.n2id_)
      {
        return true;
      }
      else if (graph_ == e.graph_ and n1id_ == e.n2id_ and n2id_ == e.n1id_)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      if (graph_ < e.graph_)
      {
        return true;
      }
      else if (graph_ == e.graph_ and eid_ < e.eid_)
      {
        return true;
      }
      else
      {
        return false;
      }
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type *graph_; //graph where it belongs to
    size_type n1id_;    //first  node
    size_type n2id_;    //second node
    size_type eid_;     //ID of the edge
    Edge(const graph_type *graph, size_type n1id, size_type n2id, size_type eid)
        : graph_(const_cast<graph_type *>(graph)), n1id_(n1id), n2id_(n2id), eid_(eid)
    {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    if (i < num_edges())
    {
      return Edge(this, edges_.at(i).n1id, edges_.at(i).n2id, i);
    }
    else
    {
      return Edge();
    }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    size_type aid = a.nid_;
    size_type bid = b.nid_;
    auto it = n2e_.find(aid);
    if (it != n2e_.end())
    {
      return (n2e_.at(aid).count(bid) != 0);
    }
    else
    {
      return false;
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
  Edge add_edge(const Node &a, const Node &b)
  {
    unsigned int ida = a.nid_;
    unsigned int idb = b.nid_;
    if (has_edge(a, b))
    {
      unsigned int eidx = n2e_.at(ida).at(idb); //find the index of the edge in the nested map
      return Edge(this, ida, idb, eidx);        //return the current edge if already exists
    }
    else
    {
      unsigned int eidx = edges_.size();
      edges_.insert({eidx, EdgeInfo{ida, idb, edge_value_type()}}); //add new edge to the edges map
      if (n2e_.count(ida) == 0)
      {
        std::unordered_map<unsigned int, unsigned int> nidmap1 = {{idb, eidx}};
        n2e_.insert({ida, nidmap1});
      }
      if (n2e_.count(idb) == 0)
      {
        std::unordered_map<unsigned int, unsigned int> nidmap2 = {{ida, eidx}};
        n2e_.insert({idb, nidmap2});
      }
      n2e_.at(ida).insert({idb, eidx});
      n2e_.at(idb).insert({ida, eidx});
      return Edge(this, ida, idb, eidx);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    nodes_.clear();
    edges_.clear();
    n2e_.clear();
    i2u_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : public std::iterator<std::forward_iterator_tag, Node>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
      graph_ = nullptr;
      currentIdx_ = 0;
    }

    /** Advancing the iterator 1 step forward and returning the updated iterator */
    NodeIterator &operator++()
    {
      currentIdx_++;
      return (*this);
    }

    /** Returning whether two NodeIterator is equivalent. */
    bool operator==(const NodeIterator &node_iter) const
    {
      return (graph_ == node_iter.graph_ and currentIdx_ == node_iter.currentIdx_);
    }

    /** Returning whether two NodeIterator is different. */
    bool operator!=(const NodeIterator &node_iter) const
    {
      return !(node_iter == (*this));
    }

    /** Return the node where the interator is pointing to. */
    Node operator*() const
    {
      return Node(graph_, graph_->i2u_[currentIdx_]);
    }

  private:
    friend class Graph;
    const graph_type *graph_;
    size_type currentIdx_;

    /** Private constructor for NodeIterator */
    NodeIterator(const graph_type *graph, size_type currentIdx)
    {
      this->graph_ = graph;
      this->currentIdx_ = currentIdx;
    }
  };

  /* Construct a NodeIterator pointing to the beginning of i2u_ vector*/
  NodeIterator node_begin() const
  {
    return NodeIterator(this, 0);
  }

  /* Construct a NodeIterator pointing to the end of i2u_ vector*/
  NodeIterator node_end() const
  {
    return NodeIterator(this, this->size());
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
    IncidentIterator()
    {
      graph_ = nullptr;
      node_ = NULL;
      cit_ = NULL;
    }

    /** Return the edge where the incident interator is pointing to. */
    Edge operator*() const
    {
      return Edge(this->graph_, node_.nid_, (*cit_).first, (*cit_).second);
    }

    /** Advancing the unordered map iterator 1 step forward and returning the updated iterator */
    IncidentIterator &operator++()
    {
      cit_++;
      return (*this);
    }

    /** Returning whether two IncidentIterators are equivalent. */
    bool operator==(const IncidentIterator &inc_iter) const
    {
      return (this->graph_ == inc_iter.graph_ and this->cit_ == inc_iter.cit_);
    }
    /** Returning whether two IncidentIterators are different. */
    bool operator!=(const IncidentIterator &inc_iter) const
    {
      return !(this->cit_ == inc_iter.cit_);
    }

  private:
    friend class Graph;
    const graph_type *graph_;
    node_type node_;
    std::unordered_map<size_type, size_type>::const_iterator cit_;

    /* Private constructor for IncidentIterator. */
    IncidentIterator(const graph_type *graph, node_type node, std::unordered_map<size_type, size_type>::const_iterator cit)
    {
      this->graph_ = graph;
      this->node_ = node;
      this->cit_ = cit;
    }
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
    EdgeIterator()
    {
      graph_ = nullptr;
      cit_ = NULL;
    }

    /** Return the edge where the interator is pointing to. */
    Edge operator*() const
    {
      return Edge(this->graph_, (*cit_).second.n1id, (*cit_).second.n2id, (*cit_).first);
    }

    /**Advancing the unordered map iterator 1 step forward and returning the updated iterator */
    EdgeIterator &operator++()
    {
      cit_++;
      return (*this);
    }

    /** Returning whether two IncidentIterators are equivalent. */
    bool operator==(const EdgeIterator &edge_iter) const
    {
      return (edge_iter.graph_ == this->graph_ and edge_iter.cit_ == this->cit_);
    }
    /** Returning whether two IncidentIterators are different. */
    bool operator!=(const EdgeIterator &edge_iter) const
    {
      return !(edge_iter == (*this));
    }

  private:
    friend class Graph;
    const graph_type *graph_;
    typename std::unordered_map<size_type, EdgeInfo>::const_iterator cit_;

    EdgeIterator(const graph_type *graph, typename std::unordered_map<size_type, EdgeInfo>::const_iterator cit)
    {
      this->graph_ = graph;
      this->cit_ = cit;
    }
  };

  /* Construct a EdgeIterator pointing to the beginning of edge_ vector*/
  edge_iterator edge_begin() const
  {
    typename std::unordered_map<size_type, EdgeInfo>::const_iterator cit = this->edges_.begin();
    return EdgeIterator(this, cit);
  }

  /* Construct a EdgeIterator pointing to the end of edge_ vector*/
  edge_iterator edge_end() const
  {
    typename std::unordered_map<size_type, EdgeInfo>::const_iterator cit = this->edges_.end();
    return EdgeIterator(this, cit);
  }
};

#endif // CME212_GRAPH_HPP
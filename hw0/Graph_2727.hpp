#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iostream>
#include <list>
#include <unordered_map>
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
    using private_size_type = unsigned;
    // Unordered maps have a short amortized search complexity of O(1).
    // Worst case, O(n)
    std::unordered_map<private_size_type, Point> graph_points_;
    private_size_type num_Nodes_;
    
    // Map of edge indexes as keys, and values being a list of the indexes
    // of the nodes that this edge connects.
    std::map<private_size_type, std::list<private_size_type> > graph_edges_;
    
    // A map of node indexes as keys, with values being a map with  node indexes
    // as the keys, with the final value being the edge index.
    // Can search by node index, then by the connecting node index, and if it
    // exists, will yield the edge index
    std::map<private_size_type, std::map<private_size_type,
            private_size_type>> edge_map_;
    private_size_type num_Edges_;
    
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
    : graph_points_(), num_Nodes_(0), num_Edges_(0){
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
    }

    /** Return this node's position. */
    const Point& position() const {
        return this->fetch();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        return this->my_index_;
    }
      /** Return this node's graph pointer*/
      Graph* graph() const {
          return this->my_graph_;
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
      return (n.graph() == this->my_graph_ && n.index() == this->my_index_);
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
        return this->index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
      Graph* my_graph_;
      size_type my_index_;
      
    /** Private constructor, allows us to pass over the graph and the index value
     */
      Node(const Graph* g, size_type i)
      : my_graph_(const_cast<Graph*>(g)), my_index_(i){
      }
      
      Point& fetch() const {
          return my_graph_->graph_points_[my_index_];
      }
      
  }; /// End Node class
///   Now in public part of graph class

    
  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
      return this->num_Nodes_;
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
      size_type old_num_Nodes = num_Nodes_;
      this->num_Nodes_ += 1;
      // assert that this index doesn't already exist in the map
      assert(this->graph_points_.find(old_num_Nodes) ==
             this->graph_points_.end());
      this->graph_points_.insert( std::pair<size_type,
                                 Point>(old_num_Nodes, position));
      return Node(this, old_num_Nodes);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
      if (n.my_graph_ == this){
          return true;
      } else {
          return false;
  }
  }
     
  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
      // Returning a proxy
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
    }

    /** Return a node of this Edge */
    Node node1() const {
        // Return a proxy node
        return Graph::Node(my_edge_graph_,
                        my_edge_graph_->graph_edges_[my_edge_index_].front());
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // Return a proxy node
        return Graph::Node(my_edge_graph_,
                        my_edge_graph_->graph_edges_[my_edge_index_].back());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
        return norm(node1().position()) == norm(e.node1().position());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
        // Sort by norm of node1, then tiebreak via norm of node2
        if (norm(node1().position()) == norm(e.node1().position())) {
            return norm(node2().position()) < norm(e.node2().position());
        } else {
        return norm(node1().position()) < norm(e.node1().position());
    }
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
      Graph* my_edge_graph_;
      size_type my_edge_index_;
      size_type node1_index_;
      size_type node2_index_;
      
      /** Private constructor
       * An edge is defined by a unique index, a graph, and two nods.
       * The nodes are identified by their indexes.
       */
      Edge(const Graph* g, size_type i, const Node& a, const Node& b)
      : my_edge_graph_(const_cast<Graph*>(g)), my_edge_index_(i),
      node1_index_(a.index()), node2_index_(b.index()){
      }
  }; // End Edge class. Now in public part of graph class

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
      return this->num_Edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      // Return a proxy Edge
      return Edge(this, i, node(graph_edges_.at(i).front()),
                  node(graph_edges_.at(i).back()));
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      size_type a_index = a.index();
      size_type b_index = b.index();
      if (edge_map_.find(a_index) != edge_map_.end()){
          // First node is in the primary map.
          // See if it is linked to the second node.
          if(edge_map_.at(a_index).find(b_index) !=
             edge_map_.at(a_index).end()){
              // Edge exists
              return true;
          }
          assert(true);
      } else if (edge_map_.find(b_index) != edge_map_.end()){
          // Second node is in the primary map. See if it is
          // linked to the first node.
          if(edge_map_.at(b_index).find(a_index) !=
                    edge_map_.at(b_index).end()){
              // Edge exists
              return true;
          }
          // Node b is in map, but is not linked to node a.
          // Nor is a linked to b.
          return false;
  }
      // Neither first node nor second node are in edge map
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
      // Check to see if the edge exists
      size_type a_index = a.index();
      size_type b_index = b.index();
      if (edge_map_.find(a_index) != edge_map_.end()){
          // First node is in the primary map.
          // See if it is linked to the second node.
          if(edge_map_.at(a_index).find(b_index) !=
             edge_map_.at(a_index).end()){
              // Edge exists
              return Edge(this, edge_map_[a_index].at(b_index),a,b);
          }
          assert(true);
      } else if (edge_map_.find(b_index) != edge_map_.end()){
          // Second node is in the primary map.
          // See if it is linked to the first node.
          if(edge_map_.at(b_index).find(a_index) !=
             edge_map_.at(b_index).end()){
              // Edge exists
              return Edge(this, edge_map_[b_index].at(a_index),a,b);
          }
          // Node b is in map, but is not linked to node a.
          // Nor is a linked to b.
          assert(true);
  }
      // Neither first node nor second node are in edge map
      assert(true);
      // Create a new edge
      size_type new_edge_index = num_Edges_;
      this->graph_edges_[new_edge_index].push_back(a.index());
      this->graph_edges_[new_edge_index].push_back(b.index());
      this->num_Edges_++;
      // Add to our edge_map_ an edge between these nodes
      this->edge_map_[a.index()][b.index()] = new_edge_index;
      
      return Edge(this, new_edge_index,a,b);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      num_Nodes_ = (size_type)(0);
      num_Edges_ = (size_type)(0);
      graph_points_.clear();
      graph_edges_.clear();
      edge_map_.clear();
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

};

#endif // CME212_GRAPH_HPP

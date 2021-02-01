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

    // Initialize p to position and initialize Graph *g
    Node(Point p, Graph *g) {
      // HW0: YOUR CODE HERE
        this->p = p;
        this->graph = g;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */

    // Identify the node in the Graph by position
    size_type index() const {
      // HW0: YOUR CODE HERE
        int index = -1;
        for(unsigned int i = 0; i < graph->size(); i++)
          if(graph->nodes[i].p == this->p)
	    {
              index = i;
              break;
	    }
      return index;
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
    // HW0: YOUR CODE HERE
    //  (void) n;          // Quiet compiler warning
      return this->index() == n.index() && this->graph == n.graph;
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
    //  (void) n;           // Quiet compiler warning      
        bool result = false;
        if(this->p[0] < n.p[0])
            result = true;
        else if(this->p[0] == n.p[0] && this->p[1] < n.p[1])
            result = true;
        else if(this->p[0] == n.p[0] && this->p[1] == n.p[1] && this->p[2] < n.p[2])
            result = true;
        return result;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Graph of which node is a part
    Point p;
    Graph *graph; 
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
      return nodes.size();
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
  
  Node add_node(const Point& position)
  {
  // HW0: YOUR CODE HERE
  //    (void) position;      // Quiet compiler warning
    
      Node n(position, this);
      if(!has_node(n))
      nodes.push_back(n);
      else 
        for(unsigned int i = 0; i < nodes.size(); i++)
            if(n==nodes[i])
                return nodes[i];
    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
  // HW0: YOUR CODE HERE
  //  (void) n;            // Quiet compiler warning

    /*if(std::find(nodes.begin(), nodes.end(), n))
      return true;
    else
      return false;
    */
    for(unsigned int i = 0; i < nodes.size(); i++)
      if(n==nodes.at(i))
	return true;
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
   // (void) i;             // Quiet compiler warning
    return nodes[i];        // Invalid node
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
    
    Edge(const Node& node1, const Node& node2) {
      // HW0: YOUR CODE HERE
        this->n1 = &node1;
        this->n2 = &node2;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return *n1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return *n2;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
    // (void) e;           // Quiet compiler warning
    //HW0: YOUR CODE HERE
        return (n1 == e.n1 && n2 == e.n2) || (n1 == e.n2 && n2 == e.n1);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const
    {
    //  (void) e;           // Quiet compiler warning
    //HW0: YOUR CODE HERE
      if(this->n1 < e.n1 || (this->n1 == e.n1 && this->n2 < e.n2))
	return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    const Node *n1;
    const Node *n2;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
  // HW0: YOUR CODE HERE
  //  (void) i;             // Quiet compiler warning
    return edges[i];        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
   // HW0: YOUR CODE HERE
   // (void) a; (void) b;   // Quiet compiler warning
    for(unsigned int i = 0; i < edges.size(); i++)
      if((edges[i].n1->p==a.p && edges[i].n2->p == b.p) || (edges[i].n1->p==b.p && edges[i].n2->p == a.p))
	return true;
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
  Edge add_edge(const Node& a, const Node& b)
  {
  // HW0: YOUR CODE HERE
  //  (void) a, (void) b;   // Quiet compiler warning
    Edge e(a, b);
    if(!has_edge(a, b))
      edges.push_back(e);
    else
      {
      for(unsigned int i = 0; i < edges.size(); i++)
	if((edges[i].n1->p==a.p && edges[i].n2->p == b.p) || (edges[i].n1->p==b.p && edges[i].n2->p == a.p))
	  return edges[i];
      }
      return e;        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
  // HW0: YOUR CODE HERE
    edges.clear();
    nodes.clear();
  }

  //
  // Node Iterator
  //

  // The Viewer would not compile without node_begin(), node_end(),
  // edge_begin(), and edge_end() so they are defined below
  // even though they are listed as HW1
  
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
    NodeIterator(const Node *node) {
        this->node = node;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
      Node operator*() const {
          return *node;
      }
      NodeIterator& operator++() {
          this->node++;
          return *this;
      }
      bool operator==(const NodeIterator&it) const {
          return this->node==it.node;
      }
      bool operator!=(const NodeIterator&it) const {
          return this->node!=it.node;
      }

   private:
    friend class Graph;
    const Node *node;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    node_iterator node_begin() const {
        NodeIterator it(&(*nodes.begin()));
        return it;
    }
    node_iterator node_end() const {
        NodeIterator it(&(*nodes.end()));
        return it;
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
    EdgeIterator(const Edge *edge) {
        this->edge = edge;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
      Edge operator*() const {
          return *edge;
      }
      EdgeIterator& operator++() {
          this->edge++;
          return *this;
      }
      bool operator==(const EdgeIterator&it) const {
          return this->edge==it.edge;
      }
      bool operator!=(const EdgeIterator&it) const {
          return this->edge!=it.edge;
      }

   private:
    friend class Graph;
      const Edge *edge;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    edge_iterator edge_begin() const {
        EdgeIterator it(&(*edges.begin()));
        return it;
    }
    
    edge_iterator edge_end() const {
        EdgeIterator it(&(*edges.end()));
        return it;
    }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    
      std::vector<Node> nodes;
      std::vector<Edge> edges;

};

#endif // CME212_GRAPH_HPP

#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <iostream>
#include <vector>
#include <map>
#include "Point.hpp"


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
    this -> nodes = std::vector< Node>(0);
    this -> edges = std::vector< Edge>(0);
    size_type z = 0;
    this -> node_size = z;
    this -> edge_size = z;
    this -> node_map = std::map< size_type,  Point>();
    this -> edge_map = std::map<std::tuple<const Node, const Node>,  size_type>();
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
    Node() { //nothing here! define a real one later in private
    }

    /** Return this node's position. */
    const Point& position() const {
      return *pos; //pos is a pointer so we need to dereference it
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      //std::cout<< "in index() function *idx is " << *idx << std::endl;
      return *idx; //idx is a pointer so we need to dereference it
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
      (void) n;          // Quiet compiler warning
      if (*idx == n.index()){ return true;
      } else {return false;}

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
      (void) n;           // Quiet compiler warning
      if(*idx < n.index()){return true;
      }  else {return false;}
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    size_type* idx; // a pointer to the node's index
    const Point* pos; // a pointer to the node's position/Point.


    Node(size_type& index, const Point& position){ //removed '&' after 'point'
      this -> idx = &index; // reference to the address of index
      this -> pos = &position; // reference to the address of the Point object
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_size;
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
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    size_type ind = node_size;
    //std::cout<< "ind is " << ind << std::endl;
    Node new_node = Node(ind, position); // make the new node
    //std::cout<< "new_node.index() is " << new_node.index() << " should be " << ind << std::endl;
    nodes.push_back(new_node); // store the node in the graph's node vector
    node_size = node_size + 1; // update the node size
    //std::cout<< " ind is " << ind << std::endl;
    //int int_ind = static_cast<int>(ind);
    //std::cout<< "nodes[ind].index() is "<< (nodes[int_ind]).index() << " and is supposed to be " << ind << std::endl;
    node_map[ind] = position; // this is where the real information is stored
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    if(node_map.count(n.index()) == 0) {
      return false;
    } else {return true;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    (void) i;             // Quiet compiler warning
    int int_i = static_cast<int>(i);
    //std::cout<<"nodes[int_i].index() is "<<(nodes[int_i]).index()<< " and should be " << int_i << std::endl;
    return nodes[int_i]; //The i index of the nodes vector (zero indexed)
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
      return *left;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return *right;      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (
          ((*left == e.node1()) or (*left == e.node2())) and
          ((*right == e.node1()) or (*right == e.node2()))
          ){ return true;
      } else {return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (*left < e.node1()){ return true; //use the node's definition of "<"
      } else {return false;}
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const Node* left; // pointer to one end of edge
    const Node* right; //pointer to second end of edge
    size_type* index;
    Edge(size_type idx, const Node a, const Node b) {
      this -> left = &a;
      this -> right = &b;
      this -> index = &idx;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    (void) i;             // Quiet compiler warning
    return edges[i];        // edge at index i of the edges vector
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node a, const Node b) const {
    (void) a; (void) b;   // Quiet compiler warning
    std::tuple<const Node, const Node> temp1 = {a, b};
    std::tuple<const Node, const Node> temp2 = {b, a};
    if(edge_map.count(temp1) + edge_map.count(temp2)== 0) {
      return false;
    } else {
      return true; }
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
  Edge add_edge(const Node a, const Node b) {
    (void) a, (void) b;   // Quiet compiler warning
    // HW0: YOUR CODE HERE
    // add the (node, node) tuple as the *key* to the map. the number is the value.
    std::tuple<const Node, const Node> pair1 = {a, b};
    std::tuple<const Node, const Node> pair2 = {b, a};
    if(has_edge(a, b)){
      return edges[ edge_map[ pair1 ] ]; // edge_map[pair1] will return the index of the edge defined by (a,b)
    } else if(has_edge(b, a)){
      return edges[ edge_map[ pair2 ] ]; // edge_map[pair1] will return the index of the edge defined by (b, a)
    } else {
      size_type ind = edge_size;
      edge_map[pair1]= ind; // add edge to edge map;
      Edge new_Edge = Edge(ind, a, b); // create the edge
      edges.push_back(new_Edge); // add edge to vec edges;
      edge_size = edge_size +1; // update the size of edges
      return new_Edge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    /**
    for(size_type i=0; i < edges.size(); i++){ //delete all the actual nodes
      delete(edges[i]);
    }
    for(size_type i=0; i < nodes.size(); i++) { //delete all the actual edges
      delete (nodes[i]);
    }
    //delete(node_map); //delete the map containing info about the nodes
    //delete(edge_map); //delete the map containing info about the edges
    //delete(nodes); //delete the vector containing the nodes
    //delete(edges); //delete the vector containing the edges
    //delete(node_size);
    //delete(edge_size);*/
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
  std::vector< Node> nodes; // a vector of references to all the nodes in the graph
  std::vector< Edge> edges; // a vector of references to all the edges in the graph
  size_type node_size; // an size_type defining the size of the graph. not const because it changes
  size_type edge_size; // an size_type defining the number of edges in the graph. not const because it changes
  std::map<size_type, Point> node_map; // a map containing all the real information of the nodes
  std::map<std::tuple<const Node, const Node>, size_type> edge_map; // a map containing all the real information of the edges

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP


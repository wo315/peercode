#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <map>
#include <tuple>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// class Graph; // forward declaration

template <typename V>
class Graph {
  /** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

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
  Graph(){
    // HW0: YOUR CODE HERE
    next_nodeID = 0;
    next_edgeID = 0;
    //std::map<size_type, unsigned> node_idx;
    //std::vector<Node*> graph_nodes;
    //std::map<size_type, unsigned> edge_idx; // indices for the edges
    //std::vector<Edge*> graph_edges;
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
    Node(): parent(nullptr), nodeID{} {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return parent->graph_nodes[index()].point; 
    }

    size_type getID() const {
      // HW0: YOUR CODE HERE
      return nodeID; 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      if (parent == nullptr) return size_type(-1);
      auto graphidx = (*parent).node_idx.find(nodeID);
      if (graphidx != (*parent).node_idx.end()) return graphidx->second;
      return size_type(-1);
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return parent->graph_nodes[index()].val;
    };

    const node_value_type& value() const {
      return parent->graph_nodes[index()].val;
    };

    bool& visited() {
      return parent->graph_nodes[index()].visited;
    };

    size_type degree() const{
      auto froma = parent->node2node.find(nodeID);
      if (froma == parent->node2node.end()){
        return 0;
      }else{
        std::vector<size_type> tonodes = *(parent->node2node[nodeID]); // returns a dereferenced pointer to vector of size_types
        return tonodes.size();
      }
    };

    incident_iterator edge_begin() const{
      return incident_iterator(parent,&this);
    };
    incident_iterator edge_end() const{
      return incident_iterator();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      //nodesinParent = *parent
      //*(n.parent).node_idx[idx]
      if (parent == n.parent && index() == n.index()){
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
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;           // Quiet compiler warning
      if (index() < n.index()) return true;
      return false;
    }

    Graph* parent;

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type nodeID; //some index

    Node(const Graph* init_parent, size_type init_idx) {
      // HW0: YOUR CODE HERE
      parent    = const_cast<Graph*>(init_parent);
      nodeID    = init_idx;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return graph_nodes.size();
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()){
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    // Point* cppoint = new Point(position);
    size_type nnodes = num_nodes();
    internal_node newnode_int(position,nnodes,next_nodeID,val); // does this work or do i have to allocate "new"
    Node newnode(this, next_nodeID); // memory allocation?
    node_idx[next_nodeID] = nnodes;
    ++ next_nodeID; // update next index
    graph_nodes.push_back(newnode_int);
    return newnode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    Node tempnode(this, graph_nodes[n.index()].node_id);
    if (n == tempnode) return true;
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
    (void) i;             // Quiet compiler warning
    size_type nnodes = num_nodes();
    if (i < nnodes) {
      return Node(this, graph_nodes[i].node_id); //create a new node
    }
    return Node();
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
    Edge(): parent(nullptr), node1ID{}, node2ID{}, edgeID{} {}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return {parent,node1ID};
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return {parent,node2ID};
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if (parent != e.parent) return false; // no comparison for graphs yet..
      if ((node1() == e.node1()) && (node2() == e.node2())) return true;
      else if ((node1() == e.node2()) && (node2() == e.node1())) return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if (edgeID < e.edgeID) return true;
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* parent;
    size_type node1ID;
    size_type node2ID;
    size_type edgeID;

    Edge(const Graph* init_parent, size_type node1ID, size_type node2ID, size_type init_idx)
    : parent(const_cast<Graph*>(init_parent)), node1ID(node1ID), node2ID(node2ID), edgeID(init_idx) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_idx.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    if (i < num_edges()) {
      internal_edge ie = graph_edges[i];
      return Edge(this, ie.node1.getID(), ie.node2.getID(), ie.edge_id); // return an Edge object with same ID
    }
    return Edge(); // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    
    auto froma = nodes2edge.find(std::make_pair(a.getID(), b.getID()));
    if (froma != nodes2edge.end()) return true;

    auto fromb = nodes2edge.find(std::make_pair(b.getID(), a.getID()));
    if (fromb != nodes2edge.end()) return true;
    
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
    // HW0: YOUR CODE HERE
    (void) a, (void) b;   // Quiet compiler warning
    if (!has_edge(a,b)){
      Edge newedge(this, a.getID(), b.getID(), next_edgeID); // use proxy
      size_type ne = num_edges();
      edge_idx[next_edgeID] = ne; // add at the end. 
      
      nodes2edge[std::make_pair(a.getID(), b.getID())] = ne; // store index f

      // update indices
      std::vector<size_type>* fromavec;
      auto froma = node2node.find(a.getID());
      if (froma == node2node.end()){
        fromavec = new std::vector<size_type>;
        node2node[a.getID()] = fromavec; //assign a pointer to a new std empty vector
      } else {
        fromavec = froma->second;
      }; 
      (*fromavec).push_back(b.getID()); //map nodes to last edge on graph_edges

      auto fromb = node2node.find(b.getID());
      std::vector<size_type>* frombvec;
      if (fromb == node2node.end()){
        frombvec = new std::vector<size_type>;
        node2node[b.getID()] = frombvec; //assign a pointer to a new std empty vector
      } else {
        frombvec = fromb->second;
      };
      (*frombvec).push_back(a.getID()); //map nodes to last edge on graph_edges

      internal_edge newedge_int(a,b,ne,next_edgeID); // create data for the edge
      graph_edges.push_back(newedge_int);
      ++next_edgeID;
      return newedge;
    }
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    graph_nodes.clear();
    graph_edges.clear();

    graph_nodes = {};
    graph_edges = {};
    node_idx = {};
    next_nodeID = 0;
    edge_idx = {};
    next_edgeID = 0;
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
    NodeIterator(): nd(Node()), parent(nullptr) {};

    // constructor
    NodeIterator(const Graph* init_parent, value_type n): nd(n), parent(const_cast<Graph*>(init_parent)) {};

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      /** Derefernce operator
       * @return position of the Node. 
      **/
      return nd;
    };

    NodeIterator& operator++() {
      /** ++ operator
       * @return next node in index.. 
      **/
      size_type idx = (nd.index()) + 1;
      size_type nnodes = (*parent).num_nodes();
      if (idx < nnodes) {
        Node nextnode(parent, parent->graph_nodes[idx].node_id);
        nd = nextnode;
      }else{
        nd = Node();
        parent = nullptr;
      };
      return *this;
    };

    bool operator==(const NodeIterator& someNI) const {
      return nd == someNI.nd;
    };

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    value_type nd;
    Graph* parent;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    /** finds the first node and @returns an iterator with the given node
     **/
    Node headnode(this, graph_nodes[0].node_id);
    return NodeIterator(this, headnode);
  }

  node_iterator node_end() const{
    return NodeIterator();
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(): parent(nullptr), incidentNode(Node()), idx{}, ed(Edge()) {};

    IncidentIterator(const Graph* init_parent, Node& initNode){
      parent        = const_cast<Graph*>(init_parent);
      incidentNode  = initNode; // const_cast<Node*>(initNode);
      idx = 0;
      if (incidentNode.degree() > 0){
        ed = findEdge(idx);
      }else{
        ed = Edge();
      };
    };

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      return ed;
    };

    Edge findEdge(size_type idx){
      // find the next node in line. 
      size_type a = incidentNode.getID();
      std::vector<size_type> tonodes = *(parent->node2node[a]); // returns a dereferenced pointer to vector of size_types
      size_type b = tonodes[idx];
      
      // find the edge_id and create an edge object.
      auto froma = parent->nodes2edge.find(std::make_pair(a, b));
      if (froma != parent->nodes2edge.end()){
        return Edge(parent, a,b, parent->graph_edges[froma->second].edge_id);
      };
      auto fromb = parent->nodes2edge.find(std::make_pair(b, a));
      if (fromb != parent->nodes2edge.end()) {
        return Edge(parent, a,b, parent->graph_edges[fromb->second].edge_id);
      };
      return Edge();
    };
    
    IncidentIterator& operator++(){
      if (ed.parent == nullptr) return *this;

      idx++;
      size_type nedges = incidentNode.degree();
      if (idx < nedges) {
        ed = findEdge(idx);
      }else{
        parent = nullptr;
        incidentNode = Node();
        idx = {};
        ed = Edge();        
      };
      return *this;
    };
    
    bool operator==(const IncidentIterator& someInc) const{
      if (parent != someInc.parent) return false;
      if (incidentNode != someInc.incidentNode) return false;
      if (idx != someInc.idx) return false;
      if (ed != someInc.ed) return false;
      return true;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* parent;
    Node incidentNode; 
    size_type idx; // index in the incidence stuff
    value_type ed;
  };

  /** return an incident iterator starting from a root
   * 
   * */
  IncidentIterator begin_incidentIter(Node& root) const{
    return IncidentIterator(this, root);
  };

  IncidentIterator end_incidentIter() const{
    return IncidentIterator();
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(): ed(Edge()), idx{}, parent(nullptr) {};
    EdgeIterator(const Graph* init_parent, value_type init_edge): ed(init_edge),idx{}, parent(const_cast<Graph*>(init_parent)) {};

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      return ed;
    };

    EdgeIterator& operator++(){
      /** ++ operator
       * @return next node in index.. 
      **/
      if (parent == nullptr) return * this;

      idx++;
      size_type nedges = parent->num_edges();
      if (idx < nedges) {
        internal_edge ie = (parent->graph_edges[idx]);
        ed = Edge(parent,ie.node1.getID(), ie.node2.getID(),ie.edge_id);
      }else{
        parent = nullptr;
        ed = Edge();
        idx = 0;
      }
      return * this;
    };

    bool operator==(const EdgeIterator& someEdge) const{
      if (parent != someEdge.parent) return false;
      if (idx != someEdge.idx) return false;
      if (ed != someEdge.ed) return false;
      return true;
    };

   private:
    friend class Graph;
    value_type ed;
    size_type idx;
    Graph* parent;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const{
    Edge headedge(this, graph_edges[0].node1.getID(), graph_edges[0].node2.getID(), graph_edges[0].edge_id);
    return edge_iterator(this,headedge);
  };

  edge_iterator edge_end() const{
    return edge_iterator();
  };

node_value_type& graphlongestpath() {
  return longestshort;
};

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    Point point;
    size_type idx;
    size_type node_id;
    node_value_type val;
    bool visited;

    internal_node(const Point& point, size_type idx, size_type node_id, node_value_type val)
    : point(point), idx(idx), node_id(node_id), val(val), visited(false) {}
  };

  struct internal_edge{
    Node node1;
    Node node2;
    size_type idx; // in graph_nodes
    size_type edge_id;  // id of edge

    internal_edge(const Node& node1, const Node& node2, size_type idx, size_type edge_id)
    : node1(node1),node2(node2), idx(idx), edge_id(edge_id){};
  };

  std::map<size_type, size_type> node_idx; // map from a nodeID to index for graph_nodes
  std::vector<internal_node> graph_nodes;
  size_type next_nodeID;
  
  std::map<std::pair<size_type,size_type>, size_type > nodes2edge; // map from a pair of nodeIDs to index for graph_edges
  std::map<size_type, std::vector<size_type>*> node2node; // map from a nodeID to vector of nodeIDs
  std::map<size_type, size_type > edge_idx; // map from a pair of edgeID to index for graph_edges
  std::vector<internal_edge> graph_edges;
  size_type next_edgeID;

  node_value_type longestshort;  
};

#endif // CME212_GRAPH_HPP

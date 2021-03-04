#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <array>
#include <unordered_map>
#include <utility>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// reference https://www.geeksforgeeks.org/how-to-create-an-unordered_map-of-pairs-in-c/
struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const std::pair<T1, T2>& p) const
    { 
        auto hash1 = std::hash<T1>{}(p.first); 
        auto hash2 = std::hash<T2>{}(p.second); 
        return hash1 ^ hash2; 
    } 
}; 

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
  public:
    using size_type = unsigned;
    using node_value_type = V;
    using edge_value_type = E;

 private:
    size_type nodes_size_,nodes_uid_, edges_size_, edges_uid_;
    std::unordered_map<size_type, Point> internal_nodes; //uid --> point
    std::unordered_map<size_type, std::array<size_type,2>> internal_edges; //edges_uid -> //(n1,n2)
    
    std::unordered_map<size_type, size_type> nodes_deg; //uid ---> number of incident edge
    std::unordered_map<size_type, std::vector<size_type>> nodes_neighbour; //nodes_uid --> edge_uid //incident neighbours
    
    //vector of active edges, uid. 
    //index is node(i) / vector(i), values is corresponding uid
    std::vector<size_type> i2u, i2eu;
    struct nodeinfo {
      Point p_;
      node_value_type v_; //
      int idx_;  //public faacing index
      //constructor
      nodeinfo(const Point& pos, node_value_type val, const size_type idx): p_(pos),v_(val), idx_(idx){}
    };
    
    struct edgeinfo {
      int e_idx_; //if negative, means invalid
      size_type  euid_;
      edge_value_type e_v_;
      //constructor
      edgeinfo(const size_type idx, const size_type uid, edge_value_type val): e_idx_(idx), euid_(uid), e_v_(val){}
    };
    //index is uid, (is vector anw)
    std::vector<nodeinfo> v_nodeinfo_;
    std::vector<edgeinfo> v_edgeinfo_;

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
  
  std::unordered_map<size_type, V> nodes_val; //uid ----> value
  std::unordered_map<size_type, std::vector<size_type>> nodes_adjacent; //nodes_uid --> nodes_uid 
  std::unordered_map<size_type, edge_value_type> edges_val; // edge_uid --> edge value
  std::unordered_map<std::pair <size_type,size_type>,edge_value_type, hash_pair> edges_val2;
  std::unordered_map<std::pair <size_type,size_type>,size_type, hash_pair> edges_nodes2euid;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
    : nodes_size_(0), nodes_uid_(0), edges_size_(0), edges_uid_(0), internal_nodes(), internal_edges(),nodes_deg(),
    nodes_neighbour(), nodes_val(),nodes_adjacent(),edges_val2(), edges_nodes2euid() {}

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
  class Node: private totally_ordered<Node>{
    private:
      Graph* graph_; //pointer to a graph
      size_type uid_;

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
      graph_ = nullptr;
    }
    Node(const Graph* g,size_type given_uid) 
      : graph_(const_cast<Graph*>(g)), uid_(given_uid) {}

    /** Return this node's position. */
    Point& position() {
      //return graph_->internal_nodes[uid_];
      return graph_->v_nodeinfo_[uid_].p_;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->v_nodeinfo_[uid_].p_;
      // return graph_->internal_nodes[uid_];
    }


    node_value_type& value() {
      //return graph_->nodes_val[uid_];
      return graph_->v_nodeinfo_[uid_].v_;
    }
    const node_value_type& value () const{
      return graph_->v_nodeinfo_[uid_].v_;
    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      return graph_->v_nodeinfo_[uid_].idx_;
    }

    size_type degree() const{
      return graph_->nodes_deg[uid_];
    }

    incident_iterator  edge_begin() const {return IncidentIterator(graph_,uid_);}
    incident_iterator  edge_end() const {return IncidentIterator(graph_,uid_,graph_->nodes_adjacent[uid_].size());}
    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value(); done
    // const node_value_type& value() const; done
    // size_type degree() const; done
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      if ( this->graph_ == n.graph_ && this->index() == n.index())
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
      // (void) n;           // Quiet compiler warning
      if (this->index() < n.index())
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_size_;
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
  Node add_node(const Point& position, const node_value_type& val= node_value_type()) {
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    internal_nodes[nodes_uid_] = position;
    nodes_val[nodes_uid_] = val;
    nodes_deg[nodes_uid_] = 0;
    nodes_neighbour[nodes_uid_] = std::vector<size_type>();
    nodes_adjacent[nodes_uid_] = std::vector<size_type>();
    i2u.push_back(nodes_uid_);
    v_nodeinfo_.emplace_back(nodeinfo(position,val,nodes_size_)); //nodes_size_ must be managed when adding removing nodes
    nodes_uid_++;     //will never go down
    nodes_size_++;      //will decrease when removing nodes, first idx is 0.
    return Node(this, nodes_uid_-1);  
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    // bool res = (internal_nodes.find(n.uid_) ==internal_nodes.end())? false:true;
    if (internal_nodes.find(n.uid_) ==internal_nodes.end()) return false;
    if(v_nodeinfo_[n.uid_].idx_>-1) return true;
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
    //(void) i;             // Quiet compiler warning
    assert(i < size());
    return Node(this,i2u[i]);    
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
      graph_=nullptr;
    }

    Edge(const Graph* g,size_type e, size_type n1, size_type n2)
      : graph_(const_cast<Graph*>(g)),euid_(e) ,uid1_(n1), uid2_(n2){}

    Edge(const Graph* g,size_type n1, size_type n2)
      : graph_(const_cast<Graph*>(g)),euid_(0), uid1_(n1), uid2_(n2){
        size_type euid = graph_->edges_nodes2euid[std::make_pair(n1,n2)];
        euid_ = euid;
      }
    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_,uid1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_,uid2_);      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool res = (uid1_ == e.uid1_ && uid2_ == e.uid2_ && euid_== e.euid_ && graph_==e.graph_) ? true:false;
      return res;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    /*
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      if (uid1_<e.uid1_)
        return true;
      if (uid1_>e.uid1_)
        return false;
      if (uid2_<e.uid2_) // it means uid1 are equal! compare uid2
        return true;
      return false;
    }*/

    bool operator<(const Edge& e) const {
        //same graph comparison
        if(this->graph_ == e.graph_){
          if (uid1_<e.uid1_)
            return true;
          if (uid1_>e.uid1_)
            return false;
          if (uid2_<e.uid2_) // it means uid1 are equal! compare uid2
            return true;
          return false;
        }
        else {
            //different graph comparison
            return std::less<Graph*>{}(this->graph_, e.graph_);
        }
    }

    // get edge value
    edge_value_type& value(){
      assert(this->graph_ != NULL);
      /*
      size_type n1 = node1().index();
      size_type n2 = node2().index();
      if (n1>n2) {
        n1=n2;
        n2 = node1().index();
      }
      return (*graph_).edges_val2.at(std::make_pair(n1,n2)); 
      */
      return graph_->v_edgeinfo_[euid_].e_v_;
    }
    const edge_value_type& value() const{
      assert(this->graph_ != NULL);
      return graph_->v_edgeinfo_[euid_].e_v_;
    }
    
    double length() const{return norm(node1().position() - node2().position());}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type euid_, uid1_, uid2_ ;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i<num_edges());
    // return Edge(this,internal_edges.at(i)[0],internal_edges.at(i)[1]);
    size_type euid = i2eu[i];
    size_type n1 = internal_edges.at(euid)[0];
    size_type n2 = internal_edges.at(euid)[1];
    return Edge(this,euid,n1,n2);
    }
  // smaller test case:  

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    size_type n1 = a.uid_;
    size_type n2 = b.uid_;
    if (n1>n2) {
      n1=n2;
      n2 = a.uid_;
    }
    
    auto pair__ = std::make_pair(n1,n2);
    if (edges_nodes2euid.find(pair__)==edges_nodes2euid.end()) {
      
      return false;
    }

    size_type euid = edges_nodes2euid.at(pair__);
    if(v_edgeinfo_[euid].e_idx_>-1) {
      return true;
    }
    return false;
    //if(edges_nodes2euid.find(pair__)==edges_nodes2euid.end()) return false;
    //return true;
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
    size_type n1 = a.index();
    size_type n2 = b.index();
    if (n1>n2) {
      n1=n2;
      n2 = a.index();
    }
    if (has_edge(a,b))
      return Edge(this,edges_nodes2euid[std::make_pair(n1,n2)] , a.index(),b.index()); // found an existing edge, will return it
    std::array<size_type,2> t= {n1,n2};
    internal_edges[edges_uid_] = t;
    nodes_deg[n1] += 1;
    nodes_deg[n2] += 1;
    nodes_neighbour[n1].push_back(edges_uid_); //incident edges, we only store in n1. 
    nodes_adjacent[n1].push_back(n2);
    nodes_adjacent[n2].push_back(n1);
    i2eu.push_back(edges_uid_);
    v_edgeinfo_.emplace_back(edgeinfo(edges_size_,edges_uid_,edge_value_type()));

    edge_value_type ev = edge_value_type();
    edges_val[edges_uid_] = ev;
    edges_val2[std::make_pair(n1,n2)] = ev;
    edges_nodes2euid[std::make_pair(n1,n2)] = edges_uid_;
    edges_uid_++;
    edges_size_++;
    return Edge(this,edges_uid_-1,n1,n2);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    internal_edges.clear();
    internal_nodes.clear();
    nodes_size_ = 0;
    nodes_uid_ = 0;
    edges_size_ = 0;
    edges_uid_ = 0;
    nodes_deg.clear();
    nodes_neighbour.clear();
    nodes_val.clear();
    nodes_adjacent.clear();
    edges_val.clear();
    edges_val2.clear();
    edges_nodes2euid.clear();
    i2u.clear();
    i2eu.clear();
  
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
      graph_ =nullptr;
    }

    Node operator*() const {return Node(graph_,uid_);}
    NodeIterator& operator++() {
      uid_++;
      return *this;
    }

    bool operator==(const NodeIterator& ni2) const {return ((graph_==ni2.graph_) && (uid_==ni2.uid_));}
    bool operator!=(const NodeIterator& ni2) const {return !(*this==ni2);}
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_;
    NodeIterator(const Graph* g) :graph_(const_cast<Graph*>(g)), uid_(0) {}
    NodeIterator(const Graph* g, size_type i) : graph_(const_cast<Graph*>(g)), uid_(i) {}
   };


  node_iterator node_begin() const {return NodeIterator(this);}
  node_iterator node_end() const {return NodeIterator(this,size());}
  
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
    IncidentIterator() { graph_=nullptr;}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const {
      //size_type euid = graph_->nodes_neighbour[uid_][curr_ind];
      //size_type n2 = graph_->internal_edges[euid][1];
      size_type n2 = graph_->nodes_adjacent[uid_][curr_ind];
      size_type euid;
      if (n2<uid_) euid = graph_->edges_nodes2euid[std::make_pair(n2,uid_)];
      else euid = graph_->edges_nodes2euid[std::make_pair(uid_,n2)];
      return Edge(graph_,euid,uid_,n2);}
    IncidentIterator& operator++() {
      curr_ind++;
      return *this;
    }

    bool operator==(const IncidentIterator& e2) const {
      return ((graph_==e2.graph_) && (uid_==e2.uid_) && (curr_ind == e2.curr_ind));}
    bool operator!=(const IncidentIterator& e2) const {return !(*this==e2);}

   private:
    friend class Graph;
    Graph* graph_;
    size_type uid_; //node uid, will not change in the iterator
    size_type curr_ind;
    IncidentIterator(const Graph* g,size_type i) :graph_(const_cast<Graph*>(g)), uid_(i), curr_ind(0) {}
    IncidentIterator(const Graph* g, size_type i,size_type j) : graph_(const_cast<Graph*>(g)), uid_(i), curr_ind(j) {}
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
    EdgeIterator() {graph_ =nullptr;}

    Edge operator*() const {return graph_->edge(edge_uid);}
    
    EdgeIterator& operator++() {
      edge_uid++;
       
      return *this;
    }

    bool operator==(const EdgeIterator& e2) const {
      return ((graph_==e2.graph_) && (edge_uid==e2.edge_uid));}
    bool operator!=(const EdgeIterator& e2) const {return !(*this==e2);}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    Graph* graph_;
    size_type edge_uid; 
    EdgeIterator(const Graph* g) :graph_(const_cast<Graph*>(g)), edge_uid(0) {}
    EdgeIterator(const Graph* g, size_type i) : graph_(const_cast<Graph*>(g)), edge_uid(i) {}
  }; 

edge_iterator edge_begin() const {return EdgeIterator(this);}
edge_iterator edge_end() const {return EdgeIterator(this,num_edges());}

size_type remove_edge(const Node& a, const Node&b) {
  // must have node
  assert(has_node(a));
  assert(has_node(b));
  assert(has_node(a) && has_node(b) && a!= b);
  if (!has_edge(a,b)) {
    return (size_type) false;
  }
  size_type n1,n2;
  n1 = a.uid_; 
  n2 = b.uid_;
  if (n2<n1) {
    n1=n2;
    n2=a.uid_;
  }

  size_type euid = edges_nodes2euid[std::make_pair(n1,n2)];
  size_type eidx = v_edgeinfo_[euid].e_idx_;
  
  // we will edit i2eu
  // swap eidx with last index
  size_type last_euid = i2eu.back();
  //size_type last_ind = i2eu.size()-1;

  //swapping & popping
  std::iter_swap(i2eu.begin() +eidx, i2eu.end()-1);
  i2eu.pop_back();

  // changing v_edgeinfo_
  v_edgeinfo_[euid].e_idx_ = -1;
  if (euid!=last_euid) v_edgeinfo_[last_euid].e_idx_ = eidx;

  //nodes_deg
  nodes_deg[n1]-=1;
  nodes_deg[n2]-=1;


  
  //nodes neighour
  /// need loop through nodes_neighbour[n1] and remove the euid
  /// http://www.cplusplus.com/reference/vector/vector/erase/
  auto& vec = nodes_neighbour[n1];
  for(auto iter = vec.begin();iter !=vec.end();++iter){
    if((*iter)==euid){
      vec.erase(iter);
      break;
    }
  }
  //nodes adjacent
  // need to run twice, for nodes_adjacent[n1], nodes_adjacent[n2]
  auto& vec1 = nodes_adjacent[n1];
  for(auto iter = vec1.begin();iter !=vec1.end();++iter){
    if((*iter)==n2){
      vec1.erase(iter);
      break;
    }
  }

  auto& vec2 = nodes_adjacent[n2];
  for(auto iter = vec2.begin();iter !=vec2.end();++iter){
    if((*iter)==n1){
      vec2.erase(iter);
      break;
    }
  }
  //rmb to decrease num count edge 
  edges_size_--;
  return (size_type) true;
}

//remove edge using edge as arguement
size_type remove_edge(const Edge& e){return remove_edge(e.node1(),e.node2());}

size_type remove_edge(edge_iterator e_it){
  remove_edge(*e_it);
  return EdgeIterator(); //not valid, but who cares, no speicification on it. except must be of that class
}

size_type remove_node(const Node& a){
  if(!has_node(a)) return (size_type) false;

  size_type nuid = a.uid_;
  size_type nidx = v_nodeinfo_[nuid].idx_;



  //changing degree, a) ownself become 0, b) all adjacent node -=1
  //ignore above comment, both scenario handled in remove edges. we just need iterate through remove edges

  while(nodes_adjacent[nuid].size()>0){
    size_type nuid2 = nodes_adjacent[nuid][0];
    size_type nuid1 =nuid;
    if(nuid1>nuid2){
      nuid1=nuid2;
      nuid2=nuid;
    }
    auto pair__ = std::make_pair(nuid1,nuid2);
    auto euid = edges_nodes2euid[pair__];
    remove_edge(Edge(this,euid,nuid1,nuid2));
  }

  size_type last_nuid = i2u.back();
  //swap and pop
  std::iter_swap(i2u.begin()+nidx, i2u.end()-1);
  i2u.pop_back();

  //changing v_nodeinfo_
  v_nodeinfo_[nuid].idx_ = -1;
  if(nuid!=last_nuid) v_nodeinfo_[last_nuid].idx_ = nidx;

  nodes_size_--;
  return (size_type) true;
}
NodeIterator remove_node(NodeIterator n_it){
  remove_node(*n_it);
  return NodeIterator();
};
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};



#endif // CME212_GRAPH_HPP

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
template <typename V, typename E>
class Graph {

 public:
    using size_type = unsigned;
    typedef V node_value_type;
    typedef E edge_value_type;
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

    /**@struct internel_element
     @brief used to store value/point of each node
     */
    struct internal_element {
      Point point_coords;   //The coordinate for an element
      size_type uid;      // The unique representation of each element
      size_type index;    // The index of each element
      node_value_type val;

      /**
       @brief internal_element Constructor
       @param[in] point  The position of the internel_element object
       @param[in] uid  The unique representation of each internel_element object
       @param[in,out] index  The index of this internel element/node correpsonding to this internel_element
       @param[in] val  The value of the internel_element object
       @pre _uid_ >= 0 and _index_ >= 0
       */
        internal_element(const Point& point, size_type uid, size_type idex, const node_value_type& val_ = node_value_type {})
        : point_coords(point), uid(uid), index(idex), val(val_){
        }
    };
    std::vector<internal_element> elements_; //elements vec stores all internal_elements
    std::vector<size_type> uid2ind; /*if the node with uid _i_ is no long in the graph,
                                      uid2ind[_i_]=0, valid ind is from 1 to ...*/
    size_type size_; // the num of nodes
    size_type size_edge; // the num of edges
    size_type next_uid_; // the next uid that should be used
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
//  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph.
   */
  Graph(): size_(0), size_edge(0), next_uid_(0) {
      //HW0 code here
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
  class Node: private totally_ordered<Node> {
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
      Node(){}

    /** Return this node's position.
       @return this node's position
       @pre _uid__ is less than the size of _uid2ind_ vector*/
    Point& position()
    {
        assert(uid_ < graph_->uid2ind.size());
        size_type curr_index = graph_->uid2ind[uid_];
        return graph_->elements_[curr_index-1].point_coords;
    }

    /** Return this node's position.
     @return this node's position
     @pre _uid__ is less than the size of _uid2ind_ vector*/
    const Point& position() const {
      // HW0: YOUR CODE HERE
        assert(uid_ < graph_->uid2ind.size());
        size_type curr_index = graph_->uid2ind[uid_];
        return graph_->elements_[curr_index-1].point_coords;
    }

    /** Return this node's index, a number in the range [0, graph_size).
      @return  this node's index, a number in the range [0, graph_size)
     @pre _uid__ is less than the size of _uid2ind_ vector*/
    size_type index() const {
      // HW0: YOUR CODE HERE
        assert(uid_ < graph_->uid2ind.size());
        size_type curr_index = graph_->uid2ind[uid_];
        return curr_index-1;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

      // node_value_type& value();
    /**
    @brief returns the value of this node
    @return  the value of this node
    @pre the node exists in the graph/the index of this node is valid
     */
    node_value_type& value()
    {
        //assert node exits in the graph
        assert(graph_->uid2ind[uid_]!=(unsigned int)0);
        return graph_->elements_[graph_->uid2ind[uid_]-1].val;
    }
      // const node_value_type& value() const;
      /**
      @brief returns the value of this node
      @return the value of this node
      @pre the node exists in the graph/the index of this node is valid
       */
    const node_value_type& value() const
    {
        //assert node exits in the graph
        assert(graph_->uid2ind[uid_]!=(unsigned int)0);
        return graph_->elements_[graph_->uid2ind[uid_]-1].val;
    }

      // size_type degree() const;
      /**
      @brief returns the degree of this node
      @return the degree of this node
       */
      size_type degree() const
      {
          return graph_->adj_edges_[uid_].size();
      }
      // incident_iterator edge_begin() const;
      /**
      @brief returns an incident_iterator that points to the first edge of this node
      @return an incident_iterator that points to the first edge of this node
       */
      incident_iterator edge_begin() const
      {
          incident_iterator e=IncidentIterator(*this,0);
          return e;
      }
      // incident_iterator edge_end() const;
      /**
      @brief returns an incident_iterator that points to one pass the last edge of this node
      @return an incident_iterator that points to one pass the last edge of this node
       */
      incident_iterator edge_end() const
      {
          return IncidentIterator(*this,degree());
      }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     *@return true if the two nodes are equal, false otherwise
     */
    bool operator==(const Node& n) const {
      if (this->graph_ == n.graph_ && this->uid_ == n.uid_)
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
//        assert(this->graph_ == n.graph_);
        if (this->index() < n.index())
            return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_; //pointer to graph
    size_type uid_; // unique representation of each node
      /**
       @brief Node constructor
       @param[in] gra  a Graph pointer
       @param[in] uid  uid of this node
       */
    Node(const Graph* gra, size_type uid):
      graph_(const_cast<Graph*>(gra)), uid_(uid) {}

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
  Node add_node(const Point& position, const node_value_type& val_= node_value_type()) {
    // HW0: YOUR CODE HERE
      internal_element ie(position, next_uid_, size_, val_);
      elements_.push_back(ie);
      uid2ind.push_back(size_+1);
      ++size_;
      ++next_uid_;
      adj_edges_.push_back(std::vector<internel_edge>());
      return Node(this, next_uid_-1);
  }
    /** Remove a node from the graph, returning an unsinged number indicating remove or not
     * @param[in] n_rm The node to be removed
     * @return 1: successfully removed the given node; 0: no node removed
     * @post new num_nodes() == old num_nodes() - 1
     * @post new num_edges()==old num_edges()-degree(n_rm)
     *
     * Complexity: O(degree(n_rm)(n+m)) amortized operations.
     */
  size_type remove_node(const Node& n_rm)
  {
      //if the node doesn't exist
       if (not has_node(n_rm))
           return 0;
      std::vector<Edge> incident_edges {};
      //loop over all incident edges of the node
      for (auto n_it = n_rm.edge_begin(); n_it!=n_rm.edge_end();++n_it)
      {
          incident_edges.push_back(*n_it);
      }
      //remove all incident edges of the node
      while(!incident_edges.empty())
      {
          Edge e_rm = incident_edges.back();
          incident_edges.pop_back();
          remove_edge(e_rm);
      }
      size_type curr_idx = uid2ind[n_rm.uid_]-1;
      elements_[curr_idx] = elements_.back();
      uid2ind[elements_[curr_idx].uid] = curr_idx+1;
      elements_[curr_idx].index = curr_idx;
      elements_.pop_back();
      uid2ind[n_rm.uid_] = 0;
      size_--;
      return 1;

  }
    /** Remove a node from the graph, returning a node_iterator object
     * @param[in] n_it a node_iterator points to the node to be removed
     * @return node_iterator points to a valid node object
     * @post new num_nodes() == old num_nodes() - 1
     * @post new num_edges()==old num_edges()-degree(n_rm)
     *
     * Complexity: O(degree(n_rm)(n+m)) amortized operations.
     */
  node_iterator remove_node(node_iterator n_it)
    {
        assert(n_it!=node_end());
        size_type rm_succ = remove_node(*n_it);
        //if no node removed, return the given node_iterator
        if (rm_succ==0)
            return n_it;
        ++n_it;
        // if n_it is not pointing one past the end node element
        if (n_it!=node_end())
            return n_it;
        //--functionality_1
        //--may return node that has already been removed
        //--START
        return node_iterator(this,elements_[0].uid);
        //--END
    }
  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (uid2ind[n.uid_] != (unsigned)0 && n.graph_ == this)
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
      assert(i>=0 && i<size_);
      return Node(this, elements_[i].uid);
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_e, uid1);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_e, uid2);      // Invalid Node
    }
      /**returns the edge's value
       @return edge_value_type: edge's value
       @complexity O(degree(node1())
       */
    edge_value_type& value()
    {
        //loop over all incident edges of node1()
        for (unsigned i=0; i<(graph_e->adj_edges_[uid1]).size(); i++)
        {
            //if find the edge
            if(graph_e->adj_edges_[uid1][i].adj_uid==uid2)
            {
                return graph_e->adj_edges_[uid1][i].edge_val;
            }
        }
        std::cout<<"cannot find this edge in the graph"<<std::endl;
        return graph_e->adj_edges_[uid1][0].edge_val;
    }
      /**returns the edge's value
       @return edge_value_type: edge's value
       @complexity O(degree(node1())
       */
    const edge_value_type& value() const
      {
          return value();
      }
      /**returns the edge's length
       @return double: edge's length
       @complexity O(1)
       */
    double length() const
    {
        return norm(node1().position()-node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (e.graph_e != graph_e)
          return false;
      if ((e.node1()==this->node2() && e.node2()==this->node1()) ||
          (e.node1()==this->node1() && e.node2()==this->node2()))
          return true;
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    //--functionality_1
    //--edge comparator should not be based on length, equal-length edges e1 and e2 have e1 < e2 and e2 < e1
    //--START
    bool operator<(const Edge& e) const {
        double l1 = length();
        double l2 = e.length();
        if (l1<=l2)
            return true;
        else
            return false;
    }
    //--END

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_e;
    size_type uid1;
    size_type uid2;
    /**
       @brief Edge constructor
       @param[in] gra  a Graph pointer
       @param[in] uid1  uid of the first node of this edge
       @param[in] uid2  uid of the first node of this edge
    */
    Edge(const Graph* gra, size_type uid1, size_type uid2):
      graph_e(const_cast<Graph*>(gra)), uid1(uid1), uid2(uid2) {};

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return size_edge;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
      assert(i >= (unsigned)0 && i<size_edge);
      return elements_e[i];
}

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
      assert(this->has_node(a) && this->has_node(b));
      //loop over all nodes adjacent to a
      for (size_type i=0; i<a.degree(); i++)
      {
          /*if the uid of b and the uid of
            the adjacent node of a equal*/
          if (b.uid_==adj_edges_[a.uid_][i].adj_uid)
              return true;
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
      //if the edge has existed, then do nothing and return this edge
      if (this->has_edge(a,b))
        return Edge(this, a.uid_, b.uid_);
      Edge n(this, a.uid_, b.uid_);
      elements_e.push_back(n);
      size_edge++;
      adj_edges_[a.uid_].push_back(internel_edge(a.uid_, b.uid_));
      adj_edges_[b.uid_].push_back(internel_edge(b.uid_, a.uid_));
      return n;        // Invalid Edge
  }
    /**
     removes the edge with node1=a and node2=b
     @param[in] a: node1 of the edge to be removed
     @param[in] b: node2 of the edge to be removed
     @return size_type: 1 indicating the edge is removed successfully; 0 indicating no edge's been removed
    @post new num_edges() == old num_edges()-1
     @complexity: O(m+degree(node1())+degree(node2())
     */
  size_type remove_edge(const Node& a, const Node& b)
    {
        Edge n(this, a.uid_, b.uid_);
        //if the edge doesn't exist
        if (not has_edge(a,b))
            return 0;
        //loop over all edges and find the given one and remove it by swapping it with
        // the last edge in the vector and pop it
        for (unsigned i=0; i<elements_e.size();i++)
        {
            //if find the given edge in the vector
            if(elements_e[i]==n)
            {
                elements_e[i]=elements_e.back();
                elements_e.pop_back();
            }
        }
        size_edge--;
        //if the edge doesn't exist
        if (adj_edges_[a.uid_].size()==0 || adj_edges_[b.uid_].size()==0)
            return 0;
        //loop over all incident edges of node1, and remove the element froom adj_edges_
        for(unsigned i=0; i<adj_edges_[a.uid_].size();i++)
        {
            //if find the edge
            if (adj_edges_[a.uid_][i].adj_uid==b.uid_)
            {
                // if it's the only edge incident to node1
                if (adj_edges_[a.uid_].size()==1)
                    adj_edges_[a.uid_].pop_back();
                else
                {
                    adj_edges_[a.uid_][i] = adj_edges_[a.uid_].back();
                    adj_edges_[a.uid_].pop_back();
                }
                break;
            }
        }
        //loop over all incident edges of node2, and remove the element froom adj_edges_
        for(unsigned i =0; i<adj_edges_[b.uid_].size();i++)
        {
            //if find the edge
            if (adj_edges_[b.uid_][i].adj_uid==a.uid_)
            {
                // if it's the only edge incident to node1
                if(adj_edges_[b.uid_].size()==1)
                    adj_edges_[b.uid_].pop_back();
                else
                {
                    adj_edges_[b.uid_][i] = adj_edges_[b.uid_].back();
                    adj_edges_[b.uid_].pop_back();
                }
                break;
            }
        }
        return 1;
    }
    /**
     removes the edge given
     @param[in] e_rm: the edge to be removed
     @return size_type: 1 indicating the edge is removed successfully; 0 indicating no edge's been removed
    @post new num_edges() == old num_edges()-1
     @complexity: O(m+degree(node1())+degree(node2())
     */
  size_type remove_edge(const Edge& e_rm)
    {
        Node node1 = e_rm.node1();
        Node node2 = e_rm.node2();
        return remove_edge(node1, node2);
    }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
      size_ = 0;
      size_edge = 0;
      elements_.clear();
      elements_e.clear();
      uid2ind.clear();
      next_uid_=0;
      adj_edges_.clear();
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
      /**
       @brief dereference the iterator/returns the object the iterator points to
       @return the node that this iterator points to
       */
      Node operator*() const
      {
          return Node(graph_itr_node, uid_itr_node);
      }
    // NodeIterator& operator++()
      /**
       @brief make the iterator points to the next element in the sequence
       @return the iterator which points to the next element in the sequence
       @post one past the last element is the furtherest the iterator can point to
       */
      NodeIterator& operator++()
      {
          unsigned int new_node_idx = (graph_itr_node->uid2ind[uid_itr_node]); //ind is from 1
          //if the iterator doesn't point at one past the last element
          if (new_node_idx < graph_itr_node->elements_.size())
              uid_itr_node = (graph_itr_node->elements_[new_node_idx]).uid;
          else
              uid_itr_node = graph_itr_node->next_uid_;
          return *this;
      }
    // bool operator==(const NodeIterator&) const
      /**
       @brief comparing two iterators, equal if they belong to the same graph and point to the same node
       @param[in] i2 the iterator to compare with
       @return true if they belong to the same graph and point to the same node
       */
      bool operator==(const NodeIterator& i2) const
      {
          return (graph_itr_node == i2.graph_itr_node
                  && uid_itr_node == i2.uid_itr_node);
      }
      /**
       @brief comparing two iterators, not equal if they don't belong to the same graph or they don't point to the same node
       @param[in] i2 the iterator to compare with
       @return true if they don't belong to the same graph or they don't point to the same node
       */
      bool operator!=(const NodeIterator& i2) const
      {
          return !(*this == i2);
      }

   private:
    friend class Graph;
    Graph* graph_itr_node;
    size_type uid_itr_node;
    /**
     @brief NodeIterator Constructor
     @param[in] gra a Graph pointer
     @param[in] uid_tmp the uid of the node the iterator points to
     */
    NodeIterator(const Graph* gra, size_type uid_tmp):
      graph_itr_node(const_cast<Graph*>(gra)),uid_itr_node(uid_tmp){}
  };


  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
    /**
     @brief returns an node iterator that points to the first node in the graph
     @return an node iterator that points to the first node in the graph
     */
    node_iterator node_begin() const
    {
            return NodeIterator(this, 0);
    }
  // node_iterator node_end() const
    /**
     @brief returns an node iterator that points to one pass the last node in the graph
     @return an node iterator that points to one pass the last node in the graph
     */
    node_iterator node_end() const
    {
        return NodeIterator(this, next_uid_);
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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
      /**
       @brief dereference the iterator/returns the edge this iterator points to
       @return the edge this iterator points to
       */
      Edge operator*() const
      {
          return Edge(graph_incid_iter, center_uid_incid_iter,\
                      graph_incid_iter->adj_edges_[center_uid_incid_iter][idx_incid_iter].adj_uid);
      }
    // IncidentIterator& operator++()
      /**
       @brief make the iterator points to the next edge in the sequence
       @return an iterator points to the next edge in the sequence
       @post one past the last element is the furtherest the iterator can point to
       */
      IncidentIterator& operator++()
      {
          idx_incid_iter++;
          return *this;
      }
    // bool operator==(const IncidentIterator&) const
      /**
       @brief comparing two iterators, equal if they belong to the same graph and point to the edge
       @param[in] inc_iter the iterator to compare with
       @return true if they belong to the same graph and point to the edge
       */
      bool operator==(const IncidentIterator& inc_iter) const
      {
          return (graph_incid_iter == inc_iter.graph_incid_iter
                  && center_uid_incid_iter == inc_iter.center_uid_incid_iter
                  && idx_incid_iter == inc_iter.idx_incid_iter);
      }
      /**
       @brief comparing two iterators, not equal if they don't belong to the same graph or don't point to the edge
       @param[in] inc_iter the iterator to compare with
       @return true if they don't belong to the same graph or don't point to the edge
       */
      bool operator!=(const IncidentIterator& inc_iter) const
      {
          return !(*this==inc_iter);
      }

   private:
    friend class Graph;
      Graph* graph_incid_iter;
      size_type center_uid_incid_iter;//uid of the center node
      size_type idx_incid_iter; //index of the current adj node of the center node
      /**
       @brief IncidentIterator Constructor
       @param[in] node_inc_iter the node the incident iterator will iterate over the edges of
       @param[in] ind the index of the edge of _node_inc_iter_ that this iterator points to
       */
      IncidentIterator(const Node& node_inc_iter, size_type ind):graph_incid_iter(const_cast<Graph*>(node_inc_iter.graph_)), center_uid_incid_iter(node_inc_iter.uid_),idx_incid_iter(ind){}


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

    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
      /**
       @brief dereference the iterator/returns the edge this iterator points to
       @return the edge this iterator points to
       */
      Edge operator*() const
      {
          return Edge(graph_edge_iter, first_node_uid,\
                      graph_edge_iter->adj_edges_[first_node_uid][second_node_idx].adj_uid);
      }
    // EdgeIterator& operator++()
      /**
       @brief make the iterator points to the next edge in the graph
       @return an iterator points to the next edge in the graph
       @post one past the last element is the furtherest the iterator can point to
       */
      EdgeIterator& operator++()
      {
          second_node_idx++;
          increment_helper();
          return *this;

      }
    // bool operator==(const EdgeIterator&) const
      /**
       @brief comparing two iterators, equal if they belong to the same graph and point to the edge
       @param[in] inc_iter the iterator to compare with
       @return true if they belong to the same graph and point to the edge
       */
      bool operator==(const EdgeIterator& ie) const
      {
          return (graph_edge_iter == ie.graph_edge_iter
                  && first_node_uid == ie.first_node_uid
                  && second_node_idx == ie.second_node_idx);

      }
      /**
       @brief comparing two iterators, not equal if they don't belong to the same graph or don't point to the edge
       @param[in] inc_iter the iterator to compare with
       @return true if they don't belong to the same graph or don't point to the edge
       */
      bool operator!=(const EdgeIterator& ie) const
      {
          return !(*this==ie);
      }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
      Graph* graph_edge_iter;
      size_type first_node_uid; //uid of the first node
      size_type second_node_idx; //index of the current adj node of the center node
      /**
      @brief helper function for operator++, takes care of the case where the \
       iterator has reached the end/the edge has been visited/need to go to the next node
      @post the edge the iterator points to is reasonable and every edge will only be visited once
       */
      void increment_helper()
      {
          //if iter already reaches the end
          if (first_node_uid >= graph_edge_iter->next_uid_)
          {
              first_node_uid = graph_edge_iter->next_uid_;
              second_node_idx = 0;
          }
          else
          {
              if (second_node_idx >= graph_edge_iter->adj_edges_[first_node_uid].size())
              {
                  first_node_uid++;
                  second_node_idx = 0;
                  increment_helper();
              }
              //else if the edge has been looped over before
              else if (graph_edge_iter->adj_edges_[first_node_uid][second_node_idx].adj_uid < first_node_uid)
              {
                  second_node_idx++;
                  increment_helper();
              }
          }
      }
      /**
       @brief EdgeIterator Constructor
       @param[in] gra a Graph pointer
       @param[in] center_uid uid of the first node of the edge that this iterator points to
       @param[in] adj_ind index of the second node in the vector of nodes that adjacent to the first node
       */
      EdgeIterator(const Graph* gra, size_type center_uid, size_type adj_ind):
              graph_edge_iter(const_cast<Graph*>(gra)), first_node_uid(center_uid), second_node_idx(adj_ind)
     {
         increment_helper();
     }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
    /**
     @brief returns an edge iterator that points to the first edge in the graph
     @return an egde iterator that points to the first edge in the graph
     */
    edge_iterator edge_begin() const
    {
        return EdgeIterator(this, 0, 0);
    }
  // edge_iterator edge_end() const
    /**
     @brief returns an node iterator that points to one pass the last edge in the graph
     @return an node iterator that points to one pass the last edge in the graph
     */
    edge_iterator edge_end() const
    {
        return EdgeIterator(this, next_uid_, 0);
    }
    /**
     @brief removes the edge pointed by _e_it_returns an edge iterator that points to a valid edge in the graph
     @param[in] e_it a edge iterator object that points to an edge to be removed
     @return returns an edge iterator that points to a valid edge in the graph
     @post new num_edges() == old num_edges()-1
      @complexity: O(m+degree(node1())+degree(node2())
     */
    edge_iterator remove_edge(edge_iterator e_it)
    {
        size_type rm_suc = remove_edge(*e_it);
        //if no edge is removed
        if (rm_suc==0)
            return e_it;
        e_it.increment_helper();
        return e_it;
    }

 private:

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    /**
     @brief struct that stores edge info of each edge
     */
    struct internel_edge
    {
        size_type center_uid; //uid of the first node
        size_type adj_uid; //uid of the second node that is adjacent to the first node
        edge_value_type edge_val;
        /**
         @brief internel_edge constructor
         @param center_uid_ uid of the first node
         @param adj_uid_ uid of the second node that is adjacent to the first node
         */
        internel_edge(size_type center_uid_, size_type adj_uid_,edge_value_type val = edge_value_type()):
                  center_uid(center_uid_), adj_uid(adj_uid_), edge_val(val)
        {
        }
    };
    std::vector<Edge> elements_e; //vec stores all edges
    std::vector<std::vector<internel_edge>> adj_edges_; /*adj_edges_[i] is a vec of edges
                                                        incident to node with uid=i, and the
                                                        edge info is stored in internel_edge*/
};

#endif // CME212_GRAPH_HPP

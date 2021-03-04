#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <map>
#include <set>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/* util func for swapping arbitrary elements (taken from lec notes)*/
template < class T>
void swap (T &a, T &b) {
    T tmp { static_cast <T&& >(a )};
    a = static_cast <T&& >(b);
    b = static_cast <T&& >(tmp );
}
/* util func for swapping arbitrary elements (taken from lec notes)*/
template < class T>
void swap_graph_hpp (T &a, T &b) {
    T tmp { static_cast <T&& >(a )};
    a = static_cast <T&& >(b);
    b = static_cast <T&& >(tmp );
}



/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {
    private:
    public:
        /** Type of templated variable */
        using node_value_type = V;
        using edge_value_type = E;

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

        /* structures for nodes and edges to hold arbitrary data*/
        struct NodeInfo {
            Point p_;
            node_value_type v_;
            size_type idx_;
        };
        struct EdgeInfo {
            edge_value_type v_;
            size_type idx_;
        };
        // utility function
        bool id_in_set(std::set<size_type> s, size_type el) const {
            for (auto it=s.begin(); it!=s.end(); ++it) {
                if (*it==el)
                    return true;
            }
            return false;
        }

        //
        // CONSTRUCTORS AND DESTRUCTOR
        //

        /** Construct an empty graph. */
        Graph() {}

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
                Node() {
                    uid_=-1;
                    graph_=nullptr;
                }

                /** Return this node's position. */
                const Point& position() const {
                    return graph_->node_infos_[uid_].p_;
                }

                Point& position() {
                    return graph_->node_infos_[uid_].p_;
                }

                /** Return this node's index, a number in the range [0, graph_size). */
                size_type index() const {
                    return graph_->node_infos_[uid_].idx_;
                }

                /** Get reference to this nodes value. The caller acn use the
                 * reference to set the value.
                 * O(1)
                 */
                node_value_type& value() {
                    return graph_->node_infos_[uid_].v_;
                }

                /** Get reference to this nodes value that is caller cannot
                 * modify.
                 * O(1)
                 */
                const node_value_type& value() const {
                    return graph_->node_infos_[uid_].idx_;
                }

                /** Get number of edges incident to this node
                 *  O(1)
                 */
                size_type degree() const {
                    return graph_->node_to_nodes_[uid_].size();
                }

                incident_iterator edge_begin() const {
                    IncidentIterator ie {graph_, uid_, 0};
                    return ie;
                }

                 incident_iterator edge_end() const {
                    size_type incident_edges = graph_->node_to_edges_[uid_].size();
                    IncidentIterator ie {graph_, uid_, incident_edges};
                    return ie;
                 }

                /** Test whether this node and @a n are equal.
                 *
                 * Equal nodes have the same graph and the same index.
                 */
                bool operator==(const Node& n) const {
                    return  (this->uid_==n.uid_ && this->graph_==n.graph_);
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
                    if (uid_<n.uid_)
                        return true;
                    return false;
                }

            private:
                // Allow Graph to access Node's private member data and functions.
                friend class Graph;
                graph_type* graph_;
                size_type uid_=-1;

                Node(graph_type* graph, size_type uid)
                    : graph_(graph), uid_(uid) {}
        };

        /** Return the number of nodes in the graph. Complexity: O(1). */
        size_type size() const {
            return node_i2u_.size();
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
        Node add_node(const Point& position,  const node_value_type& val= node_value_type()) {
            // Create new Node
            size_type new_uid = nodes_.size();
            Node n {this, new_uid};
            nodes_.push_back(n);

            // Matched NodeInfo object. Idx_ is determined by size of node_i2u
            NodeInfo node_info;
            node_info.idx_ = node_i2u_.size();
            node_info.p_ = Point(position.x, position.y, position.z);
            node_info.v_ = val;
            node_infos_.push_back(node_info);

            // update node_i2u. add to back
            node_i2u_.push_back(new_uid);
            // Instantiate an entry in the nodes adjacency list, even if it's empty
            // This allows us to call degree() on this node uid.
            node_to_nodes_[new_uid] = {};
            return n;
        }

        /** Remove one node from the graph and all its incident edges.
         *  @pre this->has_node(n)==True
         *  @post num_nodes() == old num_node()-1
         *  @post num_edges() <= old num_edges()
         *
         *
         *  Existing NodeIterators are still valid only if they were pointing
         *  at nodes having index() smaller than the removed node's index().
         *
         *  EdgeIterator's and IncidentIterators may become invalid since
         *  multiple edges may be deleted as a side effect of this function.
         *
         *  This will change the result of n.index() for existing nodes.
         *
         *  Returns 1 if node removed. 0 otherwise.
         *
         *  Complexity: O(n.degree()), since incident edges are deleted.
         */
        size_type remove_node (const Node& n) {
            if (!has_node(n))
                return 0;
            // remove incident edges
            std::vector<Edge> e;
            for (incident_iterator it=n.edge_begin(); it!=n.edge_end(); ) {
                e.push_back(*it);
                ++it;
            }
            for (auto it=e.begin(); it!=e.end(); it++) {
                remove_edge(*it);
            }
            // record uid and indx
            size_type node_uid = n.uid_;
            size_type node_idx = node_infos_[node_uid].idx_;
            // get uid and idx of the last node
            size_type swap_node_idx = node_i2u_.size()-1;
            size_type swap_node_uid = node_i2u_[swap_node_idx];
            // swap
            std::swap(node_i2u_[node_idx] , node_i2u_[swap_node_idx]);
            // update the idx info
            node_infos_[swap_node_uid].idx_ = node_idx;
            // deleted node has an invalid indx
            node_infos_[node_uid].idx_ = -1;
            // delete the node
            node_i2u_.pop_back();

            return 1;
        }

        /** Remove node referenced by n_it.
         *  Return an iterator to the next node
         *  This does not invalidate anything.
         *
         *  @pre n_it is not an end() iterator
         *  @post *n_it == old *(++n_it)
         *
         * Complexity: O((*it).degree()) due to call to remove_node(*n_it)
         */
        node_iterator remove_node(node_iterator n_it) {
            const Node& n = *n_it;
            remove_node(n);
            return n_it;
        }

        /** Determine if a Node belongs to this Graph
         * @return True if @a n is currently a Node of this Graph
         *
         * Complexity: O(1).
         * Relies on invalid nodes habing NodeInfo.idx_=-1
         */
        bool has_node(const Node& n) const {
            return  (n.graph_ == this && n.index() < num_nodes());
        }

        /** Return the node with index @a i.
         * @pre 0 <= @a i < num_nodes()
         * @post result_node.index() == i
         *
         * Complexity: O(1).
         */
        Node node(size_type i) const {
            if (i>= num_nodes())
                throw std::invalid_argument("Index too large.");
            size_type uid = node_i2u_[i];
            return nodes_[uid];
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
                Edge() {}

                /** Return Edge index */
                size_type index() {
                    return graph_->edge_infos_[uid_].idx_;
                }

                edge_value_type& value() {
                    return graph_->edge_infos_[uid_].v_;
                }

                const edge_value_type& value() const {
                    return graph_->edge_infos_[uid_].v_;
                }
                /** Acts in place */
                void swap_edge_nodes_inplace_() {
                    size_type tmp = node1_uid_;
                    node1_uid_ = node2_uid_;
                    node2_uid_ = tmp;
                    //swap_graph_hpp(this->node1_, this->node2_);
                }

                /** Return a node of this Edge */
                Node node1() const {
                    return graph_->nodes_[node1_uid_];
                }

                /** Return the other node of this Edge */
                Node node2() const {
                    return graph_->nodes_[node2_uid_];
                }

                /** Test whether this edge and @a e are equal.
                 *
                 * Equal edges represent the same undirected edge between two nodes.
                 */
                bool operator==(const Edge& e) const {
                    if (
                             (
                              e.node1() == this->node1() && e.node2() == this->node2()
                             ) ||
                             (
                              e.node1() == this->node2() && e.node2() == this->node1()
                             )
                       ) return true;
                    return false;
                }

                /** Test whether this edge is less than @a e in a global order.
                 *
                 * This ordering function is useful for STL containers such as
                 * std::map<>. It need not have any interpretive meaning.
                 */
                bool operator<(const Edge& e) const {
                    if (this->graph_ == e.graph_)
                        return (this->uid_) < e.uid_;
                    else
                        return this->graph_ < e.graph_;
                }

                double length() {
                    Point diff = node1().position() - node2().position();
                    return norm(diff);
                }

            private:
                // Allow Graph to access Edge's private member data and functions.
                friend class Graph;
                // Note that its a non-const pointer to const node
                //const Graph::Node* node1_;
                //const Graph::Node* node2_;
                size_type node1_uid_;
                size_type node2_uid_;
                size_type uid_=-1;
                graph_type* graph_;

                /* Private constructor*/
                Edge(const Node& a, const Node& b, graph_type* graph, size_type uid)
                    : node1_uid_(a.uid_), node2_uid_(b.uid_), uid_(uid), graph_(graph) {}
        };

        /** Return the total number of edges in the graph.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        size_type num_edges() const {
            return edge_i2u_.size();
        }

        /** Return the edge with index @a i.
         * @pre 0 <= @a i < num_edges()
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        Edge edge(size_type i) const {
            if (i>= num_edges())
                throw std::invalid_argument("Index too large.");
            size_type uid = edge_i2u_[i];
            return edges_[uid];
        }

        /** Test whether two nodes are connected by an edge.
         * @pre @a a and @a b are valid nodes of this graph
         * @return True if for some @a i, edge(@a i) connects @a a and @a b.
         *
         * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
         */
        bool has_edge(const Node& a, const Node& b) const {
            if (node_to_nodes_.count(a.uid_)>0
               && node_to_nodes_.at(a.uid_).count(b.uid_)>0)  {
                return true;
            }
            if (node_to_nodes_.count(b.uid_)>0
               && node_to_nodes_.at(b.uid_).count(a.uid_)>0)  {
                return true;
            }


            return false;
        }

        /* Given 2 nodes that share an edge, find the edge uid */
        size_type find_edge_uid_(const Node& a, const Node& b) const {
            std::set<size_type> a_edges = node_to_edges_.at(a.uid_);
            std::set<size_type> b_edges = node_to_edges_.at(b.uid_);

            for (auto it=a_edges.begin(); it!=a_edges.end(); it++) {
                if (b_edges.count(*it)>0)
                    return *it;
            }
            return -1;
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
        Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
            // Check this isn't a duplicate edge
            if (has_edge(a,b)) {
                size_type uid = find_edge_uid_(a,b);
                Edge e = edges_[uid];
                // Ensure that a is node1() in the returned edge
                if (e.node1() != a)
                    e.swap_edge_nodes_inplace_();
                return e;
            }

            if (!(a.graph_==this && b.graph_==this))
                throw "Nodes must be from this graph";

            // create new edge and add
            size_type new_uid = edges_.size();;
            Edge e {a, b, this, new_uid};
            edges_.push_back(e);
            // create edge info
            EdgeInfo edge_info;
            edge_info.idx_ = edge_i2u_.size();
            edge_info.v_ = val;
            edge_infos_.push_back(edge_info);

            // update node_i2u. Alway add to end
            edge_i2u_.push_back(new_uid);

            // update adjacency maps
            node_to_edges_[a.uid_].insert(e.uid_);
            node_to_edges_[b.uid_].insert(e.uid_);
            node_to_nodes_[a.uid_].insert(b.uid_);
            node_to_nodes_[b.uid_].insert(a.uid_);

            return e;
        }

        /* Remove edge e from graph.
         * @pre: this.has(e) is true
         * @post: num_edges() == old num_edges()-1
         *
         * Complexity: O(max(e.node1().degree(), e.node2().degree())
         *    This is because we need to delete elements from adjacency lists,
         *    which are implemented as a map from an int to a set. The sets have
         *    the size equal to degree of that node, and we may have to iterate
         *    to their end to perform a delete.
         *
         *  This will change the result of e.index() for existing edges.
         *
         * This invalidates edge iterators if they were pointing at elements
         * having index() larger than this nodes e.index()
         *
         * If edge is not present, return 0.
         */
        size_type remove_edge(const Edge& e) {
           if (!has_edge(e.node1(), e.node2()))  {
                return 0;

           // Ged uid and idx of the edge to delete
           size_type edge_uid = e.uid_;
           size_type edge_idx = edge_infos_[edge_uid].idx_;
           // Get uid and idx of edge that is been swapped from the back of the vector
           size_type swap_edge_idx = edge_i2u_.size()-1;
           size_type swap_edge_uid = edge_i2u_[swap_edge_idx];
           // swap the indx in edge_i2u and also in edge_infos_
           std::swap(edge_i2u_[edge_idx], edge_i2u_[swap_edge_idx]);
           // update the idx_ of the swapped in element (which was previously at the end)
           edge_infos_[swap_edge_uid].idx_ = edge_idx;
           // deleted edge has an invalid index
           edge_infos_[edge_uid].idx_ = -1;
           // Now pop delete the entry
           edge_i2u_.pop_back();

           // Now we updated the adjacency lists
            delete_edge_from_nodes_to_edges(edge_uid);
            delete_edge_from_node_to_nodes(e.node1(), e.node2());

           return 1;
        }

        /** Remove edge given a pair of nodes.
         * @pre: this.has(e) is true
         * @post: num_edges() == old num_edges()-1
         *
         * If the edge does not exists, function returns 0.
         * Else the edge is removed, returning 1.
         *
         * This invalidates edge iterators if they were pointing at elements
         * having index() larger than this nodes e.index()
         *
         * Complexity: O(max(a.degree(), b.degree()) due to the calls to
         * has_edge and to remove_edge(Edge e).
         */
        size_type remove_edge(const Node& a, const Node &b) {
            // Find the edge uuid and then call the existing delete function
            if (!has_edge(a,b))
                return 0;
            size_type edge_uid = find_edge_uid_(a, b);
            Edge& e = edges_[edge_uid];
            return remove_edge(e);
        }

        /** Remove the edge pointed to by *e_it.
         * Return an iterator to the next edge.
         *
         * @pre n_it is not an end() iterator.
         * @post *n_it == old *(++n_it)
         *
         * Complexity: O(max((*e_it).node1().degree(), (*e_it).node2.degree())
         */
        edge_iterator remove_edge (edge_iterator e_it) {
            const Edge& e = *e_it;
            remove_edge(e);
            return e_it;
        }

        /** Templated utility function for remove_edge.
         * Delete from set s, the value el.
         * @post: if el is in s, then s.size()==old s.size()-1
         *
         * Complexity: O(s.size())
         * */
        template<typename A>
        size_type erase_element_from_set(std::set<A>& s, A el) {
            size_type init_size = s.size();
            for (auto it=s.begin(); it!=s.end(); ) {
                if (*it==el)
                    it=s.erase(it);
                else
                    it++;
            }
            return init_size-s.size();
        }
        /** Utility for remove_edge function.
         * Given edge_uid, remove evidence of that edge from adjacency lists.
         * @pre: Edge with e.uid_=edge_uid must be in the graph, otherwise an
         * exception is thrown.
         *
         * Complexity: O(max degree of nodes in edge_uid).
         */
        void delete_edge_from_nodes_to_edges(size_type edge_uid) {
            const Node& a = edges_[edge_uid].node1();
            const Node& b = edges_[edge_uid].node2();

            // Get sets from adjacency list
            std::set<size_type>& a_edges = node_to_edges_[a.uid_];
            std::set<size_type>& b_edges = node_to_edges_[b.uid_];
            // Delete the elements using set methods
            erase_element_from_set(a_edges, edge_uid);
            erase_element_from_set(b_edges, edge_uid);

            return;
        }

        /** Utility for remove_edge function.
         * Given 2 nodes, remove evidence of that edge from adjacency lists.
         * @pre: Exists an edge with nodes a and b, with edge in the graph. Else
         * an exception is thrown.
         *
         * Complexity: O(max degree of nodes in edge_uid).
         */
        void delete_edge_from_node_to_nodes(const Node& a, const Node& b) {
            std::set<size_type>& a_nodes = node_to_nodes_[a.uid_];
            std::set<size_type>& b_nodes = node_to_nodes_[b.uid_];
            erase_element_from_set(a_nodes, b.uid_);
            erase_element_from_set(b_nodes, a.uid_);
        }

        /** Remove all nodes and edges from this graph.
         * @post num_nodes() == 0 && num_edges() == 0
         *
         * Invalidates all outstanding Node and Edge objects.
         */
        void clear() {
            for (size_type i=0; i<num_nodes(); i++) {
                size_type uid = node_i2u_[i];
                nodes_[uid].graph_ = nullptr;
                nodes_[uid].uid_= size_type(-1);
            }
            node_i2u_.clear();
            for (size_type i=0; i<num_edges(); i++) {
                size_type uid = edge_i2u_[i];
                edges_[uid].graph_ = nullptr;
                edges_[uid].uid_= size_type(-1);
            }
            edge_i2u_.clear();
        }

        //
        // Node Iterator
        //

        /** @class Graph::NodeIterator
         * @brief Iterator class for nodes. A forward iterator. */
        class NodeIterator: private totally_ordered<NodeIterator>  {
            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Node;                     // Element type
                using pointer           = Node*;                    // Pointers to elements
                using reference         = Node&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                //using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid NodeIterator. */
                NodeIterator() {}

                Node operator*() const{
                    size_type uid = graph_->node_i2u_[node_indx_];
                    return graph_->nodes_[uid];
                }
                NodeIterator& operator++() {
                    node_indx_++;
                    return *this;
                }
                bool operator==(const NodeIterator& n) const {
                    return (n.node_indx_ == this->node_indx_
                            && n.graph_==this->graph_);
                }

            private:
                friend class Graph;
                const graph_type* graph_;
                size_type node_indx_;
                NodeIterator(const graph_type* graph, size_type node_indx) {
                   graph_ = graph; node_indx_=node_indx;
                  }
        };

        node_iterator node_begin() const {
            NodeIterator ni(this, 0);
            return ni;
        }
        node_iterator node_end() const {
            NodeIterator ni {this, num_nodes()};
            return ni;
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
                //using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid IncidentIterator. */
                IncidentIterator() {}

                Edge operator*() const {
                    // get iterator to the start of the set of neighbouring node uids
                    auto it = graph_->node_to_edges_[central_node_uid_].begin();
                    // iterate the number of spaces based on index
                    for (size_type j=0; j<indx_; j++)
                        it++;
                    // Look up that edge. The iterator over the map is a uid
                    Edge e = graph_->edges_[*it];
                    // check that the nodes are in the right ordering, if not then swap
                    if (e.node1().uid_ != central_node_uid_)
                        e.swap_edge_nodes_inplace_();
                    return e;
                }
                IncidentIterator& operator++() {
                    indx_++;
                    return *this;
                }
                bool operator==(const IncidentIterator& ie) const {
                    return  (ie.graph_==this->graph_
                        && ie.central_node_uid_==this->central_node_uid_
                        && ie.indx_==this->indx_);
                }

            private:
                friend class Graph;
                graph_type* graph_;
                size_type central_node_uid_;
                size_type indx_;

                // Constructor
                IncidentIterator(graph_type* graph, size_type central_node_uid, size_type indx)
                    : graph_(graph), central_node_uid_(central_node_uid), indx_(indx) {};
        };

        //
        // Edge Iterator
        //

        /** @class Graph::EdgeIterator
         * @brief Iterator class for edges. A forward iterator. */
        class EdgeIterator : totally_ordered<EdgeIterator> {
            public:
                // These type definitions let us use STL's iterator_traits.
                using value_type        = Edge;                     // Element type
                using pointer           = Edge*;                    // Pointers to elements
                using reference         = Edge&;                    // Reference to elements
                using difference_type   = std::ptrdiff_t;           // Signed difference
                //using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

                /** Construct an invalid EdgeIterator. */
                EdgeIterator() {}

                Edge operator*() const {
                    size_type uid = graph_->edge_i2u_[indx_];
                    return graph_->edges_[uid];
                }
                EdgeIterator& operator++() {
                    indx_++;
                    return *this;
                }
                bool operator==(const EdgeIterator& ei) const {
                    return (ei.graph_==this->graph_
                            && ei.indx_==this->indx_);
                }

            private:
                friend class Graph;
                const graph_type* graph_;
                size_type indx_;

                // Constructor
                EdgeIterator(const graph_type* graph, size_type indx)
                    : graph_(graph), indx_(indx) {};
        };

        edge_iterator edge_begin() const {
            EdgeIterator ei {this, 0};
            return ei;
        }
        edge_iterator edge_end() const {
            EdgeIterator ei {this, num_edges()};
            return ei;
        }

    private:
        std::vector<Node> nodes_;
        std::vector<Edge> edges_;
        std::vector<NodeInfo> node_infos_;
        std::vector<EdgeInfo> edge_infos_;
    public:
        std::vector<size_type> node_i2u_;
        std::vector<size_type> edge_i2u_;

        // map of nodes and the edge number they're in. And to nodes they have an edge with
        std::map<size_type, std::set<size_type>> node_to_edges_;
        std::map<size_type, std::set<size_type>> node_to_nodes_;
};

#endif // CME212_GRAPH_HPP

//--style_-1
//--Really clean Graph methods
//--END

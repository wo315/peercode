/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}


/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point f_grav = n.value().mass * Point(0, 0, -grav);
    Point f_spr = Point(0, 0, 0);
    for(auto incid_itr = n.edge_begin(); incid_itr != n.edge_end(); ++incid_itr){
      auto edge_curr = *incid_itr;
      Point edge_vec = Point(0, 0, 0);
      if(edge_curr.node1() == n)
        edge_vec = n.position() - edge_curr.node2().position();
      else
        edge_vec = n.position() - edge_curr.node1().position();
      f_spr -= edge_curr.value().K*edge_vec*(norm(edge_vec) - edge_curr.value().length)/norm(edge_vec);
    }
    (void) t;
    return f_grav + f_spr;
  }
};

struct Problem2Force {
  /** Return the spring force applying to @a n at time @a t.
   *
   * For HW2 #2, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point Fspring = Point(0, 0, 0);
    Point xi = n.position();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      double K = (*i).value().K;
      double L = (*i).value().L;
      //std::cout << K << "  " << L << std::endl;
      Point xj = (*i).node2().position();
      double len = (*i).length();
      Fspring += (-K * (len - L) / len) * (xi - xj);
    }
    return Fspring + (n.value().mass * Point(0, 0, -grav));
  }
};


/** Force function object for HW2 #3, only the gravity force. */
struct GravityForce {
  /** Return the gravity force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;     // silence compiler warnings
    return  (n.value().mass * Point(0, 0, -grav));
  }
};

/** Force function object for HW2 #3, only the spring force. */
struct MassSpringForce {
  /** Return the mass spring force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    Point Fspring = Point(0, 0, 0);
    Point xi = n.position();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i) {
      double K = (*i).value().K;
      double L = (*i).value().L;
      //std::cout << K << "  " << L << std::endl;
      Point xj = (*i).node2().position();
      double len = (*i).length();
      Fspring += (-K * (len - L) / len) * (xi - xj);
    }
    return Fspring;
  }
};

/** Force function object for HW2 #3, only the damping force. */
struct DampingForce {
  double c_;
  DampingForce(const double c) : c_(c) {
  };
  /** Return the damping force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return -c_ * n.value().vel;
  }
};

/** Zero Force, used in make_combined when only two parameters are passed. */
struct ZeroForce {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) n;
    (void) t;
    return Point(0, 0, 0);
  }
};

template <typename f1, typename f2, typename f3>
struct make_combined_force {
//--design_1
//--These should be in a vector
//--START
  f1 f1_;
  f2 f2_;
  f3 f3_;
//--END
  /** Combine 3 forces. */
  make_combined_force(f1 force1, f2 force2, f3 force3):
    f1_(force1), f2_(force2), f3_(force3){};
  /** Combine 2 forces. */
  make_combined_force(f1 force1, f2 force2):
    f1_(force1), f2_(force2), f3_(ZeroForce()){};
  /** Return the combined force. */
  template<typename NODE>
  Point operator()(NODE n, double t) {
    return f1_(n, t) + f2_(n, t) + f3_(n, t);
  }
};

/** The constraint of a plane. */
struct PlaneConstraint {
  double l = -0.75;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(n.position().z < l){
        n.position().z = l;
        n.value().vel.z = 0;
      }
    }
    return;
  }
};

/** The constraint of sphere. */
struct SphereConstraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      auto x = n.position();
      auto R = (x - c)/norm(x - c);
      if(norm(x - c) < r){
        n.value().vel -= (n.value().vel * R) * R;
        n.position() = c + r * R;
      }
    }
    return;
  }
};

struct modifyvel{
  Point c_;
  double r2_;
  Node n_;
  modifyvel(Point c, double r2, Node n1) : c_(c), r2_(r2), n_(n1) {}
  void operator()(Node n){
    Point r = c_ - n.position();
    double l2 = normSq(r);
    if(n != n_ and l2 < r2_){
      n_.value().vel -= (dot(r, n_.value().vel) / l2) * r;
    }
  }
};

//--documentation_3
//--END

//--style_3
//--END


struct checkcollision {
  SpaceSearcher<Node> searcher_;
  checkcollision(const SpaceSearcher<Node>& searcher) : searcher_(searcher) {}
  void operator()(Node n){
    Point center = n.position();
    double radius2 = std::numeric_limits<double>::max();
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i){
      auto e = *i;
      radius2 = std::min(radius2, normSq(e.node2().position() - center));
    }
    radius2 *= 0.9;

    // Form a bounding box that encapsulates the constraints influence
    // Add some relaxing space for the box via multiplying the radius by 2
    double radius = std::sqrt(radius2);
    Point p1 = center - 2 * Point(radius, radius, radius);
    Point p2 = center + 2 * Point(radius, radius, radius);
    Box3D bb(p1, p2);
    // Using SpaceSeacher to iterate in the given box.
    // Note that the NeighborhoodIterator is not a random access iterator so we don't use parallel methods
    thrust::for_each(searcher_.begin(bb), searcher_.end(bb), modifyvel(center, radius2, n));
  }
};

struct SelfCollisionConstraint {
  SpaceSearcher<Node> searcher_;
  SelfCollisionConstraint(SpaceSearcher<Node> searcher) : searcher_(searcher) {}
  void operator()(GraphType& g){
    // Implement the first for-loop using thrust::for_each
    // Iterate through all nodes using NodeIterator
    thrust::for_each(thrust::omp::par, g.node_begin(), g.node_end(), checkcollision(searcher_));
  }
};


struct SelfCollisionTest {
  void operator()(GraphType& g) const {
    for(auto i = g.node_begin(); i != g.node_end(); ++i){
      auto n = *i;
      const Point & center = n.position();
      double radius2 = std::numeric_limits<double>::max();
      for(auto j = n.edge_begin(); j != n.edge_end(); ++j){
        auto e = *j;
        radius2 = std::min(radius2, normSq(e.node2().position() - center));
      }
      radius2 *= 0.9;
      for(auto k = g.node_begin(); k != g.node_end(); ++k){
        auto n2 = *k;
        Point r = center - n2.position();
        double l2 = normSq(r);
        if (n != n2 && l2 < radius2) {
          // Remove our velocity component in r
          n.value().vel -= (dot(r, n.value().vel) / l2) * r ;
        }
      }
    }
  }
};

/** The functor to remove nodes which violate the constraint. */
struct SphereRemove {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  template<typename GRAPH>
  void operator()(GRAPH& g){
    auto it = g.node_begin();
    while(it != g.node_end()){
      auto n = *it;
      auto x = n.position();
      if(norm(x - c) < r){
        it = g.remove_node(it);
      }
      else{
        ++it;
      }
    }
    return;
  }
};





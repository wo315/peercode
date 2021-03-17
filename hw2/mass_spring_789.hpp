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

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Edge spring constant
  double L;     //< Edge length
  EdgeData() : K(100.0), L(0.01) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
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
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)){
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)){
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}


template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constr) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }  
  
  // call constraint functor  
  constr(g);

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

  //double K_;
  //double L_;
  // Constructor to initialize K_ and L_
  //Problem1Force(const double K, const double L) : K_(K), L_(L) {}

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;  // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    } 
    else {
      Point x_i = n.position();
      Point f_spring_i = Point(0,0,0);
      for(auto init = n.edge_begin(); init != n.edge_end(); ++init) {
        Point x_j = (*init).node2().position();
        double K_ = (*init).value().K;
        double L_ = (*init).value().L;
        f_spring_i += ( -K_ * (x_i - x_j) / norm(x_i - x_j) ) * (
                                                           norm(x_i-x_j) - L_);
      }
      Point f_grav_i = n.value().mass * Point(0,0,-grav);
      return f_spring_i + f_grav_i;
    }
  }
};

/** Force class that implements operator() as a virtual function. */
struct Force {
  virtual Point operator()(Node n, double t) {
    (void) n; // silence compiler warnings
    (void) t; // silence compiler warnings
    return Point(0,0,0);
  }
};

/** GravityForce class that implements the force of gravity. */
struct GravityForce : Force {
  Point operator()(Node n, double t){
    (void) t; // silence compiler warnings
    Point f_grav = n.value().mass * Point(0,0,-grav);
    return f_grav;
  }
};

/** MassSpringForce class that implements the spring forces. */
struct MassSpringForce : Force {
  Point operator()(Node n, double t){
    (void) t; // silence compiler warnings
    Point x_i = n.position();
    Point f_spring = Point(0,0,0);
    // traverse adjacency list of node n using incident iterator & add forces
    for(auto init = n.edge_begin(); init != n.edge_end(); ++init) {
      Point x_j = (*init).node2().position();
      double K_ = (*init).value().K;
      double L_ = (*init).value().L;
      f_spring += ( -K_ * (x_i - x_j) / norm(x_i - x_j) ) * (
                                                           norm(x_i-x_j) - L_);
    }
    return f_spring;
  }
};

/** DampingForce class that implements damping. */
struct DampingForce : Force {
  double c_;
  // Constructor to initialize c_
  DampingForce(const double c = 0.01) : c_(c) {}
  // operator()
  Point operator()(Node n, double t){
    (void) t; // silence compiler warnings
    auto f_damp = -c_*n.value().vel;
    return f_damp;
  }
};

/** Functor to apply combined force. */
struct CombinedForce {
  std::vector<Force*> forces_;
  CombinedForce(std::vector<Force*> forces) : forces_(forces) {}
  Point operator()(Node n, double t){
    (void) t; // silence compiler warnings
    Point f_comb_;
    for(unsigned i = 0; i < forces_.size(); i++){
      f_comb_ += (*(forces_[i]))(n, t);
    }
    return f_comb_;
  }
};

/** Function to create CombinedForce (2) */
template <typename FORCE1, typename FORCE2>
CombinedForce make_combined_force(FORCE1 f1, FORCE2 f2) {
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  CombinedForce fcomb(forces);
  return fcomb;
}

/** Function to create CombinedForce (3) */
template <typename FORCE1, typename FORCE2, typename FORCE3>
CombinedForce make_combined_force(FORCE1 f1, FORCE2 f2, FORCE3 f3) {
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3);
  CombinedForce fcomb(forces);
  return fcomb;
}

/** Constraint class that implements operator() as a virtual function. */
struct Constraint {
  virtual void operator()(GraphType& g) {
    (void) g; // silence compiler warnings
  }
};


/** Pin contraint functor. */
struct PinConstraint : Constraint {
  Node n0;
  Node n1;
  // Constructor to find the special nodes
  PinConstraint(GraphType& g) {
    for(auto nitr = g.node_begin(); nitr != g.node_end(); ++nitr) {
      auto n = *nitr;
      if(n.position() == Point(0,0,0)) {n0 = n;}
      if(n.position() == Point(1,0,0)) {n1 = n;}
    }
  }
  void operator() (GraphType& g){
    (void) g; // silence compiler warnings
    n0.position() = Point(0,0,0);
    n0.value().vel = Point(0,0,0);
    n1.position() = Point(1,0,0);
    n1.value().vel = Point(0,0,0);    
  }
};

/** Plane contraint functor. */
struct PlaneConstraint : Constraint {
  double z_limit = -0.75;
  void operator() (GraphType& g){
    for(auto nitr = g.node_begin(); nitr != g.node_end(); ++nitr) {
      auto n = *nitr;
      if(dot(n.position(), Point(0,0,1)) < z_limit) {
        n.position().z = z_limit;
        n.value().vel.z = 0.0; 
      }
    }
  }
};

/** Sphere contraint functor. */ 
struct SphereConstraint : Constraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator() (GraphType& g){
    for(auto nitr = g.node_begin(); nitr != g.node_end(); ++nitr) {
      auto n = *nitr;
      auto x_i = n.position();
      auto norm_val_i = norm(x_i - c);
      if(norm_val_i < r) {
        auto R_i = (x_i - c) / norm_val_i;
        n.position() = c + (r * R_i);
        n.value().vel -= dot(n.value().vel, R_i) * R_i; 
      }
    }
  }
};

/** Functor to apply combined constraints. */
struct CombinedConstraints {
  std::vector<Constraint*> contraints_;
  CombinedConstraints(std::vector<Constraint*> contraints) :
                                                     contraints_(contraints) {}
  void operator()(GraphType& g){
    for(unsigned i = 0; i < contraints_.size(); i++){
      (*(contraints_[i]))(g);
    }
  }
};


/** Tear constraint functor*/
struct TearConstraint : Constraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator() (GraphType& g){
    auto nitr = g.node_begin();
    while(nitr != g.node_end()) {
      auto n = *nitr;
      auto x_i = n.position();
      auto norm_val_i = norm(x_i - c);
      if(norm_val_i < r) {
        nitr = g.remove_node(nitr); 
      } else {++nitr;}
    }
  }
};

//--functionality_1
//--The figure is wrong if TearConstraint is applied.
//--END

//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


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

// Default spring constant K
static constexpr double K = 100;

// Default rest length L
static constexpr double L = 0.1;

// Default damping constant
static constexpr double C = 0.2;

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

  // Compute the t+dt velocimty
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) 
      continue;
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
};

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and constraint.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Constraint object defining the constraint per graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is reference to the graph and @a t is the current time.
 *           @a constraint fixes node position and velocity.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint){
// Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraint
  constraint(g,t);

  // Compute the t+dt velocimty
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
};

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    // Constrain the two corners
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) 
      return Point(0,0,0);

    // first add gravity force 
    Point f_spring = n.value().mass * Point(0,0,-grav);
    // Iterate through neighbors and add the corresponding spring force
    Point x = n.position();
    for (auto i = n.edge_begin(); i!= n.edge_end();++i){
      double Li = (*i).value();
      Point xi = (*i).node2().position();
      f_spring += -K * (x-xi)/norm(x-xi) *(norm(x-xi)-Li);
    }
    return f_spring;
  }
};

/** Base force function object to be inhertied */
struct Force {
  /** Return the force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    (void) n;
    (void) t;
    return Point(0);
  }
};

/** Gravity force function object */
struct GravityForce: Force {
  /** Return the gravity force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }
};

/** Spring force function object */
struct MassSpringForce: Force {
  /** Return the spring force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    (void) t;

    Point f_spring = Point(0);
    // Iterate through neighbors and add the corresponding spring force
    Point x = n.position();
    for (auto i = n.edge_begin(); i!= n.edge_end();++i){
      double Li = (*i).value();
      Point xi = (*i).node2().position();
      f_spring += -K * (x-xi)/norm(x-xi) *(norm(x-xi)-Li);
    }
    return f_spring;
  }
};

/** Damping force function object */
struct DampingForce: Force {
  /** Return the damping force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    (void) t;
    return -C * n.value().vel;
  }
};

/** Combined force functor (3 forces)
 * @tparam F1 is a function object that is inherited from Force. It is called as @a f1(n, @a t)
 * @tparam F2 is a function object that is inherited from Force. It is called as @a f2(n, @a t)
 * @tparam F3 is a function object that is inherited from Force. It is called as @a f3(n, @a t)
 * 
 * @return Sum of the two forces.
*/
template<typename F1,typename F2,typename F3>
struct CombinedForce: Force {
//--design_1
//--Should be stored in vector
//--START
  F1 f1_;
  F2 f2_;
  F3 f3_;
//--END
  /** Constructor for combining three forces. */
  CombinedForce(F1 f1, F2 f2, F3 f3): f1_(f1), f2_(f2),f3_(f3){};

  /** Return the combined force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    return f1_(n,t) + f2_(n,t) + f3_(n,t);
  }
};

/** Combine two or three forces
 * @param[in] f1 one type of force to be combined
 * @param[in] f2 one type of force to be combined
 * @param[in] f3 one type of force to be combined, default is zero force.
 * @tparam F1 is a function object that is inherited from Force. It is called as @a f1(n, @a t)
 * @tparam F2 is a function object that is inherited from Force. It is called as @a f2(n, @a t)
 * @tparam F3 is a function object that is inherited from Force. It is called as @a f3(n, @a t)
 * 
 * @return Sum of the three forces.
 */
template<typename F1,typename F2,typename F3 = Force>
CombinedForce<F1,F2,F3> make_combined_force(F1 f1, F2 f2,F3 f3 = Force()){
  return CombinedForce<F1,F2,F3>(f1,f2,f3);
}


/** Base constraint function object */
struct Constraint {
  /** Fix node position. Do nothing in base. */
  virtual void operator()(GraphType &g, double t) {
    (void) g;
    (void) t;
    return;
  }
};

/** Pin nodes at (0,0,0) and (1,0,0) */
struct PinConstraint: Constraint {
  /** Fix node position based on constraint */
  virtual void operator()(GraphType &g, double t) {
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
//--functionality_1
//--This doesn't work, and causes your nodes to be unpinned
//--START
      Point prev_position = n.position() - n.value().vel * 0.001;
//--END
      if (prev_position == Point(0,0,0) || prev_position == Point(1,0,0)){
        n.position() = prev_position;
      }
    }
  }
};

/** Plane constraint */
struct PlaneConstraint: Constraint {
  double z = -0.75;
  /** Fix node position based on constraint */
  virtual void operator()(GraphType &g, double t) {
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if (n.position().z < z){
        n.position().z = z;
        n.value().vel.z = 0;
      }
    }
  }
};

/** Sphere constraint */
struct SphereConstraint: Constraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  /** Fix node position based on constraint */
  virtual void operator()(GraphType &g, double t) {
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      Point x = n.position();
      if (norm(x - c) < r){
        Point R = (x-c)/norm(x-c);
        n.position() = c + R * r;
        Point v = n.value().vel;
        n.value().vel -= (v*R)*R;
      }
    }
  }
};

/** Sphere constraint. Remove if vioaltes constraint */
struct SphereConstraintTear: Constraint {
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  /** Fix node position based on constraint */
  virtual void operator()(GraphType &g, double t) {
    (void) t;
    auto it = g.node_begin();
    while(it != g.node_end()){
      auto n = *it;
      Point x = n.position();
      if (norm(x - c) < r){
        it = g.remove_node(it);
      }else{
        ++it;
      }
    }
  }
};

/** Combined constraint functor (3 constraints)
 * @tparam C1 is a function object that is inherited from Constraint. It is called as @a c1(n, @a t)
 * @tparam C2 is a function object that is inherited from Constraint. It is called as @a c2(n, @a t)
 * @tparam C3 is a function object that is inherited from Constraint. It is called as @a c3(n, @a t)
 * 
 * @return Apply the constraints sequentially.
*/
template<typename C1,typename C2,typename C3>
struct CombinedConstraints: Force {
  C1 c1_;
  C2 c2_;
  C3 c3_;

  /** Constructor for combining three constraints. */
  CombinedConstraints(C1 c1, C2 c2, C3 c3): c1_(c1), c2_(c2),c3_(c3){};

  /** Return the combined constraint applying to @a n at time @a t. */
  virtual void operator()(GraphType &g, double t) {
    c1_(g,t);
    c2_(g,t);
    c3_(g,t);
  }
};

/** Combine two or three three constraints
 * @param[in] c1 one type of constraint to be combined
 * @param[in] c2 one type of constraint to be combined
 * @param[in] c3 one type of constraint to be combined. Default is no constraint
 * @tparam C1 is a function object that is called as @a c1(g, @a t)
 * @tparam C2 is a function object that is is called as @a c2(g, @a t)
 * @tparam C3 is a function object that is is called as @a c3(g, @a t)
 * 
 * @return Sequential application of constraints.
 */
template<typename C1,typename C2,typename C3 = Constraint>
CombinedConstraints<C1,C2,C3> make_combined_constraints(C1 c1, C2 c2,C3 c3 = Constraint()){
  return CombinedConstraints<C1,C2,C3>(c1,c2,c3);
}

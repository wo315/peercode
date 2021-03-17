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
  Point position_t; //< Node's position at time t
  double t; 
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K;       //< Edge spring constant
  double L;     //< Edge rest length
  EdgeData() : K(100), L(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
// using GraphType = Graph<NodeData, double>;
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

/** Step of symplectic Euler without constraints */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // hard-coded version of PinConstraint()
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      n.value().vel = Point(0, 0, 0);
    }

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

/** Step of symplectic Euler with constraints */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.value().position_t = n.position();
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

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
    (void) t;
    if (n.position()== Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    // Compute mass spring force
    Point f_spr = Point(0,0,0);
    for (auto inc_itr = n.edge_begin(); inc_itr != n.edge_end(); ++inc_itr){
      Point edge_vec = Point(0,0,0);
      if ((*inc_itr).node1() == n){
        edge_vec = n.position() - (*inc_itr).node2().position();
      }
      else{
        edge_vec = n.position() - (*inc_itr).node1().position();
      }
      f_spr += -(*inc_itr).value().K * edge_vec/norm(edge_vec) * (norm(edge_vec) - (*inc_itr).value().L);

    }
    // Compute gravitational force 
    Point f_grav = n.value().mass * Point(0,0, -grav);
    return f_spr + f_grav; 
  }
};

/** Force class that implements the operator() method. */
struct Force{
  virtual Point operator()(Node n, double t) const{
    (void) t;
    (void) n;
    return  Point(0,0,0); // return zero force
  }
};   

/** Class that inherits from Force class and implements the force of gravity. */
struct GravityForce : Force  {
  Point operator()(Node n, double t) const override{
    (void) t;
    // Compute the gravity force 
    return  n.value().mass * Point(0,0, -grav);
  }
};

/** Class that inherits from Force class and implements the spring forces. */
struct MassSpringForce  : Force {
  Point operator()(Node n, double t)  const override{
    (void) t;
    Point f_spr = Point(0,0,0);
    for(auto inc_itr = n.edge_begin(); inc_itr != n.edge_end(); ++inc_itr){
      Point edge_vec = Point(0, 0, 0);
      if((*inc_itr).node1() == n){
        edge_vec = n.position() - (*inc_itr).node2().position();
      }
      else{
        edge_vec = n.position() - (*inc_itr).node1().position();
      }
      // Compute the mass spring force 
      f_spr += -(*inc_itr).value().K * edge_vec/norm(edge_vec) * (norm(edge_vec) - (*inc_itr).value().L);
    }
    return f_spr;
  }
};

/** Class that inherits from Force class and implements the damping force. */
struct DampingForce : Force  {
  double damp_coeff_;
  // Initialize constructor with default damp coefficient of 0.2 if no coefficient is passed into the constructor
  DampingForce(double damp_coeff = 0.2): Force(), damp_coeff_(damp_coeff){
  };
  
  Point operator()(Node n, double t) const override{
    (void) t;
    // Compute the damping force 
    return  -damp_coeff_ * n.value().vel;
  }
};

/** Functor that combines 2 or 3 forces. */
struct CombinedForce{
  // TODO change to vectors?!
  const Force& force1_;
  const Force& force2_;
  const Force& force3_;
  
  // Combine 2 or 3 forces. By default force3 is Force(), which returns zero force, when no force3 is passed into functor. 
  CombinedForce(const Force& force1,  const Force& force2, const Force& force3 ):
    force1_(force1), force2_(force2), force3_(force3){};

  // Return the combined forces. 
  Point operator()(Node n, double t) {
    return force1_(n, t) + force2_(n, t) + force3_(n,t);
  }
};

/** Function that combines 2 or 3 forces in arbitrary order. */
CombinedForce make_combined_force(const Force& force1,  const Force& force2, const Force& force3= Force()){
  return CombinedForce(force1, force2, force3);
}

/** Constraint class that implements the constraint method. */
struct Constraint{
  virtual void operator()(GraphType& g, double t) const{
    (void) g;
    (void) t;
    return; // return zero force
  }
};   

/** Class that inherits from Constraint class and implements the constraint method for pin constraint. */
struct PinConstraint: Constraint{
  void operator()(GraphType& g, double t) const override{
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      // Keep nodes at (0,0,0) and (1,0,0) fixed; 
      if (n.value().position_t == Point(0, 0, 0) or n.value().position_t == Point(1, 0, 0)){
        n.value().vel = Point(0, 0, 0); // revert velocity of node at time t+1 to velocity at time t; i.e. velocity of Zero
        n.position() = n.value().position_t; // revert position of node at time t+1 to position at time t
      }
    }

    return;
  }
};   

/** Class that inherits from Constraint class and implements the constraint method for plane constraint. */
struct PlaneConstraint: Constraint{
  void operator()(GraphType& g, double t) const override{
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      //TODO: check the position reset
      if(n.position().z < 0.75){
        n.value().vel.z = 0;
        n.position().z = -0.75;
      }
    }
    return;
  }
};

/** Class that inherits from Constraint class and implements the constraint method for sphere constraint. */
struct SphereConstraint: Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& g, double t) const override{
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      auto x_i = n.position();
      if(norm(x_i - c) < r && norm(x_i - c) != 0){
        auto R_i = (x_i - c)/norm(x_i - c);
        n.value().vel -= dot(n.value().vel,R_i) * R_i;
        n.position() = c + r * R_i;
      }
      else if(norm(x_i - c) == 0){
        n.value().vel.x = 0;
        n.position() = c + Point(r,0,0);
      }
    }
    return;
  }
};

/** Class that inherits from Constraint class and implements the constraint method for sphere constraint with "tearing".
 * I.e. for each node that violates the Sphere constraint, it removes the node and all its edges. */
struct SphereTearConstraint: Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& g, double t) const override{
    (void) t;
    auto it = g.node_begin();
    while (it != g.node_end()){
      auto x_i = (*it).position();
//--functionality_1
//--The figure cannot move. There should be "else{++it;}" after the if block. 
//--START
      if(norm(x_i - c) < r){
        g.remove_node(*it);
      }
//--END
    }
    return;
  }
};

/** Functor that combines 2 or 3 constraints. */
struct CombinedConstraint{
  const Constraint& constraint1_;
  const Constraint& constraint2_;
  const Constraint& constraint3_;
  
  // Combine 2 or 3 constraints. By default constraint3 is Constraint(), which doesn't apply any constraint, when no constraint3 arg is passed into functor. 
  CombinedConstraint(const Constraint& constraint1, const Constraint& constraint2, const Constraint& constraint3):
    constraint1_(constraint1), constraint2_(constraint2), constraint3_(constraint3){};

  // Return the combined constraints.
  void operator()(GraphType& g, double t) {
    constraint1_(g,t);
    constraint2_(g,t);
    constraint3_(g,t);
  }

};

/** Function that combines 2 or 3 constraints in arbitrary order. */
CombinedConstraint make_combined_constraint(const Constraint& constraint1, const Constraint& constraint2, const Constraint& constraint3=Constraint()){
  return CombinedConstraint(constraint1, constraint2, constraint3);
}

//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


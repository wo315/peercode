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
  Point pos;       //< Initial node position
  NodeData() : vel(0), mass(1), pos(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double k;           //< Spring constant
  double length;      //< Edge length
  EdgeData() : k(100.0), length(1.0) {}
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

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and constraint.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraint
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
    (void) t; // silence warnings
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point gravity_force = n.value().mass * Point(0, 0, -grav);
    Point spring_force = Point(0, 0, 0);
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i)
    {
      auto edge_curr = *i;
      Point edge_vec = Point(0, 0, 0);
      if(edge_curr.node1() == n){
        edge_vec = n.position() - edge_curr.node2().position();
      }
      else {
        edge_vec = n.position() - edge_curr.node1().position();
      }
      spring_force -= edge_curr.value().k * edge_vec * (
        norm(edge_vec) - edge_curr.value().length) / norm(edge_vec);
    }
    return gravity_force + spring_force;
  }
};

/** Force function object parent class. */
struct Force {
  virtual Point operator()(Node n, double t) {
    (void) n; // silence warnings
    (void) t; // silence warnings
    return Point(0, 0, 0);
  }
};

/** Force function object for gravity force. */
struct GravityForce : Force {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    return Point(0, 0, -grav) * n.value().mass;
  }
};

/** Force function object for spring force. */
struct MassSpringForce : Force {
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    Point spring_force = Point(0, 0, 0);
    for(auto i = n.edge_begin(); i != n.edge_end(); ++i)
    {
      Point edge_vec = Point(0, 0, 0);
      auto edge_curr = *i;
      if(edge_curr.node1() == n){
        edge_vec = n.position() - edge_curr.node2().position();
      }
      else{
        edge_vec = n.position() - edge_curr.node1().position();
      }
      spring_force -= edge_curr.value().k * edge_vec * (
        norm(edge_vec) - edge_curr.value().length)/norm(edge_vec);
    }
    return spring_force;
  }
};

/** Force function object for damping force. */
struct DampingForce : Force {
  double damping_constant;
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; // silence warnings
    return -damping_constant * n.value().vel;
  }
};

template<typename force1_type, typename force2_type>
struct CombinedForce {
  force1_type force_1;
  force2_type force_2;
  /** Return the sum of two forces on node @a n at time @a t.
   * @param[in] n   node in graph we apply force to
   * @param[in] t   the current time
   * @return the sum of the two forces @a n at time @a t.
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    return force_1(n, t) + force_2(n, t);
  }
};

/** Compute combination of two forces.
 * @param[in] force_1   force of type force1_type
 * @param[in] force_2   force of type force2_type
 * @return a new force containing the sum of the two forces.
 */
template<typename force1_type, typename force2_type>
CombinedForce<force1_type, force2_type> make_combined_force(
  force1_type force_1, force2_type force_2){
  return {force_1, force_2};
};

/** Compute combination of three forces.
 * @param[in] force_1   force of type force1_type
 * @param[in] force_2   force of type force2_type
 * @param[in] force_3   force of type force3_type
 * @return a new force containing the sum of the three forces.
 */
template<typename force1_type, typename force2_type, typename force3_type>
CombinedForce<CombinedForce<force1_type, force2_type>, force3_type> make_combined_force(
  force1_type force_1, force2_type force_2, force3_type force_3){
  return make_combined_force(make_combined_force(force_1, force_2), force_3);
};


// CONSTRAINTS

/** Constraint function object parent class. */
struct Constraint {
  virtual void operator()(GraphType& g, double t) {
    (void) g; // silence warnings
    (void) t; // silence warnings
  }
};

/** Constraint function object for pin constraint. */
struct PinConstraint : Constraint {
  void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto curr_node = *i;
      if (curr_node.value().pos == Point(0, 0, 0)){
        curr_node.value().vel = Point(0, 0, 0);
        curr_node.position() = Point(0, 0, 0);
      } else if (curr_node.value().pos == Point(1, 0, 0)){
        curr_node.value().vel = Point(0, 0, 0);
        curr_node.position() = Point(1, 0, 0);
      }
    }
  }
};

/** Constraint function object for plane constraint. */
struct PlaneConstraint : Constraint {
  void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto curr_node = *i;
      if(curr_node.position().z < -0.75){
        curr_node.value().vel.z = 0;
        curr_node.position().z = -0.75;
      }
    }
  }
};

/** Constraint function object for sphere constraint. */
struct SphereConstraint : Constraint {
  void operator()(GraphType& g, double t){
    (void) t;
    double radius = 0.15;
    Point center = Point(0.5, 0.5, -0.5);
    for(auto i = g.node_begin(); i != g.node_end(); ++i){
      auto curr_node = *i;
      if(norm(curr_node.position() - center) < radius){
        curr_node.position() = center + (curr_node.position() - center) * radius / norm(curr_node.position() - center);
        Point R_i = (curr_node.position() - center)/norm(curr_node.position() - center);
        curr_node.value().vel -= dot(curr_node.value().vel, R_i) * R_i;
      }
      else if(norm(curr_node.position() - center) == 0){
        curr_node.position() = center + Point(0.15, 0, 0);
        curr_node.value().vel.x = 0;
      }
    }
  }
};

/** Constraint function object for sphere that tears constraint. */
struct TearConstraint : Constraint {
  void operator()(GraphType& g, double t){
    (void) t;
    double radius = 0.15;
    Point center = Point(0.5, 0.5, -0.5);
    for(auto i = g.node_begin(); i != g.node_end(); ++i){
      auto curr_node = *i;
      if(norm(curr_node.position() - center) < radius){
        g.remove_node(curr_node);
      }
    }
  }
};

template<typename constr1_type, typename constr2_type>
struct CombinedConstraints {
  constr1_type constr_1;
  constr2_type constr_2;
  /** Apply two constraints to graph g at a given time @a t.
   * @param[in] g   the graph we apply constraints to
   * @param[in] t   the current time
   */
  void operator()(GraphType& g, double t) {
    constr_1(g, t);
    constr_2(g, t);
  }
};

/** Return combination of two constraints.
 * @param[in] constr_1  constraint of type constr1_type
 * @param[in] constr_2  constraint of type constr2_type
 * @return CombinedConstraints that combines both the constraints.
 */
template<typename constr1_type, typename constr2_type>
CombinedConstraints<constr1_type, constr2_type> make_combined_constraints(
  constr1_type constr_1, constr2_type constr_2){
  return {constr_1, constr_2};
};


/** Return combination of three constraints.
 * @param[in] constr_1  constraint of type constr1_type
 * @param[in] constr_2  constraint of type constr2_type
 * @param[in] constr_3  constraint of type constr3_type
 * @return CombinedConstraints that combines all three constraints.
 */
template<typename constr1_type, typename constr2_type, typename constr3_type>
CombinedConstraints<CombinedConstraints<constr1_type, constr2_type>, constr3_type> make_combined_constraints(
  constr1_type constr_1, constr2_type constr_2, constr3_type constr_3){
  return make_combined_constraints(make_combined_constraints(constr_1, constr_2), constr_3);
};

//--functionality_0
//--Passed all tests!
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


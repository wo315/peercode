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

struct EdgeData {
  double K; 
  double length;
  EdgeData() : K(100.), length(1.) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using NodeIterator = typename GraphType::NodeIterator;


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
 *    method with the given node force.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint for certain nodes
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports certain node_value_type
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint will affect certain Node n at time @a t.
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
  // apply the constraint
  constraint(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** @class Force
 *  @brief A template for forces, a function object called as @a force(n, @a t)
 *  where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t. 
 */ 
class Force {
public:
  virtual Point operator()(Node n, double t) {
    (void) t;
    (void) n;
    return Point(0, 0, 0);
  }
};

/** Force function object for HW2 #1. */
// struct Problem1Force : Force {
struct Problem1Force : Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  // template <typename NODE>
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t; // silence any compiler warnings
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }
    Point f_ms = Point(0, 0, 0);
    Point x1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto edge = *it;
      Point x2 = edge.node2().position();
      f_ms -= (edge.value().K * (norm(x1 - x2) - edge.value().length) / norm(x1 - x2)) * (x1-x2);
    }
    return f_ms + Point(0, 0, -grav) * n.value().mass;
  }
};

/* Force function for the force of gravaity */
struct GravityForce : Force {
  // template <typename NODE>
  Point operator()(Node n, double t) {
    (void) t;
    return Point(0, 0, -grav) * n.value().mass;
  }
};

/* Force function for the spring forces */
struct MassSpringForce : Force {
  // template <typename NODE>
  Point operator()(Node n, double t) {
    (void) t;
    Point f_ms = Point(0, 0, 0);
    Point x1 = n.position();
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto edge = *it;
      Point x2 = edge.node2().position();
      double len = edge.length();
//--functionality_1
//--This should be negative
//--START
      f_ms -= (edge.value().K * (len - edge.value().length) / len ) * (x1-x2);
//--END
    }
    
    return f_ms;
  }    
};

/* Force function for the damping force */
struct DampingForce : Force {
  double damp_constant;
  // template <typename NODE>
  Point operator()(Node n, double t) {
    (void) t;
    return -damp_constant * n.value().vel;
  } 
};

/** Combined forces
 *  @tparam   force_type_1 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 *  @tparam   force_type_2 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 */
template<typename force_type_1, typename force_type_2>
struct CombinedForce : Force {
  force_type_1 force_1;
  force_type_2 force_2;
  // template <typename NODE>
  Point operator()(Node n, double t) {
    (void) t;
    return force_1(n, t) + force_2(n, t);
  }
  CombinedForce(force_type_1 f1, force_type_2 f2) : force_1(f1), force_2(f2) {}
};

/** Combine two forces and return the combined force
 *  @param[in]   force_1    force_type_1
 *  @param[in]   force_2    force_type_2
 *  @return the combined force made up of two input forces
 *
 *  @tparam   force_type_1 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 *  @tparam   force_type_2 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 */
template<typename force_type_1, typename force_type_2>
CombinedForce<force_type_1, force_type_2> \
              make_combined_force(force_type_1 force_1, force_type_2 force_2) {
  return {force_1, force_2};
};

/** Combine three forces and return the combined force
 *  @param[in]   force_1    force_type_1
 *  @param[in]   force_2    force_type_2
 *  @param[in]   force_3    force_type_3
 *  @return the combined force made up of three input forces
 *
 *  @tparam   force_type_1 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 *  @tparam   force_type_2 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 *  @tparam   force_type_3 is a function object called as @a force(n, @a t)
 *            where n is a node of the graph and @a t is the current time.
 *            @a force must return a Point representing the force vector on
 *            Node n at time @a t. 
 */
template<typename force_type_1, typename force_type_2, typename force_type_3>
CombinedForce<CombinedForce<force_type_1, force_type_2>, force_type_3> \
            make_combined_force(force_type_1 force_1, force_type_2 force_2, force_type_3 force_3) {
  return make_combined_force(make_combined_force(force_1, force_2), force_3);
};


/** @class Constraint
 *  @brief A template for contraints, a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint will affect certain Node n at time @a t.
 */ 
class Constraint {
public:
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    (void) g;
  }
};

/* constraint of a pin */
struct PinConstraint : Constraint {
  Node n1;
  Node n2;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if (*it == n1) {
        (*it).position() = Point(0, 0, 0);
      }
      if (*it == n2) {
        (*it).position() = Point(1, 0, 0);
      }
    }
  }
  PinConstraint(Node point1, Node point2) : n1(point1), n2(point2) {}
};

/* constraint of a plane */
struct PlaneConstraint : Constraint {
  double l = -0.75;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto node = *it;
      if (node.position().z < l) {
        node.position().z = l;
        node.value().vel.z = 0;
      }
    }
  }
};

/* constraint of a sphere */
struct SphereConstraint : Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto node = *it;
      auto x = node.position();
      auto R = (x - c) / norm(x - c);
      if (norm(x - c) < r) {
        node.position() = c + r * R;
        node.value().vel -= (node.value().vel * R) * R;
      }
    }
  } 
};

/* constraint of a sphere tear */
struct SphereTearConstraint : Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto node = *it;
      if ((norm(node.position() - c)) < r) {
        g.remove_node(node);
      }
    }
  } 
};

/** Combined constraints
 *  @tparam   constraint_type_1 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 *  @tparam   constraint_type_2 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 */
template<typename constraint_type_1, typename constraint_type_2>
struct CombinedConstraint : Constraint {
  constraint_type_1 constraint_1;
  constraint_type_2 constraint_2;

  void operator()(GraphType& g, double t) {
    (void) t;
    constraint_1(g, t);
    constraint_2(g, t);
  }
  CombinedConstraint(constraint_type_1 c1, constraint_type_2 c2) : constraint_1(c1), constraint_2(c2) {}
};


/** Combine two contraints and apply the combined constraint
 *  @param[in]   constraint_1    constraint_type_1
 *  @param[in]   constraint_2    constraint_type_2
 *  @return the combined constraint made up of two input constraints
 *
 *  @tparam   constraint_type_1 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 *  @tparam   constraint_type_2 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 * 
 */
template<typename constraint_type_1, typename constraint_type_2>
CombinedConstraint<constraint_type_1, constraint_type_2> \
              make_combined_constraint(constraint_type_1 constraint_1, constraint_type_2 constraint_2) {
  return {constraint_1, constraint_2};
};

/** Combine three contraints and apply the combined constraint
 *  @param[in]   constraint_1    constraint_type_1
 *  @param[in]   constraint_2    constraint_type_2
 *  @param[in]   constraint_3    constraint_type_3
 *  @return the combined constraint made up of three input constraints
 *
 *  @tparam   constraint_type_1 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 *  @tparam   constraint_type_2 is a function object called as @a constraint(g, @a t),
 *            where g is the graph and @a t is the current time.
 *            @a constraint will affect certain Node n at time @a t.
 * 
 */
template<typename constraint_type_1, typename constraint_type_2, typename constraint_type_3>
CombinedConstraint<CombinedConstraint<constraint_type_1, constraint_type_2>, constraint_type_3> \
            make_combined_constraint(constraint_type_1 constraint_1, constraint_type_2 constraint_2,\
             constraint_type_3 constraint_3) {
  return make_combined_constraint(make_combined_constraint(constraint_1, constraint_2), constraint_3);
};



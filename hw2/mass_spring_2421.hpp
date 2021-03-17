/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <vector>

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
  double k; // Spring constant
  double l; // Unstretched length
  EdgeData() : k(1), l(1) {}
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
 *    method with the given node force.
 * @param[in,out] g          Graph
 * @param[in]     t          The current time (useful for time-dependent forces)
 * @param[in]     dt         The time step
 * @param[in]     force      Function object defining the force per node
 * @param[in]     constraint Constraint on nodes
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
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
    // HW2 #1: YOUR CODE HERE
    (void) t; // silence compiler warnings
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    Point spring_forces(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      Point action;
      if (n == e.node1())
        action = (n.position() - e.node2().position())/e.length();
      else
        action = (n.position() - e.node1().position())/e.length();
      Point spring_force = -1*(e.value().k)*(e.length() - e.value().l)*action;
      spring_forces = spring_forces + spring_force;
    }
    Point gravity_force(0, 0, -1*grav*n.value().mass);
    return gravity_force + spring_forces;
  }
};

///////////////////// HW2.5 - FORCES

/** @class Force
 *
 * The superclass for all forces. All forces should overload operator()
 */
class Force {
public:
  virtual Point operator()(Node n, double t) {
    (void) n; (void) t; // silence compiler warnings
    return Point(0, 0, 0);
  }
};

/** @class GravityForce
 *
 * Force of gravity pulling on the node.
 */
class GravityForce : public Force {
public:
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warnings
    return Point(0, 0, -1*grav*n.value().mass);
  }
};

/** @class MassSpringForce
 *
 * Force of adjacent nodes pulling on each other via their edges, modeled as springs.
 */
class MassSpringForce : public Force {
public:
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warnings
    Point spring_forces(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      Point action;
      if (n == e.node1())
        action = (n.position() - e.node2().position())/e.length();
      else
        action = (n.position() - e.node1().position())/e.length();
      Point spring_force = -1*(e.value().k)*(e.length() - e.value().l)*action;
      spring_forces = spring_forces + spring_force;
    }
    return spring_forces;
  }
};

/** @class DampingForce
 *
 * Damping force on a node due to its velocity.
 */
class DampingForce : public Force {
public:
  DampingForce() : c(1) {}
  DampingForce(double dampingconstant) : c(dampingconstant) {}
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warnings
    return -c*n.value().vel;
  }
private:
  double c; // Damping constant
};

// Combined Force Functor
class CombinedForce {
public:
  CombinedForce(std::vector<Force*> fs) : my_forces(fs) {}
  Point operator()(Node n, double t) {
    Point resultant(0, 0, 0);
    for (auto it = my_forces.begin(); it != my_forces.end(); ++it) {
      resultant = resultant + (*(*it))(n, t);
    }
    return resultant;
  }
private:
  std::vector<Force*> my_forces;
};

// Make Combined Force functions
template <typename force1, typename force2>
CombinedForce make_combined_force(force1 f1, force2 f2) {
  std::vector<Force*> forces = {&f1, &f2};
  return forces;
}

// Make Combined Force functions
template <typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 f1, force2 f2, force3 f3) {
  std::vector<Force*> forces = {&f1, &f2, &f3};
  return forces;
}

///////////////////// HW2.6 - CONSTRAINTS

/** @class Constraint
 *
 * The superclass for all constraints. All constraints should overload operator()
 */
class Constraint {
public:
  virtual void operator()(GraphType& g, double t) {
    (void) g; (void) t; // silence compiler warnings
  }
};

/** @class PinConstraint
 *
 * Fixes nodes at (0, 0, 0) and (1, 0, 0) in place
 */
class PinConstraint : public Constraint {
public:
  PinConstraint(GraphType& g) {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position() == Point(0, 0, 0)) {
        p000 = *it;
      }
      if ((*it).position() == Point(1, 0, 0)) {
        p100 = *it;
      }
    }
  }
  void operator()(GraphType& g, double t) {
    (void) g; (void) t; // silence compiler warnings
    p000.position() = Point(0, 0, 0);
    p100.position() = Point(1, 0, 0);
  }
private:
  Node p000;
  Node p100;
};

/** @class PlaneConstraint
 *
 * Keeps all points above the plane z = -0.75
 */
class PlaneConstraint : public Constraint {
public:
  void operator()(GraphType& g, double t) {
    (void) t; // silence compiler warnings
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position()[2] < -0.75) {
        (*it).position()[2] = -0.75;
        Point new_vel = (*it).value().vel;
        new_vel[2] = 0;
        (*it).value().vel = new_vel;
      }
    }
  }
};

/** @class SphereConstraint
 *
 * Keeps points outside of sphere centered at (0.5, 0.5, -0.5)
*  with radius 0.15.
 */
class SphereConstraint : public Constraint {
public:
  void operator()(GraphType& g, double t) {
    (void) t; // silence compiler warnings
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (norm(n.position() - c) < r) {
        Point Ri = (n.position() - c)/norm(n.position() - c);
        Point new_pos = c + r*Ri;
        Point new_vel = n.value().vel - dot(n.value().vel, Ri)*Ri;
        n.position() = new_pos;
        n.value().vel = new_vel;
      }
    }
  }
  SphereConstraint() : c(Point(0.5, 0.5, -0.5)), r(0.15) {}
  SphereConstraint(Point center, double radius) : c(center), r(radius) {}
private:
  Point c; // Center
  double r; // Radius
};

/** @class DestructoSphere
 *
 * Destroys all nodes that enter the sphere centered at
*  (0.5, 0.5, -0.5) with radius 0.15.
 */
class DestructoSphere : public Constraint {
public:
  void operator()(GraphType& g, double t) {
    (void) t; // silence compiler warnings
    std::vector<Node> nodes_to_destroy;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (norm(n.position() - c) < r) {
        nodes_to_destroy.push_back(n);
      }
    }
    for (auto it = nodes_to_destroy.begin(); it != nodes_to_destroy.end(); ++it) {
      g.remove_node(*it);
    }
  }
  DestructoSphere() : c(Point(0.5, 0.5, -0.5)), r(0.15) {}
  DestructoSphere(Point center, double radius) : c(center), r(radius) {}
private:
  Point c; // Center
  double r; // Radius
};

// Combined Constraint Functor
class CombinedConstraint {
public:
  CombinedConstraint(std::vector<Constraint*> cs) : my_constraints(cs) {}
  void operator()(GraphType& g, double t) {
    for (auto it = my_constraints.begin(); it != my_constraints.end(); ++it) {
      (*(*it))(g, t);
    }
  }
private:
  std::vector<Constraint*> my_constraints;
};

// Make Combined Constraint functions
template <typename constraint1, typename constraint2>
CombinedConstraint make_combined_constraint(constraint1 c1, constraint2 c2) {
  std::vector<Constraint*> constraints = {&c1, &c2};
  return constraints;
}

// Make Combined Force functions
template <typename constraint1, typename constraint2, typename constraint3>
CombinedConstraint make_combined_constraint(constraint1& c1, constraint2& c2, constraint3& c3) {
  std::vector<Constraint*> constraints = {&c1, &c2, &c3};
  return constraints;
}

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


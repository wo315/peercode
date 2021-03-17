/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/Util.hpp"
#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double damping = 1.0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;    //< Node velocity
  double mass;  //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double k;  // Spring constant
  double l;  // Rest length
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

    if (n.position() == Point(0) || n.position() == Point(1, 0, 0)) {
      continue;
    }

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if (n.position() == Point(0) || n.position() == Point(1, 0, 0)) {
      continue;
    }

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force, and the given constraint.
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent
 * forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Constraint object defining the applied constraints
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

  // Reset violated constraints.
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
  /** Return the force applying to @a n at time @a t. [HW2]
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0) || n.position() == Point(1, 0, 0)) {
      return Point(0);
    }
    (void)n;
    (void)t;
    (void)grav;  // silence compiler warnings

    Point force(0, 0, -n.value().mass * grav);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      double l = e.length();
      Point dx = e.node1().position() - e.node2().position();
      force += -e.value().k * dx / l * (l - e.value().l);
    }

    return force;
  }
};

/** @class Force
 * @brief Base class for force functors.
 */
class Force {
 public:
  virtual Point operator()(Node n, double t) {
    (void)n;
    (void)t;  // silence compiler warnings
    return Point(0);
  }
};

/** Force due to gravity. */
class GravityForce : public Force {
 public:
  Point operator()(Node n, double t) {
    (void)t;  // silence compiler warnings
    return Point(0, 0, -n.value().mass * grav);
  }
};

/** Force due to spring stiffness. */
class MassSpringForce : public Force {
 public:
  Point operator()(Node n, double t) {
    (void)t;  // silence compiler warnings
    Point force(0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto e = *it;
      double l = e.length();
      Point dx = e.node1().position() - e.node2().position();
      force += -e.value().k * dx / l * (l - e.value().l);
    }
    return force;
  }
};

/** Force due to viscous damping. */
class DampingForce : public Force {
 public:
  Point operator()(Node n, double t) {
    (void)t;  // silence compiler warnings
    return -damping * n.value().vel;
  }
};

/** Force due to a combination of other forces.
 * @param[in] forces A vector of pointers to Forces.
 */
class CombinedForce : public Force {
 private:
  std::vector<Force*> forces_;

 public:
  CombinedForce(std::vector<Force*> forces) : forces_(forces) {}
  Point operator()(Node n, double t) {
    Point force(0);
    for (Force* f : forces_) {
      force += (*f)(n, t);
    }
    return force;
  }
};

/** Returns the combination of two forces. */
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2) {
  std::vector<Force*> forces = {&f1, &f2};
  return CombinedForce{forces};
}

/** Returns the combination of three forces. */
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3) {
  std::vector<Force*> forces = {&f1, &f2, &f3};
  return CombinedForce{forces};
}

/** @class Constraint
 * @brief Base class for kinematic constraints.
 */
class Constraint {
 public:
  virtual void operator()(GraphType& g, double t) {
    (void)g;
    (void)t;  // silence compiler warnings
  }
};

class PinConstraint : public Constraint {
 private:
  size_t index1_;
  size_t index2_;
  Point pos1_;
  Point pos2_;

 public:
  PinConstraint(Node node1, Node node2)
      : index1_(node1.index()),
        index2_(node2.index()),
        pos1_(node1.position()),
        pos2_(node2.position()) {}
  void operator()(GraphType& g, double t) {
    (void)t;  // silence compiler warnings
    g.node(index1_).position() = pos1_;
    g.node(index1_).value().vel = Point(0);
    g.node(index2_).position() = pos2_;
    g.node(index2_).value().vel = Point(0);
  }
};

class PlaneConstraint : public Constraint {
 private:
  double z_;

 public:
  PlaneConstraint(double z) : z_(z) {}
  void operator()(GraphType& g, double t) {
    (void)t;  // silence compiler warnings
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position().z < z_) {
        n.position().z = z_;
        n.value().vel.z = 0.0;
      }
    }
  }
};

class SphereConstraint : public Constraint {
 private:
  Point c_;
  double r_;

 public:
  SphereConstraint(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double t) {
    (void)t;  // silence compiler warnings
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      double d = norm(n.position() - c_);
      if (d < r_) {
        Point normal = (n.position() - c_) / d;
        n.position() = c_ + r_ * normal;
        n.value().vel -= dot(n.value().vel, normal) * normal;
      }
    }
  }
};

class SphereRemovalConstraint : public Constraint {
 private:
  Point c_;
  double r_;

 public:
  SphereRemovalConstraint(Point c, double r) : c_(c), r_(r) {}
  void operator()(GraphType& g, double t) {
    (void)t;  // silence compiler warnings
//--design_1
//--You are missing nodes by iterating through a for loop and incrementing every time
//--START 
   for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      double d = norm(n.position() - c_);
      if (d < r_) {
        g.remove_node(n);
        // Point normal = (n.position() - c_) / d;
        // n.position() = c_ + r_ * normal;
        // n.value().vel -= dot(n.value().vel, normal) * normal;
      }
    }
//--END
  }
};

class CombinedConstraint : public Constraint {
 private:
  std::vector<Constraint*> constraints_;

 public:
  CombinedConstraint(std::vector<Constraint*> constraints)
      : constraints_(constraints) {}
  void operator()(GraphType& g, double t) {
    for (Constraint* c : constraints_) {
      (*c)(g, t);
    }
  }
};

//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
/** Returns the combination of two constraints. */
template <typename F1, typename F2>
CombinedConstraint make_combined_constraint(F1 f1, F2 f2) {
  std::vector<Constraint*> constraints = {&f1, &f2};
  return CombinedConstraint{constraints};
}
//--END

/** Returns the combination of three constraints. */
template <typename F1, typename F2, typename F3>
CombinedConstraint make_combined_constraint(F1 f1, F2 f2, F3 f3) {
  std::vector<Constraint*> constraints = {&f1, &f2, &f3};
  return CombinedConstraint{constraints};
}

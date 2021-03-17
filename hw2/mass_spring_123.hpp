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
  Point start;     //< Node start position
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< Spring constant
  double length;  //< Edge length
  EdgeData() : K(100.0), length(1.0) {}
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
 * @tparam G::node_value_type supports NodeData
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
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Constraint object defining constraints of the graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step (G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraints to graph
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
    // Set constraints for corners to have 0 force
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    // Compute gravitational force
    Point force_gravity = n.value().mass * Point(0, 0, -grav);

    // Compute spring force
    Point force_spring = Point(0, 0, 0);

    // Iterate through all adjacent nodes to compute spring force
    for(auto incident_iter = n.edge_begin(); incident_iter != n.edge_end(); ++incident_iter) {
      auto edge = *incident_iter;
      Point edge_vec = n.position() - edge.node2().position();
      force_spring -= edge.value().K * edge_vec * (edge.length() - edge.value().length)/edge.length();
    }
    (void) t;
    return force_gravity + force_spring;
  }
};



// Parent Force class
class Force {
 public:
  Point force_ = Point(0, 0, 0);
  // Return 0 force.
  virtual Point operator()(Node n, double t) {
    return force_;
  }
};


// Spring Force class, child of Parent Force class
class MassSpringForce : public Force {
 public:
 
  // Return the spring force applying to @a n at time @a t.
  Point operator()(Node n, double t) {
    Point force_ = Point(0, 0, 0);
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    // Iterate through all adjacent nodes to compute spring force
    for(auto incident_iter = n.edge_begin(); incident_iter != n.edge_end(); ++incident_iter) {
      auto edge = *incident_iter;
      Point edge_vec = n.position() - edge.node2().position();
      force_ -= edge.value().K * edge_vec * (edge.length() - edge.value().length)/edge.length();
    }
    (void) t;
    return force_;
  }
};

// Gravitational Force class, child of Parent Force class
class GravitationalForce : public Force {
 public:
  // Return the gravitational force applying to @a n at time @a t.
  Point operator()(Node n, double t) {
    Point force_ = Point(0, 0, 0);
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    force_ = n.value().mass * Point(0, 0, -grav);
    (void) t;
    return force_;
  }
};

// Damping Force class, child of Parent Force class
class DampingForce : public Force {
 public:
  double damp_coeff;
  // Return the damping force applying to @a n at time @a t.
  Point operator()(Node n, double t) {
    Point force_ = Point(0, 0, 0);
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    force_ = -n.value().vel*damp_coeff;
    (void) t;
    return force_;
  }
};


// Combined Force class
class CombinedForce {
  std::vector<Force*> forces_;
 public:
  CombinedForce(std::vector<Force*> forces) : forces_(forces) {}

  // Return a combination of forces acting on node @a n at time @a t.
  Point operator()(Node n, double t) {
    Point total = Point(0, 0, 0);
    for(auto it = forces_.begin(); it != forces_.end(); ++it) {
      total = (*it)->operator()(n, t) + total;
    }
    return total;
  }
  
};


/** Combine two forces together.
 * @tparam F1 is a function object called as @a force( @a n, @a t)
 *         where @a n is a node of the graph and @a t is the current time.
 *         @a force must return a Point representing the force vector on
 *         Node n at time @a t.
 * @tparam F2 is a function object called as @a force( @a n, @a t)
 *         where @a n is a node of the graph and @a t is the current time.
 *         @a force must return a Point representing the force vector on
 *         Node @a n at time @a t.
 */
template<typename F1, typename F2>
 /** Return the sum of the two forces.
  * @param[in]     force1      force of type F1
  * @param[in]     force2      force of type F2
  * @return a new force containing the sum of the two forces.
  *
  */
//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
CombinedForce make_combined_force(F1 force1, F2 force2) {
  std::vector<Force*> forces_{};
  forces_.push_back(&force1);
  forces_.push_back(&force2);
  CombinedForce comb_force(forces_);
  return comb_force;
}
//--END

/** Combine three forces together.
 * @tparam F1 is a function object called as @a force( @a n, @a t)
 *         where @a n is a node of the graph and @a t is the current time.
 *         @a force must return a Point representing the force vector on
 *         Node n at time @a t.
 * @tparam F2 is a function object called as @a force( @a n, @a t)
 *         where @a n is a node of the graph and @a t is the current time.
 *         @a force must return a Point representing the force vector on
 *         Node @a n at time @a t.
 * @tparam F3 is a function object called as @a force( @a n, @a t)
 *         where @a n is a node of the graph and @a t is the current time.
 *         @a force must return a Point representing the force vector on
 *         Node @a n at time @a t.
 */
template<typename F1, typename F2, typename F3>
/** Return the sum of the three forces.
  * @param[in]     force1      force of type F1
  * @param[in]     force2      force of type F2
  * @param[in]     force3      force of type F3
  * @return a new force containing the sum of the two forces.
  *
  */
CombinedForce make_combined_force(F1 force1, F2 force2, F3 force3) {
  std::vector<Force*> forces_{};
  forces_.push_back(&force1);
  forces_.push_back(&force2);
  forces_.push_back(&force3);
  CombinedForce comb_force(forces_);
  return comb_force;
}




/* Parent Constraint class. */
class Constraint {
  public:
    virtual void operator()(GraphType& graph_, double t) = 0;
};


/** Constraint function object for plane constraint. */
class PinConstraint : public Constraint {

  /** Applying plane constraint to @a g at time @a t. */
  void operator()(GraphType& graph_, double t) {
    for (auto iter = graph_.node_begin(); iter != graph_.node_end(); ++iter) {
      if ((*iter).value().start == Point(0, 0, 0)) {
        (*iter).position() = Point(0, 0, 0);
      } else if ((*iter).value().start == Point(1, 0, 0)) {
        (*iter).position() = Point(1, 0, 0);
      }
    }

    (void) t;
  }
  
};


/** Constraint function object for plane constraint. */
class PlaneConstraint : public Constraint {
  /** Applying plane constraint to @a g at time @a t. */
  void operator()(GraphType& graph_, double t) {
    for(auto node_iter = graph_.node_begin(); node_iter != graph_.node_end(); ++node_iter) {
      auto current_node = *node_iter;

      if(current_node.position().z < -0.75) {
        current_node.value().vel.z = 0;
        current_node.position().z = -0.75;
      }
    }
    (void) t;
  }
  
};



/** Constraint function object for sphere constraint. */
class SphereConstraint : public Constraint {
  /** Applying sphere constraint to @a g at time @a t. */
  void operator()(GraphType& graph_, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    for (auto node_iter = graph_.node_begin(); node_iter != graph_.node_end(); ++node_iter) {
      auto current_node = *node_iter;
      if (norm(current_node.position() - center) < radius) {
        Point R_i = (current_node.position() - center)/norm(current_node.position() - center);
        current_node.value().vel -= dot(current_node.value().vel, R_i)*R_i;
        current_node.position() = center + (current_node.position() - center)*radius/norm(current_node.position() - center);
      }
    }

    (void) t;
  }
};


// Combined Constraints class
class CombinedConstraints {
  std::vector<Constraint*> constraints_;
 public:
  CombinedConstraints(std::vector<Constraint*>& constraints) : constraints_(constraints) {}

  // Apply a combination of constraints acting on node @a n at time @a t.
  void operator()(GraphType& graph_, double t) {
    for(auto it = constraints_.begin(); it != constraints_.end(); ++it) {
      (*it)->operator()(graph_, t);
    }
  }
  
};



/** Apply two constraints together.
 * @tparam C1 is a function object called as @a constraint( @a g, @a t)
 *         where @a t is the current time, @a g is the graph.
 *         @a constraint must apply the given constraint on
 *         all nodes at time @a t.
 * @tparam C2 is a function object called as @a constraint( @a g, @a t)
 *         where @a t is the current time, @a g is the graph.
 *         @a constraint must apply the given constraint on
 *         all nodes at time @a t.
 */
template<typename C1, typename C2>
  /** Apply the two constraints at time @a t.
   * @param[in]     c1   Constraint of type C1
   * @param[in]     c2   Constraint of type C2
   */
CombinedConstraints make_combined_constraint(C1 c1, C2 c2) {
  std::vector<Constraint*> constraints_{};
  constraints_.push_back(&c1);
  constraints_.push_back(&c2);
  CombinedConstraints comb_constraints(constraints_);
  return comb_constraints;
}


/** Apply three constraints together.
 * @tparam C1 is a function object called as @a constraint( @a g, @a t)
 *         where @a t is the current time, @a g is the graph.
 *         @a constraint must apply the given constraint on
 *         all nodes at time @a t.
 * @tparam C2 is a function object called as @a constraint( @a g, @a t)
 *         where @a t is the current time, @a g is the graph.
 *         @a constraint must apply the given constraint on
 *         all nodes at time @a t.
 * @tparam C3 is a function object called as @a constraint( @a g, @a t)
 *         where @a t is the current time, @a g is the graph.
 *         @a constraint must apply the given constraint on
 *         all nodes at time @a t.
 */
template<typename C1, typename C2, typename C3>
/** Apply the two constraints at time @a t.
   * @param[in]     c1   Constraint of type C1
   * @param[in]     c2   Constraint of type C2
   * @param[in]     c3   Constraint of type C3
   */
CombinedConstraints make_combined_constraint(C1 c1, C2 c2,C3 c3) {
  std::vector<Constraint*> constraints_{};
  constraints_.push_back(&c1);
  constraints_.push_back(&c2);
  constraints_.push_back(&c3);
  CombinedConstraints comb_constraints(constraints_);
  return comb_constraints;
}




/** Constraint function object for sphere constraint. */
class SphereConstraintTear : public Constraint {
  /** Applying sphere constraint to @a g at time @a t. */
  void operator()(GraphType& graph_, double t) {
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
//--design_1
//--You're skipping nodes by using a for loop (the one at the place you just deleted)
//--START
    for (auto node_iter = graph_.node_begin(); node_iter != graph_.node_end(); ++node_iter) {
      auto current_node = *node_iter;
      if (norm(current_node.position() - center) < radius) {
        graph_.remove_node(current_node);
      }
    }
//--END
    (void) t;
  }
};

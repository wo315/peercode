/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <fstream>
#include <functional>
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
  bool is_fixed;   //< Whether node is pinned in place
  Point orig_pos;  //< Original location for node

  NodeData() : vel(0), mass(1), is_fixed(false), orig_pos(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double spring_constant;  //< Edge spring constant
  double rest_length;  //< Edge rest length
  EdgeData() : spring_constant(100), rest_length(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Base force function object for use with mass-spring simulation **/
struct Force {
  /** Return the force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) const {
    (void) n;  // Silence compiler warnings
    (void) t;  // Silence compiler warnings

    return Point(0);
  }
};

/** Base constraint object for use with mass-spring simulation **/
struct Constraint {
  virtual void operator()(GraphType& graph, double t) const {
    (void) graph;  // Silence compiler warnings
    (void) t;  // Silence compiler warnings
  }
};

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
  return symp_euler_step(g, t, dt, force, Constraint());
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and node constraint.
 * @param[in,out] g          Graph
 * @param[in]     t          The current time (useful for time-dependent forces)
 * @param[in]     dt         The time step
 * @param[in]     force      Function object defining the force per node
 * @param[in]     constraint Function object defining the constraint per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData and G::edge_value_type supports EdgeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint modifies in-place nodes that violate the contraint.
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

  // Apply the constraints
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Gravity force function object for use with mass-spring simulation **/
struct GravityForce : public Force {
  /** Return the gravity force applying to @a n at time @a t. */
  Point operator()(Node n, double t) const {
    (void) t;  // Silence compiler warnings

    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Spring force function object for use with mass-spring simulation **/
struct MassSpringForce : public Force {
  /** Return the spring force applying to @a n at time @a t.
   * @pre node @a n has adjacent edges @a e_i s.t. @a e_i.length() != 0
  */
  Point operator()(Node n, double t) const {
    (void) t;  // Silence compiler warnings

    return std::accumulate(n.edge_begin(), n.edge_end(), Point(0),
      [&n](Point accum, const Edge& edge) -> Point {
        const EdgeData& edge_value = edge.value();

        Point diff = n.position() - edge.node2().position();
        double length = edge.length();

        Point edge_force = edge_value.spring_constant * diff
                            * (1 - edge_value.rest_length / length);

        return accum - edge_force;
      });
  }
};

/** Damping/friction force function object
 * for use with mass-spring simulation **/
struct DampingForce : public Force {
  double damping_const_;

  explicit DampingForce(double damping_const = 0.05)
    : damping_const_(damping_const) {
  }

  /** Return the damping force applying to @a n at time @a t. */
  Point operator()(Node n, double t) const {
    (void) t;  // Silence compiler warnings

    return - damping_const_ * n.value().vel;
  }
};

/**
 * @brief Force functor that combines a parameter pack of other force functors
 * 
 * @tparam TForces types of individual force objects
 */
template <typename ...TForces>
struct CombinedForce : public Force {
  using ForceFunction = std::function<Point(Node, double)>;
  std::vector<ForceFunction> forces;

  explicit CombinedForce(TForces... forces)
    : forces({ forces... }) {
  }

  /** Return the force applying to @a n at time @a t,
   * accumulated over the force functors. */
  Point operator()(Node n, double t) const {
    return std::accumulate(forces.begin(), forces.end(), Point(0),
      [&n, &t](Point accum, const ForceFunction& force_fn) {
        return accum + force_fn(n, t);
      });
  }
};

/**
 * @brief Helper function to combine multiple force functors
 * 
 * @tparam TForces types of individual force objects
 * @param forces parameter pack of force functors
 * @return CombinedForce<TForces...> force functor
 */
template <typename... TForces>
CombinedForce<TForces...> make_combined_force(TForces... forces) {
  return CombinedForce<TForces...>(forces...);
}

/** Pinning constraint object for use with mass-spring simulation **/
struct PinConstraint : Constraint {
  /** Constrains nodes at (0,0,0) and (1,0,0) to remain fixed */
  void operator()(GraphType& graph, double t) const {
    (void) t;  // Silence compiler warnings

    std::for_each(graph.node_begin(), graph.node_end(),
      [this](Node node) { constrain_node(node); });
  }

 private:
  // Node-specific function that applies pin constraint
  void constrain_node(Node& node) const {
    if (node.value().is_fixed) {
      node.position() = node.value().orig_pos;
      node.value().vel = Point(0);
    }
  }
};

/** Plane constraint object for use with mass-spring simulation **/
struct PlaneConstraint : public Constraint {
  Point support_vec = Point(0, 0, 1);
  double support_val = -0.75;

  /** Constrains nodes to be one side of the half-plane given
   * by @a support_vec and @a support_val.
   *
   * @post Graph has nodes @a n_i s.t. if
   *  1) original @a n_i.position() * @a support_vec >= @a support_val
   *       then final @a n_i is unchanged
   *  2) original @a n_i.position() * @a support_vec < @a support_val
   *       then final @a n_i.position() is nearest point on half-plane
   *        and the component of @a n_i.value().vel in direction of @a support_vec is zeroed
   */
  void operator()(GraphType& graph, double t) const {
    (void) t;  // Silence compiler warnings

    std::for_each(graph.node_begin(), graph.node_end(),
      [this](Node node) { constrain_node(node); });
  }

 private:
  // Node-specific function that applies plane constraint
  void constrain_node(Node& node) const {
    double pos_dot = dot(node.position(), support_vec);

    if (pos_dot < support_val) {
      double vel_dot = dot(node.value().vel, support_vec);
      double support_dot = dot(support_vec, support_vec);

      node.position() += (support_val - pos_dot) / support_dot * support_vec;
      node.value().vel -= vel_dot / support_dot * support_vec;
    }
  }
};

/** Sphere constraint object for use with mass-spring simulation **/
struct SphereConstraint : public Constraint {
  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  /** Constrains nodes to be outside of the sphere given
   *  by @a center and @a radius.
   *
   * @post Graph has nodes @a n_i s.t. if
   *  1) original | @a n_i.position() - @a center | >= @a radius
   *       then final @a n_i is unchanged
   *  2) original | @a n_i.position() - @a center | < @a radius
   *       then final @a n_i.position() is nearest point on surface of the sphere
   *       and the component of @a n_i.value().vel in
   *       direction normal to the sphere's surface is zeroed
   */
  void operator()(GraphType& graph, double t) const {
    (void) t;  // Silence compiler warnings

    std::for_each(graph.node_begin(), graph.node_end(),
      [this](Node node) { constrain_node(node); });
  }

 private:
  // Node-specific function that applies sphere constraint
  void constrain_node(Node& node) const {
    Point support_vec = node.position() - center;
    double support_dot = dot(support_vec, support_vec);
    double support_rad = norm(support_vec);

    if (support_rad < radius) {
      double vel_dot = dot(node.value().vel, support_vec);

      node.position() += (radius / support_rad - 1) * support_vec;
      node.value().vel -= vel_dot / support_dot * support_vec;
    }
  }
};

/** Sphere tearing constraint object for use with mass-spring simulation **/
struct SphereTearConstraint : public Constraint {
  Point center = Point(0.5, 0.5, -0.5);
  double radius = 0.15;

  /** Removes nodes inside of the sphere given
   *  by @a center and @a radius.
   *
   * @post Graph has nodes @a n_i s.t. if
   *  1) original | @a n_i.position() - @a center | >= @a radius
   *       then final @a n_i is unchanged
   *  2) original | @a n_i.position() - @a center | < @a radius
   *       then @a n_i is removed from the graph
   */
  void operator()(GraphType& graph, double t) const {
    (void) t;  // Silence compiler warnings

    auto node_it = graph.node_begin();
    while (node_it != graph.node_end()) {
      Point support_vec = (*node_it).position() - center;
      double support_rad = norm(support_vec);

      if (support_rad < radius)
        node_it = graph.remove_node(node_it);
      else
        ++node_it;
    }
  }
};

/**
 * @brief Constraint functor that combines a parameter pack of other constraint functors
 * 
 * @tparam TConstraints types of individual constraint objects
 */
template <typename ...TContraints>
struct CombinedConstraint : public Constraint {
  using ConstraintFunction = std::function<void(GraphType&, double)>;
  std::vector<ConstraintFunction> constraints;

  explicit CombinedConstraint(TContraints... constraints)
    : constraints({ constraints... }) {
  }

  /** Applies each constraint applying to @a n at time @a t sequentially. */
  void operator()(GraphType& graph, double t) const {
    std::for_each(constraints.begin(), constraints.end(),
      [&graph, &t](const ConstraintFunction& constraint_fn) {
        constraint_fn(graph, t);
      });
  }
};

/**
 * @brief Helper function to combine multiple constraint functors
 * 
 * @tparam TContraints types of individual constraint objects
 * @param constraints parameter pack of constraint functors
 * @return CombinedConstraint<TConstraints...> constraint functor
 */
template <typename... TConstraints>
CombinedConstraint<TConstraints...> make_combined_constraint(
  TConstraints... constraints
) {
  return CombinedConstraint<TConstraints...>(constraints...);
}

/** Force function object for HW2 #1. */
struct Problem1Force : public Force {
  /** Return the force applying to @a n at time @a t.
   *
   * This is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double t) const {
    (void) t;  // Silence compiler warnings

    if (n.position() == Point(0) || n.position() == Point(1, 0, 0))
      return Point(0);

    Point f_spring = MassSpringForce()(n, t);
    Point f_grav = GravityForce()(n, t);

    return f_spring + f_grav;
  }
};

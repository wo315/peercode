/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <chrono>
#include <fstream>
#include <thread>

#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  // Initial location of each node, this attribute mostly
  // helps with the constrains for pining a node
  Point initial_loc;
  Point vel;   //< Node velocity
  double mass; //< Node mass
  NodeData() : initial_loc(0), vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  // Edge spring constant
  double K;
  // Edge rest length
  double L;
  // Default constructor
  EdgeData() : K(100.), L(1.) {}
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
 * @tparam G::node_value_type supports data of type NodeData defined
 * for each node
 * @tparam G is the type of graph, also showing data type used for
 * 					 its nodes and edges
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F>
double symp_euler_step(G &g, double t, double dt, F force) {
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

/** @struct Problem1Force
 *  @brief A struct used to find the force on a node at
 *  each time step.
 *
 *	We use this struct as a functor to find the force on
 *	each node of the graph.
 */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * This is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE> Point operator()(NODE n, double t) {
    // For constrains
    if (n.position() == Point(0, 0, 0) or n.position() == Point(1, 0, 0))
      // Return a zero vector force
      return Point(0, 0, 0);

    // The location of the point that we are currently on
    const Point point_i = n.position();

    // Total force acting on the node
    Point tot_force = Point(0., 0., 0.);

    // We go over all the edges connected to this node
    for (auto first = n.edge_begin(); first != n.edge_end(); ++first) {
      Point delta_x = point_i - (*first).node2().position();
      // Norm of the vector connecting two nodes
      double l2_norm = norm(delta_x);
      // adding the spring force
      tot_force += -(*first).value().K * (delta_x) / l2_norm *
                   (l2_norm - (*first).value().L);
    }

    // Force contribution from gravity
    tot_force += n.value().mass * Point(0, 0, -grav);

    (void)t; // To silence the compiler for now

    return tot_force;
  }
};

// ########################################################
// ==> Force Part
// ########################################################

/** @class Force
 * 	@brief A parent class for all forces used.
 *
 *	We use this class as a functor to find the force on
 *	each node of the graph.
 */
class Force {
public:
  virtual Point operator()(const Node &n, double t) const {
    (void)n; // To silence the compiler
    (void)t; // To silence the compiler
    return Point(0., 0., 0.);
  }
};

/** A class for gravity force */
class GravityForce : public Force {
public:
  Point operator()(const Node &n, double t) const override {
    (void)t; // To silence the compiler
    return (n.value().mass * Point(0, 0, -grav));
  }
};

/** The class for mass-spring forces */
class MassSpringForce : public Force {
public:
  Point operator()(const Node &n, double t) const override {
    (void)t; // To silence the compiler
    // The location of the point that we are currently on
    const Point point_i = n.position();
    // Total force acting on the node
    Point tot_force = Point(0., 0., 0.);
    // We go over all the edges connected to this node
    for (auto first = n.edge_begin(); first != n.edge_end(); ++first) {
      Point delta_x = point_i - (*first).node2().position();
      // Norm of the vector connecting two nodes
      double l2_norm = norm(delta_x);
      tot_force += -(*first).value().K * (delta_x) / l2_norm *
                   (l2_norm - (*first).value().L);
    }
    return tot_force;
  }
};

/** A class for damping force */
class DampingForce : public Force {
  // Damping constant;
  double damping_cnst = .01;

public:
  // Default c-tor
  DampingForce() : damping_cnst(.01) {}
  // Ordinary c-tor
  DampingForce(double c) : damping_cnst(c) {}
  Point operator()(const Node &n, double t) const override {
    (void)t; // To silence the compiler
    return (-damping_cnst * n.value().vel);
  }
};

/** @struct CombinedForce
 *  @brief a structure to combine different forces
 *  and add their forces.
 */
struct CombinedForce {
  // A vector (of pointer) of all
  // the forces that we are going to act
  // on the system
  std::vector<const Force *> all_forces;
  // Ordinary c-tro
  CombinedForce(const std::vector<const Force *> &given_f)
      : all_forces(given_f) {}
  // Operator of this structure that calls the appropriate methods of classes
  Point operator()(const Node &n, double t) {
    Point totforce = Point(0., 0., 0.);
    for (auto &f : all_forces)
      totforce += f->operator()(n, t);
    return totforce;
  }
};

/** A function to combine multiple forces (inherited from Force)
 * and give back a CombinedForce object.
 *
 * @tparam T type of the force used (from Force family)
 * @tparam E type of the force used (from Force family)
 * @tparam G type of the force used (default Force)
 */
template <typename T, typename E, typename G = Force>
CombinedForce make_combined_force(const T &force1, const E &force2,
                                  const G &force3 = Force{}) {
  std::vector<const Force *> forces;
  forces.push_back(&force1);
  forces.push_back(&force2);
  // @note Default value of last input argument is a zero
  // force that would not add anything to the system if
  // left to its default value.
  forces.push_back(&force3);
  return CombinedForce(forces);
}

// ########################################################
// ==> Constrain Part
// ########################################################

/** @class Constraint
 *  @brief The parent class to apply constraint on the graph
 */
class Constraint {
public:
  virtual void operator()(GraphType &g, double t) const {
    (void)g; // To silence the compiler
    (void)t; // To silence the compiler
  }
};

/** A class for pining two nodes of the graph */
class PinConstraint : public Constraint {
public:
  void operator()(GraphType &g, double t) const override {
    (void)t; // To silence the compiler
    // We go over all the nodes and apply the constraint
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // We check and see if the initial location of nodes
      // are those that we are interested in, we put their
      // location back
      auto initial_loc = n.value().initial_loc;
      if (initial_loc == Point(0, 0, 0))
        n.position() = Point(0, 0, 0);
      else if (initial_loc == Point(1, 0, 0))
        n.position() = Point(1, 0, 0);
      else
        continue;
    }
  }
};

/** A class for a plane constrain */
class PlaneConstraint : public Constraint {
  // Shows the plane z
  double plane_z = -0.75;

public:
  // Default c-tor
  PlaneConstraint() : plane_z(-0.75) {}
  // Ordinary c-tor
  PlaneConstraint(double given_z) : plane_z(given_z) {}
  void operator()(GraphType &g, double t) const override {
    (void)t; // To silence the compiler
    // We go over all the nodes and apply the constraint
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if (n.position().z < plane_z) {
        // Set the position to the nearest point on the plane
        n.position().z = plane_z;
        // Set z-component of the Node velocity to zero
        n.value().vel.z = 0.;
      }
    }
  }
};

/** A class for sphere constrain */
class SphereConstraint : public Constraint {
  // Data showing the center point of a sphere
  // as well as its radius
  Point center;
  double r;

public:
  // Default c-tor
  SphereConstraint() : center(0.5, 0.5, -0.5), r(0.15) {}
  // Ordinary c-tor
  SphereConstraint(Point p, double g_r) : r(g_r) { center = p; }
  void operator()(GraphType &g, double t) const override {
    (void)t; // To silence the compiler
    // We go over all the nodes and apply the constraint
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // delta is the vector from the node to the center
      // of the sphere
      auto delta = n.position() - center;
      // If we need to apply the constrain
      if (norm(delta) < r) {
        // R is the unit vector from the node
        // to the center of the sphere
        Point R = delta / norm(delta);
        // This shows how much the point should
        // move in the direction of R to be on
        // the sphere
        // @note This is always a positive number
        double dist = r - norm(delta);
        // Set the position to the nearest point on the surface of the sphere
        n.position() += dist * R;
        n.value().vel += -dot(n.value().vel, R) * R;
      }
    }
  }
};

/** A class for sphere constrain
 * This class is like the previous class
 * In this class we remove the nodes violating the constrain
 */
class RemoveSphereConstraint : public Constraint {
  // Data showing the center point of a sphere
  // as well as its radius
  Point center;
  double r;

public:
  // Default c-tor
  RemoveSphereConstraint() : center(0.5, 0.5, -0.5), r(0.15) {}
  // Ordinary c-tor
  RemoveSphereConstraint(Point p, double g_r) : r(g_r) { center = p; }
  void operator()(GraphType &g, double t) const override {
    (void)t; // To silence the compiler
    // We go over all the nodes and apply the constraint
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // delta is the vector from the node to the center
      // of the sphere
      auto delta = n.position() - center;
      // If we need to apply the constrain
      if (norm(delta) < r) {
        // Remove this node
        g.remove_node(n);
      }
    }
  }
};

/** @struct CombinedConstraints
 * @brief a structure to combine different constraint
 * and apply them on a given graph
 */
struct CombinedConstraints {
  // A vector (of pointer) of all
  // the constrains that we are going to apply
  // to the system
  std::vector<const Constraint *> all_cons;
  // Ordinary c-tro
  CombinedConstraints(const std::vector<const Constraint *> &given_c)
      : all_cons(given_c) {}
  // Operator of this structure that calls the appropriate methods of classes
  void operator()(GraphType &g, double t) const {
    for (auto &c : all_cons)
      c->operator()(g, t);
  }
};

/** A function to combine multiple constrains (inherited from Constraint)
 * and give back a CombinedConstraints object.
 *
 * @tparam T type of the constrain used (from Constraint family)
 * @tparam E type of the constrain used (from Constraint family)
 * @tparam G type of the constrain used (default Constraint)
 */
template <typename T, typename E, typename G = Constraint>
CombinedConstraints make_combined_constraint(const T &c1, const E &c2,
                                             const G &c3 = Constraint{}) {
  std::vector<const Constraint *> all_cons;
  all_cons.push_back(&c1);
  all_cons.push_back(&c2);
  // @note Default value of this argument would not change
  // the graph
  all_cons.push_back(&c3);
  return CombinedConstraints(all_cons);
}

// ########################################################
// ==> symp euler step
// ########################################################

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constrain Function object defining the constraint on
 * 								system
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports data of type NodeData defined
 * for each node
 * @tparam G is the type of graph, also showing data type used for
 * 					 its nodes and edges
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a functor called as constraint(_g_, _dt_)
 * 					 and fixes the system based on a
 * 					 specific constrain
 */
template <typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constrains on the system
  // @note after this call x^{n+1} and v^{n} are fixed
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

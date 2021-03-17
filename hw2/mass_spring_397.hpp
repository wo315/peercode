/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <algorithm>
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
  Point ori_pos;   //< Node original position
  NodeData() : vel(Point(0)), mass(1), ori_pos(Point(0)) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double spr_const;       //< Spring constant
  double rest_len;        //< Spring rest_length
  EdgeData() : spr_const(100), rest_len(1) {}
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
    if (n.position() != Point(0) && n.position() != Point(1, 0, 0)) {
      n.position() += n.value().vel * dt;
    }
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() != Point(0) && n.position() != Point(1, 0, 0)) {
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}


struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * This is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t;
    // Constrain two corners, (0, 0, 0) and (1, 0, 0), of the cloth.
    Point force = Point(0);
    if (n.position() == Point(0) || n.position() == Point(1, 0, 0)) {
      return force;
    }
    // Compute and add mass-spring force.
    for (auto n_adj_it = n.edge_begin(); n_adj_it != n.edge_end(); ++n_adj_it) {
      auto edge_curr = *n_adj_it;
      force += -1 * edge_curr.value().spr_const
                  * ((edge_curr).node1().position() -
                     (edge_curr).node2().position()) / ((edge_curr).length())
                  * ((edge_curr).length() - edge_curr.value().rest_len);
    }
    // Compute and add gravity.
    force += n.value().mass * Point(0, 0, -grav);
    return force;
  }
};

//
// Force
//

/** @class Force
 * @brief Functor represents a force applying to a node at some time. */
class Force {
 public:
  /** Return the force applying to node @a n at time @a t.
   * @return A Point that represents the direction and magnitude of the force.
   *
   * @param[in] n  The applied node.
   * @param[in] t  The current time.
   */
  virtual Point operator()(Node n, double t) {
    (void) n; (void) t;
    return Point(0);
  }
};

/** @class GravityForce
 * @brief Functor represents gravity applying to a node at some time. */
class GravityForce : public Force {
 public:
  /** Return the gravitational force applying to node @a n at time @a t.
   * @return A Point that represents the direction and magnitude of the
   *         gravitational force on Node _n_ at time _t_.
   * @param[in] n  The applied node.
   * @param[in] t  The current time.
   */
  Point operator()(Node n, double t) {
    (void) t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** @class MassSpringForce
 * @brief Functor represents mass spring force applying to a node at some time.
 */
class MassSpringForce : public Force {
 public:
  /** Return the mass spring force applying to node @a n at time @a t.
   * @return A Point that represents the direction and magnitude of the
   *         mass spring force on Node _n_ at time _t_.
   * @param[in] n  The applied node.
   * @param[in] t  The current time.
   */
  Point operator()(Node n, double t) {
    (void) t;
    Point springForce = Point(0);
    // Sum up mass spring force from different edges incident to the node.
    for (auto n_adj_it = n.edge_begin(); n_adj_it != n.edge_end(); ++n_adj_it) {
      auto edge_curr = *n_adj_it;
      springForce += -1 * edge_curr.value().spr_const *
                    ((edge_curr).node1().position() -
                     (edge_curr).node2().position()) / ((edge_curr).length()) *
                    ((edge_curr).length() - edge_curr.value().rest_len);
    }
    return springForce;
  }
};

/** @class DampingForce
 * @brief Functor represents damping force applying to a node at some time.
 */
class DampingForce : public Force {
 public:
  DampingForce() : damp_const(1.0) {}
  /**
   * Construct a DampingForce by the damping constant.
   * @param[in] damping_const  Damping constant.
   */
  DampingForce(double damping_const) : damp_const(damping_const) {}

  /** Return the damping force applying to node @a n at time @a t.
   * @return A Point that represents the direction and magnitude of the
   *         damping force on Node _n_ at time _t_.
   * @param[in] n  The applied node.
   * @param[in] t  The current time.
   */
  Point operator()(Node n, double t) {
    (void) t;
    return -n.value().vel * damp_const;
  }
 private:
  double damp_const;
};

/** @class CombinedForce
 * @brief Functor represents a combined force sourced from different forces
 *        applying to a node at some time.
 */
class CombinedForce {
  std::vector<Force*> forces;
 public:
  /**
   * Construct a CombinedForce by a vector of pointers to Force.
   * @param[in] ext_forces  A vector of pointers to Force to different external
   *                        forces.
   */
  CombinedForce(std::vector<Force*> ext_forces) {
    forces = ext_forces;
  }

  /** Return the combined force applying to node @a n at time @a t.
   * @return A Point that represents the direction and magnitude of the
   *         combined force on Node _n_ at time _t_.
   * @param[in] n  The applied node.
   * @param[in] t  The current time.
   */
  Point operator()(Node n, double t) {
    Point force = Point();
    for (auto indi_force : forces) {
      force += indi_force->operator()(n, t);
    }
    return force;
  }
};

/** Update the vector of pointers to Force for the combined force in the base
 * case.
 *
 * @param[in,out] comb_forces  The vector of pointers to Force.
 * @param[in]     force        A force in the combined force.
 * @tparam        F            A Force type of the forces.
 */
template <typename F>
void combined_force_helper_in_place(std::vector<Force*>& comb_forces,
                                    F& force) {
  comb_forces.push_back(&force);
}

/** Update the vector of pointers to Force for the combined force in the
 * recursive case.
 *
 * @param[in,out] comb_forces  The vector of pointers to Force.
 * @param[in]     force        A force in the combined force.
 * @param[in]     forces       Many forces in the combined force.
 * @tparam        F            A Force type of the forces.
 * @tparam        Fs           A parameter pack of Force types.
 */
template <typename F, typename... Fs>
void combined_force_helper_in_place(std::vector<Force*>& comb_forces,
                                    F& force, Fs&&... forces) {
  comb_forces.push_back(&force);
  combined_force_helper_in_place(comb_forces, forces...);
}

/** A helper function to construct a CombinedForce functor by constructing
 * the vector of pointers to Force.
 *
 * @param[in]     force        A force in the combined force.
 * @param[in]     forces       Many forces in the combined force.
 * @tparam        F            A Force type of the forces.
 * @tparam        Fs           A parameter pack of Force types.
 */
template <typename F, typename... Fs>
CombinedForce make_combined_force(F force, Fs... forces) {
  std::vector<Force*> combined_forces;
  combined_force_helper_in_place(combined_forces, force, forces...);
  return CombinedForce(combined_forces);
}

//
// Constraint
//

/** @class Constraint
 * @brief Functor represents a constraint applying to a graph at some time. */
class Constraint {
 public:
  /** Update position and velocity for all nodes violating a constraint in
   * Graph @a graph at time @a t.
   * @param[in, out] graph  The applied graph.
   * @param[in]      t      The current time.
   */
  virtual void operator()(GraphType& graph, double t) {
    (void) graph; (void) t;
  }
};

/** @class PinConstraint
 * @brief Functor represents a pin constraint applying to a graph at some time,
 *        where the pinned Points (0, 0, 0) and (1, 0, 0) cannot move during
 *        time.
 */
class PinConstraint : public Constraint {
 public:
  /** Reset position for all nodes violating a pin constraint
   * (i.e., the nodes with starting position (0, 0, 0) or (1, 0, 0))
   * in Graph @a graph at time @a t.
   * @param[in, out] graph  The applied graph.
   * @param[in]      t      The current time.
   */
  void operator()(GraphType& graph, double t) {
    (void) t;
    for (auto node_iter = graph.node_begin(); node_iter != graph.node_end();
         ++node_iter) {
      if ((*node_iter).value().ori_pos == Point(0)) {
        (*node_iter).position() = (*node_iter).value().ori_pos;
      } else if ((*node_iter).value().ori_pos == Point(1, 0, 0)) {
        (*node_iter).position() = (*node_iter).value().ori_pos;
      }
    }
  }
};

/** @class PlaneConstraint
 * @brief Functor represents a plane constraint applying to a graph at some
 *        time, where the nodes cannot be below Plane z = -0.75 during
 *        time.
 */
class PlaneConstraint : public Constraint {
 public:
  /** Update position and velocity for all nodes violating a plane constraint
   * (x_i).z < -0.75 (x_i denotes the position of the node) in Graph @a graph
   * at time @a t.
   * @param[in,out] graph  The applied graph.
   * @param[in]     t      The current time.
   */
  void operator()(GraphType& graph, double t) {
    (void) t;
    for (auto node_iter = graph.node_begin(); node_iter != graph.node_end();
         ++node_iter) {
      if (dot((*node_iter).position(), Point(0, 0, 1)) < -0.75) {
        // Set the position to the nearest point on the plane.
        (*node_iter).position().z = -0.75;
        // Set the z-component of the Node velocity to zero.
        (*node_iter).value().vel.z = 0;
      }
    }
  }
};

/** @class SphereConstraint
 * @brief Functor represents a sphere constraint applying to a graph at some
 *        time, where the nodes cannot get into a sphere defined by a center c
 *        and a radius r during time.
 */
class SphereConstraint : public Constraint {
 public:
  SphereConstraint() : center(Point(0.5, 0.5, -0.5)), radius(0.15) {}
  /**
   * Construct a SphereConstraint by the center position _c_ and radius _r_.
   * @param[in] c  Center position.
   * @param[in] r  Radius of the sphere.
   */
  SphereConstraint(const Point& c, double r) : center(c), radius(r) {}

  /** Update position and velocity for all nodes violating a sphere constraint
   * |x_i - c| < r (x_i denotes the position of the node) in Graph @a graph
   * at time @a t.
   * @param[in,out] graph  The applied graph.
   * @param[in]     t      The current time.
   */
  void operator()(GraphType& graph, double t) {
    (void) t;
    for (auto node_iter = graph.node_begin(); node_iter != graph.node_end();
         ++node_iter) {
      if (norm((*node_iter).position() - center) < radius) {
        // Set the position to the nearest point on the surface of the sphere.
        auto direction = unitDirection((*node_iter).position());
        (*node_iter).position() = center + direction * radius;
        // Set the component of the velocity that is normal to the sphere’s
        // surface to zero.
        (*node_iter).value().vel -= dot((*node_iter).value().vel, direction)
                                    * direction;
      }
    }
  }
 private:
  Point  center;
  double radius;

  /** Return a normalized direction from center to the given _point_.
   * @param[in] point  The given Point position.
   */
  Point unitDirection(const Point& point) {
    return (point - center) / norm(point - center);
  }
};

/** @class DestroySphereConstraint
 * @brief Functor represents a sphere constraint applying to a graph at some
 *        time, where the nodes getting into a sphere defined by a center c
 *        and a radius r during time will be removed from the graph.
 */
class DestroySphereConstraint : public Constraint {
 public:
  DestroySphereConstraint() : center(Point(0.5, 0.5, -0.5)), radius(0.15) {}
  /**
   * Construct a SphereConstraint by the center position _c_ and radius _r_.
   * @param[in] c  Center position.
   * @param[in] r  Radius of the sphere.
   */
  DestroySphereConstraint(const Point& c, double r) : center(c), radius(r) {}

  /** Remove all nodes violating a sphere constraint |x_i - c| < r
   * (x_i denotes the position of the node) in Graph @a graph at time @a t.
   * @param[in,out] graph  The applied graph.
   * @param[in]     t      The current time.
   */
  void operator()(GraphType& graph, double t) {
    (void) t;
//--design_1
//--You are missing nodes by iterating through a for loop
//--START 
   for (auto node_iter = graph.node_begin(); node_iter != graph.node_end();
         ++node_iter) {
      if (norm((*node_iter).position() - center) < radius) {
        graph.remove_node(node_iter);
      }
    }
//--END
  }
 private:
  Point  center;
  double radius;

  /** Return a normalized direction from center to the given _point_.
   * @param[in] point  The given Point position.
   */
  Point unitDirection(const Point& point) {
    return (point - center) / norm(point - center);
  }
};

/** @class CombinedConstraints
 * @brief Functor represents combined constraints sourced from different
 *        constraints applying to a graph at some time.
 */
class CombinedConstraints {
  std::vector<Constraint*> constraints;
 public:
  /**
   * Construct a CombinedConstraints by a vector of pointers to Constraint.
   * @param[in] ext_constraints  A vector of pointers to Constraint to
   *                             different constraints.
   */
  CombinedConstraints(std::vector<Constraint*> ext_constraints) {
    constraints = ext_constraints;
  }

  /** Update position and velocity for all nodes violating the combined
   * constraint in Graph @a graph at time @a t.
   * @param[in,out] graph  The applied graph.
   * @param[in]     t      The current time.
   */
  void operator()(GraphType& graph, double t) {
    for (auto indi_constraints : constraints) {
      indi_constraints->operator()(graph, t);
    }
  }
};

/** Update the vector of pointers to Constraint for the combined constraints
 * in the base case.
 *
 * @param[in,out] comb_constraints  The vector of pointers to Constraint.
 * @param[in]     constraint        A constraint in the combined constraint.
 * @tparam        C                 A Contraint type of the constraints.
 */
template <typename C>
void combined_constraint_helper_in_place(std::vector<Constraint*>&
                                         comb_constraints, C& constraint) {
  comb_constraints.push_back(&constraint);
}

/** Update the vector of pointers to Constraint for the combined constraints
 * in the recursive case.
 *
 * @param[in,out] comb_constraints  The vector of pointers to Constraint.
 * @param[in]     constraint        A constraint in the combined constraints.
 * @param[in]     constraints       Many constraints in the combined
 *                                  constraints.
 * @tparam        C                 A Constraint type of the constraints.
 * @tparam        Cs                A parameter pack of Constraint types.
 */
template <typename C, typename... Cs>
void combined_constraint_helper_in_place(std::vector<Constraint*>&
                                         comb_constraints, C& constraint,
                                         Cs&&... constraints) {
  comb_constraints.push_back(&constraint);
  combined_constraint_helper_in_place(comb_constraints, constraints...);
}

/** A helper function to construct a CombinedConstraints functor by constructing
 * the vector of pointers to Constraint.
 *
 * @param[in]     constraint      A constraint in the combined constraints.
 * @param[in]     constraints     Many constraints in the combined constraints.
 * @tparam        C               A Constraint type of the constraints.
 * @tparam        Cs              A parameter pack of Constraint types.
 */
template <typename C, typename... Cs>
CombinedConstraints make_combined_constraints(C constraint, Cs... constraints) {
  std::vector<Constraint*> combined_constraints;
  combined_constraint_helper_in_place(combined_constraints, constraint,
                                      constraints...);
  return CombinedConstraints(combined_constraints);
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and graph constraints.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraints  Function object defining the Graph constraints
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(@a g, @a t),
 *           where @a g is the Graph and @a t is the current time.
 *           It updates position and velocity for all nodes violating
 *           constraints in Graph @a g at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Update the position and velocity of all nodes according to graph
  // constraints.
  constraints(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//--design_-1
//--Good job combining rvalues constraints and forces
//--END



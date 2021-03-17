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
  Point vel;     //< Node velocity
  double mass;   //< Node mass
  Point orig;    //< Node original position

  NodeData(Point v = Point(0), double m = 1.0, Point p0 = Point(0))
    : vel(v), mass(m), orig(p0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;      //< Edge spring constant
  double length; //< Edge length (at rest)
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
 * @return the next time step (usually _t_ + _dt_)
 *
 * @tparam G::node_value_type supports NodeData (mass, vel, and first position)
 * @tparam F is a function object called as _force(n, _t_)_,
 *           where n is a node of the graph and _t_ is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time _t_.
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
 *    method with the given node force and constraint on motion
 * @param[in,out] g          Graph
 * @param[in]     t          The current time (useful for time-dependent forces)
 * @param[in]     dt         The time step
 * @param[in]     force      Function object defining the force per node
 * @param[in]     constraint Function object defining constraints on nodes
 * @return the next time step (usually _t_ + _dt_)
 *
 * @tparam G::node_value_type supports NodeData (mass, vel, and first position)
 * @tparam F is a function object called as _force_(n, _t_),
 *           where n is a node of the graph and _t_ is the current time.
 *           _force_ must return a Point representing the force vector on
 *           Node n at time _t_.
 * @tparam C is a function object called as _constraint(_g_, _t_)_, where
 *           where _g_ is the graph and _t_ is the current time. _constraint_
 *           modifies the node positions and velocity such that the constraint
 *           is not violated when the next timestep is taken.
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

  constraint(g, t); // reset positions

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;

}


/** Force function object for HW2 Problem 1 */
struct Problem1Force {

  /** Return the force applying to @a n at time @a t */
 private:

  // private spring force function for individual edges
  Point force_spring(Edge e) {
    EdgeData e_val = e.value();
    Point diff = e.node1().position() - e.node2().position();
    double cur_dist = norm(diff);
    double rest_len = e_val.length;
    double k = e_val.K;
    return -k * diff * (1.0 - rest_len / cur_dist);
  }

 public:
  /** Operator that returns spring and gravity forces applied to
   *  node, as well as fixing the positions of nodes at (0,0,0) and (1,0,0).
   * @param[in] _n_ Node where we evaluate the force
   * @param[in] _t_ current timestep (double) - not used
   * @return Point representing 3d force applied to _n_ at time _t_
   *
   * @pre _n_.graph_ must have edges with edge_value_type containing
   *      K (spring constant) and length (length of edge at rest) as
   *      attributes.
   * @pre _n_ must be provided a node_value_type with a mass attribute.
   *
   * Complexity: O(_n_.graph_->num_edges() / _n_.graph_->num_nodes())
   *             on average
   */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    // gravitational force
    Point f_grav = n.value().mass * Point(0, 0, -grav);

    // spring force
    Point f_spring = Point(0);
    GraphType::incident_iterator ii = n.edge_begin();
    GraphType::incident_iterator ii_last = n.edge_end();

    // for each adjacent node defined by an edge, compute the spring
    // force and add to the existing force.
    for (; ii != ii_last; ++ii) f_spring += force_spring(*ii);

    (void) t; // silence compiler warning
    return f_spring + f_grav;
  }
};

/** @class Force
 *  @brief base class for Force functors that take a node and timestep
 *  and return the force acting on the node at that timestep.
 */
class Force {

 public:
  /** Return 0 as the force acting on the given node at current time
   * @param[in] _n_ the node whose force we are calculating
   * @param[in] _t_ the current time (double), for time dependent forces
   * @return Point force of 0 acting on the node
   *
   * Virtual method that simply provides a template for other types of
   * force functors.
   *
   * Complexity: O(1)
   */
  virtual Point operator()(Node n, double t) {
    (void) n; (void) t;
    return Point(0);
  }
};

/** @class GravityForce
 *  @brief functor that computes the gravitational force acting on the given
 *         node in the z direction based on mass and static constexp grav
 *         (gravitational constant)
 */
class GravityForce : public Force {
 public:

  /** Return the force of gravity acting on the given node
   * @param[in] _n_ the node whose force we are calculating
   * @param[in] _t_ the current time (double) - not used
   * @return Point for gravitational force at current node
   *         (-grav in z direction)
   *
   * @pre _n_ must be provided a node_value_type with a mass attribute.
   *
   * Complexity: O(1)
   */
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warning
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** @class MassSpringForce
 *  @brief functor that computes the spring force acting on the given
 *         node based on positions of adjacent nodes and spring
 *         constants of outbound edges
 */
class MassSpringForce : public Force {
 public:

  /** Return the spring force acting on the given node
   * @param[in] _n_ the node whose force we are calculating
   * @param[in] _t_ the current time (double) - not used
   * @return Point of total spring force applied on the current node
   *         by adjacent nodes
   *
   * @pre _n_.graph_ must have edges with edge_value_type containing
   *      K (spring constant) and length (length of edge at rest) as
   *      attributes.
   *
   * Complexity: O(_n_.graph_->num_edges() / _n_.graph_->num_nodes())
   *             on average
   */
  Point operator()(Node n, double t) {

    Point f_spring = Point(0); // spring force (running sum)
    GraphType::incident_iterator ii = n.edge_begin();
    GraphType::incident_iterator ii_last = n.edge_end();

    // for each adjacent node defined by an edge, compute the spring
    // force and add to the existing force.
    for (; ii != ii_last; ++ii) {
      EdgeData e_val = (*ii).value();
      Node n_j = (*ii).node2();
      Point diff = n.position() - n_j.position();
      double cur_dist = norm(diff);
      double rest_len = e_val.length;
      double k = e_val.K;
      f_spring += - k * diff * (1.0 - rest_len / cur_dist);
    }

    (void) t; // silence compiler warning
    return f_spring;
  }
};

/** @class DampingForce
 *  @brief functor that computes the damping force acting on the given
 *         node based on its velocity and a damping constant.
 */
class DampingForce : public Force {

 double c_; // damping constant

 public:
  /** DampingForce Constructor
   * @param _c_ damping constant (double) - defaults to 0.0 (no damping)
   */
  DampingForce(double c = 0.0) : c_(c) {}

  /** Return the damping force acting on the given node (resistance)
   * @param[in] _n_ the node whose force we are calculating
   * @param[in] _t_ the current time (double) - not used
   * @return Point of damping force applied on the current node
   *
   * @pre _n_ has a node_value_type with a vel (velocity) attribute.
   *
   * Complexity: O(1)
   */
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warning
    return -c_ * n.value().vel;
  }
};

/** @class CombinedForce
 *  @brief functor that computes a combination of various
 *         forces known to apply to a given node
 */
class CombinedForce : public Force {

  std::vector<Force*> forces_; // different force functors (as pointers)

 public:

  /** CombinedForce Constructor
   * @param[in] _forces_ vector of pointers to derived force functors
   *            to be applied to graph nodes.
   */
  CombinedForce(const std::vector<Force*>& forces) : forces_(forces) {}

  /** Returns the combined force acting on the given node
   * @param[in] _n_ the node whose force we are calculating
   * @param[in] _t_ the current time (double) - not used
   * @return Point of combined (sum) force applied on the current node
   *
   * Applies and accumulates each force in forces_ in order. Note that the
   * force must be applied in the same scope in which it is constructed due
   * to the reliance on pointers.
   *
   * Complexity: Dependent on the number and nature of force functors
   */
  Point operator()(Node n, double t) {
    Point f_combo = Point(0);
    for (unsigned int i = 0; i < forces_.size(); i++) {
      f_combo += (*(forces_[i]))(n,t);
    }
    return f_combo;
  }
};

/** Creates a CombinedForce from 2 or 3 force functors
 * @param[in] _f1_ const reference to first force functor to be applied
 * @param[in] _f2_ const reference to second force functor to be applied
 * @param[in] _f3_ const reference to third force functor to be applied (if
 *                 not provided, uses base functor that returns 0)
 * @return CombinedForce functor that applies all of the given forces
 *
 * Abstracts away the creation of a force vector from the client, such that
 * a CombinedForce can be obtained as e.g.
 * make_combined_force(GravityForce(), MassSpringForce()).
 *
 * Complexity: O(1)
 */
CombinedForce make_combined_force(const Force& f1, const Force& f2,
                                  const Force& f3 = Force()) {
  std::vector<Force*> fs;
  fs.push_back(const_cast<Force*>(&f1));
  fs.push_back(const_cast<Force*>(&f2));
  fs.push_back(const_cast<Force*>(&f3));
  return CombinedForce(fs);
}

/** @class Constraint
 *  @brief base class from which constraints on node positions and velocities
 *         within the graph are derived
 */
class Constraint {
 public:

  /** Apply no constraint to the graph
   *  @param[in] _g_ const reference to graph to which we apply the constraint
   *  @param[in] _t_ current time (double)
   *
   *  Default (null) constraint, with no effect on graph, just provides
   *  template for other constraint functors
   *
   *  Complexity: O(1)
   */
  virtual void operator()(const GraphType& g, double t) {
    (void) g; (void) t;
    return;
  }
};

/** @class PinConstraint
 *  @brief applies a constraint such that two of the graph's nodes
 *         retain a fixed position as forces are applied
 */
class PinConstraint : public Constraint {
  Point fp1_ = Point(0,0,0); // first fixed point (pin)
  Point fp2_ = Point(1,0,0); // second fixed point

 public:

  /** Apply a constraint that fixes the graph at two nodes
   * @param[in,out] _g_ const reference to graph to which we apply constraint
   * @param[in]     _t_ current time (double) - not used
   *
   * @pre  _g_ is given a node_value_type with attribute orig corresponding
   *       to the original position of each node.
   * @pre  there exist two nodes n1 and n2 in _g_ such that n1.value().orig
   *       = fp1_ || n2.value().orig = fp2_
   * @post the positions of n1 and n2 will be reset if n1.position()
   *       != fp1_ || n2.position() != fp2_.
   *
   * Complexity: O(_g_.num_nodes())
   */
  void operator()(const GraphType& g, double t) {
    GraphType::node_iterator ni = g.node_begin();
    GraphType::node_iterator ni_end = g.node_end();

    for (; ni != ni_end; ++ni) {
      Node n = *ni;
      if (n.value().orig == fp1_) n.position() = fp1_;
      if (n.value().orig == fp2_) n.position() = fp2_;
    }
    (void) t;
  }
};

/** @class PlaneConstraint
 *  @brief applies a constraint to a graph such that no nodes fall below a
 *         specified 2d plane
 */
class PlaneConstraint : public Constraint {
  Point plane_;   // normalized point that define plane at origin
  double affine_; // translation from the origin

 public:
  /** PlaneConstraint Constructor
   *  @param[in] _p_ point that defines 2d plane at origin (orthogonal vector)
   *  @param[in] _a_ affine translation of the plane from the origin
   *
   *  The 2d plane is defined as (p.x * x) + (p.y * y) + (p.z * z) = a
   *  Note that _p_ need not of unit length - it is normalized in construction.
   *  By default, the plane is horizontal, offset by -0.75
   */
  PlaneConstraint(Point p = Point(0,0,1), double a = -0.75)
    : plane_(p / norm(p)), affine_(a) {}

  /** Apply a constraint that fixes the graph above a 2d plane
   * @param[in,out] _g_ const reference to graph to which we apply constraint
   * @param[in]     _t_ current time (double) - not used
   *
   * @pre  for all nodes n in _g_, dot(n.position(), plane_) >= affine_
   *       at t < _t_ and in particular at t = 0
   * @post for all nodes n in _g_, dot(n.position(), plane_) >= affine_
   *       (constraint is maintained at the next step)
   *
   * Complexity: O(_g_.num_nodes())
   */
  void operator()(const GraphType& g, double t) {
    GraphType::node_iterator ni = g.node_begin();
    GraphType::node_iterator ni_end = g.node_end();

    for (; ni != ni_end; ++ni) {
      Node n = *ni;
      if (dot(n.position(), plane_) < affine_) {
        // nearest point on plane
        n.position() += (affine_ - dot(plane_, n.position())) * plane_;

        // zero out velocity perpendicular to plane
        n.value().vel -= dot(n.value().vel, plane_) * plane_;
      }
    }
    (void) t;
  }
};

/** @class SphereConstraint
 *  @brief apply a constraint to the graph such that no nodes
 *         enter a specified 3d sphere
 */
class SphereConstraint : public Constraint {
  Point ctr_; // center of the sphere
  double r_;  // radius of the sphere

 public:
  /** SphereConstraint Constructor
   * @param[in] _c_ point (coordinates) of sphere center
   * @param[in] _r_ radius of sphere
   *
   * @pre _r_ > 0.0
   *
   * Sphere has center at (0.5,0.5,-0.5) and radius 0.15 by default
   */
  SphereConstraint(Point c = Point(0.5, 0.5, -0.5), double r = 0.15)
    : ctr_(c), r_(r) {}

  /** Apply a constraint that fixes a graph outside of a sphere
   * @param[in,out] _g_ const reference to graph to which we apply constraint
   * @param[in]     _t_ current time (double) - not used
   *
   * @pre  for all nodes n in _g_, norm(n.position() - ctr_) >= r_
   *       at t < _t_ and in particular at t = 0
   * @post for all nodes n in _g_, norm(n.position() - ctr_) >= r_
   *       (constraint is maintained at the next step)
   *
   * Complexity: O(_g_.num_nodes())
   */
  void operator()(const GraphType& g, double t) {
    GraphType::node_iterator ni = g.node_begin();
    GraphType::node_iterator ni_end = g.node_end();
    for (; ni != ni_end; ++ni) {
      Node n = *ni;
      double r_cur = norm(n.position() - ctr_);  // distance from ctr_
      Point R_i = (n.position() - ctr_) / r_cur; // unit direction wrt ctr_
      if ( r_cur < r_ ) {
        n.position() = r_ * R_i + ctr_; // nearest point on sphere
        n.value().vel -= dot(n.value().vel, R_i) * R_i;
      }
    }
    (void) t;
  }
};

/** @class CombinedConstraints
 *  @brief Combine effects of derived Constraints into a single Constraint
 */
class CombinedConstraints : public Constraint {
  using C = Constraint;
  std::vector<C*> cons_; // vector of constraints to apply

 public:
  /** CombinedConstraints Constructor
   * @param _c1_ const reference to first constraint
   * @param _c2_ const reference to second constraint
   * @param _c3_ (opt) const reference to third constraint (default to base)
   */
  CombinedConstraints(const C& c1, const C& c2, const C& c3 = Constraint()) {
    cons_.push_back(const_cast<C*>(&c1));
    cons_.push_back(const_cast<C*>(&c2));
    cons_.push_back(const_cast<C*>(&c3));
  }

  /** Apply constraints sequentially to graph
   * @param _g_ const reference to graph to which we apply constraints
   * @param _t_ current time at which constraint is applied - not used
   */
  void operator()(const GraphType& g, double t) {
    for (unsigned int i = 0; i < cons_.size(); i++) (*(cons_[i]))(g,t);
  }

};

/** @class SphereCutConstraint
 *  @brief Eliminate nodes and corresponding edges that enter spherical region
 */
class SphereCutConstraint : public Constraint {
  Point ctr_; // center of the sphere
  double r_;  // radius of the sphere

 public:
  /** SphereCutConstraint Constructor
   * @param[in] _c_ point (coordinates) of sphere center
   * @param[in] _r_ radius of sphere
   *
   * @pre _r_ > 0.0
   *
   * Sphere has center at (0.5,0.5,-0.5) and radius 0.15 by default
   */
  SphereCutConstraint(Point c = Point(0.5, 0.5, -0.5), double r = 0.15)
    : ctr_(c), r_(r) {}

 /** Apply a constraint that cuts a graph with a spherical constraint
   * @param[in,out] _g_ const reference to graph to which we apply constraint
   * @param[in]     _t_ current time (double) - not used
   *
   * @pre  for all nodes n in _g_, norm(n.position() - ctr_) >= r_
   *       at t < _t_ and in particular at t = 0
   * @post for all nodes n in _g_, norm(n.position() - ctr_) >= r_
   *       (constraint is maintained at the next step)
   *
   * Complexity: O(_g_.num_nodes())
   */
  void operator()(const GraphType& g, double t) {
    GraphType::node_iterator ni = g.node_begin();
    GraphType::node_iterator ni_end = g.node_end();

    while(ni != ni_end) {
      Node n = *ni;
      double r_cur = norm(n.position() - ctr_);
      // remove nodes that fall inside sphere (tear)
      if ( r_cur < r_ ) ni = const_cast<GraphType&>(g).remove_node(ni);
      else ++ni;
    }
    (void) t;
  }

};

//--functionality_0
//--Passed all tests!
//--END

//--design_0
//--Well designed!
//--END

//--documentation_0
//--good
//--END


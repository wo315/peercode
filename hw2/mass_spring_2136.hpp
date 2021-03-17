/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <cmath>

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

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // Initialize zero force.
    Point f(0);
    for (const Edge e : n) {
      double len = e.length();
      double L_ij = e.value();
      double K_ij = 100;
      f += - K_ij * (e.node1().position() - e.node2().position())
           * (len - L_ij)
           / len;
    }
    // Add gravity
    f += Point(0, 0, -grav);
    return f;
  }
};

/** @brief Parent class for force functors. */
struct Force {
  virtual Point operator()(Node n, double t) const {
    // Silence compiler warning.
    (void) t;
    (void) n;
    return Point(0);
  }
};

/** @brief Derived class implementing gravity. */
struct GravityForce : public Force {
  Point operator()(Node n, double t) const {
    // Silence compiler warning.
    (void) t;
    (void) n;
//--design_1
//--Missing node mass
//--START
    return Point(0, 0, -grav);
//--END
  }
};

/** @brief Derived class implementing spring force. */
struct MassSpringForce : public Force {
  double K; // Spring constant. Default 100, but values in the range 1000~5000
            // give more realistic approximations of cloth.
  MassSpringForce(double K = 100) : K(K) {}

  Point operator()(Node n, double t) const {
    // Silence compiler warning.
    (void) t;
    Point f(0);
//--functionality_1
//--your spring force is not working behaving properly
//--(no inversion of force signs)
//--START 
    for (const Edge e : n) {
//--END
      double len = e.length();
      double L_ij = e.value();
      // NOTE: A spring constant of ~100 makes the sheet too stretchy to drape
      // over the sphere. A spring constant of ~5000 is closer to how a sheet
      // might behave.
      double K_ij = K;
      f += - K_ij * (e.node1().position() - e.node2().position())
           * (len - L_ij)
           / len;
    }
    return f;
  }
};

/** @brief Derived class implementing damping force. */
struct DampingForce : public Force {
  // Default constructor allowing for adjustment to damping force.
  DampingForce(double c = 1) : c(c) {}
  Point operator()(Node n, double t) const {
    (void) t;
    return - n.value().vel * c;
  }
 private:
  double c;
};

/** @brief Functor for applying forces of different types.
 *  Forces must be used via pointers to enable polymorphism, but must also be
 *  stored in case the forces used to initialize the CombinedForce go out of
 *  scope before it does. Unique pointers are used to automatically manage this
 *  memory.
 */
struct CombinedForces {
  /** Constructor for combined forces. Pointers to all of the forces are
   *  allocated on the heap.
   */
  CombinedForces(std::vector<std::unique_ptr<Force> >& forces)
    : forces_(std::move(forces)) {
  }

  /** Call operator. Sums all of the forces and returns the result. */
  Point operator()(Node n, double t) {
    auto accumulator = [&n, t] (Point p, const std::unique_ptr<Force>& f) {
      return p + f->operator()(n, t);
    };
    return std::accumulate(forces_.begin(), forces_.end(), Point(0),
                           accumulator);
  }

 private:

  std::vector<std::unique_ptr<Force> >  forces_;
};

template <typename F0, typename F1, typename F2 = Force>
CombinedForces make_combined_force(const F0& f0, const F1& f1,
                                   const F2& f2 = {}) {
  // Initialize vector to hold pointers to forces.
  std::vector<std::unique_ptr<Force> > forces;

  // Add pointers to forces to back of vector.
  forces.emplace_back(new F0 {f0});
  forces.emplace_back(new F1 {f1});
  forces.emplace_back(new F2 {f2});

  return CombinedForces(forces);
}

/** @brief Parent class for constraint functors. */
struct Constraint {
  virtual void operator()(GraphType& g, double t) const {
    // Silence compiler warning.
    (void) g;
    (void) t;
  }
};

/** @brief Constraint to pin certain nodes in position. */
class PinConstraint : public Constraint {
  struct set_lt_ {
    bool operator()(const Point& p0, const Point& p1) const {
      if (p0.x < p1.x) {
        return true;
      } else if (p0.x == p1.x && p0.y < p1.y) {
        return true;
      } else if (p0.x == p1.x && p0.y == p1.y && p0.z < p1.z) {
        return true;
      } else {
        return false;
      }
    }
  };

  std::set<Point, set_lt_> fixed_points_;

 public:

  /** Default constructor. */
  PinConstraint() : fixed_points_ {} {}

  /** Initializer list constructor. */
  PinConstraint(std::initializer_list<Point> l) : fixed_points_ {} {
    for (auto it = l.begin(); it != l.end(); ++it) {
      fixed_points_.insert(*it);
    }
  }
  
  void operator()(GraphType& g, double t) const {
    // Silence compiler warning.
    (void) t;
    
    for (const Point& p : fixed_points_) {
      // Find the closet point to a fixed point.
      auto dist_to_p = [&p] (const Node& n0, const Node& n1) {
        return norm(n0.position() - p) < norm(n1.position() - p);
      };
      Node n = *std::min_element(g.node_begin(), g.node_end(), dist_to_p);

      // Snap to the point and set the velocity to 0.
      n.value().vel = Point(0);
      n.position() = p;
    }
  }
};

/** @brief Constraint to prevent nodes from falling below plane. */
struct PlaneConstraint : public Constraint {
  void operator()(GraphType& g, double t) const {
    // Silence compiler warning.
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;

      // If the node is below -0.75, set its vertical velocity to 0 and snap it
      // above the plane.
      if (n.position().z < -0.75) {
        n.position().z = -0.75;
        n.value().vel.z = 0;
      }
    }
  }
};

/** @brief Constraint to retain nodes within a sphere of given radius and
 * center.
 */

struct SphereConstraint : public Constraint {
  void operator()(GraphType& g, double t) const {
    // Silence compiler warning.
    (void) t;
    static const Point center = Point(0.5, 0.5, -0.5);
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;

      // If the distance between the node and the center is less than 0.15,
      // snap back to nearest point in the sphere and make velocity tangent to
      // the sphere.
      Point diff = n.position() - center;
      double dist = norm(diff);
      Point n_diff = diff / dist;
      if (dist < 0.15) {
        n.position() -= (dist - 0.15) * n_diff;
        n.value().vel -= inner_prod(n.value().vel, n_diff) * n_diff;
      }
    }
  }
};

struct TearConstraint : public Constraint {
  void operator()(GraphType& g, double t) const {
    // Silence compiler warning.
    (void) t;
    static const Point center = Point(0.5, 0.5, -0.5);
    auto it = g.node_begin();
    while (it != g.node_end()) {
      auto n = *it;

      // If the distance between the node and the center is less than 0.15,
      // delete the node.
      double dist = norm(n.position() - center);
      if (dist < 0.15) {
        it = g.remove_node(it);
      } else {
        ++it;
      }
    }
  }
};

/** @brief Functor for applying constraints of different types.  Constraints
 *  must be used via pointers to enable polymorphism, but must also be stored in
 *  case the constraints used to initialize the CombinedConstraints go out of
 *  scope. Unique pointers are used to automatically manage this memory.
 */
struct CombinedConstraints {
  /** Constructor for combined constraints. Pointers to all of the constraints
   *  have to be allocated on the heap.
   */
  CombinedConstraints(std::vector<std::unique_ptr<Constraint> >& constraints)
    : constraints_(std::move(constraints)) {
  }

  /** Call operator. Applies all of the constraints sequentially. */
  void operator()(GraphType& g, double t) const {
    auto apply = [&g, &t] (const std::unique_ptr<Constraint>& c) {
      c->operator()(g, t);
    };

    std::for_each(constraints_.begin(), constraints_.end(), apply);
  }

 private:
  std::vector<std::unique_ptr<Constraint> > constraints_;
};

template <typename C0, typename C1, typename C2 = Constraint>
CombinedConstraints make_combined_constraints(const C0& c0, const C1& c1,
                                             const C2& c2 = {}) {
  // Initialize vector to hold pointers to forces.
  std::vector<std::unique_ptr<Constraint> > constraints;

  // Add pointers to forces to back of vector.
  constraints.emplace_back(new C0 {c0});
  constraints.emplace_back(new C1 {c1});
  constraints.emplace_back(new C2 {c2});

  return CombinedConstraints(constraints);
}

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
template <typename G, typename F, typename C = PinConstraint>
double symp_euler_step(G& g, double t, double dt, F force,
                       C constraint = {Point(0, 0, 0), Point(1, 0, 0)}) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraints.
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

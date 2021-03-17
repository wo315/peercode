/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <unordered_set>

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

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using Size = typename GraphType::size_type;


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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      continue;
    }
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      continue;
    }
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Fix violated constraints
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
  double k_;
  Problem1Force() : k_(100) {}
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0);
    }
    Point res(0, 0, -grav * (n.value().mass));
    for (auto ci = n.edge_begin(); ci != n.edge_end(); ++ci) {
      double delta_length = (*ci).length() - (*ci).value();
      res += k_ * delta_length * ((*ci).node2().position() - n.position()) / (*ci).length();
    }
    (void) t;
    return res;
  }
};

class Force {
 public:
  virtual Point operator()(Node n, double t) const = 0;
  virtual ~Force() = default;
};

class GravityForce : public Force {
 public:
  Point operator()(Node n, double t) const {
    (void) t;
    return Point(0, 0, -grav * (n.value().mass));
  }
};

class MassSpringForce : public Force {
 public:
  MassSpringForce(double k = 100) : k_(k) {}
  Point operator()(Node n, double t) const {
    Point res(0);
    for (auto ci = n.edge_begin(); ci != n.edge_end(); ++ci) {
      double delta_length = (*ci).length() - (*ci).value();
      res += k_* delta_length * ((*ci).node2().position() - n.position()) / (*ci).length();
    }
    (void) t;
    return res;
  }
 private:
  double k_;
};

class DampingForce : public Force {
 public:
  DampingForce(double c = 0.1) : c_(c) {}
  Point operator()(Node n, double t) const {
    (void) t;
    return -c_ * n.value().vel;
  }
 private:
  double c_;
};

class CombinedForce : public Force {
 public:
  CombinedForce(const std::vector<Force*>& forces) : forces_(forces) {}
  Point operator()(Node n, double t) const {
    Point combined_force(0);
    for (const auto& force : forces_) {
      combined_force += (*force)(n, t);
    }
    return combined_force;
  }
 private:
  std::vector<Force*> forces_;
};

template <class F1, class F2, class F3> inline
CombinedForce make_combined_force(const F1& f1, const F2& f2, const F3& f3) {
  std::vector<Force*> forces;
  forces.push_back(&const_cast<F1&>(f1));
  forces.push_back(&const_cast<F2&>(f2));
  forces.push_back(&const_cast<F3&>(f3));
  return CombinedForce(forces);
}

template <class F1, class F2> inline
CombinedForce make_combined_force(const F1& f1, const F2& f2) {
  std::vector<Force*> forces;
  forces.push_back(&const_cast<F1&>(f1));
  forces.push_back(&const_cast<F2&>(f2));
  return CombinedForce(forces);
}

class Constraint {
 public:
  virtual void operator()(GraphType& graph, double t) = 0;
  virtual ~Constraint() = default;
};

class PinConstraint : public Constraint {
 public:
  PinConstraint(const GraphType& graph) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      if ((*ni).position() == Point(0)) {
        constrainted0_idx_.insert((*ni).index());
      } else if ((*ni).position() == Point(1, 0, 0)) {
        constrainted1_idx_.insert((*ni).index());
      }
    }
  }
  void operator()(GraphType& graph, double t) {
    for (Size idx : constrainted0_idx_) {
      graph.node(idx).position() = Point(0);
    }
    for (Size idx : constrainted1_idx_) {
      graph.node(idx).position() = Point(1, 0, 0);
    }
    (void) t;
  }
 private:
  std::unordered_set<Size> constrainted0_idx_;
  std::unordered_set<Size> constrainted1_idx_;
};

class PlaneConstraint : public Constraint {
 public:
  PlaneConstraint(double z = -0.75) : z_(z) {}
  void operator()(GraphType& graph, double t) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      if ((*ni).position().z < z_) {
        (*ni).position().z = z_;
        (*ni).value().vel.z = 0;
      }
    }
    (void) t;
  }
 private:
  double z_;
};

class SphereConstraint : public Constraint {
 public:
  SphereConstraint(Point c = Point(0.5, 0.5, -0.5), double r = 0.15) : c_(c), r_(r) {}
  void operator()(GraphType& graph, double t) {
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      Point dis = (*ni).position() - c_;
      Point R = dis / norm_2(dis);
      if (norm_2(dis) < r_) {
        (*ni).position() = c_ + r_ * R;
        (*ni).value().vel -= (*ni).value().vel * R * R;
      }
    }
    (void) t;
  }
 private:
  Point c_;
  double r_; 
};

class CombinedConstraint : public Constraint {
 public:
  CombinedConstraint(std::vector<Constraint*> constraints) : constraints_(constraints) {}
  void operator()(GraphType& graph, double t) {
    for (const auto& constraint : constraints_) {
      (*constraint)(graph, t);
    }
  }
 private:
  std::vector<Constraint*> constraints_;
};

template <class C1, class C2, class C3> inline
CombinedConstraint make_combined_constraint(const C1& c1, const C2& c2, const C3& c3) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&const_cast<C1&>(c1));
  constraints.push_back(&const_cast<C2&>(c2));
  constraints.push_back(&const_cast<C3&>(c3));
  return CombinedConstraint(constraints);
}

template <class C1, class C2> inline
CombinedConstraint make_combined_constraint(const C1& c1, const C2& c2) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&const_cast<C1&>(c1));
  constraints.push_back(&const_cast<C2&>(c2));
  return CombinedConstraint(constraints);
}

class VanishSphereConstraint : public Constraint {
 public:
  VanishSphereConstraint(Point c = Point(0.5, 0.5, -0.5), double r = 0.15) : c_(c), r_(r) {}
  void operator()(GraphType& graph, double t) {
//--design_1
//--If *ni is removed, then after removal ni actually points to the next node.
//--You should use a while loop. Consider using if {...} else{++ni;}.
//--START
    for (auto ni = graph.node_begin(); ni != graph.node_end(); ++ni) {
      Point dis = (*ni).position() - c_;
      if (norm_2(dis) < r_) {
        graph.remove_node(*ni);
      }
    }
//--END
    (void) t;
  }
 private:
  Point c_;
  double r_; 
};



//--functionality_0
//--Passed all tests!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


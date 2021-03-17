/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <utility>

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
  double dt;
  NodeData() : vel(0), mass(1), dt(0.001) {}
};

struct EdgeData {
  double K;
  double L;
  EdgeData() : K(100.0), L(1.0) {}
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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
        continue;
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

    // Apply constraints
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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
        return Point(0, 0, 0);
    Point f_spring {0, 0, 0};
    for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        auto e = *ii;
        Point diff = n.position() - e.node2().position();
        f_spring -= e.value().K * (diff / e.length()) * (e.length() - e.value().L);
    }
    (void) t;
    return f_spring + (n.value().mass * Point(0, 0, -grav));
  }
};

class Force {
  public:
    virtual Point operator()(Node n, double t) {
        (void) t;
        (void) n;
        return {0, 0, 0};
    }

    virtual ~Force() {}
};

class GravityForce : public Force {
  public:
    Point operator()(Node n, double t) {
        (void) t;
        return n.value().mass * Point(0, 0, -grav);
    }
};

class MassSpringForce : public Force {
  public:
    Point operator()(Node n, double t) {
        Point f_spring {0, 0, 0};
        for (auto ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
            auto e = *ii;
            Point diff = n.position() - e.node2().position();
            f_spring -= e.value().K * (diff / e.length()) * (e.length() - e.value().L);
        }
        (void) t;
        return f_spring;
    }
};

class DampingForce : public Force {
    double c_ {0.0};
  public:
    Point operator()(Node n, double t) {
        (void) t;
        return -c_ * n.value().vel;
    }
};

struct CombinedForce {
    std::vector<Force*> forces_;
    explicit CombinedForce(std::vector<Force*> forces) : forces_(std::move(forces)) {}

    ~CombinedForce() {
        for (auto f : forces_)
            delete f;
    }

    Point operator()(Node n, double t) {
        Point f_total = Point(0, 0, 0);
        for (auto f : forces_)
            f_total += (*f)(n, t);
        return f_total;
    }
};

template <typename F1, typename F2, typename F3 = Force>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3 = Force()) {
    std::vector<Force*> forces;
    auto f1_ptr = new F1(f1);
    auto f2_ptr = new F2(f2);
    auto f3_ptr = new F3(f3);
    forces.push_back(f1_ptr);
    forces.push_back(f2_ptr);
    forces.push_back(f3_ptr);
    return CombinedForce(forces);
}

class Constraint {
  public:
    virtual void operator()(GraphType& g, double t) {
        (void) g;
        (void) t;
    }

    virtual ~Constraint() {}
};

class PinConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t) {
        for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
            auto n = *ni;
            if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
                n.value().vel -= (MassSpringForce()(n, t) + GravityForce()(n, t)) * (n.value().dt / n.value().mass);
            }
        }
        (void) t;
    }
};

class PlaneConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t) {
        for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
            auto n = *ni;
            if (n.position().z < -0.75) {
                n.position().z = -0.75;
                n.value().vel.z = 0;
            }
        }
        (void) t;
    }
};

class SphereConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t) {
        Point c {0.5, 0.5, -0.5};
        double r {0.15};
        for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
            auto n = *ni;
            if (norm(n.position() - c) < r) {
                auto R = (n.position() - c) / norm(n.position() - c);
                n.position() = c + r * R;
                n.value().vel -= inner_prod(n.value().vel, R) * R;
            }
        }
        (void) t;
    }
};

class TearConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t) {
        Point c {0.5, 0.5, -0.5};
        double r {0.15};
        auto n_it = g.node_begin();
        while (n_it != g.node_end()) {
            auto n = *n_it;
            if (norm(n.position() - c) < r)
                n_it = g.remove_node(n_it);
            else
                ++n_it;
        }
        (void) t;
    }
};

struct CombinedConstraints {
    std::vector<Constraint*> constraints_;
    explicit CombinedConstraints(std::vector<Constraint*> constraints) : constraints_(std::move(constraints)) {}

    ~CombinedConstraints() {
        for (auto c : constraints_)
            delete c;
    }

    void operator()(GraphType& g, double t) {
        for (auto c : constraints_) {
            (*c)(g, t);
        }
    }
};

template <typename C1, typename C2, typename C3 = Constraint>
CombinedConstraints make_combined_constraint(C1 c1, C2 c2, C3 c3 = Constraint()) {
    std::vector<Constraint*> constraints;
    auto c1_ptr = new C1(c1);
    auto c2_ptr = new C2(c2);
    auto c3_ptr = new C3(c3);
    constraints.push_back(c1_ptr);
    constraints.push_back(c2_ptr);
    constraints.push_back(c3_ptr);
    return CombinedConstraints(constraints);
}

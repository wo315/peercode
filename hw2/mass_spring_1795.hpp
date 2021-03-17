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
// spring constant
static constexpr double spring_const = 20;
// damping constant
static constexpr double damping_const = 20;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

class Force {
public:
    virtual Point operator()(Node n, double t) const{
        (void) n;
        (void) t;
        return Point(0);
    }
};

class GravityForce: public Force {
    public:
        Point operator()(Node n, double t) const{
            (void) t;
            Point p(0, 0, -grav);
            p *= n.value().mass;
            return p;
        }
};

class MassSpringForce: public Force {
    public:
        Point operator()(Node n, double t) const{
            (void) t;
            Point p(0,0,0);
            for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
                auto edge = *it;
                auto xi = n.position();
                auto xj = edge.node2().position();
                // direction of the spring force
                Point spring_force = Point(0, 0, 0);
                spring_force = (- spring_const * (edge.length() - edge.value()) / edge.length()) * (xi - xj);
                p += spring_force;
            }
            return p;
        }
};

class DampingForce: public Force {
    public:
        Point operator()(Node n, double t) const{
            (void) t;
            Point p;
            p = -damping_const * n.value().vel;
            return p;
        }
};

class ZeroForce: public Force {
    public:
        Point operator()(Node n, double t) const{
            (void) t;
            (void) n;
            return Point(0);
        }
};

class CombinedForce{
    std::vector<const Force*> vec_forces;
    public:
    CombinedForce(std::vector<const Force*> vec) : vec_forces(vec) {}
    Point operator() (Node n, double t) const {
            (void) t;
            Point p;
            for (auto const force : vec_forces) {
                p += (*force)(n, t);
            }
            return p;
    }
};

CombinedForce make_combined_force(const Force& force_a, const Force& force_b, const Force& force_c = ZeroForce()){
    std::vector<const Force*> vec;
    vec.push_back(&force_a);
    vec.push_back(&force_b);
    vec.push_back(&force_c);
    auto comb_f = CombinedForce(vec);
    return comb_f;
}

class Constraint {
public:
    virtual void operator()(GraphType& g, double t) const {
        (void) g;
        (void) t;
    }
};

class PinConstraint: public Constraint {
public:
    const std::vector<GraphType::node_type*> pin_nodes;
    PinConstraint() : pin_nodes(std::vector<GraphType::node_type*>()) {}
    PinConstraint(const std::vector<GraphType::node_type*> pin_nodes) : pin_nodes(pin_nodes) {}
    void operator()(GraphType& g, double t) const{
        (void) t;
        for(auto node_it = g.node_begin(); node_it != g.node_end(); ++node_it) {
            for (auto pin_node : pin_nodes) {
                if (*node_it == *pin_node) {
//--functionality_1
//--This won't work because all the points are getting pinned to 0,
//--instead one should be pinned at Point(1,0,0)
//--START
                    (*node_it).position() = Point(0);
//--END
                }
            }
        }
    }
};

class PlaneConstraint: public Constraint {
public:
    double plane = -0.75;
    void operator()(GraphType& g, double t) const{
        (void) t;
        for(auto node_it = g.node_begin(); node_it != g.node_end(); ++node_it) {
            auto n = *node_it;
            if (n.position()[2] < plane){
                n.position()[2] = plane;
                n.value().vel[2] = 0;
            }
        }
    }
};

class SphereConstraint: public Constraint {
public:
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    void operator()(GraphType& g, double t) const{
        (void) t;
        for(auto node_it = g.node_begin(); node_it != g.node_end(); ++node_it) {
            auto n = *node_it;
            double norm_diff = norm(n.position() - c);
            Point vector_diff = n.position() - c;
            if (norm_diff < r){
                n.position() = c + (vector_diff)/norm_diff*r;
                Point r_i = vector_diff/norm_diff;
                n.value().vel = n.value().vel - dot(n.value().vel, r_i)*r_i;
            }
        }
    }
};

class TearConstraint: public Constraint {
public:
    Point c = Point(0.5,0.5,-0.5);
    double r = 0.15;
    void operator()(GraphType& g, double t) const{
        (void) t;
        for(auto node_it = g.node_begin(); node_it != g.node_end(); ++node_it) {
            auto n = *node_it;
            double norm_diff = norm(n.position() - c);
            if (norm_diff < r){
                node_it = g.remove_node(node_it);
            }
        }
    }
};

class ZeroConstraint: public Constraint {
public:
    void operator()(const GraphType& g, double t) const {
        (void) g;
        (void) t;
    }
};

class CombinedConstraints{
    std::vector<const Constraint*> constraints;
    public:
    CombinedConstraints(std::vector<const Constraint*> cons) : constraints(cons) {}
    void operator() (GraphType& g, double t) const {
            for (auto cons : constraints)
                (*cons)(g, t);
    }
};

CombinedConstraints make_combined_constraints(const Constraint& cons_a = ZeroConstraint(), const Constraint& cons_b = ZeroConstraint(), const Constraint& cons_c = ZeroConstraint()){
    std::vector<const Constraint*> vec;
    vec.push_back(&cons_a);
    vec.push_back(&cons_b);
    vec.push_back(&cons_c);
    auto comb_cons = CombinedConstraints(vec);
    return comb_cons;
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // apply the constraints
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
    (void) t;
    if (n.position () ==  Point (0,0,0) || n.position () ==  Point (1,0,0))
        return  Point(0,0,0);
    // gravity force
    Point p(0, 0, -grav);
    p *= n.value().mass;

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        auto edge = *it;
        auto xi = n.position();
        auto xj = edge.node2().position();
        // direction of the spring force
        Point spring_force = Point(0, 0, 0);
        spring_force = (- spring_const * (edge.length() - edge.value()) / edge.length()) * (xi - xj);
        p += spring_force;
    }
    return p;
  }
};




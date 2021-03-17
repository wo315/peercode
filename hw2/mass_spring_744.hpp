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
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  double valid;    //< Node in constraint
  NodeData() : vel(0), mass(1), valid(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double spring_const;       //< Node velocity
  double rest_len;     //< Node mass
  EdgeData() : spring_const(100), rest_len(0.1) {}
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Only update velocities for points that are not pinned
    if (n.value().valid == 1) {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
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
    if (n.position() ==  Point(0,0,0) || n.position() ==  Point(1,0,0)) {
      return  Point(0,0,0);
    }
    else {
      Point f_grav = Point(0, 0, -n.value().mass*grav);
      Point f_spring = Point(0, 0, 0);
      for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
        Edge e = *iit;
        Node n2 = e.node2();
        double dist_val = norm(n.position() - n2.position());
        f_spring -= (n.position() - n2.position())/dist_val * (dist_val - e.value().rest_len)*e.value().spring_const;
      }
      Point out = f_spring + f_grav;
      return f_spring + f_grav;
    }
  }
};

/** Define parent class Force */
class Force {
  public:
    virtual Point operator()(Node n, double t){
      return Point(0, 0, 0);
    };
};

/** Gravity Force Derived class */
class GravityForce : public Force {
  public:
    Point operator()(Node n, double t){
      return Point(0, 0, -n.value().mass*grav);
    };
};

/** Damping Force Derived Class */
class DampingForce : public Force {
  public:
    Point operator()(Node n, double t){
      return -c_val*n.value().vel;
    };
  private:
    double c_val = 10;
};

/** Mass spring Force derived class */
class MassSpringForce : public Force {
  public:
    Point operator()(Node n, double t){
      Point f_spring = Point(0, 0, 0);
      for (auto iit = n.edge_begin(); iit != n.edge_end(); ++iit) {
        Edge e = *iit;
        Node n2 = e.node2();
        double dist_val = norm(n.position() - n2.position());
        f_spring -= e.value().spring_const*(n.position() - n2.position())/dist_val * (dist_val - e.value().rest_len);
      }
      return f_spring;
    };
};

/** Force function object that adds up the forces for an input vector of force pointers */
struct CombinedForce {
  CombinedForce(std::vector<Force*> _forces) : forces(_forces) {}
  Point operator()(Node n, double t) {
    Point out = Point(0, 0, 0);
    for (unsigned int i = 0; i<forces.size(); ++i) {
      out = out + (*forces[i])(n, t);
    }
    return out;
  }
  private:
    std::vector<Force*> forces;
};

/** Function that takes in 2 or 3 arguments and returns a CombinedForce object */
template <typename F1, typename F2, typename F3 = Force>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3 = Force()) {
  std::vector<Force*> a;
  F1* ptr = &f1;
  F2* ptr2 = &f2;
  F3* ptr3 = &f3;
  a.push_back(ptr);
  a.push_back(ptr2);
  a.push_back(ptr3);
  CombinedForce comb(a);
  return comb;
};

/** Define parent constraint class */
class Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
    };
};

/** Define pin Constraint derived class */
class PinConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t){
      for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
        if ((*nit).position() ==  Point(0,0,0) || (*nit).position() ==  Point(1,0,0)) {
          (*nit).value().valid = 0;
        }
      }
    };
};

/** Define plane Constraint derived class */
class PlaneConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t){
      for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
        if ((*nit).position().z < -0.75) {
          (*nit).position().z = -0.75;
          (*nit).value().vel.z = 0;
        }
      }
    };
};

/** Define sphere Constraint derived class */
class SphereConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t){
      for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
        Point center = Point(0.5, 0.5, -0.5);
        double radius = 0.15;
        if (norm((*nit).position() - center) < radius) {
          Point unit_vec = ((*nit).position() - center)/norm((*nit).position() - center);
          (*nit).position() = center + unit_vec*radius;
          (*nit).value().vel = (*nit).value().vel - dot((*nit).value().vel, unit_vec)*unit_vec;
        }
      }
    };
};

/** Define tear Constraint derived class */
class TearConstraint : public Constraint {
  public:
    void operator()(GraphType& g, double t){
      for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
        Point center = Point(0.5, 0.5, -0.5);
        double radius = 0.15;
        Node n1 = *nit;
        if (norm((*nit).position() - center) < radius) {
          nit = g.remove_node(nit);
        }
      }
    };
};

/** Force function object that evaluates the constraints for an input vector of constraint pointers */
struct CombinedConstraints {
  CombinedConstraints(std::vector<Constraint*> _constraints) : constraints(_constraints) {}
  void operator()(GraphType& g, double t) {
    for (unsigned int i = 0; i<constraints.size(); ++i) {
      (*constraints[i])(g, t);
    }
  }
  private:
    std::vector<Constraint*> constraints;
};

/** Function that takes in 2 or 3 arguments and returns a CombinedConstraints object */
template <typename C1, typename C2, typename C3 = Constraint>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2, C3 c3 = Constraint()) {
  std::vector<Constraint*> a;
  C1* ptr = &c1;
  C2* ptr2 = &c2;
  C3* ptr3 = &c3;
  a.push_back(ptr);
  a.push_back(ptr2);
  a.push_back(ptr3);
  CombinedConstraints comb(a);
  return comb;
};

//--functionality_0
//--Passed all tests!
//--END

//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


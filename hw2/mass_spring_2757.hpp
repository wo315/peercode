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
static constexpr double c    = 0.01; 

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData>;
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

// template <typename G, typename F>
// double symp_euler_step(G& g, double t, double dt, F force) {
//   // Compute the t+dt position
//   for (auto it = g.node_begin(); it != g.node_end(); ++it) {
//     auto n = *it;

//     // Update the position of the node according to its velocity
//     // x^{n+1} = x^{n} + v^{n} * dt
//     n.position() += n.value().vel * dt;
//   }

//   //constraint(g, t);

//   // Compute the t+dt velocity
//   for (auto it = g.node_begin(); it != g.node_end(); ++it) {
//     auto n = *it;

//     // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
//     n.value().vel += force(n, t) * (dt / n.value().mass);
//   }

//   return t + dt;
// }

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

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

//unsigned number_edges = graph.num_edges();

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    int K = 100;
    // HW2 #1: YOUR CODE HERE
    //(void) n; (void) t; (void) grav;    // silence compiler warnings

     (void) t; // quiet compiler warning

    // nodes at (0, 0, 0) and (1, 0, 0) never move (no force is applied)
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0, 0, 0);
    }

    // otherwise there are spring and gravitational forces acting on the point
    Point force_grav   = n.value().mass * Point(0, 0, - grav); // gravitational force
    Point force_spring(0, 0, 0);

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      Edge e = *it;
      Point p = e.node1().position() - e.node2().position();
      force_spring -= (K- 1/norm(p)) * p;
    } 
    return force_grav + force_spring;
  }
};


//FORCE

//Base force functor
class BaseForce {
public:
  virtual Point operator()(Node n, double t)
  {
    (void) t; // quiet compiler warning
    (void) n; // quiet compiler warning
    return Point(0, 0, 0);
  }
};

//GravityForce functor which implements GravityForce
class GravityForce : public BaseForce {
public:
  Point operator()(Node n, double t)
  {
    (void) t; // quiet compiler warning
    return n.value().mass * Point(0, 0, - grav);
  }
};


//MassSpringForce functor which implements MassSpringForce
class MassSpringForce : public BaseForce {
public:
  Point operator()(Node n, double t) 
  {
    (void) t; // quiet compiler warning
    Point force_spring(0, 0, 0);
    int K = 100;

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it) 
    {
      Edge e = *it;
      Point p = e.node1().position() - e.node2().position();
      force_spring -= (K- 1/norm(p)) * p;
    }
    return force_spring;
  }
};

//DampingForce functor which implements DampingForce
class DampingForce : public BaseForce {
public:
  Point operator()(Node n, double t) 
  {
    (void) t; // quiet compiler warning
    return (- c * n.value().vel); // c is the damping constant
  }
};


 // CombinedForce to combine forces 
struct CombinedForces {
  Point operator()(Node n, double t) {
    Point f(0, 0, 0);

    for(const auto& x : forces_) {
      f = f + (*x)(n, t);
    }

    return f;
  }

  CombinedForces(std::vector<BaseForce*>& forces) : forces_(forces) {}
private:
  std::vector<BaseForce*> forces_;
};


CombinedForces make_combined_forces(GravityForce f1, MassSpringForce f2, DampingForce f3) {

  std::vector<BaseForce*> v {};

  v.push_back(&f1);
  v.push_back(&f2);
  v.push_back(&f3);

  CombinedForces comb(v);

  return comb;
}

// -------------------------- Constraint Functors ------------------------------------ //


class BaseConstraint {
public:
  virtual void operator()(GraphType& g, double t) {
    (void) g; // quiet compiler warning
    (void) t; // quiet compiler warning
  }
};


class PinConstraint : public BaseConstraint {
public:
  void operator()(GraphType& g, double t) {
    (void) t; //quiet compiler warning

    Node n = g.node(index1);
    Node m = g.node(index2);

    n.position() = Point(0, 0, 0);
    m.position() = Point(1, 0, 0);
  }

  PinConstraint(unsigned i, unsigned j) : index1(i), index2(j) {}
private:
  unsigned index1; 
  unsigned index2; 
};

 
class PlaneConstraint : public BaseConstraint {
public:
  void operator()(GraphType& g, double t)
  {
    (void) t; //quiet compiler warning

    for(auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = *it;
      if(dot(n.position(), Point(0, 0, 1)) < -0.75){
        n.value().vel.z = 0;      
        n.position().z  = -0.75;  
      }
    }
  }
};


class SphereConstraint : public BaseConstraint {
public:
  void operator()(GraphType& g, double t) {
    (void) t; //quiet compiler warning

    for(auto it = g.node_begin(); it != g.node_end(); ++it) {
      Node n = *it;

      Point direction = n.position() - center;
      double distance = norm(direction);

      if(distance < radius) {
        n.position()   = center + (radius / distance) * direction;
        n.value().vel -= ((dot(n.value().vel, direction)) / normSq(direction)) * direction;
      }
    }
  }

  SphereConstraint() : center{0.5, 0.5, -0.5}, radius{0.15} {}
private:
    Point center;
    double radius;
};


struct CombinedConstraints {
  void operator()(GraphType& g, double t) {
    (void) t; //quiet compiler warning
    for(const auto & x : constraints)
    {
      (*x)(g, t);
    }
  }
  CombinedConstraints(std::vector<BaseConstraint*>& v) : constraints(v) {}
private:
  std::vector<BaseConstraint*> constraints;
};


CombinedConstraints make_combined_constraints(PinConstraint c1, PlaneConstraint c2,
                                              SphereConstraint c3) 
{
  std::vector<BaseConstraint*> v;
  v.push_back(&c1);
  v.push_back(&c2);
  v.push_back(&c3);

  CombinedConstraints con(v);
  return con;
}

//--functionality_1
//--The figure of mass_spring is not correct.
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


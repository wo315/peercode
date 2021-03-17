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
  Point initial_position;
  NodeData() : vel(0), mass(1), initial_position() {}
};

struct EdgeData {
  double L;
	double K;
  EdgeData() : L(), K() {}
  EdgeData(double L_, double K_) : L(L_), K(K_) {}
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
    // HW2 #1: YOUR CODE HERE
    (void)t;
    if(n.position()==Point(0, 0, 0) || n.position()==Point(1, 0, 0))
      return Point(0, 0, 0);
    else{
      Point p(0, 0, -grav*n.value().mass);
      for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        p += -(*it).value().K*((*it).value().L-(*it).length())/(*it).length()*((*it).node2().position()-(*it).node1().position());
      }
      return p;
    }

  }
};

class Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void)n; (void)t;
      return Point(0, 0, 0);
    }
};

class GravityForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void)t;
      return Point(0, 0, -grav*n.value().mass);
    }
};

class MassSpringForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void)t;
      Point p(0, 0, 0);
      for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
        p += -(*it).value().K*((*it).value().L-(*it).length())/(*it).length()*((*it).node2().position()-(*it).node1().position());
      }
      return p;
    }
};

class DampingForce: public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void)t;
      double c = 1.0;
      return -c*n.value().vel;
    }
};

class CombinedForces: public Force {
  std::vector<Force*> forces;
  public:
    CombinedForces(std::vector<Force*> forces_): forces(forces_) {}
    virtual Point operator()(Node n, double t) {
      Force* ptr = nullptr;
      Point p(0, 0, 0);
      for(auto force_it = forces.begin(); force_it != forces.end(); force_it++){
        ptr = *force_it;
        p += ptr->operator()(n, t);
      }
      return p;
    }
};

CombinedForces make_combined_force(Force&& f1, Force&& f2) {
  std::vector<Force*> forces{&f1, &f2};
  return CombinedForces(forces);
}

CombinedForces make_combined_force(Force&& f1, Force&& f2, Force&& f3) {
  std::vector<Force*> forces{&f1, &f2, &f3};
  return CombinedForces(forces);
}

class Constraint {
  public:
    virtual void operator()(GraphType& graph, double t) {
      (void)graph; (void)t;
    }
};

class PinConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& graph, double t) {
      (void)t;
      for(auto it=graph.node_begin(); it!=graph.node_end(); ++it){
        if((*it).value().initial_position==Point(0, 0, 0) || (*it).value().initial_position==Point(1, 0, 0)){
          (*it).position() = (*it).value().initial_position;
        }
      }
    }
};

class PlaneConstraint: public Constraint {
  public:
    virtual void operator()(GraphType& graph, double t) {
      (void)t;
      for(auto it=graph.node_begin(); it!=graph.node_end(); ++it){
        if((*it).position().z<-0.75){
          (*it).position().z = -0.75; // orthogonal projection
          (*it).value().vel.z = 0.; // Set the z-component of the velocity to 0
        }
      }
    }
};

class SphereConstraint: public Constraint {
  Point center;
  double radius;
  public:
    SphereConstraint(): center(Point(0.5, 0.5, -0.5)), radius(0.15) {}
    virtual void operator()(GraphType& graph, double t) {
      (void)t;
      for(auto it=graph.node_begin(); it!=graph.node_end(); ++it){
        Point v = (*it).position()-center; //vector from center to position
        double d = norm(v);
        if(d<radius){
          Point R = v/d; // unit direction vector
          (*it).position() = center + radius * R; // projection on the sphere
          // Applies projection formula to the velocity
          (*it).value().vel = (*it).value().vel - dot((*it).value().vel, R) * R;
        }
      }
    }
};

class TearConstraint: public Constraint {
  Point center;
  double radius;
  public:
    TearConstraint(): center(Point(0.5, 0.5, -0.5)), radius(0.15) {}
    virtual void operator()(GraphType& graph, double t) {
      (void)t;
      auto it=graph.node_begin();
      do {
        Point v = (*it).position()-center; //vector from center to position
        double d = norm(v);
        if(d<radius){
          it = graph.remove_node(it);
        }
        else{
          ++it;
        }
      } while(it!=graph.node_end());
    }
};

class CombinedConstraints: public Constraint {
  std::vector<Constraint*> constraints;
  public:
    CombinedConstraints(std::vector<Constraint*> constraints_): constraints(constraints_) {}
    virtual void operator()(GraphType& graph, double t) {
      Constraint* ptr = nullptr;
      for(auto constraint_it = constraints.begin(); constraint_it != constraints.end(); constraint_it++){
        ptr = *constraint_it;
        ptr->operator()(graph, t);
      }
    }
};

CombinedConstraints make_combined_constraints(Constraint&& c1, Constraint&& c2) {
  std::vector<Constraint*> constraints{&c1, &c2};
  return CombinedConstraints(constraints);
}

CombinedConstraints make_combined_constraints(Constraint&& c1, Constraint&& c2, Constraint&& c3) {
  std::vector<Constraint*> constraints{&c1, &c2, &c3};
  return CombinedConstraints(constraints);
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

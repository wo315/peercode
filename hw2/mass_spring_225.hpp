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

using namespace std;

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData{
  double K;
  double L;
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
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      continue;
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

template <typename G, typename F, typename C>
double symp_euler_step (G& g, double t, double dt , F force , C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  
  // Apply the constraints
  constraint(g,t);
  
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
  Point operator()(NODE n, double t){
	  
    (void) t;   // silence compiler warnings
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
  
    Point f_spring=Point(0, 0, 0);
    Point xi=n.position();
    double mi=n.value().mass;
	
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      double Kij=(*it).value().K;
      double Lij=(*it).value().L;
      Point xj=(*it).node2().position();
      f_spring+=(-Kij*(xi-xj)/norm(xi-xj)*(norm(xi-xj)-Lij)); 
    }
	
    return f_spring+mi*Point(0,0,-grav);
  }
};

/** Parent class representing a force*/
class Force{
public:
  virtual Point operator()(Node n, double t){
    (void) t; (void) n;
    return Point(0,0,0);
  }
};

/** Inherited class from Force representing GravityForce */
class GravityForce:public Force{
public:
  Point operator()(Node n, double t){
    (void) t;
    return n.value().mass*Point(0,0,-grav);
  }
}; 

/** Inherited class from Force representing MassSpringForce */
class MassSpringForce:public Force{
public:
  Point operator()(Node n, double t){
    (void) t;
    Point f_spring=Point(0, 0, 0);
    Point xi=n.position();
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      double Kij=(*it).value().K;
      double Lij=(*it).value().L;
      Point xj=(*it).node2().position();
      f_spring+=(-Kij*(xi-xj)/norm(xi-xj)*(norm(xi-xj)-Lij)); 
    }
    return f_spring;
  }
};

/** Inherited class from Force representing DampingForce */
class DampingForce: public Force{
public:
  double c;
  DampingForce(): c(0.0){}
  DampingForce(const double damping_const):c(damping_const){}
  Point operator()(Node n, double t){
    (void) t;
    return -c*n.value().vel;
  }
}; 

/** Functor that returns the combination of a vector of forces */
struct CombinedForce{
  vector<Force*> forces;
  CombinedForce(vector<Force*> f): forces(f){};
  Point operator()(Node n, double t) {
    (void) t;
    Point combined_force=Point(0,0,0);
    for(auto i=forces.begin(); i!=forces.end(); ++i){
      combined_force += (*(*i))(n, t);
    }
    return combined_force;
  }
};

/** Function that takes in 2 forces and returns the CombinedForce */
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2){
  vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  return CombinedForce(forces);
}

/** Function that takes in 3 forces and returns the CombinedForce */
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3){
  vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3);
  return CombinedForce(forces);
}

/** Parent class representing a constraint */
class Constraint{
public:
  virtual void operator()(GraphType& g, double t){
    (void) t; (void) g;
  }
};

/** Inherited class from Constraint representing a pin constraint */
class PinConstraint:public Constraint{
public:
  double dt;
  
  PinConstraint(double d_t):dt(d_t){}
  void operator()(GraphType& g, double t){
    (void) t;
    double eps = 1e-99;
    for(auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
  
      Point prev = n.position()-n.value().vel*dt;
      if(norm(prev-Point(0, 0, 0))<eps){
        n.position()=Point(0, 0, 0);
        n.value().vel=Point(0, 0, 0);
      }
      if(norm(prev-Point(1, 0, 0))<eps){
        n.position()=Point(1, 0, 0);
        n.value().vel=Point(0, 0, 0);
      }   
    }
  }
}; 

/** Inherited class from Constraint representing a plane constraint */
class PlaneConstraint:public Constraint{
public:
  double cz=-0.75;
  
  void operator()(GraphType& g, double t){
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(n.position().z<cz){
        n.position().z=cz;
        n.value().vel.z=0;
      }
    }
  }
}; 

/** Inherited class from Constraint representing a sphere constraint */
class SphereConstraint: public Constraint{
public:
  Point c=Point(0.5,0.5,-0.5);
  double r=0.15;
  
  void operator()(GraphType& g, double t){
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(norm(n.position()-c)<r){
        Point Ri=(n.position()-c)/norm(n.position()-c);
        n.position()=c+Ri*r;
        n.value().vel-=dot(n.value().vel, Ri)*Ri;
      }
    }
  }
}; 

/** Inherited class from Constraint representing tear constraint */
class TearConstraint:public Constraint{
public:
  Point c=Point(0.5,0.5,-0.5);
  double r=0.15;
  void operator()(GraphType& g, double t){
    (void) t;
//--design_1
//--if (*it) is removed, then after removal it points to the next node. "++it" should not be applied.
//--You should use while loop instead of for loop. Consider using if{...} else{++it}.
//--START
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      auto n = *it;
      if(norm(n.position()-c)<r){
        g.remove_node(n);
      }
    }
//--END
  }
}; 


/** Functor that returns the combination of a vector of constraints */
struct CombinedConstraints{
  vector<Constraint*> constraints;
  CombinedConstraints(vector<Constraint*> c): constraints(c){};
  void operator()(GraphType& g, double t) {
    (void) t;
    for(auto i=constraints.begin(); i!=constraints.end(); ++i){
      (*(*i))(g, t);
    }
  }
};

/** Function that takes in 2 constraints and returns the CombinedConstraints */
template <typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2){
  vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  return CombinedConstraints(constraints);
}

/** Function that takes in 3 constraints and returns the CombinedConstraints */
template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2, C3 c3){
  vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);
  return CombinedConstraints(constraints);
}

/** Function that takes in 4 constraints and returns the CombinedConstraints */
template <typename C1, typename C2, typename C3, typename C4>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2, C3 c3, C4 c4){
  vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);
  constraints.push_back(&c4);
  return CombinedConstraints(constraints);
}

//--functionality_0
//--Passed all tests!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good


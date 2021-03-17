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
using size_type = unsigned;


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

/*Change a graph's nodes according to a step of the symplectic Euler
*    method with the given node force under given constraints
*/
template <typename G, typename F, typename C>
double symp_euler_step (G& g, double t, double dt , F force , C constraint ) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
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
  Point operator()(NODE n, double t) {
    (void) t;
    if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
      return Point(0, 0, 0);
    Point fs=Point(0, 0, 0);
    Point xi=n.position();
    double mi=n.value().mass;
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      double Kij=(*it).value().K;
      double Lij=(*it).value().L;
      Point xj=(*it).node2().position();
      fs+=(-Kij*(xi-xj)/norm(xi-xj)*(norm(xi-xj)-Lij));
    }
    return fs+mi*Point(0,0,-grav);
  }
};

/*parent class that represents any force
*/
class Force {
public:
  virtual Point operator()(Node n, double t) {
    (void) t;
    (void) n;
    //std::cout<<"call force"<<std::endl;
    return Point(0,0,0);
  }
};

/*inherited class from Force that represents gravity force
*/
class GravityForce: public Force {
public:
  Point operator()(Node n, double t) {
    (void) t;
    //std::cout<<"call gforce"<<std::endl;
    return n.value().mass*Point(0,0,-grav);
   }
};

/*inherited class from Force that represents mass spring force
*/
class MassSpringForce: public Force {
public:
  Point operator()(Node n, double t) {
    (void) t;
    //std::cout<<"call msforce"<<std::endl;
    Point fs=Point(0, 0, 0);
    Point xi=n.position();
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it){
      double Kij=(*it).value().K;
      double Lij=(*it).value().L;
      Point xj=(*it).node2().position();
      fs+=(-Kij*(xi-xj)/norm(xi-xj)*(norm(xi-xj)-Lij));
    }
    return fs;
   }
};

/*inherited class from Force that represents damping force
*/
class DampingForce: public Force {
public:
  double c;
  DampingForce(): c(0.0) {}
  DampingForce(const double dc) : c(dc) {}
  Point operator()(Node n, double t) {
    (void) t;
    //std::cout<<"call dforce"<<std::endl;
    return -c*n.value().vel;
   }
};

/*functor that takes in a vector of Force and combine them
*/
struct CombinedForce{
  std::vector<Force*> forces;
  CombinedForce(std::vector<Force*> v): forces(v) {};
  Point operator()(Node n, double t) {
    (void) t;
    Point cf=Point(0,0,0);
    for(auto i=forces.begin(); i!=forces.end(); ++i){
      cf+=(*(*i))(n, t);
    }
    return cf;
  }
};

/*function that takes in 3 forces and returns a CombinedForce functor
*/
//--design_1
//--argument objects are destroyed after function finishes
//--START
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 force1, F2 force2, F3 force3){
  std::vector<Force*> f;
  f.push_back(&force1);
  f.push_back(&force2);
  f.push_back(&force3);
  return CombinedForce(f);
}
//--END

/*function that takes in 2 forces and returns a CombinedForce functor
*/
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 force1, F2 force2){
  std::vector<Force*> f;
  f.push_back(&force1);
  f.push_back(&force2);
  return CombinedForce(f);
}

/*parent class that represents any constraint
*/
class Constraint {
public:
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    (void) g;
  }
};

/*inherited class from Constraint that represents pin constraint
*/
class PinConstraint: public Constraint {
public:
  double dt;
  PinConstraint(double d): dt(d){}
  void operator()(GraphType& g, double t) {
    (void) t;
    for(auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      Point prev_pos=n.position()-n.value().vel*dt;
      //--design_1
      //--unmoved points are wrongly identified as pinned
      //--START
      if(norm(prev_pos-Point(0, 0, 0))<1e-9){
        n.position()=Point(0,0,0);
      }
      if(norm(prev_pos-Point(1, 0, 0))<1e-9){
        n.position()=Point(1,0,0);
      }
      //--END
    }
  }
};


/*inherited class from Constraint that represents plane constraint
*/
class PlaneConstraint: public Constraint {
public:
  double cz=-0.75;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if(n.position().z<cz){
        n.position().z=cz;
        n.value().vel.z=0;
      }
    }
  }
};

/*inherited class from Constraint that represents sphere constraint
*/
class SphereConstraint: public Constraint {
public:
  Point c=Point(0.5,0.5,-0.5);
  double r=0.15;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if(norm(n.position()-c)<r){
        Point ri=(n.position()-c)/norm(n.position()-c);
        n.position()=c+r*ri;
        n.value().vel-=dot(n.value().vel, ri)*ri;
      }
    }
  }
};

/*inherited class from Constraint that represents tear constraint
*/
class TearConstraint: public Constraint {
public:
  Point c=Point(0.5,0.5,-0.5);
  double r=0.15;
  void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if(norm(n.position()-c)<r){
        g.remove_node(n);
      }
    }
  }
};

/*functor that takes in a vector of Constraint and combine them
*/
struct CombinedConstraints{
  std::vector<Constraint*> cons;
  CombinedConstraints(std::vector<Constraint*> v): cons(v) {};
  void operator()(GraphType& g, double t) {
    (void) t;
    for(auto i=cons.begin(); i!=cons.end(); ++i){
      (*(*i))(g, t);
    }
  }
};

/*function that takes in 3 constraints and returns a CombinedConstraints functor
*/
template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 con1, C2 con2, C3 con3){
  std::vector<Constraint*> c;
  c.push_back(&con1);
  c.push_back(&con2);
  c.push_back(&con3);
  return CombinedConstraints(c);
}

/*function that takes in 2 constraints and returns a CombinedConstraints functor
*/
template <typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 con1, C2 con2){
  std::vector<Constraint*> c;
  c.push_back(&con1);
  c.push_back(&con2);
  return CombinedConstraints(c);
}

/*function that takes in 4 constraints and returns a CombinedConstraints functor
*/
template <typename C1, typename C2, typename C3, typename C4>
CombinedConstraints make_combined_constraints(C1 con1, C2 con2, C3 con3, C4 con4){
  std::vector<Constraint*> c;
  c.push_back(&con1);
  c.push_back(&con2);
  c.push_back(&con3);
  c.push_back(&con4);
  return CombinedConstraints(c);
}

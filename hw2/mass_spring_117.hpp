/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

///////////////////////////////////////////////////////////
#include <iostream>
#include <string>
#include <vector>
///////////////////////////////////////////////////////////

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
  Point init_pos;  //< Initial position of the node
  NodeData() : vel(0), mass(1) {}
};


///////////////////////////////////////////////////////////
// HW2 #1: YOUR CODE HERE
/** Custom structure of data to store with Edges */
struct EdgeData {
  double spring_constant;       //< Edge spring constant
  double rest_length;           //< Edge rest length
  EdgeData() : spring_constant(100.0), rest_length(0.01) {}
};
///////////////////////////////////////////////////////////


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>; // modified for part 4
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
    ///////////////////////////////////////////////////////////
    if ((n.position() == Point(0,0,0)) || (n.position() == Point(1,0,0))){
      // Corner points which shouldn't move
      return Point(0,0,0);
    }
    else{
      Point force = n.value().mass * Point(0, 0, -grav); // gravity force

      //mass-spring force
      double edge_length;
      double K;
      double L;
      NODE neigbour;
      Edge edge;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        edge = *it;
        edge_length = edge.length();
        K = edge.value().spring_constant;
        L = edge.value().rest_length;
        neigbour = edge.node2();
        force += (-1.0*K*(edge_length-L)/edge_length) * \
                 (n.position()-neigbour.position());
      }
      (void) t; //quiet compiler warning
      return force;
    }
    ///////////////////////////////////////////////////////////
  }

  ///////////////////////////////////////////////////////////
  // NOW USELESS AS WE HAVE EDGE VALUES SINCE PART 4
  // HW2 #1: YOUR CODE HERE
  // spring constant
  //const double K = 100;
  // initial length constant
  //const double L = 0.001;
  ///////////////////////////////////////////////////////////

};



// HW2 #5: YOUR CODE HERE
///////////////////////////////////////////////////////////

/** interface for all kind of forces */
class Force { 
 public:

  /** Return the force applying to @a n at time @a t. */
  virtual Point operator()(Node n, double t) {
    (void) t; //quiet compiler warning
    (void) n; //quiet compiler warning
    return Point(0,0,0);
  }

  /** Destructor */
  virtual ~Force() {}

};

/** Implements gravity force, derived from class Force */
class GravityForce : public Force{
 public:
  virtual Point operator()(Node n, double t) {
    Point force = n.value().mass * Point(0, 0, -grav); // gravity force
    (void) t; //quiet compiler warning
    return force;
  }
};

/** Implements mass-spring force, derived from class Force */
class MassSpringForce : public Force{
 public:
  virtual Point operator()(Node n, double t) {
    Point force = Point(0,0,0);//mass-spring force
    double edge_length;
    double K;
    double L;
    Node neigbour;
    Edge edge;
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      edge = *it;
      edge_length = edge.length();
      K = edge.value().spring_constant;
      L = edge.value().rest_length;
      neigbour = edge.node2();
      force += (-1.0*K*(edge_length-L)/edge_length) * \
               (n.position()-neigbour.position());
    }
    (void) t; //quiet compiler warning
    return force;
  }
};

/** Implements damping force, derived from class Force */
class DampingForce : public Force{
 private:
  double damping_constant;
 
 public:
 /** Constructor */
  DampingForce(const double c=0.005) : damping_constant(c) {  }

  virtual Point operator()(Node n, double t) {
    Point force = - damping_constant * n.value().vel; // damping force
    (void) t; //quiet compiler warning
    return force;
  }
};


/** Combined force functor to add forces */
class CombinedForce: public Force {
 private:
  std::vector<Force*> forces;

 public:
  /** Return the combined force applying to @a n at time @a t.
  */
  CombinedForce(std::vector<Force*> forces_):forces(forces_) {}

  virtual Point operator()(Node n, double t) {
    Point force = Point(0.0, 0.0, 0.0);
    for(auto it=forces.begin(); it!=forces.end(); ++it){
      force += (*(*it))(n,t);
    }
    return force;
  }

  ~CombinedForce () {
    forces.clear(); 
  }
};



/** Combines the force f1 and f2, which must inherit from Force */
template <typename T1, typename T2>
CombinedForce make_combined_force(T1 f1, T2 f2){
  std::vector<Force*> forces_vec;
  forces_vec.push_back(&f1);
  forces_vec.push_back(&f2);
  return CombinedForce(forces_vec);
}

/** Combines the force f1 f2, and f3, which must inherit from Force */
template <typename T1, typename T2, typename T3>
CombinedForce make_combined_force(T1 f1, T2 f2, T3 f3){
  std::vector<Force*> forces_vec;
  forces_vec.push_back(&f1);
  forces_vec.push_back(&f2);
  forces_vec.push_back(&f3);
  return CombinedForce(forces_vec);
}

///////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////
// HW2 #6: YOUR CODE HERE
/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     force  Function object defining the constraint per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a constraint must update the node n (position and value) to 
 *           respect the purpose of the constrint
 */
template  <typename G, typename F, typename C>
double  symp_euler_step(G& g, double t, double dt , F force , C constraint){
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the system constraints
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}



/** interface for all kind of constraints */
class Constraint { 
 public:

  /** Modify the Graph @a g so that it respect the defined constraint
   *  at time @a t. */
  virtual void operator()(GraphType& graph, double t) {
    (void) t; //quiet compiler warning
    (void) graph; //quiet compiler warning
    return;
  }

  /** Destructor */
  virtual ~Constraint() {}

};

/** Implements pin constrainy, derived from class Constraint */
class PinConstraint : public Constraint{
 public:
  virtual void operator()(GraphType& graph, double t) {
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      if (((*it).value().init_pos == Point(0,0,0)) \
         || ((*it).value().init_pos  == Point(1,0,0))){
        // these points shouldn't move
        (*it).position() = (*it).value().init_pos;
        (*it).value().vel = Point(0.0, 0.0, 0.0);
      }     
    }
    (void) t; //quiet compiler warning
    return;
  }
};


/** Implements plane constraint, derived from class Constraint */
class PlaneConstraint : public Constraint{
 private:
  double z_plan_value;
 
 public:
  virtual void operator()(GraphType& graph, double t) {
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      if ((*it).position().z < z_plan_value){
        // these points shouldn't move
        (*it).position().z = z_plan_value;
        (*it).value().vel.z = 0.0;
      }     
    }
    (void) t; //quiet compiler warning
    return;
  }

  /** Constructor */
  PlaneConstraint(double z_val=-0.75): z_plan_value(z_val) {}
};



/** Implements sphere constraint, derived from class Constraint */
class SphereConstraint : public Constraint{
 private:
  Point c;
  double r;
 
 public:
  virtual void operator()(GraphType& graph, double t) {
    Point R_i;
    Point old_v;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      if (norm((*it).position() - c) < r){
        // these points are in the sphere
        R_i = (1.0/norm((*it).position()-c))*((*it).position()-c);
        old_v = (*it).value().vel;
        (*it).position() = c+r*R_i;
        (*it).value().vel = old_v - dot(old_v,R_i)*R_i;
      }     
    }
    (void) t; //quiet compiler warning
    return;
  }

  /** Constructor */
  SphereConstraint(Point center=Point(0.5,0.5,-0.5), double radius=0.15):
    c(center), r(radius) {}
};


/** Combined constraint functor to add different constraints */
class CombinedConstraints: public Constraint {
 private:
  std::vector<Constraint*> constraints;

 public:
  /** Return the combined constraint that the graph should respact.
  */
  CombinedConstraints(std::vector<Constraint*> constraints_) : 
      constraints(constraints_) {}

  virtual void operator()(GraphType& graph, double t) {
    for(auto it = constraints.begin(); it!= constraints.end(); ++it){
      (*(*it))(graph, t);
    }
    return;
  }

  ~CombinedConstraints () {
    constraints.clear(); 
  }
};



/** Combines the constraints c1 and c2, which must inherit from Constraint */
template <typename T1, typename T2>
CombinedConstraints make_combined_constraint(T1 c1, T2 c2){
  std::vector<Constraint*> constraints_vec;
  constraints_vec.push_back(&c1);
  constraints_vec.push_back(&c2);
  return CombinedConstraints(constraints_vec);
}

/** Combines the constraints c1 c2 and c3, which must inherit from Constraint */
template <typename T1, typename T2, typename T3>
CombinedConstraints make_combined_constraint(T1 c1, T2 c2, T3 c3){
  std::vector<Constraint*> constraints_vec;
  constraints_vec.push_back(&c1);
  constraints_vec.push_back(&c2);
  constraints_vec.push_back(&c3);
  return CombinedConstraints(constraints_vec);
}



/** Implements sphere constraint + disaparition, derived from class Constraint*/
class TearConstraint : public Constraint{
 private:
  Point c;
  double r;
 
 public:
  virtual void operator()(GraphType& graph, double t) {
    Point R_i;
    Point old_v;
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it){
      if (norm((*it).position() - c) < r){
        // these points are in the sphere
        graph.remove_node(it);
      }     
    }
    (void) t; //quiet compiler warning
    return;
  }

  /** Constructor */
  TearConstraint(Point center=Point(0.5,0.5,-0.5), double radius=0.15):
    c(center), r(radius) {}
};


///////////////////////////////////////////////////////////

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




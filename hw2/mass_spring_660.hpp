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
  bool valid;      //< Node validity
  NodeData() : vel(0), mass(1), valid(1) {}
};

// /** Custom structure of data to store with Edges */
// struct EdgeData {
//   double k;       //< Edge spring constant
//   double l;       //< Edgle length
//   EdgeData() : k(100), l(0.001) {}
// };

// Define the Graph type
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

// Force parent class
class Force
{
public:
  Force() {}

  /** Functor to return Force as a Point
  */
  virtual Point operator()(Node n, double t)
  {
    (void) n;
    (void) t;
    return Point(0.0,0.0,0.0);
  }
};

// Constraint parent class
class Constraint
{
public:
  Constraint() {}

  /** Functor to modify a graph reference
   * applying a defined constraint
  */
  virtual void operator()(GraphType& graph, double t)
  {
    (void) graph;
    (void) t;
  }
};

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
template <typename G, typename F, typename C = Constraint>
double symp_euler_step(G& g, double t, double dt, F force, C constraint=Constraint()) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // apply constraint before force
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (double)n.value().valid * (dt / n.value().mass);
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
    //(void) n; (void) t; (void) grav;    // silence compiler warnings
    (void) t;
    //return Point(0);

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point (0 ,0 ,0);

    else {
      //double K = 100.0; // defined as a const for this problem.
      //double L = 0.001; // uniform for all edges.
      double K,L;
      Point fg = n.value().mass * Point(0.0,0.0, -grav);
      Node temp;
      Point temp_pos;
      
      Point fk = Point(0.0,0.0,0.0);
      for(GraphType::IncidentIterator ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        temp = (*ii).node2();
        K = 100.0;
        L = (*ii).value();
        temp_pos = n.position() - temp.position();
        double two_norm = norm(temp_pos);
        double scale = (K * (two_norm - L)) / two_norm;
        fk -= temp_pos * scale;
      }
      return fg + fk;
    }
  }
};

/** GravityForce derives from Force and implements
 * force due to gravity as -m*g.
 */
class GravityForce : public Force
{
public:
  GravityForce() {}


  Point operator()(Node n, double t)
  {
    (void)t;
    return n.value().mass * Point(0.0,0.0, -grav);
  }
};

/** MassSpringForce derives from Force and implements
 * spring force iterating over all the neighbors
 * and aggregating their contribution.
 */
class MassSpringForce : public Force
{
public:
  MassSpringForce() {}

  Point operator()(Node n, double t)
  {
    (void) t;
    Node temp;
    Point temp_pos;
    double K,L;
    Point fk = Point(0.0,0.0,0.0);

    for(GraphType::IncidentIterator ii = n.edge_begin(); ii != n.edge_end(); ++ii) {
        temp = (*ii).node2();
        K = 100.0;
        L = (*ii).value();
        temp_pos = n.position() - temp.position();
        double two_norm = norm(temp_pos);
        double scale = (K * (two_norm - L)) / two_norm;
        fk -= temp_pos * scale;
      }

    return fk;
  }
};

/** DampingForce derives from Force and implements
 * dampling force
 */
class DampingForce : public Force
{
public:
  DampingForce(double c = 0.001) {_c = c;}

  Point operator()(Node n, double t)
  {
    (void)t;
    return -1.0 * n.value().vel * _c;
  }

  private:
   double _c;
};

/** Functor to combine multiple forces
 *  Takes a vector of pointers to forces.
 */
struct CombinedForce {

  CombinedForce(std::vector<Force*> forces) :
    _forces(forces) {}

  Point operator()(Node n, double t)
  {
    Point f = Point(0.0);
    for(auto i=_forces.begin(); i!=_forces.end(); ++i) {
      f  = f + (*(*i))(n,t);
    }
    return f;
  }

  private:
  std::vector<Force*> _forces;
};

/** Creates a CombinedForce Object
 *  Takes 2 or 3 Force objects / children objects
 */
template <typename Force1, typename Force2, typename Force3 = Force>
CombinedForce make_combined_force(Force1 f1, Force2 f2, Force3 f3 = Force()) {

  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3);

  return CombinedForce(forces);
}

/** Pins some points. Derives from Constraint
 */
class PinConstraint : public Constraint
{
  public:
    PinConstraint() {
      // semi configurable.
      _pins.push_back(Point(0,0,0));
      _pins.push_back(Point(1,0,0));
    }

  void operator()(GraphType& graph, double t) 
  {
    (void) t;
    // iterate through nodes
    for(GraphType::NodeIterator niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
      // iterate through pin points
      for (auto iter = _pins.begin(); iter != _pins.end(); ++iter) {
        // if node's position is same as pin point
        if((*niter).position() == *iter) {
          // set its value to 0 - invalidate it so force won't act on it.
          (*niter).value().valid = 0;
        }
      }
  }
  }
  private:
   std::vector<Point> _pins;
};

/** Modifies nodes wrt a plane. Derives from Constraint
 */
class PlaneConstraint : public Constraint
{
  public:
    PlaneConstraint() {
      // semi configurable.
      z = -0.75;
      fix = Point(0,0,1);

    }

  void operator()(GraphType& graph, double t) 
  {
    (void) t;
    // iterate through nodes
    for(GraphType::NodeIterator niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
      // iterate through pin points
        // if node's position is same as pin point
        double dot_prod = dot( (*niter).position(), fix);
        if( dot_prod < z) {
        
          (*niter).position()[2] = z;
          // set velocity z to 0
          (*niter).value().vel[2] = 0.0;
        }
    }
  }
  private:
   double z;
   Point fix;
};

/** Modifies nodes wrt a sphere. Derives from Constraint
 */
class SphereConstraint : public Constraint
{
  public:
    SphereConstraint() {
      centre = Point(0.5,0.5,-0.5);
      r = 0.15;
    }

  void operator()(GraphType& graph, double t) 
  {
    (void) t;
    // iterate through nodes
    for(GraphType::NodeIterator niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
        Point diff = (*niter).position() - centre;
        Point ri = diff / norm(diff);
        if( norm(diff) < r) {

          (*niter).position() = centre + r*ri;
          Point vi = (*niter).value().vel;
          (*niter).value().vel = vi - ri*dot(vi,ri);
        }
    }
  }
  private:
   double r;
   Point centre;
};

/** Modifies nodes wrt a sphere. 
 * Removes nodes and edges.
 * Derives from Constraint
 */
class TearConstraint : public Constraint
{
  public:
    TearConstraint() {
      centre = Point(0.5,0.5,-0.5);
      r = 0.15;
    }

  void operator()(GraphType& graph, double t) 
  {
    (void) t;
    // iterate through nodes
//--design_1
//--You're missing nodes by removing with a for loop that increments each time
//--START
    for(GraphType::NodeIterator niter = graph.node_begin(); niter != graph.node_end(); ++niter) {
        Point diff = (*niter).position() - centre;
        if( norm(diff) < r) {
          niter = graph.remove_node(niter);
        }
    }
//--END
  }
  private:
   double r;
   Point centre;
};

/** Functor to combine multiple constraints
 *  Takes a vector of pointers to constraints.
 */
struct CombinedConstraint {

  CombinedConstraint(std::vector<Constraint*> constraints) :
    _constraints(constraints) {}

  void operator()(GraphType& g, double t)
  {
    for(auto i=_constraints.begin(); i!=_constraints.end(); ++i) {
      (*(*i))(g,t);
    }
  }

  private:
  std::vector<Constraint*> _constraints;
};

/** Creates a CombinedConstraint Object
 *  Takes 2 or 3 Constraint objects / children objects
 */
template <typename Con1, typename Con2, typename Con3 = Constraint>
CombinedConstraint make_combined_constraint(Con1 c1, Con2 c2, Con3 c3 = Constraint()) {

  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);

  return CombinedConstraint(constraints);
}

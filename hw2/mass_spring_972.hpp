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
static constexpr double K = 1.0;
static constexpr double C = 0.1;
double eps = 0;

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
double symp_euler_step(G& g, double t, double dt, F force, C constraint){
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // apply constaint
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel = n.value().vel +  force(n, t) * (dt / n.value().mass);
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
  template <typename N>
  Point operator()(N n, double t) {
    static_cast<void>(t);
    Point force(0,0,0);
    Edge e_ij;
    Point delta;
    double norm_delta;

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) 
      return force;

    // initialize force to (0,0,-mg)
    force = force + Point(0, 0, -grav) * n.value().mass;

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      e_ij = *it;

      // preserve the order; always subtract from x_i, current node
      delta = e_ij.node1().position() - e_ij.node2().position();
      if (e_ij.node1() != n){
        delta = delta * (-1);
      }

      norm_delta = norm(delta);
      force = force - K * delta *( norm_delta  - e_ij.value()) / norm_delta ;
    }

    return force;
  }
};


/**
  Parent class for forces
*/
struct Force{
  virtual Point operator()(Node n, double t) = 0;
};


/**
  GravityForce applied to node, inherits from Force
*/
struct GravityForce: Force
{
  virtual Point operator()(Node n, double t){
    static_cast<void>(t);
    Point force = Point(0, 0, -grav) * (n.value().mass);

    return force;

  }

  virtual ~GravityForce(){};
  
};


/**
  MassSpringForce applied to node, inherits from Force
*/

//--functionality_1
//--Your spring forces are not behaving as expected (no change of signs)
//--END
struct MassSpringForce: Force
{
  virtual Point operator()(Node n, double t){
    static_cast<void>(t);

    Point force(0,0,0);
    Edge e_ij;
    Point delta;
    double norm_delta;

//--style_1
//--This shouldn't be here now that you have a PinConstraint
//--START
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) 
      return force;
//--END

    // initialize force to (0,0,-mg)
    force = force + Point(0, 0, -grav) * n.value().mass;

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
      e_ij = *it;

      // preserve the order; always subtract from x_i, current node
      delta = (e_ij.node1()).position() - (e_ij.node2()).position();
      if (e_ij.node1() != n){
        delta = delta * (-1);
      }

      norm_delta = norm(delta);
      force = force - K * delta *( norm_delta  - e_ij.value()) / norm_delta ;
    }

    return force;

  }

  virtual ~MassSpringForce(){};
  
};


/**
  DampingForce applied to node, inherits from Force
*/
struct DampingForce: Force
{
  virtual Point operator()(Node n, double t){
    static_cast<void>(t);
    Point force = - C * n.value().vel;

    return force;

  }

  virtual ~DampingForce(){};
  
};


/**
  Functor that combines all forces
*/
struct CombinedForce
{
  Point operator()(Node n, double t){
      static_cast<void>(t);
      Point total_forces(0,0,0);
      for (unsigned i = 0; i < forces_.size(); ++i)
      {
          total_forces = total_forces + (*forces_.at(i))(n,t);
      }
      return total_forces;
  }

  CombinedForce(std::vector<Force*> forces) : forces_(forces) {}

  std::vector<Force*> forces_;
};


template<typename Force1, typename Force2>
CombinedForce make_combined_force(Force1 f1, Force2 f2) {
    std::vector<Force*> forces {&f1, &f2};
    
    return CombinedForce(forces);
};

template<typename Force1, typename Force2, typename Force3>
CombinedForce make_combined_force(Force1 f1, Force2 f2, Force3 f3) {
    std::vector<Force*> forces {&f1, &f2, &f3};
    
    return CombinedForce(forces);
};


/**
  Parent class for constraints
*/
struct Constraint{
  virtual void operator()(GraphType& g, double t) = 0;
};


/**
  PinConstraint applied to node, inherits from Constraint
*/
struct PinConstraint: Constraint
{
  virtual void operator()(GraphType& g, double t){
    Node n;
    Point pos;
    static_cast<void>(t);

    for (auto it = g.node_begin(); it != g.node_end(); ++it){

          n = *it;
          // previous position
//--design_1
//--this doesn't work, your points are only getting pinned due because of 
//--the lines you kept in the force functors
//--START
          pos = n.position() - (n.value()).vel * eps;
//--END
          if (pos == Point(0,0,0)){
              n.position() = Point(0,0,0);
          }
          else if (pos == Point(1,0,0)){
              n.position() = Point(1,0,0);
          }
      }
  }

  virtual ~PinConstraint(){};
  
};

/**
  PlaneConstraint enforces plane constraint on the graph, inherits from Constraint
*/
struct PlaneConstraint: Constraint
{
    double Z = -0.75;

    virtual void operator()(GraphType& g, double t) {
        static_cast<void>(t);
        
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {

            Node n = *it;
            if (n.position().z < Z){

                n.position().z = Z;
                n.value().vel.z = 0.0;

            }
        }
    }
    virtual ~PlaneConstraint(){}
};

/**
  SphereConstraint enforces sheric constraint on the graph, inherits from Constraint
*/
struct SphereConstraint: Constraint
{
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;

    virtual void operator()(GraphType& g, double t) {
      static_cast<void>(t);
      Node n;
      Point r_i;
      double norm_delta;
      
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {

          n = *it;
          norm_delta = norm(n.position() - c);
          if ( norm_delta < r){
              r_i = (n.position() - c) / norm_delta;
              n.position() = c + r_i * r;
              n.value().vel = n.value().vel - dot(n.value().vel, r_i)*r_i ;
          }
      }
    }
    virtual ~SphereConstraint(){}
};

/**
  TearConstraint enforces sheric constraint on the graph, inherits from Constraint
  remove edges, nodes within sphere
  @post Graph with removed nodes/edges
*/
struct TearConstraint: Constraint
{
    Point c = Point(0.5, 0.5, -0.5);
    double r = 0.15;

    virtual void operator()(GraphType& g, double t) {
      static_cast<void>(t);
      Node n;
      Point r_i;
      double norm_delta;

//--design_1
//--You're not looping over all your nodes with this for loop (skipping some nodes)
//--START      
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
          n = *it;
          norm_delta = norm(n.position() - c);
          if ( norm_delta < r){
              g.remove_node(n);
          }
      }
//--END
    }
    virtual ~TearConstraint(){}
};


/**
  Functor that combines all constraints
*/
struct CombinedConstraints
{
  void operator()(GraphType& g, double t){
      static_cast<void>(t);
      for (unsigned i = 0; i < constraints_.size(); ++i)
      {
          (*constraints_.at(i))(g,t);
      }
  }

  CombinedConstraints(std::vector<Constraint*> constraints) : constraints_(constraints) {}

  std::vector<Constraint*> constraints_;
};

//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template<typename C1, typename C2>
CombinedConstraints make_combined_constraint(C1 c1, C2 c2) {
    std::vector<Constraint*> constraints {&c1, &c2};
    
    return CombinedConstraints(constraints);
};
//--END

template<typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraint(C1 c1, C2 c2, C3 c3) {
    std::vector<Constraint*> constraints {&c1, &c2, &c3};
    
    return CombinedConstraints(constraints);
};


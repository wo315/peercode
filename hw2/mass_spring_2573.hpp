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
  Point initpos;
  NodeData() : vel(0), mass(1), initpos(0) {}
};

struct EdgeData {
  double k;     //< Node mass
  double l;     //< Node mass
  EdgeData() : k(100), l(1) {};
  EdgeData(double k,double l) : k(k), l(l) {};
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

using NodeIter  = typename GraphType::node_iterator;
using NIncIter  = typename GraphType::incident_iterator;
using EdgeIter  = typename GraphType::edge_iterator;
//using valType   = typename GraphType::node_value_type;

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData; G::edge_value_type supports EdgeData; both classes defined above. 
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t) constraint 
 *           where g is the graph and @a t is the current time.
 *           @a constraint changes the data for the nodes in the graph to satisfy given constraints
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
};

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // apply constraints
  constraint(g,t); 

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
};

// Forces

/** @class Base class for Force
*/
class Force {
  public:
    virtual Point operator() (){
      return Point(0,0,0);
    };

    virtual const Point operator()(Node n, double t) const{
      (void) n; (void) t; //silence compiler
      return Point(0,0,0);
    };

    // virtual ~Force(){};
};

/** @class Gravity Force
*/
class GravityForce : public Force {
  public:
    virtual const Point operator()(Node n, double t) const{
      (void) t; //silence compiler
      double g = -1* grav * n.value().mass;
      return Point(0,0,g);
    };
    //virtual ~GravityForce(){};
};

/** @class Mass Spring Force
*/
class MassSpringForce : public Force{
  public:
    virtual const Point operator()(Node n, double t) const {
      (void) t; //silence compiler
      Point force = Point(0,0,0);
      
      // breadth first. update direct children.
      for (auto first = n.edge_begin(); first !=  n.edge_end() ; ++first){
        auto n2 = (*first).node2(); // Get node on the other side of the edge. (root should always be node 1)
        double d = norm(n.position()- n2.position());
        double L = (*first).value().l;
        double K = (*first).value().k;
        Point spring_force = -1 * K * (n.position()- n2.position())/d * (d-L);
        force = force + spring_force;
      };
      return force; // effect of gravity. 
  }
};

/** @class Damping force.  
*/
class DampingForce : public Force{
  public:
    virtual const Point operator()(Node n, double t) const {
      (void) t; //silence compiler
      Point force = -1 * c * n.value().vel;
      return force; // effect of gravity. 
  }
  protected:
    double c; // dampling constant
};

/** @class Functor for representing set of different forces 
 * Implements operator() to apply all the forces,
*/
class CombinedForce{
  public:
    // Constructor
    CombinedForce(std::vector<Force*> forces_array): forces_array(forces_array) {};

    void add_force(Force* f){
      forces_array.push_back(f);
    };

    const Point operator()(Node n, double t) const{
      (void) t; //silence compiler
      Point force = Point(0,0,0);
      for (auto it = forces_array.begin(); it != forces_array.end(); it++){
        // auto f = (*(*it)); // becomes base class
        // Point p = f(n,t); // this uses the base operator() again..
        // Point p1 = (*(*it))(n,t); // this uses the derived operator()
        force = force + (*(*it))(n,t);
      };
      return force; // effect of gravity. 
    };

    ~CombinedForce(){
      forces_array.clear();
    };

  protected:
    std::vector<Force*> forces_array;
};

/** make_combined_force
 * @brief Function for combining 2-3 forces
 * @return CombinedForce functor
*/
template<typename force1, typename force2>
CombinedForce make_combined_force(const force1& f1, const force2& f2){
  std::vector<Force*> forces_array;
  force1* ptr1 = new force1(); //create new force1 on the heap, since we cant store references.
  *ptr1 = f1;
  force2* ptr2 = new force2();
  *ptr2 = f2;
  forces_array.push_back(ptr1);
  forces_array.push_back(ptr2);
  return CombinedForce(forces_array);
};

template<typename force1, typename force2, typename force3>
CombinedForce make_combined_force(const force1& f1, const force2& f2, const force3& f3){
  std::vector<Force*> forces_array;
  force1* ptr1 = new force1();
  *ptr1 = f1;
  force2* ptr2 = new force2();
  *ptr2 = f2;
  force3* ptr3 = new force3();
  *ptr3 = f3;
  forces_array.push_back(ptr1);
  forces_array.push_back(ptr2);
  forces_array.push_back(ptr3);
  return CombinedForce(forces_array);
};


/** @struct Problem1Force
 *  @brief Force function for problem 1
 */
struct Problem1Force {
  double K = 100;
  // double L = 1;

  /** 
   * @brief function for returning the force at a node at time t
   * @return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    (void) t; //silence compiler
    // constraint at endpoints
    if ((n.position() == Point(0,0,0)) || (n.position() == Point(1,0,0)))
      return Point(0);

    double g = -1* grav * n.value().mass;
    Point force = Point(0,0,g);

    // breadth first. update direct children.
    for (auto first = n.edge_begin() ; first != n.edge_end() ; ++first){
      auto n2 = (*first).node2(); // Get node on the other side of the edge. (root should always be node 1)
      double d = norm(n.position()- n2.position());
      double L = (*first).length();
      Point spring_force = -1 * K * (n.position()- n2.position())/d * (d-L);
      force = force + spring_force;
    };
    return force; // effect of gravity. 
  }
};

// Constraints

/** @class Base class for constraints
*/
class Constraint {
  public:
    virtual void operator() (){
      return;
    };

    virtual void operator()(GraphType& g, double t) const{
      (void) g; (void) t; //silence compiler
      return;
    };
    // virtual ~Force(){};
};

/** @class Pin Constraint
 * Takes graph, iterate through all the nodes, and if they are at the fixed positions, change their velocity to 0.
*/
class PinConstraint : public Constraint {
  public:
    virtual void operator()(GraphType& g, double t) const{
      (void) g; (void) t; //silence compiler
      for (auto first0 = g.node_begin(); first0 != g.node_end(); ++first0) {
        if (((*first0).value().initpos == Point(0,0,0)) || 
            ((*first0).value().initpos == Point(1,0,0))){
              (*first0).position() = (*first0).value().initpos;
        }; // change positions back to 1.
      }
      return;
    };
    // virtual ~PinConstraint(){};
};

/** @class Plane Constraint
 * Takes graph, iterate through all the nodes, and if they are below z = -0.75, then move the position back up the plane.
*/
class PlaneConstraint : public Constraint {
  public:
    virtual void operator()(GraphType& g, double t) const{
      (void) g; (void) t; //silence compiler
      for (auto first0 = g.node_begin(); first0 != g.node_end(); ++first0) {
        if ((*first0).position()[2]<-0.75){
          (*first0).position()[2] = -0.75; // closest position on the plane is (x,y,-0.75)
          (*first0).value().vel[2] = 0; // set z-component velocity to zero
        }; 
      }
      return;
    };
    // virtual ~PlaneConstraint(){};
};

/** @class Sphere Constraint
 * Takes graph, iterate through all the nodes, and if they are within a sphere, pop them back out.
*/
class SphereConstraint : public Constraint {
  public:
    virtual void operator()(GraphType& g, double t) const{
      (void) g; (void) t; //silence compiler
      for (auto first0 = g.node_begin() ; first0 != g.node_end(); ++first0) {
        Point d = (*first0).position()- c; // vector from center to the point
        if (norm(d) < r){
          Point Ri = d/norm(d); //normal vector
          Point proj = dot((*first0).value().vel, Ri) * Ri; // projection
          (*first0).position() = c + Ri*r; // closest position on the surface
          (*first0).value().vel = (*first0).value().vel - proj ; // set velocity normal to the surface to zero
        }; 
      }
      return;
    };
  
  protected:
    Point c = Point(0.5,0.5,-0.5) ;
    double r = 0.15;
    // virtual ~PlaneConstraint(){};
};

/** @class Sphere Remove Constraint
 * Takes graph, iterate through all the nodes, and if they are within a sphere, remove the nodes.
*/
class SphereRemoveConstraint : public Constraint {
  public:
    virtual void operator()(GraphType& g, double t) const{
      (void) g; (void) t; //silence compiler
      for (auto first0 = g.node_begin() ; first0 != g.node_end(); ++first0) {
        Point d = (*first0).position()- c; // vector from center to the point
        if (norm(d) < r){
          g.remove_node(*first0);
        }; 
      }
      return;
    };
  
  protected:
    Point c = Point(0.5,0.5,-0.5) ;
    double r = 0.15;
    // virtual ~PlaneConstraint(){};
};

/** @class Functor for representing set of different constraints
 * Implements operator() to apply all the constraints,
*/
class CombinedConstraints{
  public:
    // Constructor
    CombinedConstraints(std::vector<Constraint*> constraints_array): constraints_array(constraints_array) {};

    const void operator()(GraphType& g, double t) const {
      (void) t; //silence compiler
      for (auto it = constraints_array.begin(); it != constraints_array.end(); it++){
        (**it)(g,t); // apply constraint by calling each constraint.
      };
      return; // effect of gravity. 
    };

    ~CombinedConstraints(){
      constraints_array.clear();
    };

  protected:
    std::vector<Constraint*> constraints_array;
};

/** make_combined_constraints
 * @brief Function for combining 2-3 forces
 * @return CombinedForce functor
*/
template<typename constraint1, typename constraint2>
CombinedConstraints make_combined_constraint(const constraint1& c1, const constraint2& c2){
  std::vector<Constraint*> constraints_array;
  constraint1* ptr1 = new constraint1(); //create new force1 on the heap, since we cant store references.
  *ptr1 = c1;
  constraint2* ptr2 = new constraint2();
  *ptr2 = c2;
  constraints_array.push_back(ptr1);
  constraints_array.push_back(ptr2);
  return CombinedConstraints(constraints_array);
};

template<typename constraint1, typename constraint2, typename constraint3>
CombinedConstraints make_combined_constraint(const constraint1& c1, const constraint2& c2, const constraint3& c3){
  std::vector<Constraint*> constraints_array;
  constraint1* ptr1 = new constraint1(); //create new force1 on the heap, since we cant store references.
  *ptr1 = c1;
  constraint2* ptr2 = new constraint2();
  *ptr2 = c2;
  constraint3* ptr3 = new constraint3();
  *ptr3 = c3;
  constraints_array.push_back(ptr1);
  constraints_array.push_back(ptr2);
  constraints_array.push_back(ptr3);
  return CombinedConstraints(constraints_array);
};

//--functionality_1
//--The figure ends midway with segmentation fault.
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


/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <vector>

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
  Point initial_position; //< Node initial position
  NodeData() : vel(0), mass(1), initial_position(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K; //< Node spring constant
  double L; //< Node rest length (initial)
  EdgeData() : K(100.0), L(0.1) {}
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
 * @tparam G:: node_value_type supports NodeData
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

    if (n.position() != Point(0,0,0) and n.position() != Point(1,0,0)){
      // Update the position of the node according to its velocity
      // x^{n+1} = x^{n} + v^{n} * dt
      n.position() += n.value().vel * dt;
    }
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
    (void) t; // silence warning

    // Pinned points
    if (n.position () == Point(0 ,0 ,0) || n.position () == Point(1 ,0 ,0))
      return Point(0 ,0 ,0);

    Point force = Point(0,0,0); //Initialize force
    double distance;
    Point p2;

    // Loop through neighbors of n
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){

      p2 = (*it).node2().position();
      distance = norm(n.position() - p2);

      // Substract spring force associated to current neighbor
      force = force - ((*it).value().K) *
                      ((distance - ((*it).value().L))/distance) *
                      (n.position() - p2);
    }
    // Add gravity
    force = force + Point(0,0, -n.value().mass*grav);
    return force;
  }
};

/** Parent Force class */
class Force {

  public:
    /** Return the force applying to @a n at time @a t.
     *
     * For parent class, the force is 0
     */
    virtual Point operator() (Node n, double t) {
      (void) n; (void) t;
      return Point(0,0,0);
    }

    virtual ~Force() {

    }
};

/** Derived gravity Force class */
class GravityForce: public Force {

  public:
    /** Return the gravity force applying to @a n at time @a t. */
    Point operator() (Node n, double t) {
      (void) t;
      return Point(0,0, -n.value().mass*grav);
    }
};

/** Derived mass spring Force class */
class MassSpringForce: public Force {

  public:

    /** Return the mass spring force applying to @a n at time @a t. */
    Point operator() (Node n, double t) {
      (void) t; // silence warning
      Point force = Point(0,0,0);
      double distance;
      Point p2;

      // Loop through neighbors of n
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it){

        p2 = (*it).node2().position();
        distance = norm(n.position() - p2);

        // Substract spring force associated to current neighbor
        force = force - ((*it).value().K) *
                        ((distance - ((*it).value().L))/distance) *
                        (n.position() - p2);
      }
      return force;
    }
};

/** Derived damping Force class */
class DampingForce: public Force {

  private:
    // Damping parameter
    const double c_;

  public:

    DampingForce(const double c) : c_(c) {  }

    /** Return the damping force (parameter c_) applying to @a n at time @a t.
    */
    Point operator() (Node n, double t) {
      (void) t;
      return (-c_ * n.value().vel);
    }
};

/** Class implementing the sum of other Forces, derived from Force class */
class CombinedForce: public Force {

  private:
    // Vector of all forces applied
    std::vector<Force*> forces_;

  public:
    CombinedForce(std::vector<Force*> forces) : forces_(forces) {}

    /** Return the combined forces in forces_ applying to @a n at time @a t. */
    Point operator() (Node n, double t) {

      Point resulting_force = Point(0,0,0);

      // Loop through forces_ vector
      for (auto it = forces_.begin(); it != forces_.end(); ++it){
        resulting_force += (**it)(n,t);
      }
      return resulting_force;
    }

    ~CombinedForce (){
      forces_.clear();
    }
};

/** Return combined forces based on 2 forces
 * @param[in]     f1     First force
 * @param[in]     f2     Second force
 * @return CombinedForce combining @a f1 and @a f2
 *
 * @tparam t1:: one of the derived classes of Force
 * @tparam t2:: one of the derived classes of Force
 */
//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template <typename t1, typename t2>
CombinedForce make_combined_force(t1& f1, t2& f2){
  std::vector<Force*> f_vec;
  f_vec.push_back(&f1);
  f_vec.push_back(&f2);
  return CombinedForce(f_vec);
}
//--END




/** Return combined forces based on 3 forces
 * @param[in]     f1     First force
 * @param[in]     f2     Second force
 * @param[in]     f3     Third force
 * @return CombinedForce combining @a f1 , @a f2 and @a f3
 *
 * @tparam t1:: one of the derived classes of Force
 * @tparam t2:: one of the derived classes of Force
 * @tparam t3:: one of the derived classes of Force
 */
//--functionality_1
//--You're passing the forces by copy, and then taking a pointer
//--to that force. This is a problem as the force is deleted
//--Once the function goes out of scope, and you now have a 
//--dangling pointer
//--I corrected it to be able to test the rest of your functions.
//--START
template <typename t1, typename t2, typename t3>
CombinedForce make_combined_force(t1& f1, t2& f2, t3& f3){
  std::vector<Force*> f_vec;
  f_vec.push_back(&f1);
  f_vec.push_back(&f2);
  f_vec.push_back(&f3);
  return CombinedForce(f_vec);
}
//--END
/** Parent Constraint class: constraint applied to the nodes of the graph */
class Constraint {

  public:
    /** Apply the constraint to @a g at time @a t.
     *
     * For parent class, no constraint is applied.
     */
    virtual void operator()(GraphType& g, double t) {
      (void) g; (void) t;
    }

    virtual ~Constraint() {

    }
};

/** Derived Pin Constraint class */
class PinConstraint : public Constraint {

  public:

    /** Apply the Pin constraint to @a g at time @a t.
     * The goal is to fix (0,0,0) and (1,0,0) at
    */
    void operator()(GraphType& g, double t) {
      (void) t;

      // Loop through all nodes of g
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Node n = *it;
        // Fetch initial position
        Point init_pos = n.value().initial_position;
        if (init_pos == Point(0,0,0) or init_pos == Point(1,0,0)){
          //If initial position is one of two points, get pos back to initial
          n.position() = init_pos;
        }
      }
    }
};

/** Derived Plane Constraint class */
class PlaneConstraint : public Constraint {

  public:

    /** Apply the Plane constraint to @a g at time @a t.
     * The goal is to not let the graph go below the plane z = -0.75 (simulates
     * the ground)
    */
    void operator()(GraphType& g, double t) {
      (void) t;
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Node n = *it;
        // Constraint violation
        if (n.position().z < -0.75){
          n.position().z = -0.75; // Closest point on plane
          n.value().vel.z = 0;
        }
      }
    }
};

/** Derived Sphere Constraint class */
class SphereConstraint : public Constraint {

  public:

    Point c = Point(0.5, 0.5, -0.5); // Center of the sphere
    double r = 0.15; // Radius of the Sphere

    /** Apply the Sphere constraint to @a g at time @a t.
     * This constraint is equivalent to putting a solid sphere in the space:
     * graph can't go through it
    */
    void operator()(GraphType& g, double t) {
      (void) t;

      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        Node n = *it;
        // Constraint violation
        if (norm(n.position() - c) < r){
          Point R = (n.position() - c)/norm(n.position() - c);
          n.position() = c + r*R; // Closest point on surface of sphere
          n.value().vel -= dot(n.value().vel, R) * R;
        }
      }
    }
};

/** Derived Tear Constraint class */
class TearConstraint: public Constraint {

  public:

    Point c = Point(0.5, 0.5, -0.5); // Center of the sphere
    double r = 0.15; // Radius of the Sphere

    /** Apply the Tear constraint to @a g at time @a t.
     * This constraint is equivalent to putting a solid and sharp sphere
     * in the space: graph can't go through it and gets torn at contact
     * points
    */
    void operator()(GraphType& g, double t) {
      (void) t;

      auto it = g.node_begin();
      while (it != g.node_end()){
        Node n = *it;
        // Constraint violation
        if (norm(n.position() - c) < r){
          // Remove node via iterator: returns iterator pointing to next Node
          it = g.remove_node(it);
        }
        else{
          // Get to next node
          ++it;
        }
      }
    }
};

/** Class implementing the combination of other Constraints,
 * derived from Constraint class */
class CombinedConstraints : public Constraint {

  private:
  // Vector of all constraints applied
    std::vector<Constraint*> constraints_;

  public:
    CombinedConstraints(std::vector<Constraint*> constraints) :
    constraints_(constraints) {}

    /** Apply all the constraints in constraints_  to @a g at time @a t. */
    void operator()(GraphType& g, double t) {

      for (auto it = constraints_.begin(); it != constraints_.end(); ++it){
        (**it)(g,t);
      }
    }

    ~CombinedConstraints (){
      constraints_.clear();
    }
};

/** Return combined constraints based on 2 constraints
 * @param[in]     c1     First Constraint
 * @param[in]     c2     Second Constraint
 * @return CombinedConstraints combining @a c1 and @a c2
 *
 * @tparam T1:: one of the derived classes of Constraint
 * @tparam T2:: one of the derived classes of Constraint
 */
template <typename T1, typename T2>
CombinedConstraints make_combined_constraints(T1 c1, T2 c2){
  std::vector<Constraint*> c_vec;
  c_vec.push_back(&c1);
  c_vec.push_back(&c2);
  return CombinedConstraints(c_vec);
}

/** Return combined constraints based on 3 constraints
 * @param[in]     c1     First Constraint
 * @param[in]     c2     Second Constraint
 * @param[in]     c3     Third Constraint
 * @return CombinedConstraints combining @a c1 , @a c2 and @a c3
 *
 * @tparam T1:: one of the derived classes of Constraint
 * @tparam T2:: one of the derived classes of Constraint
 * @tparam T3:: one of the derived classes of Constraint
 */
template <typename T1, typename T2, typename T3>
CombinedConstraints make_combined_constraints(T1& c1, T2& c2, T3& c3){
  std::vector<Constraint*> c_vec;
  c_vec.push_back(&c1);
  c_vec.push_back(&c2);
  c_vec.push_back(&c3);
  return CombinedConstraints(c_vec);
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and given constraint
 * @param[in,out] g           Graph
 * @param[in]     t           The current time, useful for time-dependent forces
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint
 *
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G:: node_value_type supports NodeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint returns void and applies the constraint on
 *           Graph g at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step (G& g, double t, double dt , F force , C constraint){
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraint to the graph
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

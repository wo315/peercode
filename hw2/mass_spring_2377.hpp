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

  // Sentinel value for whether a node should be pinned
  // 0 -> Free
  // 1 -> Pinned to origin
  // 2 -> Pinned to (1,0,0)
  // Note: could do bool pinned, Point pinned_point to generalize pinning
  // Downside: Have to store the pinned point, which is larger than an int
  int pinned;
  NodeData() : vel(0), mass(1), pinned(0) {}
};

/** Custom stucture of data to store with Edges */
struct EdgeData{
  float K; // Spring Constant
  float L; // Rest Length
  EdgeData(): K(100.0), L(1.0) {}
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
 * @tparam G::node_value_type supports structs containing (at least) a vel, mass, and pinned attribute
            G::edge_value_type supports sturcts containing (at least) a K (spring constant) and L (rest length) attribute
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
    //Pinned Points
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }

    //Calculate gravity and spring separately and combine
    Point grav_effect = Point(0,0,-1*grav)*n.value().mass;
    Point spring_effect = Point(0,0,0);
    //For each adjacent node, add the spring effect
    for (auto it = n.edge_begin(); it!= n.edge_end();++it){
      //Length of spring
      Point dist_  = n.position()-(*it).node2().position();
      //f_spring = sumover adj nodes -Kij (xi-xj)/(|xi-xj|)(|xi-xj|-Lij)
      //Note that edge.length() returns a constant value for this problem
      spring_effect += -1*(*it).value().K*(dist_)/(norm(dist_))*(norm(dist_)-(*it).length());
    }
    (void) t;
    return grav_effect+spring_effect;
  }
};


/** Base Force class
 *  @tparam NODE can be arbitrary for Force, but derived classes
 *          define their own restrictions
 */
template <typename NODE>
class Force {
  public:
    /** General operator() function
     *  @param n a NODE object being acted on
     *  @param t time
     *  @return Point(0,0,0)
     */
    virtual Point operator()(NODE& n, double t) const{
      (void) n;
      (void) t;
      return Point();
  }
};

/** Gravity class, derived from Force with NODE=Node */
class GravityForce : public Force<Node>{
  public:
    /** Gravity operator() function
     *  @param n the Node being acted on
     *  @param t time
     *  @return Point represnting the force of gravity on n
     *  @pre n has value() initialized as a structure with a mass value
     */
    Point operator()(Node& n, double t) const{
       (void) t;
       //mass*(0,0,-grav)
       return Point(0,0,-1*grav)*n.value().mass;
     }
};

/** MassSpring class, derived from Force with NODE=Node */
class MassSpringForce : public Force<Node> {
  public:
    /** MassSpring operator() function
     *  @param n the Node being acted on
     *  @param t time
     *  @return Point represnting the force of spring on n
     *  @pre edges have value() initialized as a structure with a spring constant K and rest_length L
     */
    Point operator()(Node& n, double t) const{
      (void) t;
      Point rslt = Point();
      for (auto it_ = n.edge_begin(); it_!= n.edge_end();++it_){
        //Second node is always the "other" node
        Point dist_  = n.position()-(*it_).node2().position();
        //f_spring = sumover adj nodes -Kij (xi-xj)/(|xi-xj|)(|xi-xj|-Lij)
        rslt += -1*(*it_).value().K*(dist_)/(norm(dist_))*(norm(dist_)-(*it_).value().L);
      }
      return rslt;
    }
};

/** Damping class, derived from Force with NODE=Node */
class DampingForce : public Force<Node> {
  public:
    /** DampingForce operator() function
     *  @param n the Node being acted on
     *  @param t time
     *  @param c the damping constant
     *  @return Point represnting the force of spring on n
     *  @pre n has value() initialized as a structure with a velocity value vel
     */
    Point operator() (Node& n, double t, double c = 1) const{
      (void) t;
      return n.value().vel*c;
    }
};

/** CombinedForce functor that combined several sequential forces
 *  @param f: a vector pointers to Force objects
 */
class CombinedForce{
  public:
    std::vector<const Force<Node>*> f_;
    CombinedForce(std::vector< const Force<Node>*> f) : f_(f) {}
    /** CombinedForce operator() function
     *  @param n the Node being acted on
     *  @param t time
     *  @return Point representing the net effects of all
     *          forces in f_
     */
    Point operator()(Node& n, double t){
      Point rslt = Point();
      for (auto it_ = f_.begin(); it_ != f_.end(); ++it_){
        //Iterator over pointers
        rslt += (**it_)(n, t);
      }
      return rslt;

  }
};


/** function called to apply a series of forces
 *  @param f A sequence of Force-derived objects
 *  @return  A CombinedForce object containing all of the forces in f
 *  @pre  f is not empty
 *  @ref https://stackoverflow.com/questions/47126273/converting-parameter-pack-into-a-vector
 */
//--design_0
//--good job with parameter pack!
//--START
template <typename... FORCES>
CombinedForce make_combined_force(FORCES const & ... f){
  std::vector<const Force<Node>*> force_vec {{&f...}};
  return CombinedForce(force_vec);
};
//--END


/** Base constraint class*/
class Constraint{
  public:
    /** General operator()
     *  @param g the graph being constrained
     *  @param t time
     *  @result none. Graph is left as is
     */
    virtual void operator()(GraphType& g, double& t) const {
        (void) g;
        (void) t;
    }
};

/**Pin constraint class*/
class PinConstraint : public Constraint {
  public:
    /** Pin operator()
     *  @param g the graph being constrained
     *  @param t time
     *  @post Any pinned nodes are moved to their pinned location
     *          Other nodes are unaffected
     *  @pre Nodes are defined as pinned via the value().pinned attribute
     *        where:
     *          0 -> Node is free
     *          1 -> Node is pinned to (0,0,0)
     *          2 -> Node is pinned to (1,0,0)
     */
    void operator()(GraphType& g, double& t) const {
      (void) t;
      for (auto it = g.node_begin(); it!=g.node_end(); ++it){
        if ((*it).value().pinned){
          (*it).position() = (*it).value().pinned == 1 ? Point(0.0,0.0,0.0) : Point(1.0,0.0,0.0);
        }
      }
    }
};

/**Plane constraint class*/
class PlaneConstraint : public Constraint{
  public:
    /** Plane operator()
     *  @param g the graph being constrained
     *  @param t time
     *  @post Any node that would be on the wrong side of the plane are set
     *          to the nearest point on the plane
     *          Other nodes are unaffected
     *  @pre n has value() initialized as a structure with a velocity value vel
     */
    void operator()(GraphType& g, double& t) const {
      (void) t;
      const Point p = Point(0,0,1);
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        if (dot((*it).position(),p) < -0.75){
          //Nearest point is normal to the actual point, i.e. (x,y,-0.75)
          (*it).position()[2] = -0.75;
          (*it).value().vel[2] = 0.0;
        }
      }
    }

};

/**Sphere constraint class to MOVE violating nodes*/
class SphereConstraint : public Constraint{
  public:
    /** Plane operator()
     *  @param g the graph being constrained
     *  @param t time
     *  @post Any node that would be inside the sphere are set
     *          to the nearest point on the sphere
     *          Other nodes are unaffected
     *  @pre n has value() initialized as a structure with a velocity value vel
     */
    void operator()(GraphType& g, double& t) const{
      (void) t;
      const Point c = Point(0.5,0.5, -0.5);
      for (auto it = g.node_begin(); it != g.node_end(); ++it){
        //Check if the point is within one radius of the center
        Point dist = (*it).position() - c;
        if (norm(dist) < 0.15){
          //Nearest point has direction dist and is set at c+r
          //https://piazza.com/class/kjx6t0a2rv4589?cid=236
          (*it).position() = c + (0.15)/norm(dist)*dist;
          (*it).value().vel = (*it).value().vel - dot((*it).value().vel,dist/norm(dist))*(dist/norm(dist));
        }
      }
    }

};

/**Sphere constraint class to REMOVE violating nodes*/
class SphereConstraintRemove : public Constraint{
  public:
    /** Plane operator()
     *  @param g the graph being constrained
     *  @param t time
     *  @post Any node that would be inside the sphere are removed from g
     *          Other nodes are unaffected
     *  @pre n has value() initialized as a structure with a velocity value vel
     */
    void operator()(GraphType& g, double& t) const{
      (void) t;
      const Point c = Point(0.5,0.5, -0.5);
      auto it = g.node_begin();
      //Re-calculates g.node_end() each iteration since remove may invalidate an iterator
      while (it != g.node_end()){
        Point dist = (*it).position() - c;
        if (norm(dist) < 0.15){
          //Returns a valid new iterator
          it = g.remove_node(it);
        }
        else{
          //Iterate
          ++it;
        }
      }
    }
};






/** CombinedConstraints functor that combined several sequential constraints
 *  @param c: a vector pointers to Constraint objects
 */
class CombinedConstraints{
  public:
    std::vector<const Constraint*> c_;
    CombinedConstraints(std::vector< const Constraint*> c) : c_(c) {}
    /** CombinedConstraints operator() function
     *  @param g the graph being acted on
     *  @param t time
     *  @post All constraints in c are applied to g seqentially
     */
    void operator()(GraphType& g, double t){
      for (auto it_ = c_.begin(); it_ != c_.end(); ++it_){
        (**it_)(g, t);
      }

  }
};

/** function called to apply a series of constraints
 *  @param f A sequence of Constraint-derived objects
 *  @return  A CombinedConstraint object containing all of the forces in c
 *  @pre  c is not empty
 *  @ref https://stackoverflow.com/questions/47126273/converting-parameter-pack-into-a-vector
 */
template <typename... CONSTRAINTS>
CombinedConstraints make_combined_constraints(CONSTRAINTS const & ... c){
  std::vector<const Constraint*> constraint_vec {{&c...}};
  return CombinedConstraints(constraint_vec);
};



/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g             Graph
 * @param[in]     t             The current time (useful for time-dependent forces)
 * @param[in]     dt            The time step
 * @param[in]     force         Function object defining the force per node
 * @param[in]     constraints   Function object defining the constraints per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports structs containing (at least) a vel, mass, and pinned attribute
            G::edge_value_type supports sturcts containing (at least) a K (spring constant) and L (rest length) attribute
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraints(g, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a constraints must return void and update all nodes at time @a t
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraints){
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  //Apply the constraints to the graph
  constraints(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}





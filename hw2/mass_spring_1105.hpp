/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <cmath>
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
  Point original_position; 
  NodeData() : vel(0), mass(1), original_position(0) {}
};

struct EdgeData {
  double edge_k;     // Edge specific spring constant 
  double edge_rest_len;     // Edge specific resting length
  EdgeData() : edge_rest_len(0), edge_k(0) {} //, edge_value(0) {}
};


// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
//using GraphType = Graph<NodeData, double>;
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

//second euler function for 6.3
template<typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // apply connstraint after updating positionn but before calculating force
  for (auto it = g.node_begin(); it != g.node_end(); ++it){
    auto n = (*it);

    constraint(g, t);

    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}


//HW 2: 3.2: Problem1Force as is

/** Problem1Force to apply mass spring force and gravity together.
* param [in] n template node and t double to indicate the time
* param [out] total force as a point applied to n at time t
*
* post new F_total + F_grabity >= F_gravity
* post this force will be applied to the nodes it is called on and swings the fabric.
*/
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    double mass_n = n.value().mass;
    Point F_total = Point(0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = (*ei);
      Node n2 = e.node2();

      double L0 = e.edge_value().edge_rest_len; //resting length
      double L = e.length(); //current length
      double K = e.edge_value().edge_k;

      // Prevent the fall into an infinite abyss
      if (n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0))
        return Point (0 ,0 ,0);

      Point u = (n2.position() - n.position()) / norm(n2.position() - n.position());
      //std::cout << "u = " << u.x << " " << u.y << " " << u.z << std::endl;
      double diff = (K * (L - L0));
      //std::cout << "diff = " << diff << std::endl;
      Point F_spring = u * diff;
      F_total = F_total + F_spring;
    }

    Point F_gravity = (Point(0, 0, -grav) * mass_n);

    return F_total + F_gravity;
  }
};



// HW2: 4.2: Generalized Mass-Spring


/* Force Base class 
* param [in] n template node and t double to indicate the time
* param [out] total force as a point applied to n at time t
*
* This base class will have one virtual function that will be
* redefined in the derived classes. This will allow for polymorphism
*/
struct Force {
  Point F_ = Point(0);
  virtual Point operator()(Node n, double t) {
    return F_;
  }  
};


/* Derived MassSpringForce
* param [in] n template node and t double to indicate the time
* param [out] total spring force applied to n at time t
* 
* post notice that the final force is in the direction needed
* the spring force on its own does not move the graph.
* instead, it's interaction with gravity is what'll move it.
*/
struct MassSpringForce: public Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  Point operator()(Node n, double t) {
    Point F_ = Point(0);
    if (n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0))
        return Point (0 ,0 ,0);

    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei) {
      Edge e = *ei;
      Node n2 = e.node2();
      double L0 = e.edge_value().edge_rest_len; //resting length
      double L = e.length(); //current length
      double K = e.edge_value().edge_k;

      Point u = (n2.position() - n.position()) / norm(n2.position() - n.position());
      double diff = (K * (L - L0));
      
      Point F_spring = u * diff;
      F_ = F_ + F_spring;
    }
  
    return F_;
  }
};

/* Derived GravityForce
* param [in] n template node and t double to indicate the time
* param [out] total gravitational force applied to n at time t
* 
* post notice that the final force is in the direction needed
*
*/
struct GravityForce: public Force 
{
  Point F_ = Point(0);
  Point operator()(Node n, double t) {
    if (n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0))
      return Point (0 ,0 ,0);
    double mass_n = n.value().mass;
    F_ = Point(0, 0, -grav) * mass_n;
    return F_;
  }
};

/* Derived DampingForce
* param [in] n template node and t double to indicate the time
* param [out] total damping force applied to n at time t
*
*/
struct DampingForce: public Force {
  template <typename NODE>
  Point operator()(Node n, double t){
    Point vel_n = n.value().vel;
    double damp_coef = 0.1;
    F_ = - damp_coef * vel_n;
    return F_;    
  }
};

/* A functor that combines given forces
* param [in] a vector of pointers to the base force class
* param [out] a summation of the given forces
* 
* The constructor takes the vector of pointers.
*/
class CombinedForce {
  std::vector<Force*> force_;
public:
  CombinedForce(std::vector<Force*> force) 
  : force_(force) {}
  
  Point operator() (Node n, double t) {
    Point combined_force = Point(0);
    for(auto it = force_.begin(); it != force_.end(); ++it) {
      combined_force = combined_force + (*it)->operator()(n, t);
    }
     ;
    return combined_force;
  }  
};

/* A funnction that combines the given forces
* param [in] types of forces
* param [out] total force applied to n at time t
* 
* We'll allow this function to take two or three
* arguments by using templates.
*/
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 a, F2 b) {
  std::vector<Force*> forces {};
  forces.push_back(&a);
  forces.push_back(&b);
  return CombinedForce(forces);
}

/* A funnction that combines the given forces
* param [in] types of forces
* param [out] total force applied to n at time t
* 
* We'll allow this function to take two or three
* arguments by using templates.
*/
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 a, F2 b, F3 c){
  std::vector<Force*> forces {};
  forces.push_back(&a);
  forces.push_back(&b);
  forces.push_back(&c);
  return CombinedForce(forces);
}


/* Constraint Base class 
* param [in] a graph and a and a time t
*
* post: the constraints will be applied one by one
*
* This base class will have one pure virtual function that must be
* defined in the derived classes. This will allow polymorphism.
*/
struct Constraint {
  //I'll make this pure virtual since every class needs to have one
  virtual void operator()(GraphType graph, double t) = 0; 
};

/* PinConstraint Derived class 
* param [in] a graph and a and a time t
*
* post: the pin constraint will fix points 0,0,0 and 1,0,0
* this is a way to generalize the constraints inside the Problem1Force
*/
struct PinConstraint: public Constraint {
  void operator()(GraphType graph, double t) {
    // First, find the original position of each node
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        if ((*it).value().original_position == Point(0,0,0)){
          (*it).position() = Point(0,0,0);
        } else if ((*it).value().original_position == Point(1,0,0)) {
          (*it).position() = Point(1,0,0);
        }
    }
  }
};

/* PinPlaneConstraint Derived class 
* param [in] a graph and a and a time t
*
* post: the plane constraint will constraint nodes at the -0.75 line
* this is a way to generalize the constraints inside the Problem1Force
*/
struct PlaneConstraint: public Constraint {
  void operator()(GraphType graph, double t) {
    Point p = Point(0, 0, 1);
    //Iterate through each node to find those that violate
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
        if (dot((*it).position(), p) < -0.75) {
          //set the new position to the point nearest (i.e. perpendicular) to the plane with z = 0.75
          (*it).position() = Point((*it).position().x, (*it).position().y, -0.75); 
          (*it).value().vel.z = 0; 
        }
      }
    }
};

/* SphereConstraint Derived class 
* param [in] a graph and a and a time t
*
* post: the sphere constraint will make the graph look like it has
* a ball inside it
*/
struct SphereConstraint: public Constraint {
//citation: https://math.stackexchange.com/questions/831109/closest-point-on-a-sphere-to-another-point
  void operator()(GraphType graph, double t){
    Point c = Point(0.5, 0.5, -0.5);
    double rad = 0.15;

    //Find the nodes that violate
    for (auto it = graph.node_begin(); it != graph.node_end(); ++it) {
      if (norm((*it).position() - c) < rad) { //since norm is >0, abs can be removed
        Point Ri = ((*it).position() - c )/norm((*it).position() - c);
        Point q = c + (rad * Ri);
        (*it).position() = q;
        (*it).value().vel = (*it).value().vel - (Ri * (dot((*it).value().vel, Ri)));
      }
    }
  }
};

/* SphereConstraint Derived class 
* param [in] a graph and a and a time t
*
* post: the tear constraint will make the graph look like a ball is cut inside it
*/
struct TearConstraint: public Constraint {

  void operator()(GraphType graph, double t){
    std::cout << "started tear " << std::endl;
    Point c = Point(0.5, 0.5, -0.5);
    double rad = 0.15;

//--functionality_1
//--The figure cannot move if TearConstraint is applied.
//--First, the condition of while statement may always be true. it should be while(it != graph.node_end()).
//--Second, there should be else{++it;} after if block.
//--START

    while (graph.num_nodes() > 0) {
      auto it = graph.node_begin();
      if (norm((*it).position() - c) < rad) { //since norm is >0, abs can be removed
        graph.remove_node(it);
      }
    }
  }
};
//--END

/* A functor that combines given constraints
* param [in] a vector of pointers to Constraint
*
* post: the constraints will be applied one by one
*/
class CombinedConstraints {
  std::vector<Constraint*> constraint_;

public:
  CombinedConstraints(std::vector<Constraint*> constraint) 
  : constraint_(constraint) {}

  void operator()(GraphType graph, double t) {
    for (auto it = constraint_.begin(); it != constraint_.end(); ++it){
      (*it)->operator()(graph, t);
    }
  }
};

/* make_combined_constraint
* param [in] Constraint objects
* param [out] a CombinedConstraint object
*
* post: each constraint gets applied one at a time
* take a template of two onstraints
*/
template <typename C1, typename C2>
CombinedConstraints make_combined_constraint(C1 a, C2 b) {
  std::vector<Constraint*> constraints {};
  constraints.push_back(&a);
  constraints.push_back(&b);
  CombinedConstraints combined_constraint_called = CombinedConstraints(constraints);
  return combined_constraint_called;
}

/* combines given constraints
* param [in] Constraint objects
* param [out] a CombinedConstraint object
*
* post: each constraint gets applied one at a time
* take a template of three constraints
*/
template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraint(C1 a, C2 b, C3 c){
  std::vector<Constraint*> constraints {};
  constraints.push_back(&a);
  constraints.push_back(&b);
  constraints.push_back(&c);
  CombinedConstraints combinend_constraint_called = CombinedConstraints(constraints);
  return combinend_constraint_called;
}


//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


























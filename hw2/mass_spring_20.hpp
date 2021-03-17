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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply Constraint
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    // if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
    //     continue;
    // }  
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
    (void) t;
    const double K = 100.0;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }
    Point f_spring = Point(0,0,0);

    for (auto iter = n.edge_begin(); iter != n.edge_end(); ++iter){
      auto rest_len = (*iter).value();
      auto sub = n.position() - (*iter).node2().position();
      f_spring += -K * (sub / norm(sub)) * (norm(sub) - rest_len);
    }

    auto f_grav = n.value().mass * Point(0,0,-grav);

    return f_spring + f_grav;
  }
};


struct Force {
  /*
  Parent class for all the sub-forces. 
  */

  virtual Point operator()(Node n_, double t_){
    (void) t_;
    return Point(0,0,0);
  }

};


struct GravityForce : Force {
  
  /**
   * @param[in] Node n, time t
   *    Uses gravity constant and node's mass to return
        gravitational force (mass * (0,0,-gravity constant))
   * @return Gravitational force for node n.
  */

  virtual Point operator()(Node n_, double t_) {
    (void) t_;
    return n_.value().mass * Point(0,0,-grav);
  }

};

struct MassSpringForce : Force {
  /**
   * @param[in] Node n, time t
   *  Uses spring constant and node's rest length, and
      positions of all adjacent nodes to calculate the
      mass-spring force.
   * @return Mass-spring force for node n.
  */
  virtual Point operator()(Node n_, double t_) {
    (void) t_;
    const double K_ = 100.0;
    Point f_spring = Point(0,0,0);
    for (auto iter = n_.edge_begin(); iter != n_.edge_end(); ++iter){
      auto rest_len = (*iter).value();
      auto sub = n_.position() - (*iter).node2().position();
      f_spring += -K_ * (sub / norm(sub)) * (norm(sub) - rest_len);
    }
    return f_spring;
  }

};

struct DampingForce : Force {
  /**
   * @param[in] Node n, time t
   *  Uses damping constant and node n to calculate
   *  damping force value.
   * @return Damping force for node n.
  */
  virtual Point operator()(Node n_, double t_) {
    (void) t_;
    const double c_ = 1.0;
    return -c_ * n_.value().vel;
  }

};


struct CombinedForce {
  /**
   * @param[in]  n  Node to calculate force from
   * @param[in]  t  time parameter

   * Inputs vector of 2 or 3 derived classes from Force, 
   * uses operator () to find the total force based on
   * the combination of all forces inside the vector.

  * @return total_force
  */

  std::vector<Force*> f_vec_;
  CombinedForce(std::vector<Force*> f_vec) : f_vec_(f_vec) {}

  Point operator()(Node n, double t){
      Point tot_force = Point(0,0,0);
      for (auto f : f_vec_) {
        tot_force += (*f)(n, t);
      }
      return tot_force;
  }

};

template<typename F1, typename F2>
CombinedForce make_combined_force(F1 force1, F2 force2){
  /**
   * @param[in] GravityForce(), MassSpringForce() or DampingForce() object
   *  Inputs two of the sub-forces that inherit from the Force class
    in an arbitrary order, and returns a CombinedForces functor
    object that is used to calculate the overall force on a node.
   * @return CombinedForce object.
  */
  std::vector<Force*> forces_vec;

  forces_vec.push_back(&force1);
  forces_vec.push_back(&force2);
  
  return CombinedForce(forces_vec);
}

template<typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 force1, F2 force2, F3 force3){
  /**
   * @param[in] GravityForce(), MassSpringForce() or DampingForce() object
   *  Inputs three of the sub-forces that inherit from the Force class
    in an arbitrary order, and returns a CombinedForces functor
    object that is used to calculate the overall force on a node.
   * @return CombinedForce object.
  */
  std::vector<Force*> forces_vec;

  forces_vec.push_back(&force1);
  forces_vec.push_back(&force2);
  forces_vec.push_back(&force3);

  return CombinedForce(forces_vec);
}

struct Constraint {

  virtual void operator()(GraphType& g, double t){
    (void) t;
  }

};

struct PinConstraint : Constraint{
  /**
   * @param[in] Graph g, time t
   *  Iteratively pins any points at the corners to the origin.
  */
  virtual void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto pos = (*i).position();
      if (pos == Point(0,0,0) || pos == Point(1,0,0)){
        (*i).position() = Point(0,0,0);
        (*i).value().vel = Point(0,0,0);
      }
    }
  }

};

struct PlaneConstraint : Constraint{
  /**
   * @param[in] Graph g, time t
   * Iteratively adjusts the position and velociy 
   * of any point outside of a z-axis plane.
  */
  double plane_z = -0.75;
  virtual void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto pos = (*i).position();
      if (pos.z < plane_z){
        (*i).position().z = plane_z;
        (*i).value().vel.z = 0;
      }
    }
  }


};

struct SphereConstraint : Constraint{

   /**
   * @param[in] Graph g, time t
   * Iteratively adjusts the position and velocity of any points
   * that fall outside of a sphere
   */

  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;

  virtual void operator()(GraphType& g, double t){
    (void) t;
    for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto pos = (*i).position();
      auto sub = pos - c;
      auto this_vel = (*i).value().vel;
      if (norm(sub) < r){
        auto Ri = sub / norm(sub);
        auto nearest = c + r * Ri;
        (*i).position() = nearest;
        (*i).value().vel -= dot(this_vel, Ri) * Ri;
      }
    }
  }

};

struct TearConstraint : Constraint{

   /**
   * @param[in] Graph g, time t
   * Iteratively removes all nodes that do not fall within the 
   * sphere defined.
   */

  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;

  virtual void operator()(GraphType& g, double t){
    (void) t;
//--design_1
//-- You're missing nodes by using a for loop and incrementing at each step
//--START 
   for (auto i = g.node_begin(); i != g.node_end(); ++i){
      auto pos = (*i).position();
      auto sub = pos - c;
      auto this_vel = (*i).value().vel;
      if (norm(sub) < r){
        i = g.remove_node(i);
      }
    }
//--END
  }

};

struct CombinedConstraints {
  /**
   * @param[in]  g  Reference to graph where all nodes lie
   * @param[in]  t  time parameter

   * Applies all con
   * uses operator () to find the total force based on
   * the combination of all forces inside the vector.

  * @return total_force
  */

  std::vector<Constraint*> c_vec_;
  CombinedConstraints(std::vector<Constraint*> c_vec) : c_vec_(c_vec) {}

  void operator()(GraphType& g, double t){
      (void) t;
      for (auto c : c_vec_) {
        (*c)(g, t);
      }
  }

};

//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template<typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 constraint1, C2 constraint2){
  /**
   * @param[in] PinConstraint(), PlaneConstraint() or SphereConstraint() object
   *  Inputs two of the constraints that inherit from the Constraint class
      in an arbitrary order, and returns a CombinedConstraints functor
      object that applies these constraints to all nodes.
   * @return CombinedForce object.
  */
  std::vector<Constraint*> constraint_vec;

  constraint_vec.push_back(&constraint1);
  constraint_vec.push_back(&constraint2);

  return CombinedConstraints(constraint_vec);
}
//--END

template<typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 constraint1, C2 constraint2, C3 constraint3){
  /**
   * @param[in] PinConstraint(), PlaneConstraint() or SphereConstraint() object
   *  Inputs three of the constraints that inherit from the Constraint class
      in an arbitrary order, and returns a CombinedConstraints functor
      object that applies these constraints to all nodes.
   * @return CombinedForce object.
  */
  std::vector<Constraint*> constraint_vec;

  constraint_vec.push_back(&constraint1);
  constraint_vec.push_back(&constraint2);
  constraint_vec.push_back(&constraint3);

  return CombinedConstraints(constraint_vec);
}

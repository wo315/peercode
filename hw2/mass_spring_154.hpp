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

#include <iostream>

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};


/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;       //< spring constant
  double L;     //< initial length
  EdgeData() : K(0), L(1) {}
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
    //implement constraint directly for combiined force function
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      n.value().vel = Point(0,0,0);
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


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint constraint object definint the constraints 
 *@return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a cbonstraint must apply the constraint
 * 
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    //implement constraint directly for combiined force function
    //if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
    //  n.value().vel = Point(0,0,0);
    
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

    
  //apply constraints
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
    // HW2 #1: YOUR CODE HERE
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);

    //initialize f_spring
    Point f_spring = Point(0,0,0);
    Point x_i = n.position();
    for(auto e = n.edge_begin(); e!= n.edge_end(); ++e){
      Point x_j = (*e).node2().position();
      double K = (*e).value().K;
      double L = (*e).value().L;
      f_spring += -K *((x_i-x_j)/norm(x_i-x_j))*(norm(x_i-x_j)-L);
    }
    
    
    return f_spring + n.value().mass * Point(0,0,-grav);
  }
};


/**@brief Force Base function object,
 *@params[in] n  the node to apply the force to
 *@params[in] t  the current time 
 *@returns  this returns a force vector equal to zero
 */
struct Force{
  
  virtual Point const operator()(Node n, double t) const {
    (void) n;
    (void) t;
    return Point(0,0,0);
  }
};

/**@brief Gravity function object inherited from Force
 *        base function object 
 *@params[in] n  the node to apply the force to
 *@params[in] t  the current time 
 *@Returns force due to gravity 
 */
struct GravityForce:  Force{
  Point const operator()(Node n, double t) const {
    (void) t;
    return n.value().mass*Point(0,0,-grav);
  }
};

/**@brief spring force function object inherited from force
 *        base function object
 *@params[in] n  the node to apply the force to
 *@params[in] t  the current time  
 *@returns force due to spring 
 */
struct MassSpringForce:  Force{
  Point const operator()(Node n, double t) const {
    (void) t;
//--style_1
//--Redundant with pin constraint
//--START
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0))
      return Point(0,0,0);
//--END
    //initialize f_spring
    Point f_spring = Point(0,0,0);
    Point x_i = n.position();
    for(auto e = n.edge_begin(); e!= n.edge_end(); ++e){
      Point x_j = (*e).node2().position();
      double K = (*e).value().K;
      double L = (*e).value().L;
      f_spring += -K *((x_i-x_j)/norm(x_i-x_j))*(norm(x_i-x_j)-L);
    }
    return f_spring;
  }
};

/**@brief Damping force function object inherited from force
 *        base function object
 *@params[in] n  the node to apply the force to
 *@params[in] t  the current time 
 @returns force due to damping 
 */
struct DampingForce:  Force{
  double c_;

  DampingForce(): c_(1){
  };

  
  DampingForce(const double c): c_(c){
  };

  Point const operator()(Node n, double t){
    (void) t;
    return -c_*n.value().vel;
  }
};


/**@brief function object that takes an array of forces derived from the 
 *        base force class, and defines an operator that returns 
 *        the sum of the forces
 *@constructor params[in] forces an array of pointers to derivatives 
                          of the Force class
 *@operator params[in]  n the node to apply the force to
 *@operator params[in]  t the current time
 *@return the sum of the forces
 */
struct CombinedForce: Force{
  std::vector<const Force*> forces_;

  CombinedForce(std::vector<const Force*> forces):
    forces_(forces){};

  Point const operator()(Node n, double t) const {
    Point force_sum = Point(0,0,0);
    for(auto f = forces_.begin(); f != forces_.end(); ++f){
      force_sum += (*(*f))(n,t);
    }
    return force_sum;
  }
  
};

/**@brief funciton that combines forces
 *@params[in] the various forces
 *@return CombinedForce object
 */
//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
CombinedForce make_combined_force(const Force& force1, const Force& force2){  
    std::vector<const Force*> force_vector;
    force_vector.push_back(&force1);
    force_vector.push_back(&force2);
    return CombinedForce(force_vector);
};
//--END

//--design_-1
//--Good job handling combining rvalue forces and constraints 
//--END
CombinedForce make_combined_force(const Force& force1,
				  const Force& force2,
				  const Force& force3){  
    std::vector<const Force*> force_vector;
    force_vector.push_back(&force1);
    force_vector.push_back(&force2);
    force_vector.push_back(&force3);
    return CombinedForce(force_vector);
};


/**@brief constraint Base function object
 *@params[in] g   the graph
 *@params[in] t   the time
 *this functor does not apply any constraints
*/
struct Constraint{
  virtual void  operator()(GraphType& g, double t) const {
    (void) g;
    (void) t;
  }
};



/**@brief pin constraint class that derives from  Base functon object
          keeps specific nodes fixed
*constructor params[in] dt   the change in time
*cosntructor params[in] force  the force that will be applied to the node
*@operator params[in] g   the graph
*@opeartor params[in] t   the time
*this functor pins the two nodes by precalculating the value the velocity will
*           be update to, and assigning the velocity to negative of theat
            value, that way when the update occurs the vel is 0 and it never moves
*/
struct PinConstraint: Constraint{
  double dt_;
  CombinedForce force_;
  PinConstraint(double dt, CombinedForce force ):
    dt_(dt), force_(force){};
  void operator()(GraphType& g, double t) const {
    (void) t;
    for(auto n = g.node_begin(); n!= g.node_end(); ++n){
      auto node_ = *n;
      if (node_.position() == Point(0,0,0)){
	node_.value().vel = -force_(node_,t)*(dt_/node_.value().mass);
      }
      if (node_.position() == Point(1,0,0)){
	node_.value().vel = -force_(node_,t)*(dt_/node_.value().mass);
      }
    }
  }
};



/**@brief plane constraint class that derives from  Base functon object
 *fixes nodes that go into a certain plane
 *@params[in] g   the graph
 *@params[in] t   the time
 *this functor does not allow the cloth to go past a certain plane
 */
struct PlaneConstraint: Constraint{
  double z = -0.75;
  
  
  void  operator()(GraphType& g, double t) const {
    (void) t;
    for(auto n = g.node_begin(); n!= g.node_end(); ++n){
      auto node_ = *n;
      if (node_.position().z < z){
	node_.position().z=z;
	node_.value().vel.z = 0;
      }
    }
    
  }
};

/**@brief sphere constraint class that derives from  Base functon object
 *@params[in] g   the graph
 *@params[in] t   the time
 *this functor does not allow the cloth to enter the specified sphere by
 *fixing nodes at the surface of the sphere if they enter 
 */
struct SphereConstraint: Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  
  void operator()(GraphType& g, double t) const {
    (void) t;
    for(auto n = g.node_begin(); n!= g.node_end(); ++n){
      Node node_ = *n;
      Point x = node_.position();
      Point R = (x-c)/norm(x-c);
      
      if (norm(x-c) < r){
	node_.position() = c+r*R;
	node_.value().vel -= dot(node_.value().vel, R)*R;
	
      }
    }
    
  }
};

/**@brief tear constraint class that derives from  Base functon object
 *@params[in] g   the graph
 *@params[in] t   the time
 *this functor allows for the simulation of a tear by removing the nodes
 *and incident edges of nodes that enter the sphere
 */
struct TearConstraint: Constraint{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;
  
  void operator()(GraphType& g, double t) const {
    (void) t;
    auto n = g.node_begin();
    while(n!= g.node_end()){
      Node node_ = *n;
      Point x = node_.position();
      
      if (norm(x-c) < r){
	n = g.remove_node(n);
      }
      else{
	++n;
      }
    }
    
  }
};


/**@brief function object that takes an array of constraints derived from the 
 *        base constraint class, and defines an operator that applies 
 *        each constraint
 * @params[in] constraints an array of pointers to derivatives of 
 *             the Constraint class
 */
struct CombinedConstraints: Constraint{
  std::vector<const Constraint*> constraints_;

  CombinedConstraints(std::vector<const Constraint*> constraints):
    constraints_(constraints){};

  void operator()(GraphType& g, double t) const {
    for(auto c = constraints_.begin(); c != constraints_.end(); ++c){
      (*(*c))(g,t); //apply constraint
    }
  }
};




/**@brief function  that combines multiple constraints
 * @params[in] various constraints
 * @tparam C1 the type of constraint for the first input
 * @tparam C2 the type of constraint for the second input
 * @tparam C3 the type of constraint for the third input, 
 *            note if only two constraints are to be applied, 
 *            then making this 3rd template value Constraint(), 
 *            and only providing two inputs to the function will 
 *            provide the desired result
 * @return   CombinedConstraint object */
template<typename C1, typename C2, typename C3>
struct make_combined_constraint{
  const C1 cons1_;
  const C2 cons2_;
  const C3 cons3_;
  
  make_combined_constraint(const C1 cons1,
   			   const C2 cons2):
    cons1_(cons1), cons2_(cons2), cons3_(Constraint()){};

  make_combined_constraint(const C1 cons1,
  			   const C2 cons2,
  			   const C3 cons3):
  cons1_(cons1), cons2_(cons2), cons3_(cons3){};
  

  void operator()(GraphType& g, double t) const{
    std::vector<const Constraint*> cons_vector;
    cons_vector.push_back(&cons1_);
    cons_vector.push_back(&cons2_);
    cons_vector.push_back(&cons3_);
    return CombinedConstraints(cons_vector)(g,t); 

  }
};




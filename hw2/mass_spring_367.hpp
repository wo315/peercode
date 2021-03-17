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
using GraphType = Graph<NodeData,double>;
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
double symp_euler_step(G& g, double t, double dt, F force,C constraint) {
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
  constraint(g,t);

  return t + dt;
}

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
    Point gravF =  Point(0,0,-grav*n.value().mass);
    if (n.position() == Point(0,0,0)||n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }

  //compute spring force 
  Point springF = Point(0,0,0);
  for(auto negh_iter = n.edge_begin();negh_iter != n.edge_end();++negh_iter){
    auto edge_ = *negh_iter;//get edge
    Point edge_dir = Point(0,0,0);//find ouuuuuut direction 
    if(edge_.node1()==n){
      edge_dir = n.position()-edge_.node2().position();
    }
    else{
      edge_dir = n.position()-edge_.node1().position();
    }
    springF -=  100*edge_dir*(norm(edge_dir)-edge_.value())/norm(edge_dir);
  } 
  (void)t;
  return springF+gravF;

   
  }
};
















//force class used to inheritance
class Force{

  public:
  virtual Point  operator()(Node n,double  t){
    (void)n;(void)t;
    return   Point(0,0,0);
    }
};

//mass spring force
class MassSpringForce:public Force{
  public:
  template<typename NODE>
  Point operator()(NODE n,double t){
    (void) t;
    Point springF = Point(0,0,0);
  for(auto negh_iter = n.edge_begin();negh_iter != n.edge_end();++negh_iter){
    auto edge_ = *negh_iter;//get edge
    Point edge_dir = Point(0,0,0);//find ouuuuuut direction 
    if(edge_.node1()==n){
      edge_dir = n.position()-edge_.node2().position();
    }
    else{
      edge_dir = n.position()-edge_.node1().position();
    }
    springF -=  100*edge_dir*(norm(edge_dir)-edge_.value())/norm(edge_dir);
  } 
  return springF;
  }

};

//gravity force
class GravityForce:public Force{
  public:
  template<typename NODE>
  Point operator()(NODE n,double t){
    (void) t;
    return n.value().mass * Point(0,0,-grav);
  }
};

//damping force
class DampingForce:public Force{
  public:
  double damp_coef;
  template<typename NODE>
  Point operator()(NODE n,double t){
    (void) t;
    return -damp_coef * n.value().vel;
  }
};


//class used to combine forces
template<typename force_type_1, typename force_type_2>
struct AddForce{
  force_type_1 f1;
  force_type_2 f2;
  AddForce(force_type_1 f1v, force_type_2 f2v):f1(f1v),f2(f2v){}
  template<typename NODE>
  Point operator()(NODE n,double t){
    (void) t;
    return f1(n,t)+f2(n,t);
  }
};

//combine two forces
template<typename force_type_1, typename force_type_2>
AddForce<force_type_1,force_type_2> make_combined_force(force_type_1 f1,force_type_2 f2){
  return AddForce<force_type_1, force_type_2>(f1,f2);
}


//combine three forces
template<typename force_type_1, typename force_type_2,typename force_type_3>
AddForce<AddForce<force_type_1,force_type_2>,force_type_3> make_combined_force(force_type_1 f1,force_type_2 f2,force_type_3 f3){
  return AddForce<AddForce<force_type_1,force_type_2>, force_type_3>( AddForce<force_type_1, force_type_2>(f1,f2),f3);
}





//class to inheritance
class Constraint{
  public:
  virtual void operator()(GraphType& g, double t){
    (void)g;
    (void)t;
  }


};

//keep the nodes at 0,0,0 at 1,0,0
class PinConstraint:public Constraint{
  public:
  void operator()(GraphType & g, double t){
    (void) t;
    for(auto node_iter = g.node_begin();node_iter !=  g.node_end();++node_iter){
      auto n = *node_iter;
      if (n.position() == Point(0,0,0)||n.position() == Point(1,0,0)){
          n.value().vel = Point(0,0,0);
      }
    }
  }
};

//set the point to the nearest one on plane
//set z component of vel to zero
class PlaneConstraint:public Constraint{
  public:
  void operator()(GraphType & g, double t){
    (void) t;
    for(auto node_iter = g.node_begin();node_iter !=  g.node_end();++node_iter){
      auto n = *node_iter;
      if(n.position().z<-0.75){
        n.position().z=-0.75;
        n.value().vel.z = 0;
      }
      
    }
  }
};

//set the point to the nearest place on sphere
class SphereConstraint:public Constraint{
  public:
  Point center = Point(0.5,0.5,-0.5);
  double rad = 0.15;
  
  void operator()(GraphType & g, double t){
    (void) t;
    
    for(auto node_iter = g.node_begin();node_iter !=  g.node_end();++node_iter){
      auto n = *node_iter;
      if(norm(n.position()-center)<rad){
        Point center_dir = (n.position()-center)/norm(n.position()-center);
        n.value().vel -= dot(center_dir,n.value().vel)*center_dir;//project the vel
        n.position() = center+rad*center_dir;//move on to  sphere
      }
      
    }
  }
};

//remove the node in sphere (tear)
class SphereConstraint_remove:public Constraint{
  public:
   Point center = Point(0.5,0.5,-0.5);
  double rad = 0.15;
  
  void operator()(GraphType & g, double t){
    (void) t;
    
    for(auto node_iter = g.node_begin();node_iter !=  g.node_end();){
      auto n = *node_iter;
      if(norm(n.position()-center)<rad){
       g.remove_node(n);
      }
//--design_0
//--Careful I think g.node_end() becomes invalid when you delete nodes
//--START
      ++node_iter;
//--END
    }
  }

};

//class used to combine two constraints
template<typename constraint_type_1,typename constraint_type_2>
struct MergeConstraint{
  constraint_type_1 constraint1;
  constraint_type_2 constraint2;
  MergeConstraint(constraint_type_1 c1,constraint_type_2 c2):constraint1(c1),constraint2(c2){}
  void operator()(GraphType& g,double t){
    constraint1(g,t);
    constraint2(g,t);
  }

};

//combine two constraints
template<typename constraint_type_1,typename constraint_type_2>
MergeConstraint<constraint_type_1,constraint_type_2> make_combined_constraint(constraint_type_1 cons1,constraint_type_2 cons2){
  return MergeConstraint<constraint_type_1,constraint_type_2>(cons1,cons2);
}


//combine three constraints
template<typename constraint_type_1,typename constraint_type_2,typename constraint_type_3>
MergeConstraint<MergeConstraint<constraint_type_1,constraint_type_2>,constraint_type_3> make_combined_constraint(constraint_type_1 cons1,constraint_type_2 cons2,constraint_type_3 cons3){
  return make_combined_constraint(make_combined_constraint(cons1,cons2),cons3);
}






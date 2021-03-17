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

#define NO_OP(x) ((void)(x))


// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  bool pinned;     //< whether Node is pinned by constraint
  Point ref_pos;   //< reference position to which node is pinned
  NodeData() : vel(0), mass(1), pinned(0), ref_pos(0) {}
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

  // fix violating nodes
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    
  }

  return t + dt;
}


//
// FORCES
//

class Force {
  public:
    virtual Point operator()(Node n, double t){
        NO_OP(n); // silence compiler warnings
        NO_OP(t); 
        return Point(0, 0, 0);
    }
};

class GravityForce : public Force {

    Point operator()(Node n, double t) {
        NO_OP(t); 
        return Point(0, 0, -grav * n.value().mass);
    }
};

class MassSpringForce: public Force {

  double K = 100; // spring constant

  Point operator()(Node n, double t) {
    NO_OP(t); 
    Point spr_force(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
         spr_force += -K * ((*it).node1().position() - (*it).node2().position())
         / (*it).length() * ((*it).length() - (*it).value());
    }
    return spr_force;
  }
};

class DampingForce : public Force {

    double c = 0.05;  // damping constant
    Point operator()(Node n, double t) {
        NO_OP(t); 
    	return -c * n.value().vel;
    }
};

struct CombinedForce {

    CombinedForce(std::vector<Force*> force_vec_){
        force_vec = force_vec_;}
    std::vector<Force*> force_vec; // store pointers to component forces

    Point operator()(Node n, double t) {
        NO_OP(t); 
        Point force_sum(0, 0, 0);
        for (auto it = force_vec.begin(); it != force_vec.end(); ++it){
            force_sum += (**it)(n, t);
        }
        return force_sum;
    }
};



// function that combines an arbitrary number of forces
template <typename ...Force_Ensemble> // use template for parameter packing
    CombinedForce make_combined_force(Force_Ensemble... forces){
        std::vector<Force*> force_vec_ = {&forces...};
        return CombinedForce(force_vec_);
	}


//
// CONSTRAINTS
//

class Constraint {
  public:
    virtual void operator()(GraphType& g, double t){
        NO_OP(t); 
        NO_OP(g); 
    }
};

class PinConstraint : public Constraint {

    void operator()(GraphType& g, double t) {
        NO_OP(t); 
        if (t == 0){ // look for nodes to be pinned at the start
            for (auto it = g.node_begin(); it != g.node_end(); ++it){
                if ((*it).position() == Point(0) || 
                (*it).position() == Point(1, 0, 0)){
                    (*it).value().pinned = true;
                    (*it).value().ref_pos = (*it).position();
                    }
            }
         }
            else{ // search for the pinned nodes and keep them in place
                for (auto it = g.node_begin(); it != g.node_end(); ++it){
                    if ((*it).value().pinned){
                        (*it).value().vel = Point(0);
                        (*it).position() = (*it).value().ref_pos;
                    }
                }
            }
    }
};


class PlaneConstraint: public Constraint {

    void operator()(GraphType& g, double t) {
          NO_OP(t); 
          for (auto it = g.node_begin(); it != g.node_end(); ++it){
              if ((*it).position().z < -0.75){
                  (*it).position() = 
                      Point((*it).position().x, (*it).position().y, -0.75);
              }                                    
          }
    }
};

class SphereConstraint : public Constraint {

    Point center = Point(0.5, 0.5, -0.5);
    double rad = 0.15;
    void operator()(GraphType& g, double t) {
        NO_OP(t); 
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            if (norm((*it).position() - center) < rad){
                Point R_i = 
               ((*it).position() - center) / norm((*it).position() - center);
                (*it).position() = rad * R_i + center;
                (*it).value().vel - dot((*it).value().vel, R_i)*R_i;
            }
        }
    }
};


/* simple constraint that cuts off nodes whose z values sink below -0.6 */
class CutConstraint : public Constraint {
    void operator()(GraphType& g, double t) {
        NO_OP(t); 
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            if ((*it).position().z < -0.6){
                g.remove_node(it);
            }
        }
    }
};

class TearConstraint : public Constraint {

    Point center = Point(0.5, 0.5, -0.5);
    double rad = 0.15;
    void operator()(GraphType& g, double t) {
        NO_OP(t); 
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
//--design_1
//--Don't increment the iterator when you have removed a node
//--START
            if (norm((*it).position() - center) < rad){
                g.remove_node(it);
//--END
            }
        }
    }
}; 


struct CombinedConstraints {

    CombinedConstraints(std::vector<Constraint*> constr_vec_){
        constr_vec = constr_vec_;}
    std::vector<Constraint*> constr_vec;
    void operator()(GraphType& g, double t) { // deploy constraints sequentially
        for (auto it = constr_vec.begin(); it != constr_vec.end(); ++it){      
            (**it)(g, t);
        }
    }
};



// function that combines an arbitrary number of constraints	
template <typename ...Constraint_Ensemble> // use template for parameter packing
    CombinedConstraints make_combined_constr(Constraint_Ensemble... constrs){
        std::vector<Constraint*> constr_vec_ = {&constrs...};
        return CombinedConstraints(constr_vec_);
	}


//

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  double K = 100; // spring constant

  template <typename NODE>
  Point operator()(NODE n, double t) {

    // set force to zero for these points
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
        return Point(0, 0, 0);

    Point spr_force(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
         
         spr_force += -K * (((*it).node1().position() - (*it).node2().position())
         / (*it).length()) * ((*it).length() - (*it).value());

    }

    return (spr_force + Point(0, 0, -grav * n.value().mass));
  }
};




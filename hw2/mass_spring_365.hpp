/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

//--functionality_2
//--simulation failed
//--END

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr double K = 100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  Point og_pos;
  NodeData() : vel(0), mass(1), og_pos(0) {}
};





// Define the Graph type
//--design_1
//--did not implement edge-specific const K
//--START
using GraphType = Graph<NodeData, double>;
//--END
//using GraphType = Graph<NodeData>;
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
    //std::cout << n.position()[0] << " " << n.position()[1] << " " << n.position()[2] << std::endl;
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
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as C(g, t) where g is a reference to
 *          graph. Function operator() will apply the pre-set constraints to
 *          the graph.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
    //std::cout << n.position()[0] << " " << n.position()[1] << " " << n.position()[2] << std::endl;
  }

  constraint(g, t);


  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}





/** Force Parent Class
 * @brief The parent class of all forces
 *
*/
struct Force {

  //virtual method
  virtual Point operator()(Node n, double t){
    // supress compiler warning
    (void) t;
    (void) n;
    return Point(0);
  }

};

struct GravityForce : public Force {
  virtual Point operator()(Node n, double t){
    // supress compiler warning
    (void) t;
    // Get gravitational Force
    Point g_force = Point(0,0,-1);
    g_force *= n.value().mass * grav;
    return g_force;
  }


};


struct MassSpringForce : public Force {
    Point operator()(Node n, double t){
    // supress compiler warning
    (void) t;

    // Get spring force
    Point s_force = Point(0,0,0);
    Point x_i = n.position();

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      //
      Edge e_ij = (*it);
      Point x_j = e_ij.node2().position();
      double ij_norm = norm(x_i - x_j);
      double L = e_ij.value();

      s_force += -1 * K / ij_norm * (x_i - x_j) * (ij_norm - L);
    }
    return s_force;
  }

};


struct DampingForce : public Force {
    double c;

    // constructors, default c to be 1.0
    DampingForce(): c(1.0){}

    DampingForce(double c_input): c(c_input){}

    // OPerator overload
    Point operator()(Node n, double t){
    // supress compiler warning
    (void) t;

    // Get spring force
    Point v = n.value().vel;

    v *= -1.0 * c;

    return v;

  }

};


struct CombinedForce : public Force{
  std::vector<Force*> forces;

  CombinedForce(std::vector<Force*> input_forces) : forces(input_forces) {
  }

  Point operator()(Node n, double t){
    Point x = Point(0);
    for (auto it = forces.begin(); it != forces.end(); ++it){
      Point f = (*(*it))(n, t);
      x += f;
    }

    return x;
  }

};

// Comine two forces
//--design_1
//--argument object passed by value is destroyed after function finishes
//--START
template < class force_1, class force_2>
CombinedForce make_combined_force(force_1 f1, force_2 f2){
  std::vector<Force*> a {&f1, &f2};

  return CombinedForce(a);
}
//--END

// Comine three forces
template < class force_1, class force_2, class force_3>
CombinedForce make_combined_force(force_1 f1, force_2 f2, force_3 f3){
  std::vector<Force*> a {&f1, &f2, &f3};
  return CombinedForce(a);
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
    (void) t;
    // Constrain two corners of the cloth by returnning a zero force
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0);
    }

    // Get gravitational Force
    Point g_force = Point(0,0,-1);
    g_force *= n.value().mass * grav;


    // Get spring force
    Point s_force = Point(0,0,0);
    Point x_i = n.position();

    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
      //
      Edge e_ij = (*it);
      Point x_j = e_ij.node2().position();
      double ij_norm = norm(x_i - x_j);
      double L = e_ij.value();

      s_force += -1 * K / ij_norm * (x_i - x_j) * (ij_norm - L);
    }
    return s_force + g_force;
  }
};




//
/** Constraint Parent Class
 * @brief The parent class of all Constraints
 *
*/
struct Constraint {

  //virtual method
  virtual void operator()(GraphType& g, double t){
    // supress compiler warning
    (void) g;
    (void) t;
    std::cout << "Base Constraint Triggered" << std::endl;
  }

};



//
/** PinConstraint
 * @brief Fixing nodes at (0,0,0) and (1,0,0) by zeroing velocity
 *
*/
struct PinConstraint : public Constraint{

  //virtual method
  virtual void operator()(GraphType& g, double t){
    // supress compiler warning
    (void) t;

    // loop through all nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = (*it);
      // Constrain two corners of the cloth by returnning a zero velocity
      if (n.value().og_pos == Point(0,0,0) || n.value().og_pos == Point(1,0,0)){
        n.position() = n.value().og_pos;
      }
    }
  }
};


//
/** PlaneConstraint
 * @brief Keeping nodes away from plane z
 *
*/
struct PlaneConstraint : public Constraint{

  double z;

  PlaneConstraint() : z(-0.75) {}

  PlaneConstraint(double z_input) : z(z_input) {}

  //virtual method
  virtual void operator()(GraphType& g, double t){
    // supress compiler warning
    (void) t;
    // loop through all nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = (*it);
      // Check if node violate the condition
      if (inner_prod(n.position(), Point(0,0,1)) < z){
        //Setting the position to be the nearest point on the plane
        n.position()[2] = z;
        //Setting the z velocity to be 0
        n.value().vel[2] = 0;
      }

    }
  }
};




//
/** Sphere Constraint
 * @brief Sphere Constraint
 *
*/
struct SphereConstraint : public Constraint{

  Point c;
  double r;

  SphereConstraint() : c(Point(0.5,0.5,-0.5)), r(0.15) {}

  SphereConstraint(Point c_input, double r_input) : c(c_input), r(r_input) {}

  //virtual method
  virtual void operator()(GraphType& g, double t){
    // supress compiler warning
    (void) t;

    // loop through all nodes
    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = (*it);
      // Check if node violate the condition
      if (norm(n.position() - c) < r){
        //Setting the position to be the nearest point on the sphere
        n.position() = c + r / norm(n.position() - c) * (n.position() - c);

        //Setting the velocity
        Point r_i = (n.position() - c) / norm(n.position() - c);
        n.value().vel -= inner_prod(n.value().vel, r_i) * r_i;
      }

    }
  }
};


//
/** Sphere Constraint
 * @brief Sphere Constraint
 *
*/
struct SphereTearConstraint : public Constraint{

  Point c;
  double r;

  SphereTearConstraint() : c(Point(0.5,0.5,-0.5)), r(0.15) {}

  SphereTearConstraint(Point c_input, double r_input) : c(c_input), r(r_input){}

  //virtual method
  virtual void operator()(GraphType& g, double t){
    // supress compiler warning
    (void) t;

    // loop through all nodes

    for (auto it = g.node_begin(); it != g.node_end(); ++it){
      Node n = (*it);

      // Check if node violate the condition
      if (norm(n.position() - c) < r){
        g.remove_node(n);
        --it;
      }

    }
  }

};



struct CombinedConstraints : public Constraint {
  std::vector<Constraint*> constraints;

  CombinedConstraints(std::vector<Constraint*> input_constr)
    : constraints(input_constr) {
  }

  void operator()(GraphType& g, double t){
    (void) t;
    for (auto it = constraints.begin(); it != constraints.end(); ++it){
      (*(*it))(g,t);
    }
  }

};

// Comine two constraints
template < class constr_1, class constr_2>
CombinedConstraints make_combined_constraint(constr_1 c1, constr_2 c2){
  std::vector<Constraint*> a {&c1, &c2};
  return CombinedConstraints(a);
}

// Comine three constraints
template < class constr_1, class constr_2, class constr_3>
CombinedConstraints make_combined_constraint(constr_1 c1, constr_2 c2,
    constr_3 c3){
  std::vector<Constraint*> a {&c1, &c2, &c3};
  return CombinedConstraints(a);
}

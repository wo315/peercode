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
  Point init_pos;  //< Initial position
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double L;       //< Edge rest-length L_ij
  double K;     //< Edge own spring constant K_ij
  EdgeData() : L(0.5), K(100.0) {}
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
 * @param[in,out] g           Graph
 * @param[in]     t           The current time (useful for time-dependent forces)
 * @param[in]     dt          The time step
 * @param[in]     force       Function object defining the force per node
 * @param[in]     constraint  Constraint object defining the constraint for the graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a constraint object called as @a constraint(g, @a t),
 *           where g is a graph and @a t is the current time.
 *           @a constraint does not return anything but modifies the graph according
 *              to a certain constraint of being in a certain region of space (within a
 *              sphere or a semi-space).
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

  // We apply the constraint to the graph
  constraint(g, t);

  // Apply the constraint and compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}




/** Class Force for HW2 #3 */
class  Force {
  public:
    virtual Point operator()(Node n, double t){
      (void) n; (void) t; // silence compiler warnings
      return Point(0, 0, 0);
    }
};
  
/** Force function object for HW2 #1. */
struct Problem1Force {
  
  public:

    /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

    Point operator()(Node n, double t) {
      // HW2 #1: YOUR CODE HERE

      (void) t;    // silence compiler warnings

      // Defining the gravity force to which n is subject
      Point gravity = n.value().mass * Point(0, 0, -grav);

      // Contraining two corners of the cloth
      if (n.position() == Point (0 ,0 ,0) || n.position() == Point(1 ,0 ,0)){
        return Point(0, 0, 0);
      } 
      Point f_spring = Point(0, 0, 0);
      Point p = n.position();

      // Adding all the forces to which n is exposed
      double L;
      double K;
      Point p_j;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
      {

        p_j = (*it).node2().position();
        L = (*it).value().L;
        K = (*it).value().K;
        double d = norm(p - p_j);
        f_spring += -K * ((d - L)/d) * (p - p_j);
      }

      Point f_tot = f_spring + gravity;

      return f_tot;
    }
};

/** Force function object for Gravitational Force in HW2 #3. */
struct GravityForce : public Force {
  public:
    Point operator()(Node n, double t){
      (void) t; // silence compiler warnings
      return Point(0, 0, -n.value().mass * grav);
    }
};

/** Force function object for MassSpring Force in HW2 #3. */
struct MassSpringForce : public Force {
  public:
    /** Return the resulting force applied on node n */
    Point operator()(Node n, double t){

      (void) t; // silence compiler warnings

      Point f_spring = Point(0, 0, 0);
      Point p = n.position();
      
      // Contraining two corners of the cloth
//--style_1
//--This is redudant with Pin Constraint
//--START
      if (p == Point (0 ,0 ,0) || p == Point(1 ,0 ,0)){return f_spring;} 
//--END
      // Adding all the forces to which n is exposed
      double K;
      double L;
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
      {
        Point p_j = (*it).node2().position();
        double d = norm(p - p_j);
        K = (*it).value().K;
        L = (*it).value().L;
        f_spring += -K * ((d - L)/d) * (p - p_j);
      }
      return f_spring;
    }
};

/** Force function object for Damping Force in HW2 #3. */
struct DampingForce : public Force{
  public:
    const double c; // Damping constant

    DampingForce(const double C) : c(C) {}

    Point operator()(Node n, double t){
      (void) t; // silence compiler warnings
      return(-c * n.value().vel);
    }
};

/** Functor to combine forces */
struct CombinedForce{
  private:
    std::vector<Force*> F;

  public:
    CombinedForce(std::vector<Force*> Forces) : F(Forces) {}

    Point operator()(Node n, double t){
      Point res = Point(0, 0, 0);
      for (unsigned int i = 0; i < F.size(); i++)
      {
        res += (*(F[i]))(n, t);
      }
      return res;
    }
};


/** Combines three forces
 * @param[in]     F1    First Force object
 * @param[in]     F2    Second Force object
 * @param[in]     F3    Third Force object
 * 
 * @return A functor combining the effects of the three forces
 */
template <typename Force1, typename Force2, typename Force3>
CombinedForce make_combined_force(Force1 F1, Force2 F2, Force3 F3) {
  std::vector<Force*> F = std::vector<Force*>();
  F.push_back(&F1);
  F.push_back(&F2);
  F.push_back(&F3);
  return CombinedForce(F);
}


//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START

/** Combines two forces
 * @param[in]     F1    First Force object
 * @param[in]     F2    Second Force object
 * 
 * @return A functor combining the effects of the two forces
 */
template <typename Force1, typename Force2>
CombinedForce make_combined_force(Force1 F1, Force2 F2) {
  std::vector<Force*> F = std::vector<Force*>();
  F.push_back(&F1);
  F.push_back(&F2);
  return CombinedForce(F);
}

//--END

/** Class Constraint for HW2 #3 */
class Constraint{
  public:

    virtual void operator() (GraphType& g, double t) {
      (void) g; (void) t; // silence compiler warning
    };

};

/** Constraint function object for Pin Constraint in HW2 #3. */
struct PinConstraint : public Constraint{

  void operator() (GraphType& g, double t){
    (void) t;
    // Keeping fixed amounts to setting velocity to zero
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      if((*it).value().init_pos == Point(0, 0, 0)) {(*it).position() = Point(0, 0, 0);}
      if((*it).value().init_pos == Point(1, 0, 0)) {(*it).position() = Point(1, 0, 0);}
      //if((*it).position() == Point(0, 0, 0) || (*it).position() == Point(1, 0, 0)){
      //  (*it).value().vel = Point(0, 0, 0);
      //}
    }
  }
};


/** Constraint function object for Plane Constraint in HW2 #3. */
struct PlaneConstraint : public Constraint{

  public:
    const double alt = -0.75;

    void operator()(GraphType& g, double t){
      (void) t;
      for (auto it = g.node_begin(); it != g.node_end(); ++it)
      {
        Node curr_node = *it;
        if(curr_node.position().z < alt){
          // Project on the plane z = alt
          curr_node.position().z = alt;
          
          // Zero out v_z
          curr_node.value().vel.z = 0;
        } 
      }
    }
};


/** Constraint function object for Sphere Constraint in HW2 #3. */
struct SphereConstraint : public Constraint{

  public:
    const Point center = Point(0.5, 0.5, -0.5);
    const double radius = 0.15;

    void operator()(GraphType& g, double t){
      (void) t;
      Point p;
      for (auto it = g.node_begin(); it != g.node_end(); ++it)
      {
        Node curr_node = *it;
        p = curr_node.position();
        if(norm(p - center) < radius){
          // Project on the sphere
          Point dir = (p - center)/norm(p-center);
          curr_node.position() = center + radius * dir;
          
          // Zero out the component of v that is normal to the sphere
          Point v = curr_node.value().vel; 
          curr_node.value().vel -= (v.x*dir.x + v.y*dir.y + v.z*dir.z)*dir;
        } 
      }
    }  
};

/** Functor to combine constraints */
struct CombinedConstraints{
  private:
    std::vector<Constraint*> C;

  public:
    CombinedConstraints(std::vector<Constraint*> Constraints) : C(Constraints) {}

    void operator()(GraphType& g, double t){
      // Apply all your constraints successively
      for (unsigned int i = 0; i < C.size(); i++)
      {
        (*(C[i]))(g, t);
      }
    }
};


/** Combines three constraints
 * @param[in]     C1    First Constraint object
 * @param[in]     C2    Second Constraint object
 * @param[in]     C3    Third Constraint object
 * 
 * @return A functor combining the effects of the three constraints
 */
template <typename Constraint1, typename Constraint2, typename Constraint3>
CombinedConstraints make_combined_constraints(Constraint1 C1, Constraint2 C2, Constraint3 C3) {
  std::vector<Constraint*> C = std::vector<Constraint*>();
  C.push_back(&C1);
  C.push_back(&C2);
  C.push_back(&C3);
  return CombinedConstraints(C);
}

/** Combines two constraints
 * @param[in]     C1    First Constraint object
 * @param[in]     C2    Second Constraint object
 * 
 * @return A functor combining the effects of the two constraints
 */
template <typename Constraint1, typename Constraint2>
CombinedConstraints make_combined_constraints(Constraint1 C1, Constraint2 C2) {
  std::vector<Constraint*> C = std::vector<Constraint*>();
  C.push_back(&C1);
  C.push_back(&C2);
  return CombinedConstraints(C);
}



struct RemoveSphereConstraint : public Constraint{
  public:
    const Point center = Point(0.5, 0.5, -0.5);
    const double radius = 0.15;

    void operator()(GraphType& g, double t){
      (void) t;
      Point p;
//--design_1
//--You are skipping some nodes with this (the one you just swapped)
//--START
      for (auto it = g.node_begin(); it != g.node_end(); ++it)
      {
        Node curr_node = *it;
        p = curr_node.position();
        if(norm(p - center) < radius){
          // Violation of the constraint
          g.remove_node(curr_node);
        } 
      }
//--END
    }  
};





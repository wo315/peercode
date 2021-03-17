/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <stdarg.h>
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
  bool pin_bool;
  Point pin_loc;
  NodeData() : vel(0), mass(1), 
  pin_bool(false), pin_loc(0){}
};

struct EdgeData{
  double length;
  double spring_const;
  EdgeData() : length(0), spring_const(1) {}
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
    if (!(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))){
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

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint Constraint object defining constraints on the graph
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint){
    // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  g = constraint(g, t);

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
    (void) t;
    Point x_i = n.position();
    if (x_i == Point(0, 0, 0) || x_i == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }
    Point total_force = Point(0, 0, -grav * n.value().mass);

    //Iterate over adjacent edges/nodes to calculate spring force from the equation in 4.2
    for (auto edge_iter = n.edge_begin(); edge_iter != n.edge_end(); ++edge_iter){
      Point x_j = (*edge_iter).node2().position();
      double norm_ij = norm(x_i - x_j);
      total_force -= (*edge_iter).value().spring_const * ((x_i - x_j) / norm_ij)
       * (norm_ij - (*edge_iter).value().length);
    }
    return total_force;

  }
};

/** Parent force class that defaults to zero force */
class Force{
public:
  virtual Point operator()(Node n, double t){
    (void) n;
    (void) t;
    return Point(0, 0, 0);
  }
};

/** Gravity force class */
class GravityForce : public Force{
public:

  /** Computes gravity force
  *@param[in] n Node
  *@param[in] t The current time (useful for time-dependent forces)
  *@param[out] force The gravity force that node n experiences
  *@pre n.value().mass has been initialized
  *@pre grav has been defined as in static constexpr double grav = 9.81;
  *Complexity: O(1)
  */
  virtual Point operator()(Node n, double t){
    (void) t;
    return Point(0, 0, -grav * n.value().mass);
  }
};

/**Mass-spring force class*/
class MassSpringForce : public Force{
public:
  /** Computes mass-spring force
  *@param[in] n Node
  *@param[in] t The current time (useful for time-dependent forces)
  *@param[out] force The mass spring force that node n experiences
  *@pre edge.value().spring_const has been initialized
  *@pre edge.value().length has been initialized
  *@pre grav has been defined as in static constexpr double grav = 9.81;
  *Complexity: O(n.degree())
  */
  virtual Point operator()(Node n, double t){
    (void) t;
    Point mass_spring = Point(0., 0., 0.);

    Point x_i = n.position();
    for (auto edge_iter = n.edge_begin(); edge_iter != n.edge_end(); ++edge_iter){
      if ((*edge_iter).index() == -1){continue;}
      Point x_j = (*edge_iter).node2().position();
      double norm_ij = norm(x_i - x_j);
      Point update = (*edge_iter).value().spring_const * ((x_i - x_j) / norm_ij)
      * (norm_ij - (*edge_iter).value().length);
      mass_spring -= update;
    }
    return mass_spring;
  }
};

/** Damping force class*/
class DampingForce : public Force{
public:
  double damping_const_;

  /** Default initializes damping constant to 0.001*/
  DampingForce() : damping_const_ (0.001) {}

  /** Intializes when damping constant is specified*/
  DampingForce(double damping_const)
  : damping_const_(damping_const){}

  /** Computes damping force
  *@param[in] n Node
  *@param[in] t The current time (useful for time-dependent forces)
  *@param[out] force The damping force that node n experiences
  *Complexity: O(1)
  */
  virtual Point operator()(Node n, double t){
    (void) t;
    return -this->damping_const_ * n.value().vel;
  }

};

/** Struct to combine multiple forces to one force object*/
// Note: I didn't see the point in separating this into a ConmbinedForce struct
// and a separate function that just calls the initializer
struct make_combined_force{
  std::vector<Force*> forces_;

  make_combined_force(){}

  template <class ...Ts>
  make_combined_force(Ts... inputs)
    : forces_{&inputs...}{}

  template <typename NODE>
  Point operator()(NODE n, double t){
    (void) t;
    Point total_force = Point(0., 0., 0.);
    for (auto force_iter = this->forces_.begin();
     force_iter != this->forces_.end(); ++force_iter){
      total_force += (*(*force_iter))(n, t);
    }
    return total_force;
  }
};

/** Parent constraint class that imposes no constraints*/
class Constraint{
  public:
  virtual GraphType operator()(GraphType g, double t){
    (void) t;
    return g;
  }
};

/** Pin constraint class to pin corners of grid*/
class PinConstraint : public Constraint{
  public:

  /** Applies pin constraint by reseting position of pinned nodes
  *@param[in, out] g Graph
  *@param[in] t The current time (useful for time-dependent constraints)
  *Complexity: O(g.num_nodes())
  */
  virtual GraphType operator()(GraphType g, double t){
    if (t == 0.){
      for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter){
        Point x_i = (*node_iter).position();
        if (x_i == Point(0, 0, 0) || x_i == Point(1, 0, 0)){
          (*node_iter).value().pin_bool = true;
          (*node_iter).value().pin_loc = x_i;
        }
      }
    } else {
      for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter){
        if ((*node_iter).value().pin_bool == true){
          (*node_iter).position() = (*node_iter).value().pin_loc;
        }
      }
    }
    return g;
  }
};

/** Plane constraint class to ensure nodes stay above z = -.75 plane*/
class PlaneConstraint : public Constraint{
  public:

  /** Applies plane constraint by reseting position and velocity of violating nodes
  *@param[in, out] g Graph
  *@param[in] t The current time (useful for time-dependent constraints)
  *Complexity: O(g.num_nodes())
  */
  virtual GraphType operator()(GraphType g, double t){
    (void) t;
    for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter){
      Point x_i = (*node_iter).position();
      if (x_i.z < -0.75){
        (*node_iter).position().z = -0.75;
        (*node_iter).value().vel.z = 0.;
      }
    }
    return g;
  }

};

/** Sphere constraint class to ensure nodes stay out of sphere*/
class SphereConstraint : public Constraint{
  public:
  Point center_;
  double radius_;

  /**Default constructor sets center and radius as in HW2*/
  SphereConstraint() : center_(Point(0.5, 0.5, -0.5)), radius_(0.15) {}

  /**Constructor when center and radius are specified*/
  SphereConstraint(Point center, double radius) : center_(center), radius_(radius){}

  /** Applies sphere constraint by reseting position and velocity of violating nodes
  *@param[in, out] g Graph
  *@param[in] t The current time (useful for time-dependent constraints)
  *Complexity: O(g.num_nodes())
  */
  virtual GraphType operator()(GraphType g, double t){
    (void) t;
    for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter){
      Point x_i = (*node_iter).position();
      Point to_center = x_i - this->center_;
      double dist_to_center = norm(to_center);
      if (dist_to_center < this->radius_){
        Point R_i = to_center / dist_to_center;
        (*node_iter).position() = (this->radius_ * R_i) + this->center_; 
        (*node_iter).value().vel -= inner_prod((*node_iter).value().vel, R_i) * R_i;
      }
    }
    return g;
  }
};

/** Tear sphere constraint class to tear any nodes that come in contact with the sphere of death*/
class TearSphereConstraint : public Constraint{
  public:
  Point center_;
  double radius_;

  /**Default constructor sets center and radius as in HW2*/
  TearSphereConstraint() : center_(Point(0.5, 0.5, -0.5)), radius_(0.15) {}

  /**Constructor when center and radius are specified*/
  TearSphereConstraint(Point center, double radius) : center_(center), radius_(radius){}


  virtual GraphType operator()(GraphType g, double t){
    (void) t;
    for (auto node_iter = g.node_begin(); node_iter != g.node_end(); ++node_iter){
      if (int((*node_iter).index()) == -1){
        std::cout << "skipping" << std::endl;
        continue;}
      Point x_i = (*node_iter).position();
      Point to_center = x_i - this->center_;
      double dist_to_center = norm(to_center);
      if (dist_to_center < this->radius_){
        node_iter = g.remove_node(node_iter);
      }
    }
    return g;
  }
};

/** Struct to combine multiple constraints to one constraint object*/
// Note: I didn't see the point in separating this into a Conmbined Constraint struct
// and a separate function that just calls the initializer
struct make_combined_constraint{
  std::vector<Constraint*> constraints_;
  make_combined_constraint(){}

  template <class ...Ts>
  make_combined_constraint(Ts... inputs)
  : constraints_{&inputs...}{}

  GraphType operator()(GraphType g, double t){
    for (auto constr_iter = constraints_.begin(); 
      constr_iter != constraints_.end(); ++constr_iter){
      g = (*(*constr_iter))(g, t);
    }
    return g;
  }
};


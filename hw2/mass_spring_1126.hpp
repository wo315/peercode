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
  Point initial_position; //< Initial Position of each node
  NodeData() : vel(0), mass(1), initial_position(Point(1, 1, 1)) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;          //< Edge K coefficient
  double L0;     //< Edge length
  EdgeData() : K(100.0), L0(0.1) {}
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

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      continue; //fixed nodes, skipping
    }
    else{
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
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
    (void) t;    // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
      return Point(0,0,0); //fixed nodes, force is null
    }
    Point Fg = - grav * n.value().mass * Point(0, 0, 1);
    Point Fs = Point(0, 0, 0);
    for (auto itr_edge = n.edge_begin(); itr_edge != n.edge_end(); ++itr_edge){
      auto connected_edge = *itr_edge;
      Point edge_vec = Point(0, 0, 0);
      if (connected_edge.node1() == n){
        edge_vec = connected_edge.node2().position() - n.position();
      }
      else{
        std::cerr << "Target node should be node2" << std::endl;
        exit(0);
      }
      double dist = norm(edge_vec);
      Fs += connected_edge.value().K * (dist - connected_edge.value().L0) * edge_vec / dist;
    }
    return Fg + Fs;
  }
};

/** Base Class for a Force */
class Force {
  public:
    // Can't put a template here, virtual functions don't allow for templates
    virtual Point operator()(Node n, double t) {
      (void) n;
      (void) t;
      // std::cout << "why you here? " << std::endl;
      return Point(0, 0, 0);
    }

    virtual ~Force(){
    }

    virtual void identify(){
      std::cout << "Force" << std::endl;
    }
};
/**
 * @brief Gravity force inheriting from Force
 * 
 */
struct GravityForce: public Force{
  virtual Point operator()(Node n, double t){
    (void) t; // Silence compiler warning
    return -1 * grav * n.value().mass*Point(0, 0, 1);
  }

  virtual void identify(){
      std::cout << "GravityForce" << std::endl;
  }
};

/**
 * @brief MassSpringForce, Each edge act as a little spring trying pulling the 
 * edge to it's original position
 * 
 */
struct MassSpringForce: public Force{
  virtual Point operator()(Node n, double t){
    (void) t; // Silence compiler warning
    Point Fs = Point(0, 0, 0);
    for (auto edge_itr = n.edge_begin(); edge_itr != n.edge_end(); ++edge_itr){
      Edge connected_edge = *edge_itr;
      Point edge_vec = connected_edge.node2().position() - n.position();
      double dist = norm(edge_vec);
      Fs += connected_edge.value().K * (dist - connected_edge.value().L0) * edge_vec / dist;
    }
    return Fs;
  }
};

/**
 * @brief DampingForce, reduce the velocity of each nodes with a private
 * parameter c_. Default is 0.05
 * 
 */
class DampingForce: public Force{
  private:
    double c_;
  public:
    /** Constructor for valid friction coeficient */
    DampingForce(const double c) : c_(c){};
    /** Constructor for default friction coeficient */
    DampingForce() : c_(0.05){};
    /**
     * @brief Operator for DampingForce
     * return simple -c_ * vel
     * 
     * @param n : Node, node of interest where the force is to be computed
     * @param t : double, time t but unused
     * @return Point representing the Damping force
     */
    virtual Point operator()(Node n, double t){
      (void) t; // Silence compiler warning
      return -c_ * n.value().vel;
    }
};

/**
 * @brief Functor to combine forces between. is called by 
 * make_combined_force
 * 
 */
struct CombinedForce{
  CombinedForce(std::vector<Force*> F): F_(F){}
  Point operator()(Node n, double t) const {
    (void) t;
    Point res_force = Point(0, 0, 0);
    for (unsigned int i = 0; i < F_.size(); ++i){
      res_force += (*F_[i])(n, t);
    }
    return res_force;
  }

  private:
    std::vector<Force*> F_;
};

/**
 * @brief make a CombinedForce functor which defines an operator obtained by
 * the addition of the force @a f1 and @a f2 in each points
 * 
 * @tparam force1, class inheriting from Force
 * @tparam force2, class inheriting from Force
 * @param f1 Force
 * @param f2 Force
 * @return CombinedForce 
 */
//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template<typename force1, typename force2>
CombinedForce make_combined_force(force1 f1, force2 f2) {
  Force* F1 = &f1;
  Force* F2 = &f2;
  std::vector<Force*> F;
  F.push_back(F1);
  F.push_back(F2);
  return CombinedForce(F);
}
//--END
/**
 * @brief make a CombinedForce functor which defines an operator obtained by
 * the addition of the force @a f1, @a f2 and @a f2 in each points
 * 
 * @tparam force1 
 * @tparam force2 
 * @tparam force3 
 * @param f1 Force
 * @param f2 Force
 * @param f3 Force
 * @return CombinedForce functor representing the addition of the 3 forces
 */
template<typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 f1, force2 f2, force3 f3) {
  Force* F1 = &f1;
  Force* F2 = &f2;
  Force* F3 = &f3;
  std::vector<Force*> F;
  F.push_back(F1);
  F.push_back(F2);
  F.push_back(F3);
  return CombinedForce(F);
}

/**
 * @brief Base Class for Constraint
 */
class Constraint {
  public:
    virtual void operator()(GraphType& g, double t) {
      (void) g;
      (void) t;
    }

    virtual void identify(){
      std::cout << "Constraint" << std::endl;
    }

    virtual ~Constraint(){}
};

/**
 * @brief PinConstraint that is going to enforce the position of the two points
 * P(0, 0, 0) and P(1, 0, 0). 
 * @pre The two anchor points are defined as initial position in the graph. 
 * @post the two anchor point will be fixed
 * If the anchor points aren't defined, then the constraint is void
 */
class PinConstraint : public Constraint {
  virtual void operator()(GraphType& g, double t){
    (void) t;
    for(auto node_iterator = g.node_begin(); node_iterator != g.node_end(); ++node_iterator){
      if ((*node_iterator).value().initial_position == Point(0, 0, 0)){
        (*node_iterator).position() = Point(0, 0, 0);
        (*node_iterator).value().vel = Point(0, 0, 0);
      }
      if ((*node_iterator).value().initial_position == Point(1, 0, 0)){
        (*node_iterator).position() = Point(1, 0, 0);
        (*node_iterator).value().vel = Point(0, 0, 0);
      }
    }
  }
  virtual void identify(){
      std::cout << "PinConstraint" << std::endl;
  }
};

/**
 * @brief Plane Constraint, prohibits the half-space z<-0.75
 */
class PlaneConstraint : public Constraint{
  public:
  void operator()(GraphType& g, double t){
    (void) t;
    for(auto node_iterator = g.node_begin(); node_iterator != g.node_end(); ++node_iterator){
      if((*node_iterator).position().z < -0.75)
      {
        (*node_iterator).value().vel.z = 0;
        (*node_iterator).position().z = -0.75;
      }
    }
  }
  virtual void identify(){
    std::cout << "PlaneConstraint" << std::endl;
  }
};

/**
 * @brief Prohibits the sphere ((0.5, 0.5, -0.5), 0.15)
 * If a nodes are found in that sphere they are cast on the sphere and their
 * velocity is set to 0
 */
class SphereConstraint : public Constraint{
  public:
  void operator()(GraphType& g, double t){
    (void) t;
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    for(auto node_iterator = g.node_begin(); node_iterator != g.node_end(); ++node_iterator){
      double distance_to_center = norm((*node_iterator).position() - center);
      if (distance_to_center < radius){
        Point R = ((*node_iterator).position() - center)/norm((*node_iterator).position() - center);
        (*node_iterator).position() = center + radius*R;
        (*node_iterator).value().vel -= dot((*node_iterator).value().vel, R)*R;
      }
    }
  }
  virtual void identify(){
    std::cout << "SphereConstraint" << std::endl;
  }
};

/** Functor to combine constraints */
struct CombinedConstraints {
  public:
    std::vector<Constraint*> C_;
    CombinedConstraints(std::vector<Constraint*> constraints) 
      : C_(constraints) {}
    void operator()(GraphType &g, double t) {
      for (unsigned int i = 0; i < C_.size(); i++) {
        (*(C_[i]))(g, t);
      }
    }
};

/**
 * @brief Functions to combine 2 constraints
 * 
 * @tparam constraint1, a constraint inheriting from Constraint class
 * @tparam constraint2, a constraint inheriting from Constraint class
 * @param c1 Constraint
 * @param c2 Constraint
 * @return CombinedConstraints returns a functor CombinedConstraints 
 * implementing the addition of constraint
 */
template <typename constraint1, typename constraint2>
CombinedConstraints make_combined_constraints(constraint1 c1, constraint2 c2) {
  Constraint* C1 = &c1;
  Constraint* C2 = &c2;
  std::vector<Constraint*> C;
  C.push_back(C1);
  C.push_back(C2);

  return CombinedConstraints(C);
}

/**
 * @brief Functions to combine 3 constraints
 * 
 * @tparam constraint1, a constraint inheriting from Constraint class
 * @tparam constraint2, a constraint inheriting from Constraint class
 * @tparam constraint3, a constraint inheriting from Constraint class
 * @param c1 Constraint
 * @param c2 Constraint
 * @param c3 Constraint
 * @return CombinedConstraints returns a functor CombinedConstraints 
 * implementing the addition of constraint
 */
template <typename constraint1, typename constraint2, typename constraint3>
CombinedConstraints make_combined_constraints(constraint1 c1, constraint2 c2, 
                                              constraint3 c3) {
  Constraint* C1 = &c1;
  Constraint* C2 = &c2;
  Constraint* C3 = &c3;
  std::vector<Constraint*> C;
  C.push_back(C1);
  C.push_back(C2);
  C.push_back(C3);

  return CombinedConstraints(C);
}

/**
 * @brief Prohibits the sphere ((0.5, 0.5, -0.5), 0.15)
 * If a node violates this constraint it is removed from the graph
 */
class BreakSphereConstraint : public Constraint{
  public:
  void operator()(GraphType& g, double t){
    (void) t;
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    auto node_iterator = g.node_begin();
    // auto end_iterator = g.node_end();
    while (node_iterator != g.node_end()){
      auto current_node = *node_iterator;
      double distance_to_center = norm(current_node.position() - center);
      if (distance_to_center < radius){
        node_iterator = g.remove_node(node_iterator);
      }
      else{
        ++node_iterator;
      }
    }
  }
  virtual void identify(){
    std::cout << "BreakSphereConstraint" << std::endl;
  }
};

/**
 * @brief Implements an euler method step with both forces and constraints
 * 
 * @tparam G template graphe
 * @tparam F template force
 * @tparam C template constraint
 * @param g reference to the graph
 * @param t double of time
 * @param dt double of time step for simulation
 * @param force Force to be applied to the graph (Can be a Functor of summed 
 * forces with make_combined_force)
 * @param constraint Constraint to be applied to the graph (Can be a Functor of summed 
 * constraints with make_combined_force)
 * @return double t + dt
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

  constraint(g, t);
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }
  return t + dt;
}

/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <unordered_map>

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
  Point initial_position; // add
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
  double K;
  double L;
  EdgeData() : K(0), L(0) {}
};

// Define the Graph type
// using GraphType = Graph<NodeData>;
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

struct CombinedConstraints;

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
  return symp_euler_step(g, t, dt, force, CombinedConstraints());
}

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraints after updating the node positions
  // but before calculating the forces applied
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  // constraint(g, t);

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
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }

    else{
      Point f_spring = Point(0, 0, 0);
      for(auto incident_iter = n.edge_begin(); incident_iter != n.edge_end(); ++incident_iter){
        Point adj_n = (*incident_iter).node2().position();
        double K = (*incident_iter).value().K;
        double L = (*incident_iter).value().L;
        f_spring += K * (norm(adj_n-n.position())-L) * (adj_n - n.position())/norm(adj_n-n.position());
      }
      return f_spring + n.value().mass*Point(0, 0, -grav);
    }
  }
};

// FORCE
class Force{
  public:
    virtual Point operator()(Node n, double t){
      (void) n; // no-op cast to silence any compiler warnings
      (void) t; // no-op cast to silence any compiler warnings
      return Point(0);
    }
    virtual ~Force() {} // add
};

class GravityForce: public Force{
  public:
    Point operator()(Node n, double t){
      (void) t; // no-op cast to silence any compiler warnings
      return n.value().mass*Point(0, 0, -grav);
    }
};

class MassSpringForce: public Force{
  public:
    Point operator()(Node n, double t){
      (void) t; // no-op cast to silence any compiler warnings
      Point f_spring = Point(0, 0, 0);
      for(auto incident_iter = n.edge_begin(); incident_iter != n.edge_end(); ++incident_iter){
        Point adj_n = (*incident_iter).node2().position();
        double K = (*incident_iter).value().K;
        double L = (*incident_iter).value().L;
        Point f_spring_ = K * (norm(adj_n-n.position())-L) * (adj_n - n.position())/norm(adj_n-n.position());
        f_spring += f_spring_;
      }
      return f_spring;
    }
};

class DampingForce: public Force{
  public:
    Point operator()(Node n, double t){
      (void) t; // no-op cast to silence any compiler warnings
      return -c_*n.value().vel;
    }
  private:
    double c_;
};

struct CombinedForce{
  std::vector<Force*> forces_;

  Point operator()(Node n, double t){
    Point combined_force(0);
    for(unsigned int i=0; i<forces_.size(); i++){
      combined_force += (*forces_.at(i))(n, t);
    }
    return combined_force;
  }
};

struct make_combined_force{
  CombinedForce combined_force;

  template <typename force1_type, typename force2_type, typename force3_type>
  make_combined_force(force1_type force1, force2_type force2, force3_type force3){
    combined_force.forces_.push_back(&force1);
    combined_force.forces_.push_back(&force2);
    combined_force.forces_.push_back(&force3);
  }

  template <typename force1_type, typename force2_type>
  make_combined_force(force1_type force1, force2_type force2){
    combined_force.forces_.push_back(&force1);
    combined_force.forces_.push_back(&force2);
  }

  Point operator()(Node n, double t){
    return combined_force(n, t);
  }
};


// CONSTRAINTS
/**
 * Parent constraint class with a virtual operator().
 */
class Constraint{
  public:
    virtual void operator()(GraphType& g, double t){
      (void) g; // no-op cast to silence any compiler warnings
      (void) t; // no-op cast to silence any compiler warnings
    }
    virtual ~Constraint() {}
};

/**
 * Keeps the nodes at (0, 0, 0) and (1, 0, 0) fixed.
 */
class PinConstraint: public Constraint{
  public:
    void operator()(GraphType& g, double t){
      (void) t; // no-op cast to silence any compiler warnings
      for(auto it=g.node_begin(); it!=g.node_end(); ++it){
        if((*it).value().initial_position == Point(0, 0, 0) || (*it).value().initial_position == Point(1, 0, 0)){
          (*it).position() = (*it).value().initial_position;
          (*it).value().vel = Point(0);
        }
      }
    }
};

/**
 * Sets the position to the nearest point on the plane.
 * Sets the z-component of the Node velocity to zero.
 */
class PlaneConstraint: public Constraint{
  public:
    void operator()(GraphType& g, double t){
      (void) t; // no-op cast to silence any compiler warnings
      for(auto it=g.node_begin(); it!=g.node_end(); ++it){
        if((*it).position().z < -0.75){
          (*it).position().z = -0.75;
          (*it).value().vel.z = 0;
        }
      }
    }
};

/**
 * Set the position to the nearest poin on the surface of the sphere.
 * Set the component of the velocity that is normal to the sphere surface to zero.
 */
class SphereConstraint: public Constraint{
  private:
    Point center_;
    double radius_;

  public:   
    SphereConstraint() : center_({0.5, 0.5, -0.5}), radius_(0.15) {};

    void operator()(GraphType& g, double t){
      (void) t; // no-op cast to silence any compiler warnings
      for(auto it=g.node_begin(); it!=g.node_end(); ++it){
        if(norm((*it).position() - center_) < radius_){
          Point unit_vector = ((*it).position() - center_) / norm((*it).position() - center_);
          (*it).position() = center_ + radius_ * unit_vector;
          (*it).value().vel -= ((*it).value().vel * unit_vector) * unit_vector;
        }
      }
    }
};

/**
 * CombinedConstraints functor
 * Takes in 2 or 3 constraints in an arbitrary order and applies them to the graph.
 */
struct CombinedConstraints{
  std::vector<Constraint*> constraints_;

  CombinedConstraints() {};

  template <typename constraint1_type, typename constraint2_type>
  CombinedConstraints(constraint1_type constraint1, constraint2_type constraint2){
    constraints_.push_back(&constraint1);
    constraints_.push_back(&constraint2);
  }

  template <typename constraint1_type, typename constraint2_type, typename constraint3_type>
  CombinedConstraints(constraint1_type constraint1, constraint2_type constraint2, constraint3_type constraint3){
    constraints_.push_back(&constraint1);
    constraints_.push_back(&constraint2);
    constraints_.push_back(&constraint3);
  }

  template <typename G>
  void operator()(G& g, double t){
    for(unsigned int i=0; i<constraints_.size(); i++){
      // Apply the constraints sequentially
      (*constraints_.at(i))(g, t);
    }
  }
};

/**
 * Similar to SphereConstraint except that the parts of the grid
 * that touch the sphere should disappear
 */
class TearConstraint: public Constraint{
  private:
    Point center_;
    double radius_;

  public:   
    TearConstraint() : center_({0.5, 0.5, -0.5}), radius_(0.15) {};

    void operator()(GraphType& g, double t){
      (void) t; // no-op cast to silence any compiler warnings
      for(auto it=g.node_begin(); it!=g.node_end(); ++it){
        if(norm((*it).position() - center_) < radius_){
          g.remove_node(*it);
//--design_1
//--If you remove a node, the it is already updated to point to a new node
//-- so you don't need to increment it
//--END
        }
      }
    }
};

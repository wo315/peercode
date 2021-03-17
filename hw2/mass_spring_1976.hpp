/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <chrono>
#include <fstream>
#include <memory>
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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  
  // Apply the constraints
  constraints(g, t);

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

/** Force function object for HW2 #2. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  double K = 100;

  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #2: YOUR CODE HERE
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0,0,0);
    }
    Point force_spring = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto edge = *it;
      double rest_length = edge.value();
      NODE n2 = edge.node2();
      Point diff = n.position() - n2.position();
      double dist = norm(diff);
      force_spring += -1 * K * diff * (dist - rest_length) / dist;
    }

    Point force_gravity = n.value().mass * Point(0,0,-grav);
    (void) t; // silence compiler warning
    return force_spring + force_gravity;

  }
};

class Force {

 public:
  virtual Point operator()(Node n, double t) const {
    // std::cout<< "inside Force::operator()" <<std::endl;
    (void) n; (void) t; // silence compiler warning
    return Point(0,0,0);
  }
  
  virtual ~Force() {
  }

};

class GravityForce: public Force{

 public:

  Point operator()(Node n, double t) const {
    // std::cout<< "inside GravityForce::operator()" <<std::endl;
    Point force = n.value().mass * Point(0,0,-grav);
    (void) t; // silence compiler warning
    return force;
  }
   
};

class MassSpringForce: public Force {

  double K_ = 100;

 public:

  Point operator()(Node n, double t) const {
    // std::cout<< "inside MassSpringForce::operator()" <<std::endl;
    Point force = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      auto edge = *it;
      double rest_length = edge.value();
      Node n2 = edge.node2();
      Point diff = n.position() - n2.position();
      double dist = norm(diff);
      force += -1 * K_ * diff * (dist - rest_length) / dist;
    }
    (void) t; // silence compiler warning
    return force;
  }

};

class DampingForce: public Force {

  double c_ = 0.05;  

 public:

  Point operator()(Node n, double t) const {
    // std::cout<< "inside DampingForce::operator()" <<std::endl;
    Point force = -1 * c_ * n.value().vel;
    (void) t; // silence compiler warning
    return force;
  }

};

class CombinedForce {
  
  std::vector<const Force*> forces_{};

 public:

  CombinedForce(const Force& force_1, const Force& force_2, const Force& force_3) {
    forces_.push_back(&force_1);
    forces_.push_back(&force_2);
    forces_.push_back(&force_3);
  }

  CombinedForce(const Force& force_1, const Force& force_2) {
    forces_.push_back(&force_1);
    forces_.push_back(&force_2);
  }  

  Point operator()(Node n, double t) {
    Point force = Point(0,0,0);
    for (unsigned int i = 0; i < forces_.size(); ++i) {
      force += (*forces_.at(i))(n, t);
    }
    return force;
  }

};

CombinedForce make_combined_force(const Force& force_1, const Force& force_2, const Force& force_3) {
  CombinedForce force = CombinedForce(force_1, force_2, force_3);
  return force;
}

CombinedForce make_combined_force(const Force& force_1, const Force& force_2) {
  CombinedForce force = CombinedForce(force_1, force_2);
  return force;
}


class Constraint {
 public:
  virtual void operator()(GraphType& g, double t) const = 0;
  virtual ~Constraint() {}
};


class PinConstraint: public Constraint {
  
  unsigned i1_, i2_;
  
 public:
  
  PinConstraint(unsigned i1, unsigned i2): i1_(i1), i2_(i2) {
  }

  void operator()(GraphType& g, double t) const {
    // fetch the node objects
    auto n1 = g.node(i1_);
    auto n2 = g.node(i2_);
    // set the velocity to 0
    n1.value().vel = Point(0,0,0);
    n2.value().vel = Point(0,0,0);
    // set their position back to their original position
    n1.position()  = Point(0,0,0);
    n2.position()  = Point(1,0,0);
   (void) t; // silence compiler warning
  }

};

class PlaneConstraint: public Constraint {

  Point plane_vector_ = Point(0,0,1);
  double thresh_      = -0.75; 

 public:

  void operator()(GraphType& g, double t) const {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // if the plane constraint is violated
      if (dot(n.position(), plane_vector_) < thresh_) {
        // set the position of the node to the closest point on the plane
        Point old_pos = n.position();
        n.position()  = Point(old_pos[0], old_pos[1], thresh_);
        
        // set the component of the velocity perpendicular to the plane
        // to 0
        Point old_vel = n.value().vel;
        n.value().vel = Point(old_vel[0], old_vel[1], 0);
      }
    }
    (void) t; // silence compiler warning
  }

};

class SphereConstraint: public Constraint {

  Point center_  = Point(0.5,0.5,-0.5);
  double radius_ = 0.15;

 public:

  void operator()(GraphType& g, double t) const {
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      // if the sphere constraint is violated
      if (norm(n.position() - center_) < radius_) {
        // set the position of the node to the closest point on the sphere
        Point old_pos = n.position();
        n.position()  = (old_pos - center_) / norm(old_pos - center_)  * radius_ + center_;
        
        // set the component of the velocity perpendicular to the sphere
        // to 0
        Point old_vel  = n.value().vel;
        Point new_pos  = n.position();
        Point norm_vec = (new_pos - center_) / norm(new_pos - center_);
        n.value().vel  = old_vel - (dot(old_vel, norm_vec) * norm_vec);
      }
    }
    (void) t; // silence compiler warning 
  }

};


class TearConstraint: public Constraint {

  Point center_  = Point(0.5,0.5,-0.5);
  double radius_ = 0.15;

 public:

  void operator()(GraphType& g, double t) const {
    bool has_removed = true;
    while (has_removed) {
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        has_removed = false;
        auto n = *it;
        // if the sphere constraint is violated
        if (norm(n.position() - center_) < radius_) {
          // remove the nodes
          g.remove_node(n);
          has_removed = true;
          break;
        }
      }
    }
    (void) t; // silence compiler warning 
  }

};



class CombinedConstraints {
  
  std::vector<const Constraint*> constraints_ {}; 
 
 public:
  
  CombinedConstraints(const Constraint& c1, const Constraint& c2, const Constraint& c3) {
    constraints_.push_back(&c1);
    constraints_.push_back(&c2);
    constraints_.push_back(&c3);
  }

  CombinedConstraints(const Constraint& c1, const Constraint& c2) {
    constraints_.push_back(&c1);
    constraints_.push_back(&c2); 
  }

  virtual void operator()(GraphType& g, double t) {
    for (unsigned i = 0; i < constraints_.size(); ++i) {
      (*constraints_.at(i))(g, t);
    }
    (void) t; // silence compiler warning
  }

};

CombinedConstraints make_combined_constraints(const Constraint& c1, const Constraint& c2, const Constraint& c3) {
  CombinedConstraints constraints = CombinedConstraints(c1, c2, c3);
  return constraints;
}

CombinedConstraints make_combined_constraints(const Constraint& c1, const Constraint& c2) {
  CombinedConstraints constraints = CombinedConstraints(c1, c2);
  return constraints;
}

//--functionality_1
//--Your code does ot work if Tearconstraint is applied.
//--END

//--design_0
//--Well designed!
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END



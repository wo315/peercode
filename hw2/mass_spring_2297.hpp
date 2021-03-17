/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>

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
  Point initial_position;
  NodeData() : vel(0), mass(1), initial_position(Point(0, 0, 0)) {}
};

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double K;       //< Spring constant
  double L;       //< Rest length
  EdgeData() : K(100), L(1) {}
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

    // handle fixed points
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      n.value().vel += Point(0, 0, 0);
    }
    else {
      // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
      n.value().vel += force(n, t) * (dt / n.value().mass);
    }
  }

  return t + dt;
}

/** Computes node position and velocity update for next time step */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
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

/** Force function object */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  public:
    template <typename NODE>
    Point operator()(NODE n, double t) {
      (void) t; // silence warnings
      // prevent cloth from falling to infinity
      if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
        return Point(0, 0, 0);
      }

      // get spring force
      Point SpringForce = Point(0, 0, 0);
      Point xi = n.position();
      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Point xj = (*it).node2().position();
        double K = (*it).value().K;
        double L = (*it).value().L;

        SpringForce += - K * (xi - xj) / norm(xi - xj) * (norm(xi - xj) - L);
      }
      // get gravitational force
      double mass = n.value().mass;
      Point GravForce = mass * Point(0, 0, -grav);

      // return sum of two forces
      return SpringForce + GravForce;
    }
};

/** @class Force
 * @brief A parent class for various forces on nodes
 *
 * Users can use this class as a parent class from which to 
 * inherit a () operator for forces.
 * 
 * This force object will return a 0 force vector
 * 
 */
class Force {
  public:

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty force object. */
    Force() {}

    /** Default destructor */
    virtual ~Force() {}

    //
    // METHODS
    //

    // Method to return 0 force
    virtual Point operator()(Node &n, double t) {
      (void) n; // silence compiler warnings
      (void) t; // silence compiler warnings
      return Point(0, 0, 0);
    }
};

/** @class GravityForce
 * @brief Represents a gravity force
 * 
 * This force object will return the force due to gravity
 * on a node in the negative z direction using the () operator.
 * 
 */
class GravityForce: public Force {
  public:

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty force object. */
    GravityForce() {}

    /** Default destructor */
    virtual ~GravityForce() {}

    //
    // METHODS
    //

    // Method to return 0 force
    virtual Point operator()(Node &n, double t) {
      (void) t; // silence warnings
      double mass = n.value().mass;
      return mass * Point(0, 0, -grav);
    }
};

/** @class MassSpringForce
 * @brief Represents the force applied by a spring
 *
 * This force object will return the spring force applied 
 * to a node by an adjacent edge using the () operator.
 * 
 */
class MassSpringForce: public Force {
  public:

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty force object. */
    MassSpringForce() {}

    /** Default destructor */
    virtual ~MassSpringForce() {}

    //
    // METHODS
    //

    // Method to return 0 force
    virtual Point operator()(Node &n, double t) {
      (void) t; // silence warnings
      Point SpringForce = Point(0, 0, 0);
      Point xi = n.position();

      for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        Point xj = (*it).node2().position();
        double K = (*it).value().K;
        double L = (*it).value().L;

        SpringForce += - K * (xi - xj) / norm(xi - xj) * (norm(xi - xj) - L);
      }

      return SpringForce;
    }
};

/** @class DampingForce
 * @brief Represents the force applied by a spring
 *
 * This force object will dampen the velocity of 
 * a node by a factor of c using the () operator.
 * 
 */
class DampingForce: public Force {
  public:

    //
    // CONSTRUCTORS AND DESTRUCTOR
    //

    /** Construct an empty force object. */
    DampingForce(const double c = 0.01) : c_(c) {}

    /** Default destructor */
    virtual ~DampingForce() {}

    //
    // METHODS
    //

    // Method to return 0 force
    virtual Point operator()(Node &n, double t) {
      (void) t; // silence warnings
      Point velocity = n.value().vel;
      return this->c_ * velocity * -1;
    }

  private:
    const double c_;
};

/** Functor to combine forces */
struct CombinedForce {
  public:
    // Constructor and Attributes
    std::vector<Force*> forces_;

    CombinedForce(std::vector<Force*> forces) : forces_(forces) {}

    // Methods
    Point operator()(Node &n, double t) {
      Point total_force = Point(0, 0, 0);
      for (unsigned int i = 0; i < this->forces_.size(); i++) {
        total_force += (*(forces_[i]))(n, t);
      }

      return total_force;
    }
};


//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
/** Function to combine 2 forces */
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2) {
  Force* force1 = &f1;
  Force* force2 = &f2;
  std::vector<Force*> force_vector;
  force_vector.push_back(force1);
  force_vector.push_back(force2);

  return CombinedForce(force_vector);
}
//--END

/** Function to combine 3 forces */
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, 
                                  F3 f3) {
  Force* force1 = &f1;
  Force* force2 = &f2;
  Force* force3 = &f3;
  std::vector<Force*> force_vector;
  force_vector.push_back(force1);
  force_vector.push_back(force2);
  force_vector.push_back(force3);

  return CombinedForce(force_vector);
}

/** @class Constraint
 * @brief A parent class for various constraints on nodes
 *
 * Users can use this class as a parent class from which to 
 * inherit a () operator for constraints.
 * 
 * This constraint object applies no constraint.
 * 
 */
class Constraint {
  public:

    /** Construct an empty force object. */
    Constraint() {}

    /** Default destructor */
    virtual ~Constraint() {}

    // Method to return 0 force
    virtual void operator()(GraphType &g, double t) {
      (void) g; // silence compiler warnings
      (void) t; // silence compiler warnings
      return;
    }
};

/** @class PinConstraint
 * @brief Holds two nodes in place
 *
 * This constraint object will hold nodes at points (0, 0, 0) and
 * (1, 0, 0) fixed using the () operator.
 * 
 */
class PinConstraint: public Constraint {
  public:

    /** Construct an empty force object. */
    PinConstraint() {}

    /** Default destructor */
    virtual ~PinConstraint() {}

    // Method to return 0 force
    virtual void operator()(GraphType &g, double t) {
      (void) t; // silence compiler warnings
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        Point initial_position = (*it).value().initial_position;
        if (initial_position == Point(0, 0, 0) || 
            initial_position == Point(1, 0, 0)) {

          (*it).position() = initial_position;
          (*it).value().vel = Point(0, 0, 0);
        }
      }
      
      return;
    }
};

/** @class PlaneConstraint
 * @brief Stops nodes from moving beyond a plane
 *
 * This constraint object will keep the z dimension of a node's
 * position above -0.75 using the () operator.
 * 
 */
class PlaneConstraint: public Constraint {
  public:

    /** Attributes and Constructor*/
    PlaneConstraint() : z_(-0.75) {}

    /** Default destructor */
    virtual ~PlaneConstraint() {}

    // Method to return 0 force
    virtual void operator()(GraphType &g, double t) {
      (void) t; // silence compiler warnings
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        if ((*it).position().z < this->z_) {

          (*it).position().z = this->z_;
          (*it).value().vel.z = 0;
        }
      }

      return;
    }
  
  private:
    double z_;
};

/** @class SphereConstraint
 * @brief Prevents nodes from entering a sphere
 *
 * This constraint object will push any node entering the sphere
 * with center (0.5, 0.5, -0.5) and radius 0.15 to the sphere's 
 * edge using the () operator.
 * 
 */
class SphereConstraint: public Constraint {
  public:

    /** Attributes and Constructor */
    SphereConstraint() : center_(Point(0.5, 0.5, -0.5)), r_(.15) {}

    /** Default destructor */
    virtual ~SphereConstraint() {}

    // Method to return 0 force
    virtual void operator()(GraphType &g, double t) {
      (void) t; // silence compiler warnings
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        Point position = (*it).position();
        Point velocity = (*it).value().vel;
        Point R = (position - center_) / norm(position - center_);
        if (norm(position - center_) < r_) {

          (*it).position() = r_ * R + center_;
          (*it).value().vel -= dot(velocity, R) * R;
        }
      }

      return;
    }
  
  private:
    Point center_;
    double r_;
};

/** @class DeleteSphereConstraint
 * @brief Deletes nodes that touch a sphere
 *
 * This constraint object will delete any node entering the sphere
 * with center (0.5, 0.5, -0.5) and radius 0.15 using the () operator.
 * 
 */
class DeleteSphereConstraint: public Constraint {
  public:

    /** Attributes and Constructor */
    DeleteSphereConstraint() : center_(Point(0.5, 0.5, -0.5)), r_(.15) {}

    /** Default destructor */
    virtual ~DeleteSphereConstraint() {}

    // Method to return 0 force
    virtual void operator()(GraphType &g, double t) {
      (void) t; // silence compiler warnings
      std::vector<Node> nodes_to_remove;
//--design_1
//--You should delete nodes on the spot through a while loop.
//--Storing then deleting is dangerous due to the reindexing that
//--is not happening in your container
//--START
      for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        Point position = (*it).position();
        if (norm(position - center_) < r_) {
          nodes_to_remove.push_back(*it);
        }
      }
      for (unsigned int i = 0; i < nodes_to_remove.size(); i++) {
        g.remove_node(nodes_to_remove[i]);
      }
//--END
      return;
    }
  
  private:
    Point center_;
    double r_;
};

/** Functor to combine constraints */
struct CombinedConstraints {
  public:
    // Attributes and Constructor
    std::vector<Constraint*> constraints_;

    CombinedConstraints(std::vector<Constraint*> constraints) 
      : constraints_(constraints) {}

    // Methods
    void operator()(GraphType &g, double t) {
      for (unsigned int i = 0; i < constraints_.size(); i++) {
        (*(constraints_[i]))(g, t);
      }

      return;
    }
};

/** Function to combine 2 constraints */
template <typename C1, typename C2>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2) {
  Constraint* constraint1 = &c1;
  Constraint* constraint2 = &c2;
  std::vector<Constraint*> constraint_vector;
  constraint_vector.push_back(constraint1);
  constraint_vector.push_back(constraint2);

  return CombinedConstraints(constraint_vector);
}

/** Function to combine 3 constraints */
template <typename C1, typename C2, typename C3>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2, 
                                              C3 c3) {
  Constraint* constraint1 = &c1;
  Constraint* constraint2 = &c2;
  Constraint* constraint3 = &c3;
  std::vector<Constraint*> constraint_vector;
  constraint_vector.push_back(constraint1);
  constraint_vector.push_back(constraint2);
  constraint_vector.push_back(constraint3);

  return CombinedConstraints(constraint_vector);
}

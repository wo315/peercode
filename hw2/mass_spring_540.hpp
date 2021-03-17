/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <vector>
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

/** Custom data structure to store with Edges */
struct EdgeData {
  double k;
  double l;
  EdgeData(): k(0), l(1) {}
  EdgeData(double k_value, double l_value): k(k_value), l(l_value){}
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
 * @tparam G::node_value_type supports numerical values with velocity and mass
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
};

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g          Graph
 * @param[in]     t          The current time (useful for time-dependent forces)
 * @param[in]     dt         The time step
 * @param[in]     force      Function object defining the force per node
 * @param[in]     constraint Function object defining the constraints on nodes
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports numerical values with velocity and mass
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
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
};


/** Force function object for HW2 #1. */
struct Problem1Force {
  double spring_constant = 100;
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) return Point(0,0,0);
    Point spring_force = Point(0,0,0);
    for (auto itr = n.edge_begin(); itr != n.edge_end(); ++itr) {
      spring_force -= (n.position() - (*itr).node2().position()) * (spring_constant / (norm(n.position() - (*itr).node2().position())+1e-10) * (norm(n.position() - (*itr).node2().position()) - (*itr).value()));
    }
    (void) t; // silence compiler warnings
    return spring_force + Point(0,0,-n.value().mass*grav);
  }
};

/**
 * A Functor base class for forces.
 */
class Force {
  
  public:
    // Default Constructor
    Force() {}
    /**
     * Returns a force of 0.
     * @return @a Point of (0,0,0) representing a zero vector
     */
    virtual Point operator()(Node n, double t) {
      (void) n;
      (void) t; // silence compiler warnings
      return Point(0,0,0);
    };
};

/**
 * A functor implementing the force of gravity.
 */
class GravityForce : public Force {
  public:
    /**
     * Returns the force of gravity (default 9.81) felt by the node's mass. 
     * Assumes that gravity is in the z-direction.
     * @return @a Point of (0,0,-mg) 
     */
    Point operator()(Node n, double t) override {
      (void) t;
      return Point(0,0, -n.value().mass*grav);
    }
};

/**
 * A functor implementing the force felt by a spring.
 */
class MassSpringForce : public Force {
  private:
    double spring_constant = 100; // for now, spring constant is hard set
  public:
    /**
     * Returns the force of felt by the node based on its position with respect to its edges. 
     * Assumes that the edge value parameter is a double with the edge's rest length
     * @return @a Point representing the force vector
     */
    Point operator()(Node n, double t) override {
      (void) t;
      Point spring_force = Point(0,0,0);
      for (auto itr = n.edge_begin(); itr != n.edge_end(); ++itr) {
        spring_force -= (n.position() - (*itr).node2().position()) * (spring_constant / (norm(n.position() - (*itr).node2().position())+1e-10) * (norm(n.position() - (*itr).node2().position()) - (*itr).value()));
      }
      return spring_force;
    }
};

/**
 * A functor implementing a damping force.
 */
class DampingForce : public Force {
  private:
    double damping_constant = 0.1; // hard set damping constant
  public:
    /**
     * Returns the force of felt by the node based on its velocity. 
     * Assumes that the node value contains a velocity field.
     * @return @a Point representing the force vector
     */
    Point operator()(Node n, double t) override{
      (void) t;
      return n.value().vel * -1.0 * damping_constant;
    }
};

/**
 * A functor combining a vector of forces into a net force on a node.
 */
class CombinedForce : public Force {
  private:
    /** A vector representing all of the forces acting on a node*/
    typename std::vector<Force*> forces_;
  public:
    /** Constructor taking advantage of all of the unique forces */
    CombinedForce(std::vector<Force*> forces) {
      forces_ = forces;
    }
    /**
     * Returns the net force of felt by the node based on the sum of its forces. 
     * @return @a Point representing the force vector
     */
    Point operator()(Node n, double t) {
      (void) t;
      Point total_force = Point(0,0,0);
      for (auto itr = forces_.begin(); itr!=forces_.end();++itr) {
        total_force += (*(*itr))(n,t);
      }
      return total_force;
    }
};

/**
 * Method to generate a single combined force out of several force functors
 * @param[in]     force      Function object defining the force per node
 * @return a @a CombinedForce functor representing total force on a node
 * @tparam F is a function object, inheriting from Force, called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 */
template<typename... F>
CombinedForce make_combined_force(F... force) {
  std::vector<Force*> forces= {{&force...}};
  return CombinedForce(forces);
}

/**
 * A functor base class for graph constraints.
 */
class Constraint {
  public: 
    /** Default Constructor */
    Constraint() {}
    /**
     * Functor allowing us to maintain constraints.
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  constraints are enforced on the graph
     */
    virtual void operator()(GraphType& g, double t) {
      (void) g;
      (void) t;
      return;
    }
};

/**
 * A functor base class for graph constraints.
 */
class  PinConstraint : public Constraint {
  private:
    Node origin_node;
    Node one_zero_node;
  public:
    PinConstraint(Node origin, Node one_zero) {
      origin_node = origin;
      one_zero_node = one_zero;
    }
    /**
     * Functor to maintain the fixed endpoints constraint.
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  @a origin_node is held at the origin, @a one_zero_node held in place
     */
    void operator()(GraphType& g, double t) override {
      origin_node.position() = Point(0,0,0);
      origin_node.value().vel = Point(0);
      one_zero_node.position() = Point(1,0,0);
      one_zero_node.value().vel = Point(0);
      (void) g;
      (void) t;
    }
};

/**
 * A constraint functor enforcing nodes to remain above a plane.
 */
class PlaneConstraint : public Constraint {
  public:
    /**
     * Functor to maintain nodes above a certain plane.
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  nodes remain above the plane z=-0.75, and the velocity is 0 for points on the plane.
     */
    void operator()(GraphType& g, double t) override {
      for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
        if(dot((*itr).position(),Point(0,0,1)) < -0.75) {
          (*itr).value().vel = Point(0);
          (*itr).position().z = -0.75;
        }
      }
      (void) t;
    }
};

/**
 * A constraint functor enforcing nodes to remain out of a sphere of some size and position
 */
class SphereConstraint : public Constraint {
  private:
    /** Location of the sphere */
    Point sphere_center;
    /** Radius of the sphere */
    double radius;
  public:
    /** Constructor for spherical constraint
      * @param[in] center   @a Point representing the center of the sphere
      * @param[in] r        radius of the sphere
      */
    SphereConstraint(Point center, double r) {
      sphere_center = center;
      radius =r;
    }
    
    /**
     * Functor to maintain nodes outside of some sphere.
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  nodes remain out of defined sphere, and the velocity is chosen as to bounce away from sphere
     */
    void operator()(GraphType& g, double t) override {
      for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
        if(norm((*itr).position()-sphere_center) < radius) {
          (*itr).value().vel -= ((*itr).position()-sphere_center) / norm((*itr).position()-sphere_center) * dot((*itr).value().vel, ((*itr).position()-sphere_center) / norm((*itr).position()-sphere_center));
          (*itr).position() = (((*itr).position()-sphere_center) * radius / norm((*itr).position()-sphere_center)) + sphere_center;
        }
      }
      (void) t;
    }
};

/**
 * A constraint functor class that allows multiple constraints to be enforced.
 */
class CombinedConstraint : public Constraint {
  private:
    /** Container of constraints */
    std::vector<Constraint*> constraints;
  public:
    /** Constructor taking advantage of all of the unique @a Constraint functors */
    CombinedConstraint(std::vector<Constraint*> constr) {
      constraints = constr;
    }
    /**
     * Functor to enforce multiple constraints on some graph
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  nodes follow given constraints
     */
    void operator()(GraphType& g, double t) override {
      for (auto itr = constraints.begin(); itr!=constraints.end();++itr) {
        (*(*itr))(g,t);
      }
    }
};

/**
 * A constraint functor that removes nodes upon entry of some sphere region.
 */
class SphericalScissorConstraint : public Constraint {
  private:
    /** Location of the sphere */
    Point sphere_center;
    /** Radius of the sphere */
    double radius;
  public:
    /** Constructor for spherical constraint
      * @param[in] center   @a Point representing the center of the sphere
      * @param[in] r        radius of the sphere
      */
    SphericalScissorConstraint(Point center, double r) {
      sphere_center = center;
      radius =r;
    }
    /**
     * Functor to maintain nodes outside of some sphere.
     * @param[in] g   @a GraphType reference to apply constraints to
     * @param[in] t   time reference for enforcing constraints
     * @post  nodes and associated edges are removed whenever they enter the sphere
     */
    void operator()(GraphType& g, double t) override {
      for (auto itr = g.node_begin(); itr != g.node_end(); ++itr) {
        if(norm((*itr).position()-sphere_center) < radius) {
          g.remove_node(*itr);
        }
      }
      (void) t;
    }
};
/**
 * Method to generate a single combined force out of several constraint functors
 * @param[in]     constr      Function object defining the constraint to enforce
 * @return a @a CombinedConstraint functor enforcing all given constraints
 * @tparam C is a function object, inheriting from Constraint, called as @a constraint(g, @a t),
 *           where g is a reference to the graph and @a t is the current time.
 */
template<typename... C>
CombinedConstraint make_combined_constraint(C... constr) {
  std::vector<Constraint*> constraints= {{&constr...}};
  return CombinedConstraint(constraints);
}

//--functionality_0
//--Passed all tests!
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


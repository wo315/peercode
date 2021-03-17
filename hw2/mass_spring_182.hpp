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

struct EdgeData {
  double K; //< Edge spring constant
  double Lij;     //< Edge length
  EdgeData() : K(0), Lij(0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
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

template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // std::cout << "calling symp euler with constraint" << std::endl;
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  
  // Apply the constraints
  constraint(g, dt);

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
    (void) t;    // silence compiler warnings
    if (n.position () == Point (0 ,0 ,0) || n.position () == Point (1 ,0 ,0)) {
      return Point (0 ,0 ,0);
    }
    Point force = Point(0,0,0);
    for (auto ni = n.edge_begin(); ni != n.edge_end(); ++ni) {
      Point diff =  n.position() - (*ni).node2().position();
      const double K = (*ni).value().K;
      double Lij = (*ni).value().Lij;
      force += - K * diff / norm(diff) * (norm(diff) - Lij);
    }
    force += Point(0,0,-1*grav) * n.value().mass;
    return force;
  }
};


// GENERALIZED FORCES

// template <typename NODE>
struct Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void) n, (void) t;
      return Point(0,0,0);
    }
    Force() {}
    virtual ~Force() {}
};

class GravityForce : public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void) t;
      return Point(0,0, -1*grav) * n.value().mass;
    }
    GravityForce() {}
    virtual ~GravityForce() {}
};

// using force_type_ref = Force&;

class MassSpringForce : public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void) t;
      if (n.position () == Point (0 ,0 ,0) || n.position () == Point (1 ,0 ,0)) {
        return Point (0 ,0 ,0);
      }
      Point spring_force = Point(0,0,0);
      for (auto ni = n.edge_begin(); ni != n.edge_end(); ++ni) {
        Point diff =  n.position() - (*ni).node2().position();
        const double K = (*ni).value().K;
        double Lij = (*ni).value().Lij;
        spring_force += - K * diff / norm(diff) * (norm(diff) - Lij);
      }
      return spring_force;
    }
    MassSpringForce() {}
    virtual ~MassSpringForce() {}
};

class DampingForce : public Force {
  public:
    virtual Point operator()(Node n, double t) {
      (void) t;
      double c = 1; // damping constant
      return -1*c * n.value().vel;
    }
    DampingForce() {}
    virtual ~DampingForce() {}
};

struct CombinedForce {
  CombinedForce(std::vector<Force*> forces) : forces_(forces) {
  }
  Point operator() (Node n, double t) {
    Point p = Point(0,0,0);
    for (auto it = forces_.begin(); it != forces_.end(); ++it) {
      p += (*(*it))(n, t);
    }
    return p;
  }
  private:
    std::vector<Force*> forces_;
};

CombinedForce make_combined_force(Force&& f1, Force&& f2) {
  std::vector<Force*> forces;
  forces.push_back(&f1);
  forces.push_back(&f2);
  return CombinedForce(forces);
}

// GENERALIZED CONSTRAINTS

class Constraint {
  public:

  virtual void operator() (GraphType& g, double dt) {
    (void) g; (void) dt;
  }
  Constraint() {}
  virtual ~Constraint() {}
};

class PinConstraint : public Constraint {
  public:

  PinConstraint() {}
  virtual ~PinConstraint() {}
  virtual void operator() (GraphType& g, double dt) {
    // Keep nodes at (0,0,0) and (1,0,0) fixed
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      Point prev_pos = (*ni).position() - ((*ni).value().vel * dt);
      if (prev_pos == Point(0,0,0)) {
        (*ni).position() = Point(0,0,0);
      }
      else if (prev_pos == Point(1,0,0)) {
        (*ni).position() = Point(1,0,0);
      }
    }
  }
};

class PlaneConstraint : public Constraint {
  public:

  PlaneConstraint() {}
  virtual ~PlaneConstraint() {}

  virtual void operator() (GraphType& g, double dt) {
    (void) dt;
    double z = -0.75;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      
      // Check if the point violates the constraint
      if ( (*ni).position()[2] < z ) {
        
        // Set position to nearest point on the plane
        Point pt = (*ni).position();
        (*ni).position() = Point(pt[0], pt[1], z);

        // Set z-component of the Node velocity to 0
        Point old_vel = (*ni).value().vel;
        (*ni).value().vel = Point(old_vel[0], old_vel[1], 0);
      }
    }
  }
};

class SphereConstraint : public Constraint {
  public:

  SphereConstraint() {}
  virtual ~SphereConstraint() {}
  virtual void operator() (GraphType& g, double dt) {
    (void) dt;
    Point c = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      
      // Check if the point violates the constraint (i.e. is inside the sphere)
      if ( norm((*ni).position() - c) < radius) {
        
        // Set position to nearest point on sphere by:
        //    1. Construct a unit vector from the point to the centers
        Point Ri = ((*ni).position() - c) / norm((*ni).position() - c);

        //    2. Scale the unit vector by the radius
        (*ni).position() = c + Ri * radius;

        // Set component of velocity normal to sphere's surface = 0
        Point old_vel = (*ni).value().vel;
        (*ni).value().vel = old_vel - (Ri * dot(old_vel, Ri));
      }
    }
  }
};

class SphereConstraint2 : public Constraint {
  public:

  SphereConstraint2() {}
  virtual ~SphereConstraint2() {}
  virtual void operator() (GraphType& g, double) {
    Point c = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    for (auto ni = g.node_begin(); ni != g.node_end(); ++ni) {
      
      // Check if the point violates the constraint (i.e. is inside the sphere)
      if ( norm((*ni).position() - c) < radius) {  
        // Remove the node
        Node n = (*ni);
        g.remove_node(n);
        std::cout << "removed node" << std::endl;
      }
    }
  }
};

struct CombinedConstraints {
  // Constructor
  CombinedConstraints(std::vector<Constraint*> constraints) : constraints_(constraints) {
  }
  void operator() (GraphType& g, double dt) {
    // For each constraint in the vector, apply it to the Graph g
    for (auto it = constraints_.begin(); it != constraints_.end(); ++it) {
      (*(*it))(g, dt);
    }
  }
  private:
    // A vector of Constraint objects
    std::vector<Constraint*> constraints_;
};

CombinedConstraints make_combined_constraint(Constraint&& c1, Constraint&& c2) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  return CombinedConstraints(constraints);
}

CombinedConstraints make_combined_constraint(Constraint& c1, Constraint& c2, Constraint& c3) {
  std::vector<Constraint*> constraints;
  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);
  return CombinedConstraints(constraints);
}
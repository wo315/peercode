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

/** Custom structure of data to store with Nodes */
struct EdgeData {
  double l0;     //< edge rest length
  double K;      //< edge spring constant
  EdgeData() : l0(0), K(0) {}
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

    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0) )
      continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0) )
      continue;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force and graph constraints.
 * @param[in,out] g            Graph
 * @param[in]     t            The current time (useful for time-dependent forces)
 * @param[in]     dt           The time step
 * @param[in]     force        Function object defining the force per node
 * @param[in]     constraint   Constraint object defining the graph's constraints
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(g, @a t),
 *           where g is the graph and @a t is the current time.
 *           @a constraint has no return value.
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0) )
      continue;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraint
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

//--style_1
//--This is redundant now that we have a pin contraint
//--START
    if ( n.position() == Point(0,0,0) || n.position() == Point(1,0,0) )
      continue;
//--END
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
    (void) t; // silence compiler warnings
    // zero force on points with the position
    if (n. position () == Point (0 ,0 ,0) || n. position () == Point (1 ,0 ,0))
      return Point(0 ,0 ,0);
    double l0 = 0.001;
    double K = 100.;
    Point f {0,0,0};

    // spring
    for (auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it) {
      Edge e_i = (*e_it);
      NODE n_j = e_i.node2(); // adjacent node
      Point sgn_dist = n.position() - n_j.position();
      double dist = norm(sgn_dist);
      f -= K * sgn_dist/dist*(dist-l0);
    }

    // gravity
    f += (n.value().mass * Point(0,0,-grav));

    return f;
  }
};

struct Force {
  /** a general Force structor for 3 dimensional force acting on the input node
  * @param[in] n node which force acts on
  * @param[in] t time
  * @returns a zero 3 dimensional force
  */
  virtual Point operator()(Node n, double t) {
    (void) n; (void) t;
    return Point {0,0,0};
  }
  virtual ~Force() {}
};

struct GravityForce : public Force {
  /** a Gravity Force structor for 3 dimensional gravity acting on input node
  * @param[in] n node which force acts on
  * @param[in] t time
  * @returns gravity force in z direcation
  */
  virtual Point operator()(Node n, double t) {
    (void) t;
    double f_z = -n.value().mass*grav;
    return {0,0,f_z};
  }
};


struct MassSpringForce : public Force {
  /** a Spring Force structor for 3 dimensional spring force acting on input node
  * depending on distance of adjecent nodes and stiffness of edge between them.
  * @param[in] n node which force acts on
  * @param[in] t time
  * @returns total spring force acting on node
  */
  virtual Point operator()(Node n, double t) {
    (void) t;
    Point f {0,0,0};

    for (auto e_it = n.edge_begin(); e_it != n.edge_end(); ++e_it) {
      Edge e_i = (*e_it);
      Node n_j = e_i.node2(); // adjacent node
      Point sgn_dist = n.position() - n_j.position();
      double dist = norm(sgn_dist);
      if (dist == 0.0)
        throw std::runtime_error("0 division: Non unique positions of nodes.");
      f -= e_i.value().K * sgn_dist/dist*(dist-e_i.value().l0);
    }
  
    return f;
  }
};

/** a Damping Force structor for 3 dimensional damping acting on input node
 * @param[in] n node which force acts on
 * @param[in] t time
 * @returns total damping force acting on node
 */
struct DampingForce : public Force {
  
  double c_; // damping constant
  
  DampingForce() {};
  DampingForce(double c) : c_(c) {}

  virtual Point operator()(Node n, double t) {
    (void) t;
    return -c_*n.value().vel;
  }
  
};


struct CombinedForce : public Force {
  std::vector<Force*> forces_;
  
  CombinedForce(std::vector<Force*> forces) :
   forces_(forces) {}

  ~CombinedForce() {
    for (Force* f : forces_)
      delete f;
    forces_.clear();
  }

  /** a Combined Force structor for 3 dimensional forces acting on input node
   * @param[in] n node which force acts on
   * @param[in] t time
   * @returns total force acting on node
   */
  Point operator()(Node n, double t) {
    Point f_tot {0,0,0};
    for (Force* f : forces_)
      f_tot += (*f)(n, t);
    return f_tot;
  }
};
//--design_-1
//--Good job handling combining rvalue forces and constraints
//--END
template <typename F1, typename F2, typename F3 = Force>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3 = Force()) {
  std::vector<Force*> forces;
  forces.push_back(new F1(f1));
  forces.push_back(new F2(f2));
  forces.push_back(new F3(f3));
  return CombinedForce(forces);
}

struct Constraint {
  /** a general constraint
  * @param[in] g graph s.t. to the constraint
  * @param[in] t time
  * @returns nothing
  */
  virtual void operator()(GraphType& g, double t) {
    (void) g; (void) t;
    return;
  }

  virtual ~Constraint() {}
};

struct PinConstraint : public Constraint {

  std::vector<Point> pin_pos_ {{0,0,0},{1,0,0}};

  /** a pin constraint 
  * @param[in] g graph s.t. to pin constraint
  * @param[in] t time
  * @returns nothing
  * @post nodes with pinned positions (pin_pos_) have 0 velocity.
  * (and are therefore pinned there)
  */
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    for (auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it) {
      Node n = *n_it;
      for (Point pin : pin_pos_) {
        if (n.position() == pin)
          n.value().vel = Point(0,0,0);
      }
    }
    return;
  }

};

struct PlaneConstraint : public Constraint {

  double z_ = -0.75;

  /** a plane constraint
  * @param[in] g graph s.t. to pin constraint
  * @param[in] t time
  * @returns nothing
  * @post no nodes in g are past the z boundary 
  */
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    for (auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it) {
      Node n = *n_it;
      if (n.position().z < z_) {
        // constraint violated, set to nearest point on plane and z-vel to zero
        n.position().z = z_;
        n.value().vel.z = 0.;
      }
    }
    return;
  }

};

struct SphereConstraint : public Constraint {

  Point c_ = {0.5, 0.5, -0.5};
  double r_ = 0.15;

  /** a sphere constraint
  * @param[in] g graph s.t. to pin constraint
  * @param[in] t time
  * @returns nothing
  * @post no nodes inside sphere and velocity of boundary nodes is reduced
  */
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    for (auto n_it = g.node_begin(); n_it != g.node_end(); ++n_it) {
      Node n = *n_it;
      double r_i = norm(n.position() - c_);
      if (r_i < r_) {
        // constraint violated, set to nearest point on sphere + reduce velocity
        Point R_i = (n.position() - c_)/r_i;
        n.position() = c_ + r_*R_i;
        n.value().vel -= R_i*dot(n.value().vel, R_i);
      }
    }
    return;
  }

};

struct TearConstraint : public Constraint {

  Point c_ = {0.5, 0.5, -0.5};
  double r_ = 0.15;

  /** a tear constraint
  * @param[in] g graph s.t. to pin constraint
  * @param[in] t time
  * @returns nothing
  * @post no nodes inside sphere
  * @post nodes found in sphere removed along with their edges (teared)
  * @post new graph num_nodes <= old graph num_nodes
  * @post new graph num_edges <= old graph num_edges
  */
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    for (auto n_it = g.node_begin(); n_it != g.node_end();) {
      Node n = *n_it;
      double r_i = norm(n.position() - c_);
      // if constraint violated, remove node. then other node has same n_it
      if (r_i < r_)
        g.remove_node(n_it);
      else
        ++n_it;
    }
    return;
  }

};

struct CombinedConstraints : public Constraint {
  std::vector<Constraint*> constr_;

  CombinedConstraints(std::vector<Constraint*> constraints) :
    constr_(constraints) {}

  ~CombinedConstraints() {
    for (Constraint* c : constr_)
      delete c;
    constr_.clear();
  }

  /** a Combined Constraint
  * @param[in] g graph s.t. to pin constraint
  * @param[in] t time
  * @returns nothing
  * @post depends on member constraints
   */
  void operator()(GraphType& g, double t) {
    // apply constraints
    for (Constraint* c : constr_)
      (*c)(g, t);
    return;
  }
};

template <typename C1, typename C2, typename C3 = Constraint>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2, C3 c3 = Constraint()) {
  std::vector<Constraint*> constr;
  constr.push_back(new C1(c1));
  constr.push_back(new C2(c2));
  constr.push_back(new C3(c3));
  return CombinedConstraints(constr);
}

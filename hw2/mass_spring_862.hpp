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
  Point x0;        //< Node initial position 
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData
{
  double K; // Edge spring constant
  double L; // Edge rest length 
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
 * @tparam G::node_value_type contains the node's position, velocity, and mass
 * @tparam G::edge_value_type contains the edge's rest length and spring constant 
 * @tparam F is a function object called as @a force(n, @a t), 
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(n, @a t), 
 *           where n is a node of the graph and @a t is the current time.
 *           @a c modifies the graph directly and has no return value. 
 */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {

  // Store the initial positions of nodes 
  if (t==0){
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      n.value().x0 = n.position(); 
    }
  }

  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply constraints in succession 
  constraint(g, t);

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
    (void) t; // silence compiler warning

    // Compute spring force by iterating over the edges this node belongs to 
    Point fspring = Point(0,0,0);
    for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei)
    {
      Edge e = *ei;

      // node_disp = x_i - x_j where x_i is position of the current node @a n 
      Point node_disp = e.node1().position() - e.node2().position();
      double edge_length = norm(node_disp); 

      // Hooke's law applied to node @a n
      double L = e.value().L; // rest length Lij for the edge connecting nodes i and j 
      double K = e.value().K; // spring constant Kij for the edge connecting nodes i and j 
      fspring -= K * node_disp / edge_length * (edge_length - L);
    }

    // To prevent the cloth from falling to infinity, we constrain two of its corners
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)){
      return Point(0, 0, 0);
    }else{
      // force = spring force + force due to gravity
      return fspring - n.value().mass * Point(0, 0, grav);
      // return fspring;
    }
  }
};

/** Parent force class for making specific force function objects */
class Force{

  public:

    /** Return zero force for any provided node @a n or time @a t 
     *  (since the derived force classes should be called instead). 
     *  This is achieved by specifying the method as virtual. */
    virtual Point operator()(Node n, double t){
      (void)t; // silence compiler warning
      (void)n; // silence compiler warning

      return Point(0,0,0); 
    };
};

/** Force function object for modeling the force due to gravity. */
class GravityForce : public Force{

  public:
    /** Return the force due to gravity for node @a n at time @a t. */
    Point operator()(Node n, double t){
      (void)t; // silence compiler warning

      return -n.value().mass * Point(0, 0, grav);
    };
};

/** Force function object for modeling a simple Hookean spring force. */
class MassSpringForce : public Force{

  public:
    /** Return the spring force acting on node @a n at time @a t. */
    Point operator()(Node n, double t)
    {
      (void)t; // silence compiler warning

      // Compute spring force by iterating over the edges this node belongs to
      Point fspring = Point(0, 0, 0);
      for (auto ei = n.edge_begin(); ei != n.edge_end(); ++ei)
      {
        Edge e = *ei;

        // node_disp = x_i - x_j where x_i is position of the current node @a n
        Point node_disp = e.node1().position() - e.node2().position();
        double edge_length = norm(node_disp);

        // Hooke's law applied to node @a n
        double L = e.value().L; // rest length Lij for the edge connecting nodes i and j
        double K = e.value().K; // spring constant Kij for the edge connecting nodes i and j
        fspring -= K * node_disp / edge_length * (edge_length - L);
      }

      return fspring;
    };
};

/** Force function object for modeling a damping force. */
class DampingForce : public Force
{
  public:
    /** Return the damping force acting on node @a n at time @a t. */
    Point operator()(Node n, double t)
    {
      (void)t; // silence compiler warning

      return -c_*n.value().vel;
    };

    // Simple constructor that initializes the damping constant c 
    DampingForce(double c = 0.0) : c_(c) {}

  private:
    double c_; // damping constant c 
};

/** Force function object for combining several forces. */
struct CombinedForce
{
  private: 
    std::vector<Force*> forces_;

  public:
    /* capture function - used to capture forces we wish to combine */
    CombinedForce(const std::vector<Force*> &forces) : forces_(forces) {}

    /* Access operator for our functor */
    Point operator()(Node n, double t)
    {
      Point total_force = Point(0,0,0); 

      for (auto it = forces_.begin(); it != forces_.end(); ++it){
        Force& f = *(*it);
        total_force += f(n,t);
      }
      return total_force; 
    }
};

/** Template function that allows one to combine
 *  two force objects into a single combined 
 *  force (i.e. object of type @a CombinedForce ) */
template <typename force1, typename force2>
CombinedForce make_combined_force(force1 f1, force2 f2)
{
  std::vector<Force*> forces;
  
  forces.push_back(&f1);
  forces.push_back(&f2); 

  return CombinedForce(forces);
}

/** Template function that allows one to combine
 *  three force objects into a single combined 
 *  force (i.e. object of type @a CombinedForce ) */
template <typename force1, typename force2, typename force3>
CombinedForce make_combined_force(force1 f1, force2 f2, force3 f3)
{
  std::vector<Force *> forces;

  forces.push_back(&f1);
  forces.push_back(&f2);
  forces.push_back(&f3); 

  return CombinedForce(forces);
}

/** Parent constraint class for making specific constraint function objects */
class Constraint
{
public:
  /** Apply no constraint for any provided node @a n or time @a t 
   *  (since the derived constraint classes should always be called instead). 
   *  This is achieved by specifying the method as virtual. */
  virtual void operator()(GraphType &g, double t)
  {
    (void)g; // silence compiler warning
    (void)t; // silence compiler warning
  };
};

/** Constraint function object for constraining two nodes of the 
 *  mass-spring system (specifically the points at locations
 *  (0,0,0) and (1,0,0)). */
class PinConstraint : public Constraint
{
public:
  /** Fix the position of nodes located at (0,0,0) and (1,0,0)
   *  so that their position is "pinned" down throughout
   *  the simulation. */
  void operator()(GraphType &g, double t)
  {
    (void)t; // silence compiler warning

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;

      // Pin node at (0,0,0)
      if (n.value().x0 ==Point(0,0,0)){
        n.position() = Point(0,0,0);
      }

      // Pin node at (1,0,0)
      if (n.value().x0 ==Point(1,0,0)) n.position() = Point(1,0,0);
    }
  };
};

/** Constraint function object for preventing nodes of the
 *  mass-spring system from passing a plane located at a 
 *  constant z-value (specifically z = -0.75). */
class PlaneConstraint : public Constraint
{
public:
  /** For nodes passing through the plane located at z = -0.75,
   *  set their position to the nearest point on this plane (i.e. 
   *  set their z-position to -0.75) and remove their z-velocity. 
   *  These two modifications prevent nodes from passing through
   *  this plane. */
  void operator()(GraphType &g, double t)
  {
    (void)t; // silence compiler warning

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;

      // check if z coordinate is less than -0.75 
      if (n.position()[2] < -0.75){
        // set position to the nearest point on the plane z = -0.75. 
        // i.e. the point with the same x,y coordiantes but z = -0.75. 
        n.position()[2] = -0.75; 

        // set the z-component of the node velocity to zero 
        n.value().vel[2] = 0.0; 
      }
    }
  };
};

/** Constraint function object for preventing nodes of the
 *  mass-spring system from passing through a sphere
 *  of radius 0.15 located at (0.5, 0.5, -0.5). */
class SphereConstraint : public Constraint
{
public:
  /** For nodes passing through the sphere located at (0.5,0.5,-0.5) 
   *  having a radius 0.15, set their position to the nearest point 
   *  on this sphere and set their velocity normal to the sphere
   *  equal to 0. These two modifications prevent nodes from passing through
   *  this sphere. */
  void operator()(GraphType &g, double t)
  {
    (void)t; // silence compiler warning
    Point c = Point(0.5, 0.5, -0.5); // center of the sphere
    double r = 0.15;  // radius of the sphere 

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;

      // check if node falls within this sphere's radius 
      Point disp = n.position()-c; // displacement of the node from the center
      if (norm(disp) < r)
      {
        // set position to the nearest point on the sphere's surface.
        // this point is given by c + r/|p-c|*(p-c) where (p-c)
        // is the displacement of the point from the center @ c
        n.position() = c + r*disp/norm(disp);

        // set the component of velocity normal to the sphere's surface equal to 0
        Point unit_normal = disp/norm(disp); // sphere unit normal vector 
        n.value().vel = n.value().vel - dot(n.value().vel, unit_normal)*unit_normal; 
      }
    }
  };
};

/** Constraint function object for removing nodes of the
 *  mass-spring system that pass through a sphere
 *  of radius 0.15 located at (0.5, 0.5, -0.5). This in effect
 *  creates a "tear" in the mass-spring system for the application
 *  of modeling sheet dynamics. */
class TearConstraint : public Constraint
{
public:
  /** Remove nodes passing through the sphere located at (0.5,0.5,-0.5) 
   *  that has a radius of 0.15, effectively creating a "hole" or "tear"
   *  in the mass-spring system as viewed as a sheet. */
  void operator()(GraphType &g, double t)
  {
    (void)t;                         // silence compiler warning
    Point c = Point(0.5, 0.5, -0.5); // center of the sphere
    double r = 0.15;                 // radius of the sphere

    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;

      // check if node falls within this sphere's radius
      Point disp = n.position() - c; // displacement of the node from the center
      if (norm(disp) < r)
      {
        // delete this node
        g.remove_node(n); 
      }
//--design_1
//--Don't increment the iterator if you delete a node
//--END
    }
  };
};

/** Function object for combining several constraints. */
struct CombinedConstraints
{
  private:
    std::vector<Constraint*> constraints_;

  public:
    /* capture function - used to capture the constraints we'd like to combine */
    CombinedConstraints(const std::vector<Constraint*> &constraints) : constraints_(constraints) {}

    /* Access operator for our functor */
    void operator()(GraphType &g, double t)
    {
      for (auto it = constraints_.begin(); it != constraints_.end(); ++it)
      {
        Constraint &c = *(*it);
        c(g,t);
      }
    }
};

/** Template function that allows one to combine
 *  two constraint objects into a single combined 
 *  constraint (i.e. object of type @a CombinedConstraint ) */
template <typename constraint1, typename constraint2>
CombinedConstraints make_combined_constraints(constraint1 c1, constraint2 c2)
{
  std::vector<Constraint*> constraints;

  constraints.push_back(&c1);
  constraints.push_back(&c2);

  return CombinedConstraints(constraints);
}

/** Template function that allows one to combine
 *  three constraint objects into a single combined 
 *  constraint (i.e. object of type @a CombinedConstraint ) */
template <typename constraint1, typename constraint2, typename constraint3>
CombinedConstraints make_combined_constraints(constraint1 c1, constraint2 c2, constraint3 c3)
{
  std::vector<Constraint *> constraints;

  constraints.push_back(&c1);
  constraints.push_back(&c2);
  constraints.push_back(&c3);

  return CombinedConstraints(constraints);
}

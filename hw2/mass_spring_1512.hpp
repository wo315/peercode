/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

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
struct NodeData
{
  Point vel;   //< Node velocity
  double mass; //< Node mass
  bool is_pinned; //< Whether the point should be pinned
  NodeData() : vel(Point(0, 0, 0)), mass(1), is_pinned(false){}
};

/** Custom structure of data to store with Edges */
struct EdgeData
{
  double K;
  double L;
  EdgeData() : K(100.0), L(1.0) {}
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
template <typename G, typename F, typename C>
double symp_euler_step(G &g, double t, double dt, F force, C constraint)
{
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Applying the contraints
  constraint(g, t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;
    
    Point fc;
    if (n.value().is_pinned) {
      fc = Point(0,0,0);
    }else{
      fc = force(n, t);}

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += fc * (dt / n.value().mass);
  }
  return t + dt;
}

/** Force function object for HW2 #1. */
struct Problem1Force
{
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */

  template <typename NODE>
  Point operator()(NODE n, double t)
  {
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    {
      return Point(0, 0, 0);
    }

    Point force = Point(0, 0, 0);

    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      Edge e = *it;
      Point diff = n.position() - e.node2().position();
      Point spring = -e.value().K * (diff) / norm(diff) * (norm(diff) - e.value().L);
      force += spring;
    }

    Point gravity = n.value().mass * Point(0, 0, -grav);
    force += gravity;

    (void)t;

    return force;
  }
};

/** @class Force
 * @brief Parent class for different types of forces
 * All forces can be added.
 */
class Force
{
public:
  Force()
  {
  }

  virtual ~Force()
  {
  }

  virtual Point operator()(Node n, double t)
  {
    (void)t; (void)n;
    return Point(0, 0, 0);
  };
};

/** @class GravityForce
 * @brief Derived Force class to compute gravity given a node and time
 */
class GravityForce : public Force
{
public:
  GravityForce()
  {
  }

  virtual ~GravityForce()
  {
  }

  virtual Point operator()(Node n, double t)
  {
    (void)t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** @class MassSpringForce
 * @brief Derived Force class to compute spring force given a node and time
 */
class MassSpringForce : public Force
{
public:
  MassSpringForce() {}

  virtual ~MassSpringForce() {}

  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point force=Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      Edge e = *it;
      Point diff = n.position() - e.node2().position();
      Point spring = -e.value().K * (diff) / norm(diff) * (norm(diff) - e.value().L);
      force += spring;
    }
    return force;
  }
};

/** @class DampingForce
 * @brief A force that implements damping force
 */
class DampingForce : public Force
{
public:
  DampingForce()
  {
  }

  DampingForce(double c): c_(c){}

  virtual ~DampingForce()
  {
  }

  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point force=Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it)
    {
      Edge e = *it;
      force += -c_ * e.node2().value().vel;
    }
    return force;
  }
  
  private:
    double c_ = 0.0;
};

/** Functor to combine different forces with a vector of Force pointers. */
struct CombinedForce
{
  std::vector<Force *> fc_;
  CombinedForce(std::vector<Force *> &fc) : fc_(fc) {}
  Point operator()(Node n, double t)
  {
    (void)t;
    Point result = Point(0, 0, 0);
    for (auto force_it = fc_.begin(); force_it != fc_.end(); force_it++)
    {
      result += (*force_it)->operator()(n, t);
    }
    return result;
  }
};

/** Function to make a vector of forces and call CombinedForce functor */
template <typename fc1, typename fc2>
CombinedForce make_combined_force(fc1 f1, fc2 f2)
{
  std::vector<Force *> fc;
  fc.push_back(&f1);
  fc.push_back(&f2);
  return CombinedForce(fc);
}


template <typename fc1, typename fc2, typename fc3>
CombinedForce make_combined_force(fc1 f1, fc2 f2, fc3 f3)
{
  std::vector<Force *> fc;
  fc.push_back(&f1);
  fc.push_back(&f2);
  fc.push_back(&f3);
  return CombinedForce(fc);
}

/** @class Constraints
 * @brief Parent class for different types of constraints
 * Different constraints can be enforced together
 */
class Constraint
{
public:
  Constraint()
  {
  }

  virtual ~Constraint()
  {
  }

  virtual void operator()(GraphType &g, double t)
  {
    (void)g;
    (void)t;
  }
};

/** @class PinConstraints
 * @brief Constraint that keep the nodes at (0,0,0) and (1,0,0) fixed
 */
class PinConstraint : public Constraint
{
public:
  PinConstraint()
  {
  }

  virtual ~PinConstraint()
  {
  }

  virtual void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = (*it);
      if (n.position() == p1 || n.position() == p2)
      {
        n.value().is_pinned = true;
      }
    }
  }

private:
  Point p1 = Point(0, 0, 0);
  Point p2 = Point(1, 0, 0);
};

/** @class PlaneConstraints
 * @brief Constraint that keeps the nodes above a plane. 
 */
class PlaneConstraint : public Constraint
{
public:
  PlaneConstraint()
  {
  }

  virtual ~PlaneConstraint()
  {
  }

  virtual void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = (*it);
      if (n.position().z < Z)
      {
        n.position().z = Z;
        n.value().vel.z=0;
      }
    }
  }

private:
  double Z = -0.75;
};

/** @class SphereConstraints
 * @brief Constraint that keeps the nodes within a sphere. 
 */
class SphereConstraint : public Constraint
{
public:
  SphereConstraint()
  {
  }

  virtual ~SphereConstraint()
  {
  }

  virtual void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = (*it);
      double distance = norm(n.position()-c);
      if (distance<r)
      {
        Point np = (n.position()-c)/distance*r+c;
        n.position() = np;
        Point Ri=(n.position()-c)/distance;
        Point vi=n.value().vel;
        n.value().vel =vi-dot(vi,Ri)*Ri;
      }
    }
  }

private:
  Point c = Point(0.5,0.5,-0.5);
  double r = 0.15;
};

/** @class TearConstraints
 * @brief Constraint that removes nodes that touchs a sphere. 
 */
class TearConstraint : public Constraint
{
public:
  TearConstraint()
  {
  }

  virtual ~TearConstraint()
  {
  }

  virtual void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      Node n = (*it);
      double distance = norm(n.position()-c);
      if (distance<r)
      { 
        g.remove_node(n);
      }
    }
  }

private:
  Point c = Point(0.5,0.5,-0.5);
  double r = 0.15;
};

/** Functor to combine different constraints with a vector of Constraint pointers */
struct CombinedConstraints
{
  std::vector<Constraint *> cst_;
  CombinedConstraints(std::vector<Constraint *> &cst):cst_(cst) {}
  void operator()(GraphType &g, double t)
  {
    (void)t;
    for (auto cons_it = cst_.begin(); cons_it != cst_.end(); cons_it++)
    {
      (*cons_it)->operator()(g, t);
    }
  }
};

/** Function to combine constraints into a vector and call CombinedConstraints.*/
template <typename cons1, typename cons2>
CombinedConstraints make_combined_constraint(cons1 c1, cons2 c2)
{
  std::vector<Constraint *> cons;
  cons.push_back(&c1);
  cons.push_back(&c2);
  return CombinedConstraints(cons);
}

/** Function to combine constraints into a vector and call CombinedConstraints.*/
template <typename cons1, typename cons2, typename cons3>
CombinedConstraints make_combined_constraint(cons1 c1, cons2 c2, cons3 c3)
{
  std::vector<Constraint *> cons;
  cons.push_back(&c1);
  cons.push_back(&c2);
  cons.push_back(&c3);
  return CombinedConstraints(cons);
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


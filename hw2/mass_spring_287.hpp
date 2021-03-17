/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <typeinfo>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** @struct NodeData
 * @brief Custom structure of data to store with Nodes. 
 */
struct NodeData
{
  Point vel;   //< Node velocity
  double mass; //< Node mass
  double t;
  Point init_pos;
  NodeData() : vel(0), mass(1), t(0) {}
};

/** @struct EdgeData
 * @brief Custom structure of data to store with Edges. 
 */
struct EdgeData
{
  double rest_length_; //< Edge rest_length
  double K_;           //< Edge K_
  EdgeData() : rest_length_(0), K_(1) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/** @struct Force.
 * @brief Functor to apply a force to the nodes of the graph.
 */
struct Force
{
  /**
   * @brief Applies force to a node.
   * @param[in] @a n, a node.
   * @param[in] @t, the time instant.
   * 
   * @return Force to the node, in a Point object.
   */
  virtual Point operator()(Node n, double t)
  {
    (void)n;
    (void)t;
    return Point(0, 0, 0);
  }
};

/** @struct Constraint.
 * @brief Functor to apply a constraint to the graph.
 */
struct Constraint
{
  /**
   * @brief Applies constraint to a graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  virtual void operator()(GraphType &g, double t)
  {
    (void)g;
    (void)t;
    return;
  };
};

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint for the graph
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
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  constraint(g, t);
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it)
  {
    auto n = *it;
    //if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    //{
    //  continue;
    //}

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
double symp_euler_step(G &g, double t, double dt, F force)
{
  Constraint c = Constraint();
  return symp_euler_step(g, t, dt, force, c);
}

/** @struct GravityForce.
 * @brief Functor to apply a gravity force to the nodes of the graph.
 */
struct GravityForce : Force
{
  /**
   * @brief Applies gravity force to a node.
   * @param[in] @a n, a node.
   * @param[in] @t, the time instant.
   * 
   * @return Force to the node, in a Point object.
   */
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** @struct MassSpringForce.
 * @brief Functor to apply a mass spring force to the nodes of the graph.
 */
struct MassSpringForce : Force
{
  /**
   * @brief Applies mass spring force to a node.
   * @param[in] @a n, a node.
   * @param[in] @t, the time instant.
   * 
   * @return Force to the node, in a Point object.
   */
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    Point force = Point(0, 0, 0);
    for (auto edge_iterator = n.edge_begin(); edge_iterator != n.edge_end(); ++edge_iterator)
    {
      Edge adjacent_edge = *edge_iterator;
      Node adjacent_node = adjacent_edge.node2();
      double distance = norm(adjacent_node.position() - n.position());
      force -= adjacent_edge.value().K_ * (n.position() - adjacent_node.position()) * (distance - adjacent_edge.value().rest_length_) / distance;
    }
    return force;
  }
};

/** @struct DampingForce.
 * @brief Functor to apply a damping force to the nodes of the graph.
 */
struct DampingForce : Force
{
  double c_;

  /**
   * @brief Constructor
   * @param[in] @a c, parameter of the damping force.
   */
  DampingForce(double c = 1) : c_(c) {}

  /**
   * @brief Applies damping force to a node.
   * @param[in] @a n, a node.
   * @param[in] @t, the time instant.
   * 
   * @return Force to the node, in a Point object.
   */
  virtual Point operator()(Node n, double t)
  {
    (void)t;
    return -c_ * n.value().vel;
  }
};

/** @struct CombinedForce.
 * @brief Functor to apply a combined force to the nodes of the graph.
 */
struct CombinedForce : Force
{
  std::vector<Force *> force_vector_;

  /**
   * @brief Constructor
   * @param[in] @a force_vector, a vector of forces to combine.
   */
  CombinedForce(std::vector<Force *> &force_vector) : force_vector_(force_vector)
  {
  }

  /**
   * @brief Applies combined force to a node.
   * @param[in] @a n, a node.
   * @param[in] @t, the time instant.
   * 
   * @return Force to the node, in a Point object.
   */
  virtual Point operator()(Node n, double t)
  {
    Point force = Point(0, 0, 0);
    for (auto force_it = force_vector_.begin(); force_it != force_vector_.end(); ++force_it)
    {
      force += (**force_it)(n, t);
    }
    return force;
  }
};

/**
 * @brief creates a CombinedForce as a combination of forces.
 * @param[in] @a force1, first force to combine.
 * @param[in] @a force2, second force to combine.
 * @return CombinedForce of the input forces.
 */
template <typename T1, typename T2>
CombinedForce make_combined_force(T1 force1, T2 force2)
{
  std::vector<Force *> force_vector;
  force_vector.push_back(&force1);
  force_vector.push_back(&force2);
  return CombinedForce(force_vector);
};

/**
 * @brief creates a CombinedForce as a combination of forces.
 * @param[in] @a force1, first force to combine.
 * @param[in] @a force2, second force to combine.
 * @param[in] @a force3, third force to combine.
 * @return CombinedForce of the input forces.
 */
template <typename T1, typename T2, typename T3>
CombinedForce make_combined_force(T1 force1, T2 force2, T3 force3)
{
  std::vector<Force *> force_vector;
  force_vector.push_back(&force1);
  force_vector.push_back(&force2);
  force_vector.push_back(&force3);
  return CombinedForce(force_vector);
};

/** @struct Problem1Force.
 * @brief Functor to apply the force of Problem 1 to the nodes.
 */
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
    (void)t;
    // HW2 #1: YOUR CODE HERE
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
    {
      return Point(0, 0, 0);
    }

    Point force = Point(0, 0, 0);
    for (auto edge_iterator = n.edge_begin(); edge_iterator != n.edge_end(); ++edge_iterator)
    {
      Edge adjacent_edge = *edge_iterator;
      Node adjacent_node = adjacent_edge.node2();
      double distance = norm(adjacent_node.position() - n.position());
      force -= adjacent_edge.value().K_ * (n.position() - adjacent_node.position()) * (distance - adjacent_edge.value().rest_length_) / distance;
    }
    force += n.value().mass * Point(0, 0, -grav);
    return force;
  }
};

/** @struct PinConstraint.
 * @brief Functor to apply a pin constraint to the graph.
 */
struct PinConstraint : Constraint
{
  /**
   * @brief Applies pin constraint to the graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  void operator()(GraphType &g, double t)
  {
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      if (n.value().init_pos == Point(0, 0, 0) || n.value().init_pos == Point(1, 0, 0))
      {
        n.position() = n.value().init_pos;
        n.value().vel = Point(0, 0, 0);
      }
      n.value().t = t;
    }
    return;
  }
};

/** @struct PlaneConstraint.
 * @brief Functor to apply a plane constraint to the graph.
 */
struct PlaneConstraint : Constraint
{
  double z = -0.75;

  /**
   * @brief Applies plane constraint to the graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  void operator()(GraphType &g, double t)
  {
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      Point &position = n.position();
      Point &vel = n.value().vel;
      if (dot(position, Point(0, 0, 1)) < z)
      {
        position = position * Point(1, 1, 0) + Point(0, 0, z);
        vel = vel * Point(1, 1, 0);
      }

      n.value().t = t;
    }
    return;
  }
};

/** @struct SphereConstraint.
 * @brief Functor to apply a sphere constraint to the graph.
 */
struct SphereConstraint : Constraint
{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;

  /**
   * @brief Applies sphere constraint to the graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  void operator()(GraphType &g, double t)
  {
    for (auto it = g.node_begin(); it != g.node_end(); ++it)
    {
      auto n = *it;
      Point &position = n.position();
      Point &vel = n.value().vel;
      if (norm(position - c) < r)
      {
        Point R = (position - c) / norm(position - c);
        position = r * R + c;
        vel = vel - dot(vel, R) * R;
      }

      n.value().t = t;
    }
    return;
  }
};

/** @struct TearConstraint.
 * @brief Functor to apply a tear constraint to the graph.
 */
struct TearConstraint : Constraint
{
  Point c = Point(0.5, 0.5, -0.5);
  double r = 0.15;

  /**
   * @brief Applies tear constraint to the graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  void operator()(GraphType &g, double t)
  {
    (void)t;
    auto it = g.node_begin();
    while (it != g.node_end())
    {
      auto n = *it;
      Point &position = n.position();
      if (norm(position - c) < r)
      {
        g.remove_node(n);
        it = g.node_begin();
      }
      else
      {
        ++it;
      }
    }
    return;
  }
};

/** @struct CombinedConstraints.
 * @brief Functor to apply combined constraints to the graph.
 */
struct CombinedConstraints : Constraint
{
  std::vector<Constraint *> constraint_vector_;

  /**
   * @brief Constructor.
   * @param[in] @a constraint_vector, vector of constraints.
   */
  CombinedConstraints(std::vector<Constraint *> &constraint_vector) : constraint_vector_(constraint_vector)
  {
  }

  /**
   * @brief Applies combined constraints to a graph.
   * @param[in,out] @a g, a graph.
   * @param[in] @t, the time instant.
   */
  void operator()(GraphType &g, double t)
  {
    for (auto constraint_it = constraint_vector_.begin(); constraint_it != constraint_vector_.end(); ++constraint_it)
    {
      (**constraint_it)(g, t);
    }
    return;
  }
};

/**
 * @brief creates a CombinedConstraints as a combination of constraints.
 * @param[in] @a constraint1, first constraint to combine.
 * @param[in] @a constraint2, second constraint to combine.
 * @return CombinedConstraints of the input constraints.
 */
template <typename T1, typename T2>
CombinedConstraints make_combined_constraints(T1 constraint1, T2 constraint2)
{
  std::vector<Constraint *> constraint_vector;
  constraint_vector.push_back(&constraint1);
  constraint_vector.push_back(&constraint2);
  return CombinedConstraints(constraint_vector);
};

/**
 * @brief creates a CombinedConstraints as a combination of constraints.
 * @param[in] @a constraint1, first constraint to combine.
 * @param[in] @a constraint2, second constraint to combine.
 * @param[in] @a constraint3, third constraint to combine.
 * @return CombinedConstraints of the input constraints.
 */
template <typename T1, typename T2, typename T3>
CombinedConstraints make_combined_constraints(T1 constraint1, T2 constraint2, T3 constraint3)
{
  std::vector<Constraint *> constraint_vector;
  constraint_vector.push_back(&constraint1);
  constraint_vector.push_back(&constraint2);
  constraint_vector.push_back(&constraint3);
  return CombinedConstraints(constraint_vector);
};

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


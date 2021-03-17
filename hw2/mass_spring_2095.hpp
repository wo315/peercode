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
static constexpr double K = 100;
double damp_coef = 0;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double K;          //< Spring constant
  double L;          //< The initial edge length
  EdgeData(): K(0), L(0) {}
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
    //Skip update srep for boundary nodes
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
    	continue;
    }
    
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
    // HW2 #1: YOUR CODE HERE
    (void) t;    // silence compiler warnings
    
    //Prevents cloth from falling to infinity by returning a zero force at end points
    if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
	return Point(0,0,0);
      }
    //Zero initialize spring force
    Point spring_force = Point(0, 0, 0);
    for(auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      auto edge_ = *iter;
      Point diff = Point(0, 0, 0);
      if(edge_.node1() == n)
	diff = n.position() - edge_.node2().position();
      else
	diff = edge_.node1().position() - n.position();
      spring_force += -edge_.value().K * (diff / norm(diff)) * (norm (diff) - edge_.value().L);
    }
    //Return spring_force + gravitational force
    return spring_force + (n.value().mass * Point(0, 0, -grav));
  }
};

/** Class Force, returns zero force */
template <typename NODE>
class Force {
public:
  //Virtual Destructor
  virtual ~Force() {}
  
  //Virtual method
  virtual Point operator()(NODE n, double t) const {
    (void) n;   // silence compiler warnings
    (void) t;   // silence compiler warnings
    return Point(0, 0, 0);
  }
};

/** Force function object for spring gravitational force. */
class  GravityForce: public Force<Node> {
public:
  template <typename NODE>
  Point operator()(NODE n, double t) const {
    (void) t; // silence compiler warnings
    return n.value().mass * Point(0, 0, -grav);
  }
};

/** Force function object for spring force. */
class MassSpringForce: public Force<Node> {
public:
  template <typename NODE>
  Point operator()(NODE n, double t) const {
    (void) t;     // silence compiler warnings
    Point spring_force = Point(0, 0, 0);
    for(auto iter = n.edge_begin(); iter != n.edge_end(); ++iter) {
      auto edge_ = *iter;
      Point diff = Point(0, 0, 0);
      if(edge_.node1() == n)
	diff = n.position() - edge_.node2().position();
      else
	diff = edge_.node1().position() - n.position();
      spring_force += -edge_.value().K * (diff / norm(diff)) * (norm (diff) - edge_.value().L);
  }
    return spring_force;
  }
};

/** Force function object for damping force. */
class DampingForce: public Force<Node> {
public:
  template <typename NODE>
  Point operator()(NODE n, double t) const {
    (void) t; // silence warnings
    return -damp_coef * n.value().vel;
  }
};

/** 
 * @brief Combines two forces together by adding them.
 * 
 * @tparam Force1 A function object which is called by @a f1(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 * @tparam Force2 A function object which is called by @a f2(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 */
template<typename Force1, typename Force2>
class CombinedForce: public Force<Node> {
private:
  Force1 f1;
  Force2 f2;
public:
  CombinedForce(Force1 f1_, Force2 f2_) :
    f1(f1_),
    f2(f2_) {}

  //Destructor
  ~CombinedForce() {}

  /** 
   * @brief Returns the sum of the forces @a n at time @a t.\
   * 
   * @param[in] n      node in the graph
   * @param[in] t      the current time
   * @return the sum of the two forces @a n
   *
   * @tparam NODE is a class representing the graph nodes.
   */  
  template<typename NODE>
  Point operator()(NODE  n, double  t) const {
    return f1(n, t) + f2(n, t);
  }
};

/** 
 * @brief Combines two forces together.
 * 
 * @tparam Force1 A function object which is called by @a f1(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 * @tparam Force2 A function object which is called by @a f2(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 */
template<typename Force1, typename Force2>
CombinedForce<Force1, Force2> make_combined_force(Force1 f1, Force2 f2) {
  return {f1, f2};
}

/** 
 * @brief Combines three forces together.
 * 
 * @tparam Force1 A function object which is called by @a f1(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 * @tparam Force2 A function object which is called by @a f2(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 * @tparam Force3 A function object which is called by @a f3(n, @a t)
 *                Here, @a n is a node of the graph and @a t is the current time
 */
template<typename Force1, typename Force2, typename Force3>
CombinedForce<CombinedForce<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3) {
 return make_combined_force(make_combined_force(f1,f2), f3);
}


/** Class Constraint */
template <typename G>
class Constraint {
public:
  //Virtual Destructor
  virtual ~Constraint() {}
  
  //Virtual method
  virtual void operator()(G& graph, double t) {
    (void) graph;   // silence compiler warnings
    (void) t;       // silence compiler warnings
  }
};


//template<typename G>
class PinConstraint: public Constraint<GraphType> {
public:
  template<typename G>
  void operator()(G& graph, double t) {
    (void) t;       // silence compiler warnings
    
    for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
      auto n = *iter;
      if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
	n.value().vel = Point(0, 0, 0);
	}
    }
  }
};


class PlaneConstraint: public Constraint<GraphType> {
public:
  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;       // silence compiler warnings
    
    for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
      auto n = *iter;
      if(n.position().z < -0.75) {
	n.position().z = -0.75;
	n.value().vel.z = 0;
	}
    }
  }
};


class SphereConstraint: public Constraint<GraphType> {
public:
  const Point center = Point(0.5, 0.5, -0.5);
  const double radius = 0.15;
  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;       // silence compiler warnings
    
    for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
      auto n = *iter;
      Point diff = n.position() - center;
      if(norm(diff) < radius) {
	Point R = (diff/norm(diff));
	n.position() = center + R * radius;
	n.value().vel -= inner_prod(n.value().vel, R) * R;
      }
    }
  }
};


class SphereConstraint2: public Constraint<GraphType> {
public:
  const Point center = Point(0.5, 0.5, -0.5);
  const double radius = 0.15;
  template <typename G>
  void operator()(G& graph, double t) {
    (void) t;       // silence compiler warnings

    for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter) {
      auto n = *iter;
      Point diff = n.position() - center;
//--design_1
//--You should not increment the iterator when you actually remove a node
//--START
      if(norm(diff) < radius) {
	graph.remove_node(n);
//--END
      }
    }
  }
};


template<typename constraint_type1, typename constraint_type2>
class CombinedConstraints: public Constraint<GraphType> {
private:
  constraint_type1 c1;
  constraint_type2 c2;
public:
  CombinedConstraints(constraint_type1 c1_, constraint_type2 c2_) :
    c1(c1_),
    c2(c2_) {}
  
  //Destructor
  ~CombinedConstraints() {}

  template<typename G>
  void operator()(G& graph, double t) {
    c1(graph, t);
    c2(graph, t);
  }
};

template<typename constraint_type1, typename constraint_type2>
CombinedConstraints<constraint_type1, constraint_type2> make_combined_constraints(constraint_type1 c1, constraint_type2 c2) {
  return {c1, c2};
}

/** 
 * @brief Combines three constraints together
 * 
 * @tparam constraint_type1 is a function object called as @a c1(@a, t)
 *                          Here , @a t is the current time and @a a must be applied on all nodes at time @a t.
 * @tparam constraint_type2 is a function object called as @a c2(@a, t)
 *                          Here , @a t is the current time and @a a must be applied on all nodes at time @a t
 * @tparam constraint_type3 is a function object called as @a c3(@a, t)
 *                          Here , @a t is the current time and @a a must be applied on all nodes at time @a t
 */
template<typename constraint_type1, typename constraint_type2, typename constraint_type3>
CombinedConstraints<CombinedConstraints<constraint_type1, constraint_type2>, constraint_type3> make_combined_constraints(constraint_type1 c1, constraint_type2 c2, constraint_type3 c3) {
  return make_combined_constraints(make_combined_constraints(c1, c2), c3);
}

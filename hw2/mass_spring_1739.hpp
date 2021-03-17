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
//#include "shortest_path.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
static constexpr int k = 100;



/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double length;
  double constant;
  EdgeData() : length(0), constant(100) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using node_iterator = typename GraphType::NodeIterator;
using NodeIter = typename GraphType::NodeIterator;

// copy and pasted from shortest_path.hpp to avoid type redefinitions

/**
 * @brief Custom modification of std::min_element to override forwarditerator
 * @param first a NodeIterator to start at
 * @param last an end NodeIterator
 * @param comp a comparison functor
 * @return An iterator to the minimum, as determined by comp functor
 *           last if we reach end without finding a minimum
 */
template<class Compare>
  NodeIter custom_min(NodeIter first, NodeIter last, Compare comp) {
    if (first == last) return last;

    NodeIter smallest = first;
    ++first;
    for (; first != last; ++first) {
        if (comp(*first, *smallest)) {
            smallest = first;
        }
    }
    return smallest;
  }

// copy and pasted froms shortest_path.hpp to avoid type redefinitions

/** Find the node with the minimum euclidean distance to a point.
 * @param g  The graph of nodes to search.
 * @param point  The point to use as the query.
 * @return An iterator to the node of @a g with the minimun Eucliean
 *           distance to @a point.
 *           graph.node_end() if graph.num_nodes() == 0.
 *
 * @post For all i, 0 <= i < graph.num_nodes(),
 *          norm(point - *result) <= norm(point - g.node(i).position())
 */
node_iterator nearest_node (const GraphType& g, const Point& point) {
  // HW1 #3: YOUR CODE HERE

  /**
   * @brief Functor that stores x, y, z coordinates of a point to compare all
   * nodes' points against
   */
  class Eucd {
    double x;
    double y;
    double z;
   public:

    /**
     * @brief constructor for Eucd functor
     * @param[in] point a point object to compare all nodes against
     * @post functor is constructed and ready to call ()
     */
    Eucd(const Point& point) : x{point.x}, y{point.y}, z{point.z} {}
    /**
     * @brief overloaded call operator to compare two nodes to return if one
     * closer to the point used to construct functor
     * @param[in] node1 a node with a point associated
     * @param[in] node2 a second node with a point associated
     * @return returns true is node1 is closer to construction point than
     * node2
     */
    bool operator() (Node node1, Node node2) const {
      return sqrt(pow(node1.position().x-x,2.0) + pow(node1.position().y-y,2)\
        + pow(node1.position().z-z,2))\
        < sqrt(pow(node2.position().x-x,2) + pow(node2.position().y-y,2) +\
        pow(node2.position().z-z,2));
    }
  };
  return custom_min(g.node_begin(), g.node_end(), Eucd(point));



}

//forward declare constraint class
class Constraint;

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
template <typename G, typename F, typename C = Constraint>
double symp_euler_step(G& g, double t, double dt, F force, C constraint = C()) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    n.position() += n.value().vel * dt;
  }
  // Apply all constraints
  constraint(g,t);

  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

/**
 * @brief Parent functor class that takes node and double t and
 * returns a point at 0,0,0 to be overloaded
 * by child classes
 */
class Force{
  public:
  /**
   * @brief () operator that returns point at 0,0,0 to be
   * overloaded
   * @param[in] Node n
   * @param[in] double t (not currently used)
   * @return point at 0,0,0
   */
  virtual Point operator() (Node n, double t) {
  (void) n, (void) t; //silence compiler
  // https://stackoverflow.com/questions/7978620/whats-a-portable-way-to-implement-no-op-statement-in-c?rq=1
  return Point(0,0,0);
  }
};

/**
 * @brief Functor class that takes node and double t and
 * returns a gravity force in negative z direction
 */
class GravityForce: public Force {
 public:
 /**
   * @brief () operator that returns point
   * that acts as vector with -g in z
   * direction
   * @param[in] Node n
   * @param[in] double t (not currently used)
   * @return point/vector with -g*mass
   */
  Point operator()(Node n, double t) {
    (void) t;
    return Point (0,0,-grav*n.value().mass);
  }
};

/**
 * @brief Functor class that takes node and double t and
 * returns a spring force based on its edges
 */
class MassSpringForce: public Force {
 public:
 /**
   * @brief () operator that returns point
   * that acts as vector with forces based
   * on node's edges
   * @param[in] Node n
   * @param[in] double t (not currently used)
   * @return point/vector with sum of forces
   * based on mass-spring equation
   */
  Point operator()(Node n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;
    Point spring = Point(0,0,0);
    // iterate over all incident nodes, deref as an edge, and calculate
    // mass_spring force that edge, sum and return over all edges
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
      auto this_edge = *it;
        Point temp_pos = this_edge.node2().position();
        double distance = norm(n.position()-temp_pos);
        spring += -(this_edge.value().constant) * (n.position()-temp_pos)\
          *(distance-(this_edge.value().length))/(distance);
    }
    return spring;
  }
};

/**
 * @brief Functor class that takes node and double t and
 * returns a damping force vector
 */
class DampingForce: public Force {
  double constant_;
 public:
   /**
   * @brief default constructor with
   * creates damping force with constant
   * of 0
   */
  DampingForce() : constant_{0.0} {}
  /**
   * @brief constructor for functor that stores input as
   * as force constant
   * @param[in] constant to act as force constant
   */
  DampingForce(double& constant) : constant_{constant} {}
 /**
   * @brief () operator that returns point
   * that acts as vector with forces based
   * on damping constant
   * @param[in] Node n
   * @param[in] double t (not currently used)
   * @return point/vector with forces based
   * on damping constant
   */
  Point operator()(Node n, double t) {
    (void) t;
    return n.value().vel*constant_*-1;
  }
};

/**
 * @brief Functor class that stores a vector of pointers
 * to inherited force functors to sum their forces
 */
class CombinedForce {
  std::vector<Force*> forces_;
 public:
  /**
   * @brief constructor for functor that stores vector
   * of pointers to inherited force instances
   * @param[in] vector of pointers to forces
   */
  CombinedForce(std::vector<Force*> forces) : forces_{forces} {}
 /**
   * @brief () operator that returns point
   * that acts as vector that sums all forces
   * in forces vector
   * @param[in] Node n to apply all forces to
   * @param[in] double t (not currently used)
   * @return point/vector with sum of forces applied
   */
  Point operator()(Node n, double t) {
    Point combined = Point(0,0,0);
    for(auto it=forces_.begin(); it!=forces_.end(); ++it) {
      combined += (**it)(n, t);
    }
  return combined;
  }
};

/**
 * @brief Parent functor class that takes node and double t and
 * with virtual operator to be overloaded
 * by child classes
 */
class Constraint {
 public:
  /**
   * @brief () operator to be
   * overloaded
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */
  virtual void operator() (GraphType& g, double t) {
    (void) g, (void) t;
  }
};

/**
 * @brief Functor class inheriting from Constraint
 * that pins a node if it was the nearest node to
 * 0,0,0 or 1,0,0 at t=0
 */
class PinConstraint: public Constraint {
 Node pin1_;
 Node pin2_;
 public:
  /**
   * @brief constructor for functor that stores two node
   * uid's based on what two nodes were closed to 0,0,0
   * and 1,0,0 at t=0
   * @param[in] g reference to a graph
   */
  PinConstraint(GraphType& g) {
    Point pin_point = Point(0,0,0);
    pin1_ = *nearest_node(g, pin_point);
    pin_point = Point(1,0,0);
    pin2_ = *nearest_node(g,pin_point);
  }
  /**
   * @brief () operator that keeps pinned
   * nodes whose uid matches those that were closest
   * to 0,0,0 or 1,0,0 at t=0
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */
  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if(n==pin1_) {
        n.position()=Point(0,0,0);
      }
      if(n==pin2_) {
        n.position()=Point(1,0,0);
      }
    }
  }
};

/**
 * @brief Functor class inheriting from Constraint
 * that keeps a node from crossing a plane
 */
class PlaneConstraint: public Constraint {
  double z_;
 public:
 /**
 * @brief constructor
 * that keeps a node from crossing a plane at -.75
 */
  PlaneConstraint() : z_{-.75} {}
  /**
   * @brief () operator that keeps nodes from crossing
   * plane at z=-.75 and sets any such nodes to have
   * z-direction forces to 0
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */
  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;
      if(n.position().z<z_) {
        n.position().z=z_;
        n.value().vel.z=0;
      }
    }
  }
};

/**
 * @brief Functor class inheriting from Constraint
 * that keeps a node from crossing a sphere
 */
class SphereConstraint: public Constraint {
  Point center_;
  double radius_;
 public:
 /**
 * @brief constructor
 * that keeps a node from crossing a sphere
 */
  SphereConstraint() : center_{Point(.5,.5,-.5)}, radius_{.15} {}
  /**
   * @brief () operator that keeps nodes from crossing
   * sphere with center (.5,.5,-.5) and radius .15
   * and sets their normal velocity to 0
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */

  void operator() (GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      auto n = *it;

      if(norm(n.position()-center_)<radius_) {
      // https://piazza.com/class/kjx6t0a2rv4589?cid=236
        n.position()=center_+radius_*(n.position()-center_)/norm(n.position()\
          -center_);
        Point Ri=(n.position()-center_)/norm(n.position()-center_);
        n.value().vel=n.value().vel-(dot(n.value().vel,Ri)*Ri);
      }
    }
  }
};

/**
 * @brief Functor class inheriting from Constraint
 * that removes a node and its edges if it crosses a sphere
 */
class TearConstraint: public Constraint {
  Point center_;
  double radius_;
 public:
  TearConstraint() : center_{Point(.5,.5,-.5)}, radius_{.15} {}
  /**
   * @brief () operator that removes nodes and
   * their edges if they cross
   * sphere with center (.5,.5,-.5) and radius .15
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */
  void operator() (GraphType& g, double t) {
    (void) t;
    auto it = g.node_begin();
    // loops over all nodes, re-trying a node if it violates the constraint
    // because it swaps with a new valid node
    while(it!=g.node_end()) {
      auto n = *it;

      // similar to Ex2, will keep applying constraint on it until it no
      // longer violates constraint, as we are swapping in new nodes
      if(norm(n.position()-center_)<radius_) {
      // https://piazza.com/class/kjx6t0a2rv4589?cid=236
        g.remove_node(n);
      }
      else{
        ++it;
      }
    }
  }
};

/**
 * @brief Functor class that stores a vector of pointers
 * to inherited constraints to apply each consecutively
 */

class CombinedConstraint {
  std::vector<Constraint*> constraints_;
 public:
  /**
   * @brief constructor for functor that stores vector
   * of pointers to inherited constraint instances
   * @param[in] constraints vector of pointers
   */
  CombinedConstraint(std::vector<Constraint*> constraints) :\
    constraints_{constraints} {}
 /**
   * @brief () operator that iterates over vector of pointers to constraints
   * to apply each constraint consecutively
   * calling each constraint's call operator
   * @param[in] g reference to a graph
   * @param[in] double t (not currently used)
   */
  void operator()(GraphType& g, double t) {
    for(auto it=constraints_.begin(); it!=constraints_.end(); ++it) {
      (**it)(g, t);
    }
  }
};

/** Create a CombinedConstraint class instance from 2 or 3 constraints
 * that inherit from Constraint class
 * @param[in] constraint_a first constraint
 * @param[in] constraint_b second constraint
 * @param[in] constraint_c optional constraint that defaults to constraint
 * instance with operator that runs without doing anything
 * @return CombinedConstraint instance that will apply each constraint
 * to graph consecutively
 *
 * @tparam A class that inherits from constraint
 * @tparam B class that inherits from constraint
 * @tparam C optional class that inherits from constraint
 */

//--design_1
//--argument objects passed by value will be destroyed after function finishes
//--START
template<typename A,typename B,typename C=Constraint>
CombinedConstraint make_combined_constraint(A constraint_a, B constraint_b,\
  C constraint_c=C()) {
  Constraint* p1 = &constraint_a;
  std::vector<Constraint*> internal_vector;
  internal_vector.push_back(p1);
  Constraint* p2 = &constraint_b;
  internal_vector.push_back(p2);
  Constraint* p3 = &constraint_c;
  internal_vector.push_back(p3);
  return CombinedConstraint{internal_vector};
}
//--END

/** Create a CombinedForce class instance from 2 or 3 forces
 * that inherit from Force class
 * @param[in] force_a first force
 * @param[in] force_b second force
 * @param[in] force_c optional force that defaults to base force
 * instance that effectively adds 0 force
 * @return CombinedForce instance that will apply each constraint
 * to node
 *
 * @tparam A class that inherits from force
 * @tparam B class that inherits from force
 * @tparam C optional class that inherits from force
 */
template<typename A,typename B,typename C=Force>
CombinedForce make_combined_force(A force_a, B force_b, C force_c =C()) {
  Force* p1 = &force_a;
  std::vector<Force*> internal_vector;
  internal_vector.push_back(p1);
  Force* p2 = &force_b;
  internal_vector.push_back(p2);
  Force* p3 = &force_c;
  internal_vector.push_back(p3);
  return CombinedForce{internal_vector};
}

/** Force function object for HW2 #1. */
struct Problem1Force {
  /** Return the force applying to @a n at time @a t.
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
 /**
   * @brief () operator that returns point
   * that acts as vector with forces based
   * on gravity and spring on each node
   * @param[in] Node n
   * @param[in] double t (not currently used)
   * @return point/vector with forces based
   * on gravity and spring of each node
   */
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t; // silence compiler
    Point spring = Point(0,0,0);

    // iterate over all incident edges to sum forces for each according to
    // mass_spring formula
    for(auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
      Point temp_pos = (*it).node2().position();
      double distance = norm(n.position()-temp_pos);
      spring += -k * (n.position()-temp_pos)*(distance-(*it)\
        .value().length)/(distance);
    }
    spring+=Point(0,0,-grav*n.value().mass);
    if (n.position () ==  Point (0,0,0) || n.position () ==  Point (1,0,0)) {
      return  Point (0,0,0);
    }
    return spring;
  }
};




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
  Point initpos; //< Node initial position
  NodeData() : vel(0), mass(1), initpos(0) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
  double springconst;  //< K spring constant
  double restlen;      //< Spring rest length
  EdgeData() : springconst(100.0), restlen(1.0) {}
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
 * @param[in]     constraint Function object defining the constraints enforced
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData struct YOU CHOOSE
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

  // apply constraints
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
    (void) t; // silence compiler warning, t unused in HW2
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
      return Point(0, 0, 0);
    }

    // generate a 0 force to add to
    Point forcespring = Point(0.0);
    // iterate over incident edges and add spring forces
    for (GraphType::IncidentIterator incit = n.edge_begin(); \
         incit != n.edge_end(); ++incit) {
      Edge e = (*incit);
      Node adj = e.node2();
      Point diff = n.position() - adj.position();
      forcespring += -e.value().springconst*(norm(diff) - e.value().restlen)\
                     *diff/(norm(diff));
    }

    double mass = n.value().mass;
    return forcespring + Point(0.0,0.0,-mass*grav);
  }
};


/** @class Force
 * @brief The base force class that other forces are derived from
 */
class Force {

  public:
    /**
     * @brief overloaded () operator to return a force vector
     * @param[in] n The node which force is being calculated for
     * @param[in] t Double of time that is an input to force
     */
    virtual Point operator()(Node n, double t) {
      (void) t; (void) n; // silence compiler warning, t unused in HW2
      return Point(0.0, 0.0, 0.0);
    }

    virtual ~Force() {
    }

};

/** @class GravityForce
 * @brief Derived force class representing gravity along z-axis
 */
class GravityForce: public Force {
  // moved gravity constant in here
  static constexpr double grav = 9.81;
  public:
    /**
     * @brief overloaded () operator to return a force vector
     * @param[in] n The node which force is being calculated for
     * @param[in] t Double of time that is an input to force
     */
    virtual Point operator()(Node n, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      return Point(0.0, 0.0, -n.value().mass*grav);
    }

    virtual ~GravityForce() {
    }
};

/** @class MassSpringForce
 * @brief Derived force class representing spring force attached to node
 */
class MassSpringForce: public Force {

  public:
    /**
     * @brief overloaded () operator to return a force vector
     * calculates spring force on a per edge basis
     * @param[in] n The node which force is being calculated for
     * @param[in] t Double of time that is an input to force
     */
    virtual Point operator()(Node n, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      Point forcespring = Point(0.0);
      // iterate over incident edges and add spring forces
      for (GraphType::IncidentIterator incit = n.edge_begin(); \
           incit != n.edge_end(); ++incit) {
        Edge e = (*incit);
        Node adj = e.node2();
        Point diff = n.position() - adj.position();
        forcespring += -e.value().springconst*(norm(diff) - e.value().restlen)\
                       *diff/(norm(diff));
      }
      return forcespring;
    }

    virtual ~MassSpringForce() {
    }
};

/** @class DampingForce
 * @brief Derived force class representing damping force, currently unused
 */
class DampingForce: public Force {
  private:
    double dconstant_;

  public:
    // Constructor
    DampingForce(double dampingconstant): dconstant_(dampingconstant) {}

    /** @brief overloaded () operator to return a force vector
     * based on damping constant and velocity
     * @param[in] n The node which force is being calculated for
     * @param[in] t Double of time that is an input to force
     */
    virtual Point operator()(Node n, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      return -dconstant_*n.value().vel;
    }

    virtual ~DampingForce() {
    }
};

/** @struct CombinedForce
 * @brief Non-derived functor which combines force vectors of different
 * types into one returned force
 */
struct CombinedForce {
  std::vector<Force*> forces_;

  // Constructor
  CombinedForce(std::vector<Force*> forces) : forces_(std::move(forces)) {}

  /**
   * @brief Overloaded operator to add forces together
   * @param[in] n Node to calculate total force on
   * @param[in] t current time-step of simulation
   */
  Point operator()(Node n, double t) {
    (void) t; // silence compiler warning, t unused in HW2
    Point totalforce = Point(0.0, 0.0, 0.0);
    for (std::vector<Force*>::iterator fit = forces_.begin(); \
         fit != forces_.end(); ++fit) {
      //deref twice and apply () operator
      totalforce += (*(*fit))(n, t);
    }
    return totalforce;
  }
};

/**
 * @brief Base Case helper function for variadic template to make
 * combined force
 * @param[in] F vector of force pointers to add forces of possibly
 * different types
 * @tparam f1type Type of force to add to _F_
 * @param[in] f1 first force to add to _F_
 */
template <typename f1type>
void createforcevec(std::vector<Force*>& F, f1type& f1) {
  // break recursion here
  Force* fptr = &f1;
  F.push_back(fptr);
}

/**
 * @brief Recursive Case helper function for variadic template to
 * make combined force
 * @param[in] F vector of force pointers to add forces of possibly
 * different types
 * @tparam f1type Type of first force to add to _F_
 * @tparam ftypes Type of other forces to add to _F_
 * @param[in] f1 first force to add to _F_
 * @param[in] fi other forces to add to _F_
 */
template<typename f1type, typename ... ftypes>
void createforcevec(std::vector<Force*>& F, f1type& f1, ftypes&...fi) {
  Force* fptr = &f1;
  F.push_back(fptr);
  createforcevec(F, fi...);
}

/**
 * @brief Wrapper for recursive functions to create CombinedForce
 * from various types of forces
 * @tparam f1type Type of first force to add to _F_
 * @tparam ftypes Type of other forces to add to _F_
 * @param f1 first force to add to _F_
 * @param fi other forces to add to _F_
 * @return CombinedForce
 */
//--design_1
//--argument object passed by value is destroyed after function finishes
//--START
template<typename f1type, typename ...ftypes>
CombinedForce make_combined_force(f1type f1, ftypes...fi) {
  std::vector<Force*> F;
  //call to recursive function
  createforcevec(F, f1, fi...);
  return CombinedForce(F);
}
//--END

/**
 * @brief Base Constraint Class, defines a constraint on a graph
 * and hopefully fixes it
 */
class Constraint {

  public:
    virtual void operator()(GraphType& g, double t) {
      (void) t; (void) g; // silence compiler warning, t unused in HW2
    }

    virtual ~Constraint() {
    }

};

/**
 * @brief Derived PinConstraint pins points at (0,0,0) and (1,0,0)
 *
 */
class PinConstraint: public Constraint {

  public:
  /**
   * @brief overloaded operator fixing pinconstraint
   *
   * @param g Graph to apply PinConstraint to
   * @param t time-step of simulation
   */
    virtual void operator()(GraphType& g, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      for (GraphType::NodeIterator it = g.node_begin(); \
           it != g.node_end(); ++it) {
        Node n = (*it);
        if (n.value().initpos == Point(0,0,0) || \
            n.value().initpos == Point(1,0,0)) {
          n.position() = n.value().initpos;
        }
      }
    }

    virtual ~PinConstraint() {
    }

};

/**
 * @brief Derived PLaneConstraint fixes points to not pass z = -0.75
 *
 */
class PlaneConstraint: public Constraint {

  public:
    /**
     * @brief Overloaded operator to enforce plane constraint on
     * a graph.
     * @param g Graph to apply Plane Constraint to
     * @param t time-step of the simulation
     */
    virtual void operator()(GraphType& g, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      for (GraphType::NodeIterator it = g.node_begin(); it != g.node_end(); \
           ++it) {
        Node n = (*it);
        if (dot(n.position(), Point(0.0, 0.0, 1.0)) < -0.75) {
          //closest point on the plane
          n.position().z = -0.75;
          n.value().vel.z = 0.0;
        }
      }
    }

    virtual ~PlaneConstraint() {
    }

};

/**
 * @brief Derived Sphere Constraint fixes points to not fall within a sphere
 * with center (0.5,0.5,-0.5) and radius 0.15
 */
class SphereConstraint: public Constraint {
  Point center = Point(0.5,0.5, -0.5);
  double radius = 0.15;

  public:
    /**
     * @brief overloaded operator applying Sphere Constraint
     *
     * @param g Graph to apply Sphere constraint to
     * @param t time-step of simulation
     */
    virtual void operator()(GraphType& g, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      for (GraphType::NodeIterator it = g.node_begin(); it != g.node_end(); \
           ++it) {
        Node n = (*it);
        Point disttocenter = n.position() - center;
        Point R = disttocenter / norm(disttocenter);
        // find the closest point on the sphere
        if (norm(disttocenter) < radius) {
          n.position() = center + \
                        (radius * (disttocenter) / norm(disttocenter));
          n.value().vel -= inner_prod(n.value().vel, R)*R;
        }
      }
    }

    virtual ~SphereConstraint() {
    }

};

/**
 * @brief Derived TearConstraint class that creates a tear in the fabric
 * if it touches the Sphere with center (0.5,0.5, -0.5) and radius 0.15
 */
class TearConstraint: public Constraint {
  Point center = Point(0.5,0.5, -0.5);
  double radius = 0.15;

  public:
    /**
     * @brief Applies sphere tear to a graph at time t
     *
     * @param g Graph to apply sphere constraint to
     * @param t time-step of simulation
     */
    virtual void operator()(GraphType& g, double t) {
      (void) t; // silence compiler warning, t unused in HW2
      for (GraphType::NodeIterator it = g.node_begin(); \
           it != g.node_end(); ++it) {
        Node n = (*it);
        Point disttocenter = n.position() - center;
        if (norm(disttocenter) < radius) {
          g.remove_node(n);
        }
      }
    }

    virtual ~TearConstraint() {
    }

};

/**
 * @brief Functor to combine constraints together from a vector of
 * pointers to constraints
 */
struct CombinedConstraints {
  std::vector<Constraint*> constraints_;

  // Constructor
  CombinedConstraints(std::vector<Constraint*> constraints) : \
                      constraints_(std::move(constraints)) {}

  void operator()(GraphType& g, double t) {
    (void) t; // silence compiler warning, t unused in HW2
    for (std::vector<Constraint*>::iterator cit = constraints_.begin(); \
         cit != constraints_.end(); ++cit) {
      //deref twice and apply () operator
      (*(*cit))(g, t);
    }
  }
};

/**
 * @brief Base case to create constraint vector, see recursive case
 * for details
 */
template <typename c1type>
void createconstraintvec(std::vector<Constraint*>& C, c1type& c1) {
  // break recursion here
  Constraint* cptr = &c1;
  C.push_back(cptr);
}

/**
 * @brief Recursive case to add to constraint vector with variable number
 * of arguments
 * @tparam c1type First constraint's type
 * @tparam ctypes Other constraint's types
 * @param C Vector of pointer to constraints to add to
 * @param c1 First constraints to add
 * @param ci Other constraints to add
 */
template<typename c1type, typename ... ctypes>
void createconstraintvec(std::vector<Constraint*>& C, \
                         c1type& c1, ctypes&...ci) {
  Constraint* cptr = &c1;
  C.push_back(cptr);
  createconstraintvec(C, ci...);
}

/**
 * @brief Wrapper function that makes the CombinedCombined constraints from
 * variable number of constraints
 * @tparam c1type Type of constraint 1
 * @tparam ctypes Type of other constraints
 * @param c1 Constraint 1
 * @param ci Other constraints
 * @return CombinedConstraints
 */
template<typename c1type, typename ...ctypes>
CombinedConstraints make_combined_constraints(c1type c1, ctypes...ci) {
  std::vector<Constraint*> C;
  createconstraintvec(C, c1, ci...);
  return CombinedConstraints(C);
}




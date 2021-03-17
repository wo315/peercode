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
  Point init_position; // <initial position
  NodeData() : vel(0), mass(1), init_position(0) {}
};

struct EdgeData {
  double K;  //< Edge spring constant
  double L;  //< Edge rest length
  EdgeData() : K(1), L(1) {}
  EdgeData(double K, double L) : K(K), L(L) {}
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
 * @param[in]     constraint  Function object defining the constraints
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint( @a g, @a t) which
 *           reset the location and velocity of nodes that violates the constraints
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
  // Enfoce the constraint
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
    // HW2 #1: YOUR CODE HERE
    (void)t;
    if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
      return Point(0);
    }

    Point result(0, 0, -n.value().mass * grav);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double len = (*it).length();
      EdgeData edgeinfo = (*it).value();
      double L = edgeinfo.L;
      double K = edgeinfo.K;
      result += (K * (L - len)/len) * (n.position() - (*it).node2().position());
    }
    return result;
  }
};

/**
 * @brief A parent force class that different forces will inherit from.
 *
 */
class Force {
 public:
  virtual Point operator()(Node n, double t) {
    (void)t;
    (void)n;
    return Point(0);
  }
  virtual ~Force() {}
  virtual Force* clone() const {return new Force;}
};

/**
 * @brief Class that returns the gravity force.
 *
 */
class GravityForce : public Force {
 public:
  virtual Point operator()(Node n, double t) {
    (void)t;
    return Point(0, 0, -n.value().mass * grav);
  }

  virtual GravityForce* clone() const {return new GravityForce;}
};

/**
 * @brief Class that returns the mass spring force of the node.
 *
 */
class MassSpringForce : public Force {
 public:
  virtual Point operator()(Node n, double t) {
    (void)t;
    Point result(0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
      double len = (*it).length();
      EdgeData edgeinfo = (*it).value();
      double L = edgeinfo.L;
      double K = edgeinfo.K;
      result +=
          (K * (L - len) / len) * (n.position() - (*it).node2().position());
    }
    return result;
  }

  virtual MassSpringForce* clone() const {return new MassSpringForce;}
};

/**
 * @brief Class that implements a damping force.
 *
 */
class DampingForce : public Force {
  double c;

 public:
  virtual Point operator()(Node n, double t) {
    (void)t;
    return -c * n.value().vel;
  }
  DampingForce(double c) : c(c) {}
  DampingForce() : c(1) {}

  virtual DampingForce* clone() const {return new DampingForce(*this);}
};

/**
 * @brief A functor that combines arbitrary number of forces
 *
//  */
class CombinedForce {
  std::vector<Force*> force_vec;

 public:
  Point operator()(Node n, double t) {
    (void)t;
    Point result(0);
    for (unsigned i = 0; i < force_vec.size(); ++i) {
      result += force_vec[i]->operator()(n, t);
    }
    return result;
  }

  template <typename F>
  void add_force(F f) {
    //--design_1
    //--add force should add f, but not construct a new force
    //--START
    (void) f;
    force_vec.push_back(new F);
    //--END
  }

  // Redefine copy constructor
  CombinedForce(const CombinedForce& other) {
    for (unsigned i = 0; i < other.force_vec.size(); ++i) {
      force_vec.push_back(other.force_vec[i]->clone());
    }
  }
  CombinedForce() {}

  // Copy assignment
  CombinedForce& operator=(const CombinedForce& other){
    for (unsigned i = 0; i < force_vec.size(); ++i) {
      delete force_vec[i];
    }
    force_vec.resize(0);
    for(unsigned i = 0; i < other.force_vec.size(); ++i){
      force_vec.push_back(other.force_vec[i]->clone());
    }
    return *this;
  }


  ~CombinedForce() {
    for (unsigned i = 0; i < force_vec.size(); ++i) {
      delete force_vec[i];
    }
  }
};

/**
 * @brief get a CombinedForce functor with three forces
 *
 * @tparam F1 A child class of Force
 * @tparam F2 A child class of Force
 * @tparam F3 A child class of Force
 * @param f1 First force to apply
 * @param f2 Second force to apply
 * @param f3 Third force to apply
 * @return CombinedForce CombinedForce the combined force of f1, f2 and f3
 */
template <typename F1, typename F2, typename F3=Force>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3=Force()) {
  CombinedForce result;
  result.add_force(f1);
  result.add_force(f2);
  result.add_force(f3);
  return result;
}

/**
 * @brief A base constraint functor that resets the position of the nodes that
 * violate the constraint
 *
 */
class Constraint {
 public:
  virtual void operator()(GraphType& g, double t) {
    (void)g;
    (void)t;
  }
  virtual Constraint* clone() { return new Constraint; }
  virtual ~Constraint(){}
};

/**
 * @brief A constraint that fix points with initial position (0,0,0) or (1,0,0)
 * at their initial position
 *
 */
class PinConstraint : public Constraint {
 public:
  virtual void operator()(GraphType& g, double t) {
    (void) t;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      Point init_position = (*it).value().init_position;
      if (init_position == Point(0, 0, 0)) {
        (*it).position() = Point(0, 0, 0);
        continue;
      }

      if (init_position == Point(1, 0, 0)) {
        (*it).position() = Point(1, 0, 0);
      }
    }
  }

  virtual PinConstraint* clone() { return new PinConstraint; }
};

/**
 * @brief A constraint that make sure all nodes are above z=-0.75, and set
 * velocity in the z direction to 0 if it's violated
 *
 */
class PlaneConstraint : public Constraint {
 public:
  virtual void operator()(GraphType& g, double t) {
    (void)t;

    const double lower_bound = -0.75;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      if ((*it).position()[2] < lower_bound) {
        (*it).position()[2] = lower_bound;
        (*it).value().vel[2] = 0;
      }
    }
  }

  virtual PlaneConstraint* clone() { return new PlaneConstraint; }
};

/**
 * @brief A constraint that makes sure the nodes are outside of a sphere
 *
 */
class SphereConstraint : public Constraint {

 public:
  virtual void operator()(GraphType& g, double t) {
    (void)t;

    const Point c(0.5, 0.5, -0.5);
    const double r = 0.15;
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
      double dist = norm(c - (*it).position());
      if (dist < r) {
        Point R = ((*it).position() - c) / dist;
        (*it).position() = c + R * r;
        (*it).value().vel = (*it).value().vel - dot((*it).value().vel, R) * R;
      }
    }
  }

  virtual SphereConstraint* clone() { return new SphereConstraint; }
};

/**
 * @brief A constraint that removes the nodes in the sphere
 *
 */
class SphereConstraintRemove : public Constraint {

 public:
  virtual void operator()(GraphType& g, double t) {
    (void)t;

    const Point c(0.5, 0.5, -0.5);
    const double r = 0.15;
    auto it = g.node_begin();
    while (it!= g.node_end()) {
      double dist = norm(c - (*it).position());
      if (dist < r) {
        it = g.remove_node(it);
      } else {
        ++it;
      }
    }
  }

  virtual SphereConstraintRemove* clone() { return new SphereConstraintRemove; }
};

/**
 * @brief A combination of three constraints. Object of this class should not be
 * copied.
 *
 */
class CombinedConstraints {
  std::vector<Constraint*> constraint_vec;

 public:
  CombinedConstraints() {}

  CombinedConstraints(const CombinedConstraints& other) {
    for (unsigned i = 0; i < other.constraint_vec.size(); ++i) {
      constraint_vec.push_back(other.constraint_vec[i]->clone());
    }
  }

  ~CombinedConstraints() {
    for (unsigned i = 0; i < constraint_vec.size(); ++i) {
      delete constraint_vec[i];
    }
  }

  void operator()(GraphType& g, double t) {
    for (unsigned i = 0; i < constraint_vec.size(); ++i) {
      constraint_vec[i]->operator()(g, t);
    }
  }

  template<typename C>
  void add_constraint(C c){
    (void)c;
    constraint_vec.push_back(new C);
  }
};

/**
 * @brief A wrapper that create a combined constraint object
 *
 * @tparam C1 First constraint type, a subclass of Constraint
 * @tparam C2 Second constraint type, a subclass of Constraint
 * @tparam C3 Third constraint type, a subclass of Constraint
 * @param c1 First constraint
 * @param c2 Second constraint
 * @param c3 Third constraint
 * @return CombinedConstraints, constraint corresponding to apply these
 * constraints sequentially
 */
template <typename C1, typename C2, typename C3 = Constraint>
CombinedConstraints make_combined_constraints(C1 c1, C2 c2,
                                              C3 c3 = Constraint()) {
  CombinedConstraints result;
  result.add_constraint(c1);
  result.add_constraint(c2);
  result.add_constraint(c3);
  return result;
}

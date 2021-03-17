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
  bool pinnedleft;
  bool pinnedright;
  NodeData() : vel(0), mass(1), pinnedleft(false), pinnedright(false){}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double L;       //< Edge rest length
    double K;       //< Edge spring force
    EdgeData() : L(0), K(100) {}
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
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }

  // Apply the constraints
  constraint(g,t);

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
    (void) t;
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }
    Point gravForce = Point(0,0,n.value().mass * grav * -1);
    Point springForce = Point(0,0,0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
        NODE other = (*it).node2();
        Point nodeForce = Point(0,0,0);
        double dist = norm(n.position() - other.position());
        nodeForce += -1 * (*it).value().K * ((n.position() - other.position()) / dist) *
                     (dist - (*it).value().L);
        springForce += nodeForce;
    }
    return gravForce + springForce;
  }
};

/** Generalized Force class to enable polymorphism. */
class Force {
    public:
    /** Return the force applying to @a n at time @a t. For the base class, this
     * always returns a null force (0,0,0). */
    virtual Point operator()(Node n, double t) {
        (void) n;
        (void) t;
        return Point(0,0,0);
    }
};

/** Represents the force of gravity on a node.*/
class GravityForce : public Force {
    public:
    /** Return the force applying to @a n at time @a t. In this case, -9.81 m/s in
     * the negative z direction. */
    Point operator()(Node n, double t) {
        (void) t;
        return Point(0,0,n.value().mass * grav * -1);
    }
};

/** Represents the force of all incident edges on a node, acting as springs.*/
class MassSpringForce : public Force {
    public:
    /** Return the force applying to @a n at time @a t. In this case, calculated
     * by applying the spring equations to all of the incident edges and summing
     * their impacts. */
    Point operator()(Node n, double t) {
        (void) t;
        Point springForce = Point(0,0,0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Node other = (*it).node2();
            Point nodeForce = Point(0,0,0);
            double dist = norm(n.position() - other.position());
            nodeForce += -1 * (*it).value().K * ((n.position() - other.position()) / dist) *
                         (dist - (*it).value().L);
            springForce += nodeForce;
        }
        return springForce;
    }
};

/** Represents a frictional damping force on a node.*/
class DampingForce : public Force {
    public:
    /** Return the force applying to @a n at time @a t. In this case, some fraction
     * of the node's current velocity c, applied in the opposite direction of travel. */
    Point operator()(Node n, double t) {
        (void) t;
        double c = .1;
        return -1 * c * n.value().vel;
    }
};

/** Struct for handling a combination of forces, all applied at once.*/
struct CombinedForce {
    std::vector<Force*> forces;
    /** Constructor. Takes in a vector of Force pointers to apply when called.*/
    CombinedForce(std::vector<Force*> f) : forces(f) {};

    /** Return the force applying to @a n at time @a t. In this case, combining
     * the effects of every force passed in to the functor. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        Point sumForces = Point();
        for(unsigned i = 0; i < forces.size(); i++) {
            sumForces += (*(forces[i]))(n,t);
        }
        return sumForces;
    }
};

/** Function to combine a set of forces, returning a CombinedForce object
 * which can be called to apply the effects of every force at once. */
template <typename... Forces>
CombinedForce make_combined_force(Forces... forces) {
     std::vector<Force*> fs = {&forces...};
     return CombinedForce(fs);
 }

/** Generalized Constraint class to enable polymorphism. */
class Constraint {
    public:
    /** Reset the position and velocity of all nodes of @a g violating constraint
     * @a t. In this case, does nothing. */
    virtual void operator()(GraphType& g, double t) {
        (void) t;
        (void) g;
    }
};

/** Represents the constraint of points being pinned. */
class PinConstraint : public Constraint {
    public:
    /** Reset the position and velocity of all nodes of @a g violating constraint
     * @a t. In this case, keeps points listed as pinned in place. */
    void operator()(GraphType& g, double t) {
        (void) t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            if ((*it).value().pinnedleft) {
                (*it).position() = Point(0,0,0);
                (*it).value().vel = Point(0,0,0);
            } else if ((*it).value().pinnedright) {
                (*it).position() = Point(1,0,0);
                (*it).value().vel = Point(0,0,0);
            }
        }
    }
};

/** Represents the constraint preventing points from moving through a plane. */
class PlaneConstraint : public Constraint {
    public:
    /** Reset the position and velocity of all nodes of @a g violating constraint
     * @a t. In this case, keeps the points from falling below the plane at z = -0.75. */
    void operator()(GraphType& g, double t) {
        (void) t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            if ((*it).position()[2] <= -.75) {
                (*it).position()[2] = -.75;
                (*it).value().vel[2] = 0;
            }
        }
    }
};

/** Represents the constraint preventing points from moving through a sphere. */
class SphereConstraint : public Constraint {
    public:
    /** Reset the position and velocity of all nodes of @a g violating constraint
     * @a t. In this case, keeps the points from falling into the sphere at (0.5,0.5,-0.5)
     * with radius 0.15. */
    void operator()(GraphType& g, double t) {
        (void) t;
        Point center = Point(0.5,0.5,-0.5);
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            if (norm((*it).position()-center) <= 0.15) {
                Point ri = ((*it).position() - center)/norm((*it).position() - center);
                (*it).position() = center + 0.15 * ri;
                (*it).value().vel -= dot((*it).value().vel,ri)*ri;
            }
        }
    }
};

/** Represents the constraint deleting points that move through a sphere. */
class TearConstraint : public Constraint {
    public:
    /** Delete all nodes of @a g violating constraint @a t. In this case, deletes
     * the points falling into the sphere at (0.5,0.5,-0.5) with radius 0.15. */
    void operator()(GraphType& g, double t) {
        (void) t;
        Point center = Point(0.5,0.5,-0.5);
        auto it = g.node_begin();
        while (it != g.node_end()) {
            if (norm((*it).position()-center) <= 0.15) {
                g.remove_node(it);
                it = g.node_begin();
            } else {
                ++it;
            }
        }
    }
};

/** Struct for handling a combination of constraints, all applied in order.*/
struct CombinedConstraint {
    std::vector<Constraint*> consts;
    /** Constructor. Takes in a vector of Constraint pointers to apply when called.*/
    CombinedConstraint(std::vector<Constraint*> c) : consts(c) {};

    /** Reset the position and velocity of all nodes of @a g violating constraint
     * @a t. In this case, applying each constraint in the list. */
    void operator()(GraphType& g, double t) {
        for(unsigned i = 0; i < consts.size(); i++) {
            (*(consts[i]))(g,t);
        }
    }
};

/** Function to combine a set of constraints, returning a CombinedConstraint object
 * which can be called to apply the effects of every constraint in order. */
template <typename... Constraints>
CombinedConstraint make_combined_constraint(Constraints... consts) {
     std::vector<Constraint*> cs = {&consts...};
     return CombinedConstraint(cs);
 }

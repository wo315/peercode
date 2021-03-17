/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <chrono>
#include <thread>
#include <algorithm>
#include <memory>
#include <type_traits>

#include "CME212/Util.hpp"
#include "CME212/Color.hpp"
#include "CME212/Point.hpp"

#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;
// For part 3.2
static constexpr double L = 0.010101;
static constexpr double K=  100;

/** Custom structure of data to store with Nodes */
struct NodeData {
  Point vel;       //< Node velocity
  Point lastPnt;
  double mass;     //< Node mass
  double c;         // Node damping constant
  NodeData() : vel(0), mass(1), c(0.1) {}
  NodeData(Point vel, double mass, double c)
      : vel(vel), mass(mass), c(c) {}
};
/** Custom structure of data to store with Nodes */
struct EdgeData {
   double K; // spring constant
   double L; // natural length
  EdgeData() : K(100), L(0.1) {}
  EdgeData(double K, double L)
      : K(K), L(L) {}
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

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    // Code added to skip the update if the node has a particular position
    if (n.position()==Point(0,0,0) || n.position()==Point(1,0,0))
        (void)dt;
    else
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

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint Constraint object resetting nodes after move action.
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports velocity, mass, node last position
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
    // Save previous position (useful for constraints)
    n.value().lastPnt = n.position();
    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // apply constraints
  constraint(&g, t);
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
    Problem1Force() {}
    /* This constructor was used for initial versions only.
      Problem1Force(double K, double L) {
          L_=L;
          K_=K;
      }
      */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    //
        (void)t;
      if (n.position()==Point (0,0,0) || n.position()==Point(1,0,0))
          return Point(0,0,0);

    // compute spring force by adding
    Point f {0,0,0};
    for (auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
        Edge e = *it;
        double K = e.value().K;
        double L = e.value().L;
        Point diff = e.node2().position() - e.node1().position();
        double dis = std::max(0.0, e.length() - L);
        f += K *dis* diff;
    }
    // compute the gravitational force
    f += n.value().mass * Point(0,0, -grav);

    return f;
  }
  private:
  /* Used for initial version only
    double L_;
    double K_;
    */
};


/** Base class for all forces functors
 * Returns zero force when called, regardless of input.
 * O(1)
 */
class Force {
    public:
        /** Virtual class for all force objects
        */
        Force() {}
        virtual ~Force() {}
        virtual Point operator()(Node n, double t) {
            // std::cout<<"Base"<<std::endl;
            (void)t;
            (void)n;
            return Point(0.0);
            }
};

/** Force functor returning gravity force, regardless of input
 * O(1)
 */
class GravityForce : public Force {
    public:
    GravityForce() {}
    virtual Point operator()(Node n, double t) {
        (void)t;
       return n.value().mass * Point(0,0, -grav);
    }
};

/** Mass Spring functor, returning the net spring force on Node n passed to
 * operator().
 *
 * @pre: For edges e, incident to Node n passed to operator(), there is a
 * spring constant in e.value().K, a slack length at e.value().L.
 *
 * O(n.degree())
 */
class MassSpringForce : public Force {
    public:
    MassSpringForce() {}
    Point operator()(Node n, double t) {
        (void)t;
        Point f {0,0,0};
        //iterate over incident edges, adding spring forces together.
        for (auto it=n.edge_begin(); it!=n.edge_end(); ++it) {
            Edge e = *it;
            double K = e.value().K;
            double L = e.value().L;
            Point diff = e.node2().position() - e.node1().position();
            double dist = norm(diff);
            f += K*diff/dist*(dist-L);
        }
    return f;
    }
};

/** Force functor implementing damping force on Node n passed to operator(),
 * @pre: the damping coefficient is stored in n.value().c.
 * O(1)
 */
class DampingForce : public Force {
    public:
    DampingForce() {}
    Point operator()(Node n, double t) {
        (void)t;
        return n.value().c*n.value().vel;
    }
};

/** Functor for producing a force from a vector of forces
 * @pre: each force passed in vector forces is a pointer to a valid force opject
 * @post: the forces object are not modified by this functor.
 * O(forces.size())
 */
class CombinedForce {
    public:
    std::vector<Force*> forces;
    CombinedForce(std::vector<Force*> forces) : forces(forces)  {}
    Point operator()(Node n, double t) {
         Point f {0,0,0};
         for (auto& force  : forces) {
             f += (*force)(n, t);
         }
         return f;
    }
};

/**
 * Function combining forces. Note that the third argument is optional and
 * has a default value. The default Force() in the base class has a 0 force
 * so we can always add it without changing anything.
 * O(1), since the number of inputs is fixed.
 * @pre f1,f2,f3 are references to valid force objects.
 */
CombinedForce make_combined_force(Force&& f1, Force&& f2, Force&& f3=Force()) {
    std::vector<Force*> forces {&f1, &f2, &f3};
    return CombinedForce(forces);
}

/** Base class for all constraints. Does not modify anything.
 * O(1)
 */
class Constraint {
    public:
        Constraint(){}
        virtual void operator()(GraphType* g, double t) {(void)g; (void)t;}
};

/** Uses the Node value n.value().lastPnt, which must be set in symp_euler.
 * Fixes nodes having initial position (1,0,0) and (0,0,0)
 *
 * O(g.num_nodes())
 */
class PinConstraint : public Constraint {
    public:
        // Take a graph and vector of points that must be pinned. The points
        // to be pinned are determined before simulatino begins, so it's based
        // on initial point
        PinConstraint() {}
        // if pinConstraint flag is true, reset the value
        void operator()(GraphType* g, double t) {
            (void)t;
            for (auto it=g->node_begin(); it!=g->node_end(); ++it) {
                if ((*it).value().lastPnt==Point(0,0,0))
                    (*it).position()=Point(0,0,0);
                else if ((*it).value().lastPnt==Point(1,0,0))
                    (*it).position()=Point(1,0,0);
                else
                    ;
            }
        }
};

/* Constraint functor. If nodes move below the plane z=-0.75, then that node
 * is projected to the nearest point on the plane.
 *
 * O(g.num_nodes())
 */
class PlaneConstraint : public Constraint {
    public:
        double z_plane=-0.75;
        PlaneConstraint() {}
        void operator()(GraphType* g, double t) {
            (void)t;
            for (auto it=g->node_begin(); it!=g->node_end(); ++it) {
               if ((*it).position().z < z_plane) {
                   (*it).position().z = z_plane;
                   (*it).value().vel.z=0;
               }
            }
        }
};

/** Constraint functor. If nodes move inside sphere centred at (0.5,.5,-.5) with
 * radius 0.15, then project nodes back to the nearest point on the sphere
 * surface.
 *
 * O(g.num_nodes())
 */
class SphereConstraint : public Constraint {
    public:
        Point c=Point(0.5,0.5,-0.5);
        double r=0.15;
        SphereConstraint() {}
        void operator()(GraphType* g, double t) {
            (void)t;
            for (auto it=g->node_begin(); it!=g->node_end(); ++it) {
                Point pos = (*it).position(); // current position
                double dist_to_c = norm( pos-c ); // distance to centre
                if (dist_to_c <  r) {
                    // set new position in direction of `normal` (normalize it) and
                    (*it).position() += (pos-c)/norm(pos-c) * (r-dist_to_c);
                    // correct the velocity inplace
                    Point R_i = (pos-c) / norm(pos-c);
                    Point& v_i = (*it).value().vel;
                    v_i = v_i - dot(v_i,R_i)*R_i;
                }
            }
        }
};

/** Constraint functor. If nodes move inside sphere centred at (0.5,.5,-.5) with
 * radius 0.15, then remove the node from the graph.
 *
 * O(g.num_nodes())
 */
class TearConstraint : public Constraint {
    public:
        Point c=Point(0.5,0.5,-0.5);
        double r=0.15;
        TearConstraint() {}
        void operator()(GraphType* g, double t) {
            (void)t;
            for (auto it=g->node_begin(); it!=g->node_end(); ++it) {
                Point pos = (*it).position(); // current position
                double dist_to_c = norm( pos-c ); // distance to centre
                if (dist_to_c <  r) {
                    it=g->remove_node(it);
                }
            }
        }
};


/** Functor for producing a Constraint from a vector of forces
 * @pre the contraints vector are pointers to valid constraints
 *
 * O(constraints.size())
 */
class CombinedConstraint {
    public:
    std::vector<Constraint*> constraints;
    CombinedConstraint(std::vector<Constraint*> constraints) : constraints(constraints)  {}
    void operator()(GraphType* g, double t) {
         for (auto& constraint: constraints) {
             (*constraint)(g, t);
         }
    }
};

/** Function combining constraints. Note that the third argument is optional and
 * has  default value. The default constraint() in the base class has a 0 force
 * so we can always add it without changing anything.
 *
 * O(1)
 */
CombinedConstraint make_combined_constraints(Constraint&& c1, Constraint&& c2, Constraint&& c3=Constraint()) {
//--functionality_2
//--This causes a segfault
//--Indeed you are creating poiting to r-value references, but those r-value references fall out of range and get deleted
//--the pointers are pointing at unitialized memory
//--START
    std::vector<Constraint*> constraints{&c1, &c2, &c3};
//--END
    return CombinedConstraint(constraints);
}

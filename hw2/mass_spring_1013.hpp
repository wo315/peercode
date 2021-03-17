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
struct NodeData {
  Point vel;       //< Node velocity
  double mass;     //< Node mass
  unsigned pin;       //< 0 if PinConstraint does not fix the node
                      //< 1 if PinConstraint fixes Node position at (0,0,0)
                      //< 2 if PinConstraint fixes Node position at (1,0,0)
  NodeData() : vel(0), mass(1), pin(0) {};
  NodeData(Point v, double m, unsigned p) : vel(v), mass(m), pin(p) {};

public:
  void update_pin(bool p) {
      pin = p;
  }
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double k;            //< Spring constant
    double rest_len;     //< Rest length
    EdgeData() : k(100), rest_len(0.1) {};
    EdgeData(double k, double rl) : k(k), rest_len(rl) {};
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
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam  C is a function opject called as @a constraint(g, @a t),
 *            where g is a graph reference and @a t is the current time.
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

    // Apply constraints
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
   * @pre Gravity constant _grav_ is defined
   * @pre Node mass n.value().mass is defined
   * @pre For every Edge e of Node n, spring constant e.value().k is defined,
   *      rest length e.value().rest_len is defined, and e.length() > 0
   *
   * For HW2 #1, this is a combination of mass-spring force and gravity,
   * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
   * model that by returning a zero-valued force. */
  template <typename NODE>
  Point operator()(NODE n, double t) {
    // HW2 #1: YOUR CODE HERE
    (void) t;    // silence compiler warnings
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }

    Point Fg = Point(0, 0, -grav * n.value().mass);
    Point Fs = Point(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
        Edge e = *it;
        Point dF = - e.value().k * (e.node1().position() - e.node2().position())
                     * (e.length() - e.value().rest_len) / e.length();
        if (n == e.node1()) {
            Fs += dF;
        } else {
            Fs -= dF;
        }
    }
    return Fg + Fs;
  }
};

class Force {
public:
    /** Virtual operator() */
    virtual Point operator()(Node n, double t) {
        (void) n, (void) t;
        return Point(0,0,0);
    }
    /** Virtual deconstructor */
    virtual ~Force(){};    // prevent memory leak
};

class GravityForce : public Force {
public:
    /** Return the gravity force on a node given the current time.
     *
     * @pre Gravity constant _grav_ is defined
     * @pre Node mass n.value().mass is defined
     * @param n Node to compute the force
     * @param t Current time
     * @return point of gravity force on _n_ at _t_
     */
    Point operator()(Node n, double t) {
        (void) t;    // silence compiler warnings
        return Point(0, 0, -grav * n.value().mass);
    }
};

class MassSpringForce : public Force {
public:
    /** Return the mass spring force on a node given the current time.
     *
     * @pre For every Edge e of Node n, spring constant e.value().k is defined,
     *      rest length e.value().rest_len is defined, and e.length() > 0
     * @param n Node to compute the force
     * @param t Current time
     * @return point of mass spring force on _n_ at _t_
     */
    Point operator()(Node n, double t) {
        (void) t;    // silence compiler warnings
        Point Fs = Point(0, 0, 0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
            Edge e = *it;
            assert(e.length() != 0);
            Point dF = - e.value().k * (e.node1().position() - e.node2().position())
                       * (e.length() - e.value().rest_len) / e.length();
            if (n == e.node1()) {
                Fs += dF;
            } else {
                Fs -= dF;
            }
        }
        return Fs;
    }
};

class DampingForce : public Force {
public:
    /** Construct a DampingForce object.
     *  By default, the damping constant is 0.1.
     */
    DampingForce(): damp_const(0.1) {};

    /** Construct a DampingForce object given a damping constant. */
    DampingForce(double dc): damp_const(dc) {};

    /** Return the damping force (friction) on a node given the current time.
     *
     * @param n Node to compute the force
     * @param t Current time
     * @return point of damping force on _n_ at _t_
     */
    Point operator()(Node n, double t) {
        (void) t;    // silence compiler warnings
        return -damp_const * n.value().vel;
    }

private:
    double damp_const;
};


class CombinedForce : public Force {
public:
    /** Make a CombinedForce object from a vector of Force pointers.
     *
     * @param f Vector of Force pointers to be combined
     */
    CombinedForce(std::vector<Force*> f) : forces(f) {};

    /** Return the combined force to a node given the current time.
     *
     * @param n NODE to compute the force
     * @param t Current time
     * @return point of combined force on _n_ at _t_
     */
    template <typename NODE>
    Point operator() (NODE n, double t) const {
        Point sum_f = Point(0,0,0);
        for (Force* f : forces) sum_f += (*f)(n, t);
        return sum_f;
    }
private:
    std::vector<Force*> forces;
};

/** Return a CombinedForce from two given forces
 *
 * @tparam F1 Class of the first force
 * @tparam F2 Class of the second force
 * @param f1 First force
 * @param f2 Second force
 * @return The combined force
 */
template <typename F1, typename F2>
CombinedForce make_combined_force(F1 f1, F2 f2) {
    std::vector<Force*> f {};
    f.push_back(&f1);
    f.push_back(&f2);
    return CombinedForce(f);
}

/** Return a CombinedForce from three given forces
 *
 * @tparam F1 Class of the first force
 * @tparam F2 Class of the second force
 * @tparam F3 Class of the third force
 * @param f1 First force
 * @param f2 Second force
 * @param f3 Third force
 * @return The combined force
 */
template <typename F1, typename F2, typename F3>
CombinedForce make_combined_force(F1 f1, F2 f2, F3 f3) {
    std::vector<Force*> f {};
    f.push_back(&f1);
    f.push_back(&f2);
    f.push_back(&f3);
    return CombinedForce(f);
}


/**
 *
 */
class Constraint {
public:
    /** Virtual operator() */
    virtual void operator()(GraphType& g, double t) {
        (void) g, (void) t;
    }

    /** Virtual deconstructor */
    virtual ~Constraint(){};    // prevents memory leak
};

class PinConstraint : public Constraint {
public:
    /** Construct a PinConstraint object.
     *  By default, fix points (0,0,0) and (1,0,0).
     */
    PinConstraint() : p1(Point(0,0,0)), p2(Point(1,0,0)) {};

    /** Apply the pin constraint to a graph given the current time.
     *
     * @param g Graph (by reference) to apply the constraint
     * @param t Current time
     * @pre For each active Node n of g, if n.position() == p1 or
     *      n.position() == p2, then n.value().pin == true
     * @post For each active Node n of g, if n.value().pin == true, then
     *       n.value().vel == Point(0,0,0)
     */
    void operator()(GraphType& g, double t) {
        (void) t;    // silence compiler warnings
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (n.value().pin == 1) {
                n.position() = Point(0,0,0);
            } else if (n.value().pin == 2) {
                n.position() = Point(1,0,0);
            }

        }
    }

private:
    Point p1;
    Point p2;
};

class PlaneConstraint : public Constraint {
public:
    /** Construct a PlaneConstraint object.
     *  By default, the plane is z=-0.75.
     */
    PlaneConstraint() : z(-0.75) {};

    /** Apply the plane constraint to a graph given the current time.
     *
     * If any Node of g is below the plane z=this->z, set its position
     * to nearest point on the plane and set its z-velocity to zero
     *
     * @param g Graph (by reference) to apply the constraint
     * @param t Current time
     */
    void operator()(GraphType& g, double t) {
        (void) t;    // silence compiler warnings
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (inner_prod( n.position(), Point(0,0,1)) < z) {
                n.position()[2] = z;
                n.value().vel[2] = 0;
            }
        }
    }

private:
    double z;
};

class SphereConstraint : public Constraint {
public:
    /** Construct a SphereConstraint object.
     *  By default, the center is (0.5, 0.5, 0.5) and radius is 0.15.
     */
    SphereConstraint() : c(Point(0.5, 0.5, -0.5)), r(0.15) {};

    /** Apply the sphere constraint to a graph given the current time.
     *
     * If any Node of g is inside the sphere with center=this->z, radius=this->r,
     * set its position to nearest point on the sphere's surface and set its
     * velocity component normal to the sphere's surface to zero
     *
     * @param g Graph (by reference) to apply the constraint
     * @param t Current time
     */
    void operator()(GraphType& g, double t) {
        (void) t;    // silence compiler warnings
        for (auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (norm( n.position()-c) < r) {
                Point R = (n.position() - c) / norm(n.position() - c);
                n.position() = c + r * R;
                n.value().vel -= inner_prod(n.value().vel, R) * R;
            }
        }
    }

private:
    Point c;
    double r;
};

class SphereRemoveConstraint : public Constraint {
public:
    /** Construct a SphereRemoveConstraint object.
     *  By default, the center is (0.5, 0.5, 0.5) and radius is 0.15.
     */
    SphereRemoveConstraint() : c(Point(0.5, 0.5, -0.5)), r(0.15) {};

    /** Apply the sphere constraint to a graph given the current time.
     *
     * If any Node of g is inside the sphere with center=this->z, radius=this->r,
     * remove the Node
     *
     * @param g Graph (by reference) to apply the constraint
     * @param t Current time
     */
    void operator()(GraphType& g, double t) {
        (void) t;    // silence compiler warnings
        for (auto it = g.node_begin(); it != g.node_end();){
            auto n = *it;
            if (norm( n.position()-c) < r) {
                it = g.remove_node(it);
            } else {
                ++it;
            }
        }
    }

private:
    Point c;
    double r;
};



class CombinedConstraint : public Constraint {
public:
    /** Construct a CombinedConstraint object from a vector of Contraint pointers.
     *
     * @param c Vector of Constraint pointers to be combined
     */
    CombinedConstraint(std::vector<Constraint*> c) : constraints(c) {};

    /** Apply the combined constraint to a graph given the current time.
     *
     * @param g Graph (by reference) to apply the constraint
     * @param t Current time
     */
    void operator() (GraphType& g, double t) const {
        for (Constraint* c : constraints) (*c)(g, t);
    }
private:
    std::vector<Constraint*> constraints;
};

/** Return a CombinedConstraint from two given constraints
 *
 * @tparam C1 Class of the first constraint
 * @tparam C2 Class of the second constraint
 * @param c1 First constraint
 * @param c2 Second constraint
 * @return The combined constraint
 */
template <typename C1, typename C2>
CombinedConstraint make_combined_constraint(C1 c1, C2 c2) {
    std::vector<Constraint*> c {};
    c.push_back(&c1);
    c.push_back(&c2);
    return CombinedConstraint(c);
}

/** Return a CombinedConstraint from three given constraints
 *
 * @tparam C1 Class of the first constraint
 * @tparam C2 Class of the second constraint
 * @tparam C3 Class of the third constraint
 * @param c1 First constraint
 * @param c2 Second constraint
 * @param c3 Third constraint
 * @return The combined constraint
 */
template <typename C1, typename C2, typename C3>
CombinedConstraint make_combined_constraint(C1 c1, C2 c2, C3 c3) {
    std::vector<Constraint*> c {};
    c.push_back(&c1);
    c.push_back(&c2);
    c.push_back(&c3);
    return CombinedConstraint(c);
}

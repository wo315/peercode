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
  Point initial;
  NodeData() : vel(0), mass(1), initial(Point(1, 1, 1)) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData{
    double K;
    double length;
    EdgeData() : K(100.0), length(1.0) {}
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
 * @param[in]     constraint  Function object defining the constraint per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports ???????? YOU CHOOSE
 * @tparam F is a function object called as @a force(n, @a t),
 *           where n is a node of the graph and @a t is the current time.
 *           @a force must return a Point representing the force vector on
 *           Node n at time @a t.
 * @tparam C is a function object called as @a constraint(&g, t) where g is
 *           the reference to the graph and @a t is the current time.
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

    // Apply the constraint
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
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        return Point(0,0,0);
    }
    Point F_spring = Point(0, 0, 0);
    Point F_grav = n.value().mass * Point(0, 0, -grav);
    for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
        auto e = *it;
        Point e_vec;
        if(e.node1() == n){
            e_vec = n.position() - e.node2().position();
        }
        else{
            e_vec = n.position() - e.node1().position();
        }
        F_spring -= e.value().K * e_vec / norm(e_vec) * (norm(e_vec) - e.value().length);
      }
    (void) t;
    return F_spring + F_grav;
  }
};

/** The base Force class that has the virtual operator().*/
class Force {
public:
    virtual Point operator() () {
        return Point(0, 0, 0);
    }
};

/** The GravityForce class that inherits the Force class and has implementation of the operator() to return
 * a point representing the displacement by the gravity force. */
class GravityForce: public Force {
public:
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return n.value().mass * Point(0, 0, -grav);
    }
};

/** The MassSpringForce class that inherits the Force class and has implementation of the operator() to return
 * a point representing the displacement by the mass spring force. */
class MassSpringForce: public Force {
public:
    template <typename NODE>
    Point operator()(NODE n, double t) {
        Point F_spring = Point(0, 0, 0);
        for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
            auto e = *it;
            Point e_vec;
            if(e.node1() == n){
                e_vec = n.position() - e.node2().position();
            }
            else{
                e_vec = n.position() - e.node1().position();
            }
            F_spring -= e.value().K * e_vec / norm(e_vec) * (norm(e_vec) - e.value().length);
        }
        (void) t;
        return F_spring;
    }
};

/** The DampingForce class that inherits the Force class and has implementation of the operator() to return
 * a point representing the displacement by the damping force. */
class DampingForce: public Force {
public:
    DampingForce(): c_(1.0) {};
    DampingForce(double c) : c_(c) {};
    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        return -c_ * n.value().vel;
    }
private:
    double c_;
};

/** A functor that applies the operator() to two forces and returns the combined displacement resulted
 * by the two forces. */
//--design_0
//--Combined Force should store a vector of forces to avoid creating one CombinedForce functor for each number of forces
//--START
template<typename Force1, typename Force2>
class CombinedForce {
public:
    Force1 f1;
    Force2 f2;
    template <typename NODE>
    Point operator()(NODE n, double t) {
        return f1(n, t) + f2(n, t);
    }
};
//--END

/** Make the combination of the two forces.
 *
 * @tparam Force1 The first type of force
 * @tparam Force2 The second type of force
 * @param f1 The first force
 * @param f2 The second force
 * @return The functor CombinedForce
 */
template<typename Force1, typename Force2>
CombinedForce<Force1, Force2> make_combined_force(Force1 f1, Force2 f2){
    return {f1, f2};
}

/** Make the combination of the three forces.
 *
 * @tparam Force1 The first type of force
 * @tparam Force2 The second type of force
 * @tparam Force3 The third type of force
 * @param f1 The first force
 * @param f2 The second force
 * @param f3 The third force
 * @return The functor CombinedForce with the CombinedForce of f1 and f2 as the new f1 and f3 as the new f2
 */
template<typename Force1, typename Force2, typename Force3>
CombinedForce<CombinedForce<Force1, Force2>, Force3> make_combined_force(Force1 f1, Force2 f2, Force3 f3){
    return make_combined_force(make_combined_force(f1, f2), f3);
}

/** The base Constraint class that has the virtual operator().*/
class Constraint {
public:
    virtual void operator() () {}
};

/** The PinConstraint class that inherits the Constraint class and pin two nodes. */
class PinConstraint: public Constraint {
public:
    void operator() (GraphType& g, double t) {
        (void) t;
        for(auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (n.value().initial == Point(0, 0, 0) || n.value().initial == Point(1, 0, 0)){
                n.position() = n.value().initial;
                n.value().vel = Point(0, 0, 0);
            }
        }
    }
};

/** The PlaneConstraint class that inherits the Constraint class and impose plane constraint. */
class PlaneConstraint: public Constraint {
public:
    void operator() (GraphType& g, double t) {
        (void) t;
        double z_con = -0.75;
        for(auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            if (n.position().z < z_con) {
                n.position().z = z_con;
                n.value().vel.z = 0;
            }
        }
    }
};

/** The SphereConstraint class that inherits the Constraint class imposes sphere constraint. */
class SphereConstraint: public Constraint {
public:
    void operator() (GraphType& g, double t) {
        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;
        for(auto it = g.node_begin(); it != g.node_end(); ++it){
            auto n = *it;
            Point x = n.position();
            if (norm(x - c) < r) {
                Point R = (x - c) / norm(x - c);
                n.position() = c + r * R;
                n.value().vel -= (n.value().vel * R) * R;
            }
        }
        (void) t;
    }
};

/** The SphereConstraintRemove class that inherits the Constraint class removes nodes by constraint */
class SphereConstraintRemove: public Constraint {
public:
    void operator() (GraphType& g, double t){
        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;
        for(auto it = g.node_begin(); it != g.node_end(); ++it) {
            auto n = *it;
            if (norm(n.position() - c) < r)
                g.remove_node(n);
        }
        (void) t;
    }
};

/** A functor that applies the operator() to two constraints. */
template<typename Constraint1, typename Constraint2>
class CombinedConstraints {
public:
    Constraint1 c1;
    Constraint2 c2;
    void operator()(GraphType& g, double t) {
        c1(g, t);
        c2(g, t);
    }
};

/** Make the combination of the two constraints.
 *
 * @tparam Constraint1 The first type of constraint
 * @tparam Constraint2 The second type of constraint
 * @param c1 The first constraint
 * @param c2 The second constraint
 * @return The functor CombinedConstraints
 */
template<typename Constraint1, typename Constraint2>
CombinedConstraints<Constraint1, Constraint2> make_combined_constraints(Constraint1 c1, Constraint2 c2){
    return {c1, c2};
}

/** Make the combination of the three constraints.
 *
 * @tparam Constraint1 The first type of constraint
 * @tparam Constraint2 The second type of constraint
 * @tparam Constraint3 The third type of constraint
 * @param c1 The first constraint
 * @param c2 The second constraint
 * @param c3 The third constraint
 * @return The functor CombinedConstraints with CombinedConstraints of c1 and c2 as the new c1 and c3 as the new c2
 */
template<typename Constraint1, typename Constraint2, typename Constraint3>
CombinedForce<CombinedConstraints<Constraint1, Constraint2>, Constraint3> make_combined_constraints(Constraint1 c1, Constraint2 c2, Constraint3 c3){
    return make_combined_constraints(make_combined_force(c1, c2), c3);
}



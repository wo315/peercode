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
// Spring constant
static constexpr double K = 100;
// Gravitational constant
static constexpr double g_constant = 9.81;
// Spring rest-length
static constexpr double L = 0;
// Damping constant
static constexpr double damp_constant = 0.01;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edge */
struct EdgeData {
    // spring rest-length
    double rest_length;
    // spring constant
    double spring_constant;
    EdgeData() : rest_length(L), spring_constant(K) {}
    void set_values(double L, double K) {
        rest_length = L;
        spring_constant = K;
    }
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

/** New function definition that add constraints arguments in symp_euler step above. */
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }
    //apply constraints
    constraint(g, t);
    
    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }
    return t + dt;
}



/* Force function object for HW2 #1.*/
struct Problem1Force {
    /** Return the force applying to @a n at time @a t.
     *
     * For HW2 #1, this is a combination of mass-spring force and gravity,
     * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
     * model that by returning a zero-valued force. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        // HW2 #1: YOUR CODE HERE
        (void) n; (void) t;    // silence compiler warnings
        
        // constain two corners of cloth (never move)
        if(n.position()==Point(0,0,0) || n.position() == Point(1,0,0)){
            return Point(0,0,0);
        }
        
        // Force due to gravity
        Point grav_f = n.value().mass * Point(0,0, - g_constant);
        // By 4.2 spring force formula
        Point spring_f;
        GraphType::IncidentIterator iter = n.edge_begin();
        while(!(iter ==  n.edge_end())){
            Edge e = *iter;
            Point dx = e.node1().position() - e.node2().position();
            spring_f -= e.value().spring_constant * dx *(norm(dx) - e.value().rest_length)/ norm(dx);
            ++iter;
        }
        // Add spring force and gravitational force together
        return spring_f + grav_f;
    }
};

/******************** Force ***********************************************/

/** @class Base class Force
 *  @param[in] node n where force is added, double t to represent time *
 *  @return Point(0,0,0), no force
 */
class Force{
public:
    virtual Point operator()(Node n, double t){
        (void) n; (void) t;    // silence compiler warnings
        return Point(0,0,0);
    }
};

/** @class GravityForce
 * @param[in] node n where force is added, double t to represent time
 * @return Gravitational force on node n with magnitude m*g, m = mass, g = 9.81
 */
class GravityForce : public Force{
public:
    Point operator()(Node n, double t){
        (void) t;    // silence compiler warnings
        Point grav_force = n.value().mass * Point(0,0, - g_constant);
        return grav_force;
    }
};

/** @class MassSpringForce
 * @param[in] node n where force is added, double t to represent time
 * @return MassSpring force on node n
 */
class MassSpringForce : public Force{
public:
    Point operator()(Node n, double t){
        (void) t;    // silence compiler warnings
        Point spring_force;
        GraphType::IncidentIterator iter = n.edge_begin();
        while(!(iter ==  n.edge_end())){
            Edge e = *iter;
            Point dx = e.node1().position() - e.node2().position();
            spring_force -= e.value().spring_constant * dx *(norm(dx) - e.value().rest_length)/ norm(dx);
            ++iter;
        }
        return spring_force;
    }
};


/** @class DampingForce
 * @param[in] node n where force is added, double t to represent time
 * @return Damping force, on node n, with magnitude = damping constant * node's velocity
 */
class DampingForce : public Force{
    Point operator()(Node n, double t){
        (void) t;    // silence compiler warnings
        Point damp_force = - damp_constant * n.value().vel;
        return damp_force;
    }
};

/******************** Combined Force Functor ***********************************************/
/** @class CombinedForce
 * @param[in] node n where the combined force is added, double t to represent time
 * @return The combined force on node n
 *
 * Combine all forces (gravitational force, spring mass force, damping force etc.) together.
 */
struct CombinedForce{
public:
    Point operator()(Node n, double t){
        Point combined_force;
        for(auto& f : forces) {
            combined_force = combined_force + (*f)(n, t);
        }
        return combined_force;
    }
    CombinedForce(std::vector<Force*>& f) : forces(f) {}
private:
    std::vector<Force*> forces;
};

/** Combine  f1, f2 together. */
template <typename F, typename E>
CombinedForce make_combined_force(F f1,  E f2) {
    std::vector<Force*> v {&f1, &f2};
    return CombinedForce(v);
}

/** Combine  f1, f2, f3 together. */
template <typename F, typename E, typename D>
CombinedForce make_combined_force(F f1,  E f2, D f3) {
    std::vector<Force*> v {&f1, &f2, &f3};
    return CombinedForce(v);
}



/******************** Constraints ***********************************************/

/** @class Base class Constraint
 *  @param[in] graph g, time t
 *
 *  No constraint in this functor.
 */
class Constraint {
public:
    virtual void operator()(GraphType& g, double t) {
        (void) g; // silence compiler warnings
        (void) t; // silence compiler warnings
    }
};

/** @class PinConstraint
 * @param[in] graph g, time t
 * @post Keep the nodes at (0, 0, 0) and (1, 0, 0) fixed
 */
class PinConstraint : public Constraint{
public:
    PinConstraint(unsigned idx1, unsigned idx2) : pin_idx1(idx1), pin_idx2(idx2) {}
    void operator()(GraphType& g, double t) {
        (void) t; // silence compiler warnings
        g.node(pin_idx1).position() = Point(0, 0, 0);
        g.node(pin_idx2).position() = Point(1, 0, 0);
    }
private:
    unsigned pin_idx1;
    unsigned pin_idx2;
};

/** @class PlaneConstraint
 * @param[in] graph g, time t
 * @post if violates, z-component of the Node velocity = 0
 */
class PlaneConstraint : public Constraint {
public:
    void operator()(GraphType& g, double t) {
        (void) t;
        double plane_constr = -0.75;
        for(auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
            //violation
            if(dot((*iter).position(), Point(0, 0, 1)) <  plane_constr) {
                // set the position to the nearest point on the plane
                (*iter).position().z  =  plane_constr;
                // set the z-component of Node velocity to 0
                (*iter).value().vel.z = 0;
            }
        }
    }
};

/** @class SphereConstraint
 * @param[in] graph g, time t
 * @post if violates, set position to nearest point and set velocity component
 */
class SphereConstraint : public Constraint {
public:
    void operator()(GraphType& g, double t) {
        (void) t;
        Point sphere_center = Point(0.5, 0.5, -0.5);
        double sphere_radius = 0.15;
        for(auto iter = g.node_begin(); iter != g.node_end(); ++iter) {
            Node n = *iter;
            double ctr = norm(n.position() - sphere_center);
            // if violates the constraint
            if(ctr < sphere_radius) {
                // set position to nearest point
                n.position() = sphere_center + (sphere_radius / ctr) * (n.position() - sphere_center);
                // set velcocity
                auto R = (n.position() - sphere_center) / ctr;
                n.value().vel -= (dot(n.value().vel, R)) * R;
            }
        }
    }
};

/** @class TearConstraint: remove nodes violating the sphere constraint
 * @param[in] graph g, time t
 * @post nodes and edges that inside the sphere are all removed
 */
class TearConstraint : public Constraint {
public:
    void operator()(GraphType& g, double t) {
        (void) t;
        Point sphere_center = Point(0.5, 0.5, -0.5);
        double sphere_radius = 0.15;
        for(auto it = g.node_begin(); it != g.node_end(); ) {
            double ctr = norm((*it).position() - sphere_center);
            if(ctr < sphere_radius) {
                it = g.remove_node(it);
            }
            else {
                ++it;
            }
        }
    }
};


/******************** Combined Constraints Functor ****************************************/
/** @class CombinedConstraints
 * @param[in] graph g, time t
 * @return The combined constaints on node n at time t
 *
 * Combined constraints together.
 */
struct CombinedConstraints {
    void operator()(GraphType& g, double t) {
        (void) t;
        //apply each constraint
        for(const auto & c : constraints) {
            (*c)(g, t);
        }
    }
    CombinedConstraints(std::vector<Constraint*>& v) : constraints(v) {}
private:
    std::vector<Constraint*> constraints;
};

/** Combine  c1, c2, c3 together. */
template <typename F, typename E, typename D>
CombinedConstraints make_combined_constraints(F c1,  E c2, D c3) {
    std::vector<Constraint*> v {&c1, &c2, &c3};
    return CombinedConstraints(v);
}

/** Combine  c1, c2 together. */
template <typename F, typename E>
CombinedConstraints make_combined_constraints(F c1,  E c2) {
    std::vector<Constraint*> v{&c1, &c2};
    return CombinedConstraints(v);
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



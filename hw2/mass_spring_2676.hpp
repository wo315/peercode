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
    Point init_pos;  //< Initial position
    NodeData() : vel(0), mass(1), init_pos(0) {}
};

// Define the Graph type
//--design_1
//--did not implement per-edge constant K
//--START
using GraphType = Graph<NodeData, double>;
//--END
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

    // Impose constraint
    constraint(g, dt);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}

/** Parent class for generalized forces. */
struct Force{
    // operator()
    virtual Point operator()(Node n, double t) {
        (void) n, (void) t;
        return Point(0);
    }
    // Used to destruct objects (potentially of child classes) pointed by Force*
    virtual ~Force() = default;
};

/** Functor to implement gravity force. */
struct GravityForce : public Force{
    // operator()
    virtual Point operator()(Node n, double t) {
        (void) t;
        return Point(0, 0, -n.value().mass * grav);
    }
};

/** Functor to implement mass string force. */
struct MassSpringForce : public Force{
    // spring constant
    double K;
    MassSpringForce() : K(100) {}
    // operator()
    virtual Point operator()(Node n, double t) {
        (void) t;
        Point force(0);
        for(GraphType::incident_iterator iter = n.edge_begin(); iter != n.edge_end(); ++iter){
            Edge e = *iter;
            force += ((-K) * (e.node1().position() - e.node2().position()) * (e.length() - e.value()) / e.length());
        }
        return force;
    }
};

/** Functor to implement damping force. */
struct DampingForce : public Force{
    // damping constant
    double c = 1e-4;
    // operator()
    virtual Point operator()(Node n, double t) {
        (void) t;
        return -c * n.value().vel;
    }
};

/** functor to combine different forces. */
struct CombinedForce{
    // pointers to force functors
    std::vector<Force*> forces;
    CombinedForce(std::vector<Force*>& f) : forces(f){}
    // operator()
    Point operator()(Node n, double t){
        Point res(0);
        for(const auto& d:forces)
            res += (*d)(n, t);
        return res;
    }
};

/** function to create CombinedForce functor. */
//--functionality_1
//--should take argument forces instead of creating new ones
//--START
template <typename T1 = Force, typename T2 = Force, typename T3 = Force>
CombinedForce make_combined_force(T1 f1 = Force(), T2 f2 = Force(), T3 f3 = Force()){
    (void) f1, (void) f2, (void) f3;
    std::vector<Force*> forces;
    forces.push_back(new T1());
    forces.push_back(new T2());
    forces.push_back(new T3());
    return CombinedForce(forces);
}
//--END

/** Base class for generalized constraints. */
struct Constraint{
    virtual void operator()(GraphType& graph, double t){
        (void) graph, (void) t;
    }
    // virtual destructor
    virtual ~Constraint() = default;
};

/** functor to implement pin constraint. */
struct PinConstraint : public Constraint{
    virtual void operator()(GraphType& graph, double t){
        (void) t;
        for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter){
            Node n = *iter;
            if(n.value().init_pos == Point(0, 0, 0)){
                n.value().vel = Point(0, 0, 0);
                n.position() = Point(0, 0, 0);
            }
            else if(n.value().init_pos == Point(1, 0, 0)){
                n.value().vel = Point(0, 0, 0);
                n.position() = Point(1, 0, 0);
            }
        }
    }
};

/** functor to implement plane constraint. */
struct PlaneConstraint : public Constraint{
    // plane constant
    double z = -0.75;
    // operator()
    virtual void operator()(GraphType& graph, double t){
        (void) t;
        for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter){
            Node n = *iter;
            if(dot(n.position(), Point(0, 0, 1)) < z){
                n.position().z = z;
                n.value().vel.z = 0;
            }
        }
    }
};

/** functor to implement sphere constraint. */
struct SphereConstraint : public Constraint{
    // center point
    Point c = Point(0.5, 0.5, -0.5);
    // radius
    double r = 0.15;
    // operator()
    virtual void operator()(GraphType& graph, double t){
        (void) t;
        for(auto iter = graph.node_begin(); iter != graph.node_end(); ++iter){
            Node n = *iter;
            Point vec = n.position() - c;
            double d = norm(vec);
            vec = vec / d;
            // too close to the center
            if(d < r){
                n.position() = c + r * vec;
                n.value().vel = n.value().vel - dot(n.value().vel, vec) * vec;
            }
        }
    }
};

/** functor to implement sphere constraint with removal. */
struct RemoveSphereConstraint : public Constraint{
    // center point
    Point c = Point(0.5, 0.5, -0.5);
    // radius
    double r = 0.15;
    // operator()
    virtual void operator()(GraphType& graph, double t){
        (void) t;
        auto iter = graph.node_begin();
        while(iter != graph.node_end()){
            Node n = *iter;
            Point vec = n.position() - c;
            double d = norm(vec);
            vec = vec / d;
            if(d < r){
                graph.remove_node(n);
            }
            else{
                ++iter;
            }
        }
    }
};

/** functor to combine constraints. */
struct CombinedConstraints{
    std::vector<Constraint*> constraints;
    CombinedConstraints(std::vector<Constraint*>& c) : constraints(c) {}
    void operator()(GraphType& graph, double t){
        for(const auto& d:constraints){
            (*d)(graph, t);
        }
    }
};

/** function to create CombinedConstraints functor. */
template<typename T1 = Constraint, typename T2 = Constraint, typename T3 = Constraint>
CombinedConstraints make_combined_constraints(T1 c1 = Constraint(), T2 c2 = Constraint(), T3 c3 = Constraint()){
    (void) c1, (void) c2, (void) c3;
    std::vector<Constraint*> constraints;
    constraints.push_back(new T1());
    constraints.push_back(new T2());
    constraints.push_back(new T3());
    return CombinedConstraints(constraints);
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
        (void) t;
        // no force for (0, 0, 0) and (1, 0, 0)
        if(n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0))
            return Point(0);
        // gravity force
        Point force(0, 0, -n.value().mass * grav);
        double K = 100;
        // mass spring force
        for(GraphType::incident_iterator iter = n.edge_begin(); iter != n.edge_end(); ++iter){
            Edge e = *iter;
            force += ((-K) * (e.node1().position() - e.node2().position()) * (e.length() - e.value()) / e.length());
        }
        return force;
    }
};

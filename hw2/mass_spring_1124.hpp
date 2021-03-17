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

//damping constant
static constexpr double c_damping = 0.1;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point vel;       //< Node velocity
    double mass;     //< Node mass
    Point prev_position;
    NodeData() : vel(0), mass(1) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;       //the spring constant
    double L;     // the rest length
    EdgeData() : K(100),  L(1) {}
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
//template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        n.value().prev_position = n.position();
        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    //apply the constraints
    constraint(&g, t);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }
    return t + dt;
}

/** Parent Force class*/
class Force {

 public:
    virtual Point operator()(Node n, double t){
        (void) t; (void) n;
        return Point(0);
    }
    virtual ~Force(){};
};

/** Functor that returns gravity force on node @a n */
class GravityForce : public Force{
 public:
    ~GravityForce(){};
    virtual Point operator()(Node n, double t){
        (void) t;
        return(Point(0, 0, -grav * n.value().mass));
    }
};

/** Functor that returns mass spring force on node @a n */
class MassSpringForce : public Force{
 public:
    ~MassSpringForce(){};
    virtual Point operator()(Node n, double t){
        (void) t;

        //init force
        Point spring_force{0,0,0};

        //sum over adjacent nodes
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
            Node n2 = (*it).node2();
            Point diff = n.position() - n2.position();

            spring_force += -(*it).value().K * (diff / (norm(diff))) * (norm(diff) - (*it).value().L);
        }
        return spring_force;
    }
};

/** Functor that returns damping force on node @a n */
class DampingForce : public Force{
 public:
    ~DampingForce(){};
    virtual Point operator()(Node n, double t){
        (void) t;
        return -c_damping * n.value().vel;
    }
};

/** functor that accepts a vector of forces and combines them */
struct CombinedForce{
    std::vector<Force*> force_vector_;

    Point operator()(Node n, double t){
        Point result(0);
        for (auto it = force_vector_.begin(); it!=force_vector_.end(); it++){
            result += (*(*it))(n, t);
        }
        return result;
    }

    CombinedForce(std::vector<Force*> force_vector): force_vector_(force_vector){}
};

/** function that combines different forces of Force parent class to return a combined force */
template <typename F1, typename F2, typename F3 = Force>
CombinedForce make_combined_force(F1 force1, F2 force2, F3 force3 = Force()){
    std::vector<Force*> forces;

    forces.push_back(&force1);
    forces.push_back(&force2);
    forces.push_back(&force3);

    CombinedForce combined_force(forces);

    return combined_force;
};

/** Parent constraint class */
class Constraint{
public:
    virtual void operator()(GraphType* g, double t){(void) g; (void) t;}
    virtual ~Constraint(){};
};

/** constraint that pins the node at (0,0,0) and (1,0,0) */
class PinConstraint: public Constraint{
 public:
    void operator()(GraphType* g, double t){
        (void) t;
        for (auto it = g->node_begin(); it!=g->node_end(); ++it){
            if ((*it).value().prev_position == Point(0) ){
                (*it).position() = Point(0);
                (*it).value().vel = Point(0);
            }
            else if ((*it).value().prev_position == Point(1,0,0) ){
                (*it).position() = Point(1,0,0);
                (*it).value().vel = Point(0);
            }
        }
    }
};

/** constraint for nodes not to go through a plane at z level z_. Takes double @a z as input */
class PlaneConstraint: public Constraint{
 public:
    double z_;
    void operator()(GraphType* g, double t){
        (void) t;
        for (auto it = g->node_begin(); it!=g->node_end(); ++it){
            if ((*it).position()[2] < z_){
                (*it).position()[2] = z_;
                (*it).value().vel[2] = 0;
            }
        }
    }
    PlaneConstraint(double z): z_(z){};
};

/** constraint for nodes not to go through a sphere with center c and radius r.
 * Takes double @a c as input
 * Takes double @a r as input*/
class SphereConstraint: public Constraint{
public:
    Point c_;
    double r_;

    void operator()(GraphType* g, double t){
        (void) t;
        for (auto it = g->node_begin(); it!=g->node_end(); ++it){
            if (norm((*it).position() - c_) < r_){
                Point prev_pos = (*it).position();

                //Update the position of the node
                (*it).position() = (r_ * ((prev_pos - c_)/norm(prev_pos-c_))) + c_;

                //update the velocity of the node
                Point R = ((*it).position()-c_) / norm((*it).position()-c_);
                (*it).value().vel -= inner_prod((*it).value().vel, R) * R;

            }
        }
    }
    SphereConstraint(Point c, double r): c_(c), r_(r){};
};

/** constraint for nodes not to be deleted when hitting a sphere with center c and radius r.
 * Takes double @a c as input
 * Takes double @a r as input*/
class SphereRemovalConstraint : public Constraint {
public:
    Point c_;
    double r_;

    void operator()(GraphType* g, double t) {
        (void) t;
        for (auto it = g->node_begin(); it != g->node_end(); ++it) {
            if (norm((*it).position() - c_) < r_) {
                g->remove_node(it);
            }
        }
    }
    SphereRemovalConstraint(Point c, double r): c_(c), r_(r){};
};

/** functor that accepts a vector of forces and combines them */
struct CombinedConstraints{
    std::vector<Constraint*> constraint_vector_;

    void operator()(GraphType* g, double t){
        for (auto it = constraint_vector_.begin(); it!=constraint_vector_.end(); it++){
            (*(*it))(g, t);
        }
    }

    CombinedConstraints(std::vector<Constraint*> constraint_vector): constraint_vector_(constraint_vector){}
};

/** Combining constraints of constraint parent class */
template <typename C1, typename C2, typename C3 = Constraint>
CombinedConstraints make_combined_constraint(C1 constraint1, C2 constraint2, C3 constraint3 = Constraint()){
    std::vector<Constraint*> constraints;

    constraints.push_back(&constraint1);
    constraints.push_back(&constraint2);
    constraints.push_back(&constraint3);

    CombinedConstraints combined_constraints(constraints);
    return combined_constraints;
};

/** Force function object for HW2 #1. */
struct Problem1Force {
    /** Return the force applying to @a n at time @a t.*/

    template <typename NODE>
    Point operator()(NODE n, double t) {
        (void) t;
        if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
            return Point(0,0,0);
        }

        //init force
        Point spring_force{0,0,0};

        //sum over adjacent nodes
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it){
            NODE n2 = (*it).node2();
            Point diff = n.position() - n2.position();

            spring_force += -(*it).value().K * (diff / (norm(diff))) * (norm(diff) - (*it).value().L);
            }

        Point resulting_force = spring_force + Point(0, 0, -grav * n.value().mass);
        return resulting_force;
    }
};




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
    bool origin;
    bool x_axis;
    NodeData() : vel(0), mass(1),origin(false),x_axis(false) {}
};
/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;       //< Node velocity
    double L;     //< Node mass
    EdgeData() : K(100.0), L(1.0) {}
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using incident_iterator = typename GraphType::IncidentIterator;
using NodeIter  = typename GraphType::node_iterator;


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
 */
template <typename G, typename F>
double symp_euler_step(G& g, double t, double dt, F force) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        if(n.position() != Point(0 ,0 ,0) && n.position() != Point(1 ,0 ,0)){
            // Update the position of the node according to its velocity
            // x^{n+1} = x^{n} + v^{n} * dt
            n.position() += n.value().vel * dt;
        }
    }

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        if(n.position() != Point(0 ,0 ,0) && n.position() != Point(1 ,0 ,0)){
            // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
            n.value().vel += force(n, t) * (dt / n.value().mass);
        }
    }

    return t + dt;
}



/** Similar to previous function but with added constraint
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint  Function object defining the constraint for the system
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

    constraint();

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
        if(n.position() == Point(0 ,0 ,0) || n.position() == Point(1 ,0 ,0)){
            return Point(0 ,0 ,0);
        }else{
            Point gravity_force = n.value().mass * Point(0, 0, -grav);
            Point spring_force = Point(0, 0, 0);
            for(incident_iterator it_inc = n.edge_begin(); it_inc != n.edge_end(); ++it_inc){
                Edge cur_edge = *it_inc;
                double K = cur_edge.value().K;
                double rest_len = cur_edge.value().L;
                Point coor_edge = n.position() - cur_edge.node2().position();
                double cur_len = norm(coor_edge);
                spring_force -= K*coor_edge/cur_len*(cur_len - rest_len);
            }
            return gravity_force + spring_force;
        }

    }
};


/**
 * Parent force class: construct force to apply to the nodes
 */
class Force{
public:
    /** Virtual function which given a node and a time
     * output the current force apply to this node
     */
    virtual Point operator()(Node n,double t) const{
        (void) n; (void) t;
        return Point(0,0,0);
    }
};


/**
 * First inherited Force which is the gravity Force
 */
class GravityForce : public Force{
public:
    Point operator()(Node n, double t) const{
        (void)t;
        Point gravity_force = n.value().mass * Point(0, 0, -grav);
        return gravity_force;
    }
};


/**
 * Second inherited Force which is the mass spring Force
 */
class MassSpringForce : public Force{
public:
    Point operator()(Node n, double t) const{
        (void)t;
        Point spring_force = Point(0, 0, 0);
        for(incident_iterator it_inc = n.edge_begin(); it_inc != n.edge_end(); ++it_inc){
            Edge cur_edge = *it_inc;
            double K = cur_edge.value().K;
            double rest_len = cur_edge.value().L;
            Point coor_edge = n.position() - cur_edge.node2().position();
            double cur_len = norm(coor_edge);
            spring_force -= K*coor_edge/cur_len*(cur_len - rest_len);
        }
        return spring_force;
    }
};


/**
 * Third inherited Force which is the damping Force
 */
class DampingForce : public Force{
public:
    Point operator()(Node n, double t) const{
        (void)t;
        // no dampling constant for the moment
        double c = 1.0;
        Point damping_force = n.value().vel*c;
        return damping_force;
    }
};

/**
 * Functor to combine a vector of Forces
 */
struct CombinedForce{
    CombinedForce(std::vector<const Force*> forces): forces_(forces) {}
    Point operator()(Node n, double t) const{
        Point total_force = Point(0,0,0);
        unsigned int end_it = forces_.size();
        for (unsigned int it = 0; it<end_it; ++it){
            total_force += forces_[it]->operator()(n,t);
        }
        return total_force;
    }
private:
    std::vector<const Force*> forces_;
};


/**
 * Combine the input parameters forces to get the total force
 * @param force1 : the first type of force we need to consider, must be inherited from Force
 * @param force2 : the second type of force we need to consider, must be inherited from Force
 * @return Total force of the forces considered
 */
//--style_0
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
template <typename Force1, typename Force2>
CombinedForce make_combined_force(const Force1& force1, const Force2& force2){
    std::vector<const Force*> forces{&force1,&force2};
    return CombinedForce(forces);
}
//--END

/**
 * Same function as above but accepts 3 forces
 * @param force1 : the first type of force we need to consider, must be inherited from Force
 * @param force2 : the second type of force we need to consider, must be inherited from Force
 * @param force3 : the third type of force we need to consider, must be inherited from Force
 * @return Total force of the forces considered
 */
template <typename Force1, typename Force2, typename Force3>
CombinedForce make_combined_force(const Force1& force1, const Force2& force2, const Force3& force3){
    std::vector<const Force*> forces{&force1, &force2, &force3};
    return CombinedForce(forces);
}


/**
 * Constraint parent class to apply constraints to the nodes
 */
class Constraint{
public:
    /** Virtual function which given a node and a time
     * output the current force apply to this node
     */
    virtual void operator()() const = 0;
protected:
    GraphType * graph_;
    double t_;
    /** Private Constructor */
    Constraint(GraphType* graph, double t)
            : graph_(graph), t_(t){
    };
};


/**
 * First inherited Constraint which keeps the nodes
 * (0,0,0) and (1,0,0) fix
 */
class PinConstraint : public Constraint{
public:
    void operator()() const{
        for(auto it = graph_->node_begin(); it!=graph_->node_end();++it){
            Node cur_node = (*it);
            if(cur_node.value().origin){
                cur_node.position() = Point(0,0,0);
            }else if(cur_node.value().x_axis){
                cur_node.position() = Point(1,0,0);
            }
        }
    }
    PinConstraint(GraphType *graph, double t) : Constraint(graph, t) {};
};


/**
 * Second inherited Constraint which create a wall defined
 * by the plane z=-0.75
 */
class PlaneConstraint : public Constraint{
public:
    void operator()() const{
        double z = -0.75;
        for(auto it = graph_->node_begin(); it!=graph_->node_end();++it){
            Node cur_node = (*it);
            if(std::isless(cur_node.position().z,z)){
                cur_node.position().z = z;
                cur_node.value().vel.z = 0;
            }
        }
    }
    PlaneConstraint(GraphType *graph, double t) : Constraint(graph, t) {};
};


/**
 * Third inherited Constraint which restrict the nodes
 * outside the sphere of center (0.5,0.5,-0.5) and radius r=0.15
 */
class SphereConstraint : public Constraint{
public:
    void operator()() const{
        Point center = Point(0.5,0.5,-0.5);
        double radius = 0.15;
        for(auto it = graph_->node_begin(); it!=graph_->node_end();++it){
            Node cur_node = (*it);
            Point diff = cur_node.position()-center;
            double euc_norm = norm(diff);
            if(std::isless(0.0,euc_norm) && std::isless(euc_norm,radius)){
                Point normal = diff / euc_norm;
                cur_node.position() = center + radius * normal;
                cur_node.value().vel -= dot(cur_node.value().vel,normal)*normal;
            }else if(!std::isless(0.0,euc_norm)){
                // any point on the sphere is valid
                cur_node.position() = center + Point(radius,0,0);
                cur_node.value().vel = Point(0);
            }

        }
    }
    SphereConstraint(GraphType *graph, double t) : Constraint(graph, t) {};
};



/** Test predicate : node is in the sphere */
struct TearPredicate {
    Point center_;
    double radius_;
    template <typename NODE>
    bool operator()(const NODE& n) {
        Point diff = n.position()-center_;
        double euc_norm = norm(diff);
        return (std::isless(euc_norm,radius_));
    }
    TearPredicate(Point center, double radius) :center_(center),radius_(radius){};
};


/**
 * Sphere Constraint which delete the nodes
 * outside the sphere of center_ and radius_ (from the predicate)
 */
class TearConstraint : public Constraint{
public:
    void operator()() const{
        auto it_beg = graph_->node_begin();
        auto it_end = graph_->node_end();
        auto it_erase = std::find_if(it_beg, it_end, pred_);
        while(it_erase!=it_end){
            it_erase = graph_->remove_node(it_erase);
            it_end = graph_->node_end();
            it_erase = std::find_if(it_erase, it_end, pred_);
        }
    }
    TearConstraint(GraphType *graph, double t, Point center, double radius)
            :Constraint(graph, t), pred_(TearPredicate(center,radius)){};

private:
    TearPredicate pred_;

};


/**
 * Functor to combine a vector of Forces
 */
struct CombinedConstraint{
    CombinedConstraint(std::vector<const Constraint*> constraints): constraints_(constraints) {}
    void operator()(){
        unsigned int end_it = constraints_.size();
        for (unsigned int it = 0; it<end_it; ++it){
            constraints_[it]->operator()();
        }
    }
private:
    std::vector<const Constraint*> constraints_;
};


/**
 * Combine the input parameters constraints to get the total Constraints
 * @param Constraint1 : the first type of constraint we need to consider, must be inherited from constraint
 * @param Constraint2 : the second type of constraint we need to consider, must be inherited from constraint
 * @return Total constraint of the constraints considered
 */
template <typename Constraint1, typename Constraint2>
CombinedConstraint make_combined_constraint(const Constraint1& constraint1, const Constraint2& constraint2){
    std::vector<const Constraint*> constraints{&constraint1,&constraint2};
    return CombinedConstraint(constraints);
}

/**
 * Same function as above but accepts 3 constraints
 * @param constraint1 : the first type of constraint we need to consider, must be inherited from constraint
 * @param constraint2 : the second type of constraint we need to consider, must be inherited from constraint
 * @param constraint2 : the third type of constraint we need to consider, must be inherited from constraint
 * @return Total constraint of the constraints considered
 */
template <typename Constraint1, typename Constraint2, typename Constraint3>
CombinedConstraint make_combined_constraint(const Constraint1& constraint1, const Constraint2& constraint2,
                                            const Constraint3& constraint3){
    std::vector<const Constraint*> constraints{&constraint1, &constraint2, &constraint3};
    return CombinedConstraint(constraints);
}


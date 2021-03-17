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
  NodeData() : vel(0), mass(1) {};
};

struct EdgeData {
    double length;       //< Node velocity
    double spring_constant;     //< Node mass
    NodeData() : length(0), spring_constant(double(100)) {};
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
double symp_euler_step(G& g, double t, double dt, F force, C constraint) {
  // Compute the t+dt position
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // Update the position of the node according to its velocity
    // x^{n+1} = x^{n} + v^{n} * dt
    n.position() += n.value().vel * dt;
  }
  // Apply the constraints
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;
    constraint(n, t); //this will update the nodes affected by the constraints
  }
  // Compute the t+dt velocity
  for (auto it = g.node_begin(); it != g.node_end(); ++it) {
    auto n = *it;

    // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
    n.value().vel += force(n, t) * (dt / n.value().mass);
  }

  return t + dt;
}

EdgeData initialize(Node node_a, Node node_b){
    Point diff_point = node_a.position() - node_b.position(); //subtraction is defined in Point.hpp
    EdgeData init = EdgeData();
    init.length = normSq(diff_point); //initial length
    init.spring_constant = double(100); //spring constant is always 100
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
    //need to iterate over neighbours of n
    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
        return Point(0,0,0);
    }
    Point f_spring = Point(0,0,0);
    Point grav = Point(0, 0, -9.8);
    std::vector<Node> nbr = n.neighbours();
    for (auto it = nbr.begin(); it != nbr.end(); ++it){
        Edge e = n.find_edge(*it); //returns the edge connecting node n and neighbour *it
        double len = e.value().length;
        double K = e.value().spring_constant;
        Point diff_point = n.position() - *it->position(); //subtraction is defined in Point.hpp
        double dist = normSq(diff_point);
        f_spring += -K*(n.position()-*it->position())/dist *(dist - len);
    }
    return f_spring + n.value().mass*grav;
    (void) n; (void) t; (void) grav;    // silence compiler warnings
    return Point(0);
  }
};

class Force{
    virtual Point operator()(Node n, double t){
        return Point(0,0,0);
    }
};

class GravityForce: public Force{
    virtual Point operator()(Node n, double t){
        return n.value().mass*Point(0,0,-9.8);
    }
};

class MassSpringForce: public Force{
    virtual Point operator()(Node n, double t){
        Point f_spring = Point(0,0,0);
        std::vector<Node> nbr = n.neighbours();
        for (auto it = nbr.begin(); it != nbr.end(); ++it){
            Edge e = n.find_edge(*it); //returns the edge connecting node n and neighbour *it
            double len = e.value().length;
            double K = e.value().spring_constant;
            Point diff_point = n.position() - *it->position(); //subtraction is defined in Point.hpp
            double dist = normSq(diff_point);
            f_spring += -K*(n.position()-*it->position())/dist *(dist - len);
        }
        return f_spring;
    }
};

class DampingForce: public Force{
    virtual Point operator()(Node n, double t){
        double c = double(10);
        return -c*n.value().vel;
    }
};

class CombinedForce{
    std::vector<Force*> forces;
    CombinedForce(std::vector<Force*> F){
        this ->forces = F;
    }
    virtual Point operator()(Node n, double t){
        Point comb_force = Point(0);
        for(auto it = forces.begin(); it != forces.end(); ++it){
            comb_force += *it(n, t);
        }
        return comb_force;
    }
};

template<typename Force1, typename Force2>
CombinedForce make_combined_force(Force1 f1, Force2 f2){
    std::vector<Force*> F;
    F.push_back(f1);
    F.push_back(f2);
    return CombinedForce(F);
}

template<typename Force1, typename Force2, typename Force3>
CombinedForce make_combined_force(Force1 f1, Force2 f2, Force3 f3){
    std::vector<Force*> F;
    F.push_back(f1);
    F.push_back(f2);
    F.push_back(f3);
    return CombinedForce(F);
}

class Constraint{
    virtual operator(Node n, double t){
    }
};

class PinConstraint: public Constraint{
    virtual operator(Node n, double t){
        if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
            n.value().vel = double(0);
        }
    }
};

class PlainConstraint: public Constraint{
    double z = -0.75;
    virtual operator (Node n, double t){
        if (dot(n.position(), Point(0, 0, 1)) < z){
            n.position()->elem[2] = z;
            n.value().vel->elem[2] = double(0);
        }
    }
};

class SphereConstraint: public Constraint{
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    virtual operator (Node n, double t){
        Point diff_point = n.position() - center; //subtraction is defined in Point.hpp
        double dist = normSq(diff_point);
        if (dist < radius){
            n.position() = center + (radius/dist)*diff_point;
            n.value().vel -= dot(n.value().vel, diff_point/dist)*diff_point/dist;
        }
    }
};

class TearConstraint: public Constraint{
    Point center = Point(0.5, 0.5, -0.5);
    double radius = 0.15;
    virtual operator (Node n, double t){
        Point diff_point = n.position() - center; //subtraction is defined in Point.hpp
        double dist = normSq(diff_point);
        if (dist < radius){
            remove_node(n);
        }
    }
};

class CombinedConstraints{
    std::vector<Constraint*> constraints;
    CombinedConstraints(std::vector<Constraint*> C){
        this->constraints = C;
    }
    virtual operator()(Node n, double t){
        for(auto it = constraints.begin(); it != constraints.end(); ++it){
            *it(n, t);
        }
    }
};

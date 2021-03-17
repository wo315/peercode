/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */


#include <fstream>
#include <iostream>
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
  
  bool pin1 = false;
  bool pin2 = false;

  NodeData() : vel(0), mass(1) {}
};

// Define the Graph type

//using GraphType = Graph<NodeData>;
using GraphType = Graph<NodeData, double>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

/*
 * @brief Base force class
 *
 * each force child class has the () operation overloaded
 * the () takes in a node _n_ and a double _t_ and computes
 * the force on _n_
 */ 
class Force {
    public:

    virtual Point operator() (Node n, double t) {
        (void) n; (void) t;

        return Point(0, 0, 0);
    }
};

/* 
 * @brief Force class that implements force from gravity
 *
 * @param[in,out] n node to apply force of gravity to
 * @param[in]     t time
 * @return force of gravity
 */
class GravityForce : public Force {
    public:

    virtual Point operator() (Node n, double t) {
        (void) t;

        return Point(0, 0, -n.value().mass*grav);
    }
};

/*
 * @brief implements force from springs
 *
 * @param[in,out] n node to apply spring force to
 * @param[in]     t time
 * @return spring force on _n_
 */ 
class MassSpringForce : public Force {
    public:

    virtual Point operator() (Node n, double t) {
        (void) t;

        Point force(0, 0, 0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Edge ed = *it;
            Point vect = n.position() - ed.node2().position();

            double length = norm(vect);
            vect /= length;
            double L = ed.value();
            vect *= (length - L);

            force += vect;
        }

        force *= -100;
        return force;
    }

};

/*
 * @brief implementation of damping force
 *
 * @param[in,out] n node to apply damping force to
 * @param[in      t time
 * @return -_c_*velocity, a damping force
 */ 
class DampingForce : public Force {
    public:

    virtual Point operator() (Node n, double t) {
        (void) t;

        return -c_*n.value().vel;
    }

    private:
    double c_=0.00005;
};

/*
 * @brief combines a set of forces by summing them
 *
 * @pre @forces_ is a vector of Force objects, each of which has an overloaded ()
 *      operator that returns a force
 * @param[in,out] n node to compute forces with respect to
 * @param[in]     t time
 * @return sum of values of forces from forces_ vector
 */
struct CombinedForce {

    CombinedForce(std::vector<Force*> forces) : forces_{forces} {};

    Point operator() (Node n, double t) {
        Point total_force(0, 0, 0);    
        for (auto f : forces_) {
            total_force += (*f)(n, t);
        }
        return total_force;
    }

    private:
        std::vector<Force*> forces_;
};

/*
 * @brief combines given forces into one to compute the total force
 *
 * @param f1, f2 Force objects to sum up
 * @return functor that computes the sum of _f1_ and _f2_ when called
 */
CombinedForce make_combined_force(Force&& f1, Force&& f2) {
    std::vector<Force*> forces;
    forces.push_back(&f1);
    forces.push_back(&f2);
    return CombinedForce(forces);
}

/*
 * @brief combines given forces into one to compute the total force
 *
 * @param f1, f2, f3 Force objects to sum up
 * @return functor that returns the sum of _f1_ _f2_ and _f2_ when called
 */
CombinedForce make_combined_force(Force&& f1, Force&& f2, Force&& f3) {
    std::vector<Force*> forces;
    forces.push_back(&f1);
    forces.push_back(&f2);
    forces.push_back(&f3);
    return CombinedForce(forces);
}

/*
 * @brief base class representing constraints on a graph
 *
 * each child class has the () operator overloaded
 * the () operator takes a graph _g_ and a double _t_ and applies the 
 * constraints to all nodes in the graph
 */ 
class Constraint {
    public:
    
    virtual void operator() (GraphType& g, double t) = 0;
};

/*
 * @brief class that keeps corners of a grid fixed
 *
 * @param[in,out] g input graph
 * @param[in]     t time
 * @post resets the positions of the nodes that began with positions 
 *       (0,0,0) and (1,0,0), all other nodes left unchanged
 */
class PinConstraint : public Constraint {
    public:
    
    virtual void operator() (GraphType& g, double t) {
        (void) t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if (n.value().pin1) {
                n.position() = Point(0, 0, 0);
            }
            if (n.value().pin2) {
                n.position() = Point(1, 0, 0);
            }    
        }
    }
};

/*
 * @brief Constraint that keeps nodes from passing given plane
 *
 * @param[in,out] g input graph
 * @param[in]     t time
 * @post for all nodes n in graph, if n's z coordinate was < -0.75, 
 *       n's new z coordinate == 0.75. Else, n's old position == n's new position
 */
class PlaneConstraint : public Constraint {
    public:
    
    PlaneConstraint() {}

    virtual void operator() (GraphType& g, double t) {
        (void) t;
        Point p(0, 0, 1);
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            Point pos = n.position();
            
            if (dot(pos, p) < plane) {
                pos = Point(pos[0], pos[1], plane);
                Point vel(n.value().vel[0], n.value().vel[1], 0);
                n.value().vel = vel;
            }
        }
    }

    double plane = -0.75;
};

/*
 * @brief constraint that acts as an invisible sphere during simulations
 *
 * @param[in,out] g input graph
 * @param[in]     t time 
 * @post for all nodes n, if |n.position()-_center_| < _r_, n's new position becomes the 
 *       closest point on the surface of the sphere and n's velocity is set to be normal 
 *       to the surface of the sphere. Else, n remains unchanged
 */
class SphereConstraint : public Constraint {
    public:
    
    virtual void operator() (GraphType& g, double t) {
        (void) t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            Point pos = n.position();
            
            if (norm(pos - center) < r) {
                Point diff = pos - center;
                diff /= norm(diff);
                diff *= r;
                Point new_pos = center + diff;
                n.position() = new_pos;

                Point R = (pos - center);
                R /= norm(R);
                n.value().vel -= dot(n.value().vel, R)*R;
            }
        }
    }

    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;
};

/*
 * @brief constraint that simulates a tear in the graph
 *
 * @param[in,out] g input graph
 * @param[in]     t time
 * @post for all nodes n in graph _g_, if norm(n's position-_center_) < r, 
                                          node is removed from graph
 */
class TearConstraint : public Constraint {
    public:
    
    virtual void operator() (GraphType& g, double t) {
        (void) t;
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            Point pos = n.position();

            if (norm(pos - center) < r) {
                g.remove_node(it);
            }
        }
    }

    Point center = Point(0.5, 0.5, -0.5);
    double r = 0.15;
};

/*
 * @brief struct for applying a vector of constraints to a graph
 */
struct CombinedConstraints {
   
    CombinedConstraints(std::vector<Constraint*> constraints) : constraints_{constraints} {};
 
    void operator() (GraphType& g, double t) {
        
        for (auto c : constraints_) {
            (*c)(g, t);
        }
    }
 
    private:
        std::vector<Constraint*> constraints_;
};

/*
 * @brief function for applying 2 given constraints to a graph
 * @param c1, c2 Constraint objects to be used in graph
 *
 * @return CombinedConstraints functor that applies the given constraints when invoked
 */
CombinedConstraints make_combined_constraints(Constraint&& c1, Constraint&& c2) {
    std::vector<Constraint*> cons {};
    cons.push_back(&c1);
    cons.push_back(&c2);
    return CombinedConstraints(cons);
}

/*
 * @brief function for applying 3 given constraints to a graph
 * @param c1, c2, c3 Constraint objects to be used in graph
 *
 * @return CombinedConstraints functor that applies the given constraints when invoked
 */
CombinedConstraints make_combined_constraints(Constraint&& c1, Constraint&& c2, Constraint&& c3) {
    std::vector<Constraint*> cons {};
    cons.push_back(&c1);
    cons.push_back(&c2);
    cons.push_back(&c3);
    return CombinedConstraints(cons);
}

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
    
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
        n.position() += n.value().vel * dt;
    }

    //apply constraints
    constraint(g, t);

    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;
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
    (void) n; (void) t; (void) grav;    // silence compiler warnings

    if (n.position() == Point(0,0,0) || n.position() == Point(1,0,0)) {
        return Point(0,0,0);
    }

    Point force(0, 0, 0);
    for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {

        Edge ed = *it;

        Point vect = n.position() - ed.node2().position();

        double length = norm(vect);
        vect /= length;

        vect *= (length - ed.value());

        force += vect;
    }

    force *= -100;

    force += n.value().mass*gravity;
    return force;
  }

  Point gravity = Point(0, 0, -grav);
};

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


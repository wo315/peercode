/**
 * @file mass_spring.hpp
 * Implementation of mass-spring system using Graph
 */

#include <chrono>
#include <fstream>
#include <thread>

#include "CME212/Color.hpp"
#include "CME212/Point.hpp"
#include "CME212/Util.hpp"
#include "Graph.hpp"

// Gravity in meters/sec^2
static constexpr double grav = 9.81;

/** Custom structure of data to store with Nodes */
struct NodeData {
    Point x0;     //< Store starting position
    Point vel;    //< Node velocity
    double mass;  //< Node mass
    NodeData() : x0(0), vel(0), mass(1) {}
    NodeData(Point x, Point v, double m) : x0(x), vel(v), mass(m) {}
};

// Define the Graph type
using GraphType = Graph<NodeData, double>;
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
double symp_euler_step(G& g, double t, double dt, F force, C constraints) {
    // Compute the t+dt position
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // Update the position of the node according to its velocity
        // x^{n+1} = x^{n} + v^{n} * dt
        n.position() += n.value().vel * dt;
    }

    // apply all constraints
    constraints(g, t);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }

    return t + dt;
}

/**
 * Base class for all forces.
 */
class Force {
   public:
    virtual Point operator()(Node node, double t) const = 0;
};

class GravityForce : public Force {
   public:
    /**
     * The node must have a mass attribute in its value() attribute.
     */
    Point operator()(Node node, double t) const {
        (void)t;  // no-op
        return node.value().mass * Point(0, 0, -grav);
    }
};

class MassSpringForce : public Force {
   public:
   /**
    * `L` in the spring force should be defined in the 
    * incident edge's `value()` attribute.
    */
    Point operator()(Node node, double t) const {
        (void)t;  // no-op
        short K = 100;
        Point force = Point(0);
        for (auto it = node.edge_begin(); it != node.edge_end(); ++it) {
            Edge edge = *it;
            Node one = edge.node1();
            Node two = edge.node2();

            const Point dx = one.position() - two.position();
            const double scale = -K * (norm(dx) - edge.value()) / norm(dx);
            force += (scale * dx);
        }
        return force;
    }
};

class DampingForce : public Force {
   public:
   /**
     * The node must have a vel attribute in its value() attribute.
     */
    Point operator()(Node node, double t) const {
        (void)t;
        short c = 1;
        return c * node.value().vel;
    }
};

class CombinedForce : public Force {
   private:
    std::vector<const Force*> forces;

   public:
    CombinedForce(std::vector<const Force*> given) : forces(given) {}

    Point operator()(Node node, double t) const {
        (void)t;
        Point total = Point(0);
        for (auto it = forces.begin(); it != forces.end(); ++it) {
            const Force* force = *it;
            total += (*force)(node, t);
        }
        return total;
    }
};


//--style_1
//--This should be encompassed with a default argument
//--in the 3 arguments method
//--START
CombinedForce make_combined_force(const Force& first, const Force& second) {
    return CombinedForce(std::vector<const Force*>({&first, &second}));
};
//--END

//--design_-1
//--Good Job dealing witht the rvalues
//--START
CombinedForce make_combined_force(const Force& first, const Force& second,
                                  const Force& third) {
    return CombinedForce(std::vector<const Force*>({&first, &second, &third}));
};
//--END

class Constraint {
   public:
    virtual void operator()(GraphType& g, double t) const = 0;
};

class PinConstraint : public Constraint {
   public:
    void operator()(GraphType& g, double t) const {
        (void)t;

        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if (n.value().x0 == Point(0, 0, 0)) {
                n.position() = Point(0, 0, 0);
            }
            if (n.value().x0 == Point(1, 0, 0)) {
                n.position() = Point(1, 0, 0);
            }
        }
    }
};

class PlaneConstraint : public Constraint {
   public:
    void operator()(GraphType& g, double t) const {
        (void)t;

        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if (dot(n.position(), Point(0, 0, 1)) <= -0.75) {
                // Set the z component of velocity to zero
                n.value().vel = n.value().vel * Point(1, 1, 0);

                // Force point to closest position on z = -0.75 plane
                n.position() =
                    n.position() * Point(1, 1, 0) + Point(0, 0, -0.75);
            }
        }
    }
};

class SphereConstraint : public Constraint {
   public:
    void operator()(GraphType& g, double t) const {
        (void)t;

        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;

        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if (norm(n.position() - c) < r) {
                // Find unit vector pointing in the direction from
                // the center of the ball to the node's position
                Point R = (n.position() - c) / (norm(n.position() - c));

                // Scale vector by ball's radius, and add to center coordinate
                n.position() = c + R * r;

                Point velocity = n.value().vel;
                velocity -= dot(velocity, R) * R;
                n.value().vel = velocity;
            }
        }
    }
};

class DestructiveSphereConstraint : public Constraint {
   public:
    void operator()(GraphType& g, double t) const {
        (void)t;

        Point c = Point(0.5, 0.5, -0.5);
        double r = 0.15;

//--design_1
//--You're skipping nodes with this (the one at the taking the place you 
//--just deleted)
//--START
        for (auto it = g.node_begin(); it != g.node_end(); ++it) {
            Node n = *it;
            if (norm(n.position() - c) < r) {
                // DESTROY!
                g.remove_node(n);
            }
        }
//--END
    }
};

class CombinedConstraint : public Constraint {
   private:
    std::vector<const Constraint*> _constr;

   public:
    CombinedConstraint(std::vector<const Constraint*> given) : _constr(given) {}

    void operator()(GraphType& g, double t) const {
        for (auto c_it = _constr.begin(); c_it != _constr.end(); ++c_it) {
            const Constraint* c = *c_it;
            (*c)(g, t);
        }
    }
};

CombinedConstraint make_combined_constraint(const Constraint& first,
                                            const Constraint& second) {
    return CombinedConstraint(
        std::vector<const Constraint*>({&first, &second}));
};

CombinedConstraint make_combined_constraint(const Constraint& first,
                                            const Constraint& second,
                                            const Constraint& third) {
    return CombinedConstraint(
        std::vector<const Constraint*>({&first, &second, &third}));
};

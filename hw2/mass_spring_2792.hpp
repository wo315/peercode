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
    bool pinned_point_0;  //< Point at P(0, 0, 0)
    bool pinned_point_1;  //< Point at P(1, 0, 0)
    NodeData() : vel(0), mass(1), pinned_point_0(false), pinned_point_1(false) {}
};

/** Custom structure of data to store with Edges */
struct EdgeData {
    double K;       // Spring constant
    double L;     //< Length
    EdgeData() : K(100), L(0) {}
};

// Define the Graph type and objects
using GraphType = Graph<NodeData, EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;
using size_type = unsigned;


/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData (a struct comprising the
 * velocity and mass of a Node, and whether it is pinned to P(0, 0, 0) or
 * P(1, 0, 0))
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
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
            n.value().vel += Point(0, 0, 0) * (dt / n.value().mass);
        }
        else {
            n.value().vel += force(n, t) * (dt / n.value().mass);
        }
    }

    return t + dt;
}


/** Force function object */
struct Problem1Force {
    /** Return the force applying to @a n at time @a t.
     *
     * This is a combination of mass-spring force and gravity,
     * except that points at (0, 0, 0) and (1, 0, 0) never move. We can
     * model that by returning a zero-valued force. */
    template <typename NODE>
    Point operator()(NODE n, double t) {
        if (n.position() == Point(0, 0, 0) || n.position() == Point(1, 0, 0)) {
            return Point(0, 0, 0);
        }
        else {
            // Compute spring force
            Point spring_force = Point(0);
            for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
                Edge e_ij = *it;
                double K_ij = e_ij.value().K;
                double L_ij = e_ij.value().L;
                Point dif_ij = n.position() - e_ij.node2().position();
                double norm_dif_ij = norm(dif_ij);
                spring_force += -K_ij * (dif_ij / norm_dif_ij) * (norm_dif_ij - L_ij);
            }

            Point grav_force = n.value().mass * Point(0, 0, -grav);
            return spring_force + grav_force;
        }
    }
};

/* ------------------------- Forces --------------------------- */

/**
 * Parent Force class
 * Returns a zero force.
 */
class Force {
public:
    virtual Point operator()(Node n, double t) {
        (void)n;
        (void)t;
        return Point(0, 0, 0);
    }
};

/**
 * Derived GravityForce class
 * Implements the force of gravity
 */
class GravityForce : public Force {
public:
    virtual Point operator()(Node n, double t) {
        (void)t;
        return n.value().mass * Point(0, 0, -grav);
    }
};

/**
 * Derived MassSpringForce class
 * Implements the spring force
 */
class MassSpringForce : public Force {
public:
    virtual Point operator()(Node n, double t) {
        Point spring_force = Point(0);
        for (auto it = n.edge_begin(); it != n.edge_end(); ++it) {
            Edge e_ij = *it;
            double K_ij = e_ij.value().K;
            double L_ij = e_ij.value().L;
            Point dif_ij = n.position() - e_ij.node2().position();
            double norm_dif_ij = norm(dif_ij);
            spring_force += -K_ij * (dif_ij / norm_dif_ij) * (norm_dif_ij - L_ij);
        }
        (void)t;
        return spring_force;
    }
};

/**
 * Derived DampingForce class
 * Implements damping (a form of friction) according to the velocity of the
 * node and a damping constant
 */
class DampingForce : public Force {
public:
    virtual Point operator()(Node n, double c) {
        return -n.value().vel * c;
    }
};

/** 
 * CombinedForce Functor 
 * @brief Combines a series of forces by adding each derived force class
 * @param[in] @n a Node on which the force will act
 * @param[in] @t a double representing time
 * @return a Point representing the combined force
 */
struct CombinedForce {
    std::vector<Force*> force_vec_;
    CombinedForce(std::vector<Force*> force_vec) : force_vec_{ force_vec } {}

    Point operator()(Node n, double t) {
        Point combined_force = Point(0, 0, 0);
        for (unsigned int i = 0; i < force_vec_.size(); i++) {
            combined_force += (*force_vec_[i])(n, t);
        }
        return combined_force;
    }
};

/**
 * Make a combination of forces
 * @brief Combines different force types passed as arguments
 * @param[in] The derived force classes we wish to act upon the graph nodes
 * @return a CombinedForce struct comprising the different forces
 */
template<class...Fs>
CombinedForce make_combined_force(Fs...args) {
    std::vector<Force*> forces_ = { &args... };
    return CombinedForce(forces_);
}


/* ------------------------- Constraints --------------------------- */

/**
 * Parent Constraint
 * Implements a virtual operator() 
 */
class Constraint {
public:
    virtual void operator()(GraphType& graph_, double t) {
        (void)graph_;
        (void)t;
    }
};

/**
 * Pin Constraint
 * Keeps the nodes at P(0, 0, 0) and P(1, 0, 0) fixed
 */
class PinConstraint : public Constraint {
public:
    /**
     * PinConstraint operator
	 * @param[in, out]: @graph_ the graph to be modified with the constraints
     * @param[in]: @t a time double
     */
    virtual void operator()(GraphType& graph_, double t) {
        // Find pinned points and fix their positions
        for (auto it = graph_.node_begin(); it != graph_.node_end(); ++it) {
            Node n = *it;
            if (n.value().pinned_point_0) {
                n.position() = Point(0, 0, 0);
            }
            if (n.value().pinned_point_1) {
                n.position() = Point(1, 0, 0);
            }
        }
        (void)t;
    }
};

/**
 * Plane Constraint
 * The constraint is violated for any node falling below the plane
 * z = -0.75. Nodes in violation are fixed by re-positioning them
 * to the closest point on this plane and adjusting their velocity.
 */
class PlaneConstraint : public Constraint {
public:
    /**
     * PlaneConstraint operator
     * @param[in, out]: @graph_ the graph to be modified with the constraints
     * @param[in]: @t a time double
     */
    virtual void operator()(GraphType& graph_, double t) {
        // Define constraint
        double z = -0.75;

        // Find nodes that violate the constraint
        for (auto it = graph_.node_begin(); it != graph_.node_end(); ++it) {
            // Check predicate and fix
            Node n = *it;
            if (dot(n.position(), Point(0, 0, 0)) < z) {
                // Set the position to the nearest point on the plane
                n.position().z = -0.75;

                // Set the z-component of Node velocity to zero
                n.value().vel.z = 0;

            }
        }
        (void)t;
    }
};

/**
 * Sphere Constraint
 * The constraint is violated by any node falling within a sphere 
 * with center c and radius r. Nodes are fixed by setting their 
 * position to the nearest point on the surface of the sphere and
 * setting the component of the velocity that is normal to the sphere's
 * surface to zero.
 */
class SphereConstraint : public Constraint {
public:
    /**
     * SphereConstraint operator
     * @param[in, out]: @graph_ the graph to be modified with the constraints
     * @param[in]: @t a time double
     */
    virtual void operator()(GraphType& graph_, double t) {
        // Define constraint
        Point sphere_center = Point(0.5, 0.5, -0.5);
        double sphere_radius = 0.15;

        // Find nodes violating the constraint and fix
        for (auto it = graph_.node_begin(); it != graph_.node_end(); ++it) {
            Node n = *it;
            if (norm(n.position() - sphere_center) < sphere_radius) {
                // Set the position to the nearest point on the surface of the sphere
                Point R_i = (n.position() - sphere_center) / norm(n.position() - sphere_center);
                n.position() = sphere_center + sphere_radius * R_i;

                // Set the component of the velocity that is normal to the sphere's surface to zero 
                n.value().vel -= dot(n.value().vel, R_i) * R_i;
            }
        }
        (void)t;
    }
};


/**
 * Sphere Constraint with Node Removal
 * The constraint is violated by any node falling within a sphere
 * with center c and radius r. Nodes are fixed by removal.
 */
class SphereConstraintRemoval : public Constraint {
public:
    /**
     * SphereConstraintRemoval operator
     * @param[in, out] @graph_ the graph to be modified with the constraints
     * @param[in] @t a time double
     */
    virtual void operator()(GraphType& graph_, double t) {
        // Define constraint
        Point sphere_center = Point(0.5, 0.5, -0.5);
        double sphere_radius = 0.15;

        // Remove nodes violating the constraint
        for (auto it = graph_.node_begin(); it != graph_.node_end(); ) {
            Node n = *it;
            if (norm(n.position() - sphere_center) < sphere_radius) {
                it = graph_.remove_node(it);
            }
            else {
                ++it;
            }
        }

        (void)t;
    }
};


/**
 * CombinedConstraints Functor
 * @brief Combines a vector of derived constraints
 * @param[in, out] @graph_ the graph to be modified with the constraints
 * @param[in] @t a time double
 */
template<class...Cs>
struct CombinedConstraints {
    /**
     * Combines the derived constraint classes
     */
    std::vector<Constraint*> constraint_vec_;
    CombinedConstraints(std::vector<Constraint*> force_vec) : constraint_vec_{ force_vec } {}

    void operator()(GraphType& graph_, double t) {
        for (unsigned int i = 0; i < constraint_vec_.size(); i++) {
            (*constraint_vec_[i])(graph_, t);
        }
    }
};

/**
 * @brief Make a combination of constraints
 * @param[in] A series of derived constraints
 * @return a CombinedConstraints functor
 */
template<class...Cs>
CombinedConstraints<Cs...> make_combined_constraint(Cs...args) {
    std::vector<Constraint*> constraints_ = { &args... };
    return CombinedConstraints<Cs...>(constraints_);
}



/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g      Graph
 * @param[in]     t      The current time (useful for time-dependent forces)
 * @param[in]     dt     The time step
 * @param[in]     force  Function object defining the force per node
 * @param[in]     constraint Function object defining the graph constraints
 * @return the next time step (usually @a t + @a dt)
 *
 * @tparam G::node_value_type supports NodeData (a struct comprising the
 * velocity and mass of a Node, and whether it is pinned to P(0, 0, 0) or
 * P(1, 0, 0))
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
    constraint(g, t);

    // Compute the t+dt velocity
    for (auto it = g.node_begin(); it != g.node_end(); ++it) {
        auto n = *it;

        // v^{n+1} = v^{n} + F(x^{n+1},t) * dt / m
        n.value().vel += force(n, t) * (dt / n.value().mass);
    }
    return t + dt;
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


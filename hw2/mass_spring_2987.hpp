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
  NodeData() : vel(0), mass(1) {}
};

struct EdgeData {
    double K;     // Edge spring constant
    double L;     // Edge spring L
    EdgeData() : K(100), L(0) {}
};

struct Pin {
  Point &pos;
  Point ref_pos;
  Point &vel;

  Pin(Point &pos_, Point ref_pos_, Point &vel_)
    : pos(pos_), ref_pos(ref_pos_), vel(vel_) {
    }
};

// Define the Graph type
using GraphType = Graph<NodeData,EdgeData>;
using Node = typename GraphType::node_type;
using Edge = typename GraphType::edge_type;

struct ZeroForce {

    /** Invalid default constructor, destructor and operato to be later implemented */
    ZeroForce(){}
    virtual Point operator()(Node n, double t) = 0;
    virtual ~ZeroForce(){}

};

struct DampingForce: public ZeroForce {

    /** constructors */
    DampingForce()
        : damp_const_(0) {}
    DampingForce(double c)
        : damp_const_(c) {}

    virtual Point operator()(Node n, double t) {

        (void) t;
        Point f_damp = n.value().vel;
        return f_damp *= -damp_const_;
    }

    virtual ~DampingForce(){}

    private:
    double damp_const_;

};

struct MassSpringForce: public ZeroForce {

    MassSpringForce(){}

    virtual Point operator()(Node n, double t) {

        (void) t;
        Point f_spring = Point(0,0,0);
        for(auto i = n.edge_begin(); i!=n.edge_end(); ++i){
            auto e = *i;
            Point diff_pos = n.position() - e.node2().position();
            double norm_diff = norm(n.position() - e.node2().position());
            double K = e.value().K;
            double L = e.value().L;
            diff_pos *= -K*(norm_diff - L)/norm_diff;
            f_spring += diff_pos;
        }
        return f_spring;
    }

    virtual ~MassSpringForce(){}
};


struct GravityForce: public ZeroForce {

    GravityForce(){}

    virtual Point operator()(Node n, double t) {
        (void) t;
        return n.value().mass * Point(0,0,-grav);
    }

    virtual ~GravityForce(){}

};


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
struct CombinedForces {

    //constructor
    CombinedForces(std::vector<ZeroForce*> v)
        :input_force_vec_(v){
        };

    //operator
    Point operator()(Node n, double t) {

        //original position
        Point pt = Point(0,0,0);
        //iterate through all forces and update with forces
        for (auto it = input_force_vec_.begin(); it != input_force_vec_.end(); ++it){
            pt += (*(*it))(n,t);
        }
        return pt;
    }

    //private attribute
    private:
    std::vector<ZeroForce*> input_force_vec_;

};//end CombinedForce

template <typename F1, typename F2>
CombinedForces make_combined_force(F1 f1, F2 f2){
    std::vector<ZeroForce*> v = {&f1, &f2};
    return CombinedForces(v);
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
    (void) n; 

     if(n.position() == Point(0,0,0) || n.position() == Point(1,0,0)){
                return Point(0,0,0);
            }

      MassSpringForce msf_obj = MassSpringForce();
      GravityForce gf_obj = GravityForce();
      Point f_spring = msf_obj(n,t);
      Point f_grav = gf_obj(n,t);

      f_spring += f_grav;

      return f_spring;

  }
};



//--functionality_2
//--No constraints. No remove methods.
//--END

//--design_2
//--No constraints. No remove methods.
//--END

//--style_0
//--Good coding style!
//--END

//--documentation_0
//--good
//--END


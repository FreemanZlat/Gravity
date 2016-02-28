#pragma once

#include <SFML/Graphics.hpp>
#include <glm/glm.hpp>

#include <functional>
#include <vector>

class Gravity
{
 public:
    enum Integration
    {
        INTEGRATION_EULER1 = 0,
        INTEGRATION_EULER2,
        INTEGRATION_RUNGE_KUTTA4
    };

    enum Optimisation
    {
        OPTIMISATION_NONE = 0,
        OPTIMISATION_TREE
    };

    Gravity(Integration integration, Optimisation optimisation);
    ~Gravity();

    typedef glm::dvec2 vec;

    void Add(int count);
    void DoGravity(double dt);
    void DoCollisions();
    void DoClean();
    void Draw(sf::RenderWindow &window);
    void GetCenter(float &x, float &y);

    int GetCount();

 private:
    struct Info
    {
        double mass;
        double radius;
    };

    struct Particles
    {
        std::vector<vec> pos;
        std::vector<vec> vel;
        std::vector<Info> info;

        std::vector<sf::CircleShape> circles;

        int count;
        vec center;
        double mass_all;

        void add(vec pos, vec vel, double mass);
        void del(int idx);
    };

    struct TreeNode
    {
        int particle_idx;

        vec bbox;
        vec bbox_size;

        vec mass_center;
        double mass;

        TreeNode* nodes[4];
    };

    struct NodesPool
    {
        std::vector<TreeNode*> pool;
        int idx;

        TreeNode* get_node(int particle_idx);
        void reset();
    };

    Integration integration;
    Optimisation optimisation;

    Particles particles;

    TreeNode *root;
    NodesPool nodes_pool;

    void calc_gravity(std::vector<vec> &pos, std::vector<vec> &accel);
    void calc_gravity_all(int particle_idx, std::vector<vec> &pos, std::vector<vec> &accel);
    void calc_gravity_tree(TreeNode *node, int particle_idx, std::vector<vec> &pos, std::vector<vec> &accel);

    void calc_gravity_euler1(double dt);
    void calc_gravity_euler2(double dt);
    void calc_gravity_runge_kutta4(double dt);

    void calc_tree();
    void calc_bbox();
    void add_node1(TreeNode *node, int partice_idx);
    void add_node2(TreeNode *node, int partice_idx);
    void fill_tree(TreeNode *node);

    static double Random(double q);
};

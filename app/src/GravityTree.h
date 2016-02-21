#pragma once

#include <SFML/Graphics.hpp>

#include <vector>

class GravityTree {
 public:
    GravityTree();
    ~GravityTree();

    void Add(int count);
    void Calc(double dt);
    void Draw(sf::RenderWindow &window);
    void GetCenter(float &x, float &y);

    int GetCount();

 private:
    struct Particle {
        double x, y;
        double vx, vy;
        double mass;
        double radius;
        double ax, ay;
        bool killed;
    };

    struct TreeNode {
        Particle *particle;

        double bbox_x, bbox_y;
        double bbox_size;

        double mx, my;
        double mass;

        TreeNode* nodes[4];
    };

    int count;
    std::vector<Particle> particles;

    double bbox_x, bbox_y;
    double bbox_size;

    TreeNode *root;
    std::vector<TreeNode*> free_nodes;
    int node_idx;

    void BuildTree();
    void FillTree(TreeNode *node);
    void CalcTree(double dt);

    TreeNode* GetNode(Particle *particle);
    void TreeAdd1(TreeNode *node, Particle *particle);
    void TreeAdd2(TreeNode *node, Particle *particle);
    void CalcForces(TreeNode *node, Particle *particle);
};

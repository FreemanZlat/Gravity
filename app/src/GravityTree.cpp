#include "GravityTree.h"

#include <cmath>
#include <cstdlib>
#include <ctime>

static const double START_X = 640.0;
static const double START_Y = 360.0;
static const double START_SIZE = 400.0;

static const double GRAVITY_G = 6.673848e-11;
static const double GRAVITY_E = 0.1;
static const double GRAVITY_MASS = 100000000.0;
static const double GRAVITY_MASS_BIG = 2000000000000.0;

static const double MAX_DISTANCE = 4096;

static const double THETA = 0.5;

static double Random(double q) {
    return q * static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) + 1.0);
}

GravityTree::GravityTree()
        : count(1),
          bbox_x(START_X),
          bbox_y(START_Y),
          bbox_size(START_SIZE),
          root(nullptr),
          node_idx(0) {
    srand((unsigned) time(nullptr));

    Particle p;
    p.x = START_X;
    p.y = START_Y;
    p.vx = 0.0;
    p.vy = 0.0;
    p.mass = GRAVITY_MASS_BIG;
    p.radius = 0.007 * sqrt(sqrt(p.mass));
    p.killed = false;
    this->particles.push_back(p);
}

GravityTree::~GravityTree() {
}

void GravityTree::Add(int count) {
    for (int i = 0; i < count; ++i) {
        Particle p;

        double x = Random(2.0 * START_SIZE) - START_SIZE;
        double y = Random(2.0 * START_SIZE) - START_SIZE;
        double r = sqrt(x * x + y * y);

        p.x = START_X + x;
        p.y = START_Y + y;

        p.vx = y * 0.02 / sqrt(r);
        p.vy = -x * 0.02 / sqrt(r);

        p.mass = GRAVITY_MASS;
        p.radius = 0.007 * sqrt(sqrt(p.mass));
        p.killed = false;

        particles.push_back(p);
    }

    this->count = this->particles.size();
}

void GravityTree::Calc(double dt) {
    double min_x = this->bbox_x;
    double max_x = this->bbox_x;
    double min_y = this->bbox_y;
    double max_y = this->bbox_y;

    for (auto &particle : this->particles) {
        particle.ax = 0.0;
        particle.ay = 0.0;

        if (particle.x < min_x)
            min_x = particle.x;
        if (particle.x > max_x)
            max_x = particle.x;
        if (particle.y < min_y)
            min_y = particle.y;
        if (particle.y > max_y)
            max_y = particle.y;
    }

    double size_x = max_x - min_x;
    double size_y = max_y - min_y;
    this->bbox_size = (size_x > size_y) ? size_x : size_y;
    this->bbox_x = 0.5 * (min_x + max_x);
    this->bbox_y = 0.5 * (min_y + max_y);

    this->BuildTree();
    this->FillTree(this->root);
    this->CalcTree(dt);

    for (int i = 0; i < this->count; ++i) {
        double dx = this->particles[i].x - this->root->mx;
        double dy = this->particles[i].y - this->root->my;

        if (this->particles[i].killed || dx * dx + dy * dy > MAX_DISTANCE * MAX_DISTANCE) {
            this->count--;
            this->particles[i] = this->particles[this->count];
            this->particles.erase(this->particles.begin() + this->count);
            i--;
        }
    }
}

void GravityTree::BuildTree() {
    this->node_idx = 0;
    this->root = this->GetNode(&this->particles[0]);
    this->root->bbox_x = this->bbox_x;
    this->root->bbox_y = this->bbox_y;
    this->root->bbox_size = this->bbox_size;
    for (int i = 1; i < this->count; ++i)
        this->TreeAdd1(this->root, &this->particles[i]);
}

void GravityTree::FillTree(TreeNode *node) {
    if (node->particle != nullptr) {
        node->mx = node->particle->x;
        node->my = node->particle->y;
        node->mass = node->particle->mass;
        return;
    }

    node->mx = 0.0;
    node->my = 0.0;
    node->mass = 0.0;

    for (int i = 0; i < 4; ++i) {
        if (node->nodes[i] == nullptr)
            continue;
        this->FillTree(node->nodes[i]);
        node->mx += node->nodes[i]->mx * node->nodes[i]->mass;
        node->my += node->nodes[i]->my * node->nodes[i]->mass;
        node->mass += node->nodes[i]->mass;
    }

    node->mx /= node->mass;
    node->my /= node->mass;
}

void GravityTree::CalcTree(double dt) {
    for (auto &particle : this->particles) {
        if (particle.killed)
            continue;

        this->CalcForces(this->root, &particle);

        particle.vx += GRAVITY_G * particle.ax * dt;
        particle.vy += GRAVITY_G * particle.ay * dt;

        particle.x += particle.vx * dt;
        particle.y += particle.vy * dt;
    }
}

GravityTree::TreeNode* GravityTree::GetNode(Particle * particle) {
    if (this->node_idx == this->free_nodes.size()) {
        TreeNode *new_node = new TreeNode;
        this->free_nodes.push_back(new_node);
    }
    TreeNode *node = this->free_nodes[this->node_idx++];
    node->particle = particle;
    for (int i = 0; i < 4; ++i)
        node->nodes[i] = nullptr;
    return node;
}

void GravityTree::TreeAdd1(TreeNode *node, Particle *particle) {
    if (node->particle != nullptr) {
        this->TreeAdd2(node, node->particle);
        node->particle = nullptr;
    }
    this->TreeAdd2(node, particle);
}

void GravityTree::TreeAdd2(TreeNode *node, Particle *particle) {
    int idx = 0;
    if (particle->x > node->bbox_x)
        idx += 1;
    if (particle->y > node->bbox_y)
        idx += 2;
    if (node->nodes[idx] != nullptr) {
        this->TreeAdd1(node->nodes[idx], particle);
        return;
    }
    node->nodes[idx] = this->GetNode(particle);
    double new_size = 0.5 * node->bbox_size;
    node->nodes[idx]->bbox_x = node->bbox_x + ((particle->x > node->bbox_x) ? new_size : -new_size);
    node->nodes[idx]->bbox_y = node->bbox_y + ((particle->y > node->bbox_y) ? new_size : -new_size);
    node->nodes[idx]->bbox_size = new_size;
}

void GravityTree::CalcForces(TreeNode *node, Particle *particle) {
    if (node->particle == particle)
        return;

    double dx = node->mx - particle->x;
    double dy = node->my - particle->y;
    double dist = sqrt(dx * dx + dy * dy);

    if (node->particle == nullptr
            && (dist < 2.0 * particle->radius || node->bbox_size / (dist - particle->radius) > THETA)) {
        for (int i = 0; i < 4; ++i)
            if (node->nodes[i] != nullptr)
                this->CalcForces(node->nodes[i], particle);
        return;
    }

    if (node->particle != nullptr) {
        if (node->particle->killed)
            return;
        if (dist <= particle->radius + node->particle->radius) {
            particle->x = particle->x * particle->mass + node->particle->x * node->particle->mass;
            particle->y = particle->y * particle->mass + node->particle->y * node->particle->mass;
            particle->vx = particle->vx * particle->mass + node->particle->vx * node->particle->mass;
            particle->vy = particle->vy * particle->mass + node->particle->vy * node->particle->mass;

            particle->mass += node->particle->mass;
            particle->radius = 0.007 * sqrt(sqrt(particle->mass));

            particle->x /= particle->mass;
            particle->y /= particle->mass;
            particle->vx /= particle->mass;
            particle->vy /= particle->mass;

            node->particle->killed = true;
            return;
        }
    }
    double k = 1.0 / (dist * dist * dist + GRAVITY_E);

    particle->ax += node->mass * k * dx;
    particle->ay += node->mass * k * dy;
}

void GravityTree::Draw(sf::RenderWindow &window) {
    sf::CircleShape circle(1, 16);
    for (auto &p : this->particles) {
        circle.setPosition(p.x - p.radius, p.y - p.radius);
        circle.setRadius(p.radius);
        window.draw(circle);
    }
}

void GravityTree::GetCenter(float &x, float &y) {
    x = this->root->mx;
    y = this->root->my;
}

int GravityTree::GetCount() {
    return this->count;
}

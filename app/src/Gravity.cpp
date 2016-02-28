#include "Gravity.h"

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

Gravity::Gravity(Integration integration, Optimisation optimisation) :
        integration(integration),
        optimisation(optimisation),
        root(nullptr)
{
    srand((unsigned) time(nullptr));

    this->particles.count = 0;
    this->particles.center = vec(START_X, START_Y);
    this->particles.mass_all = 0.0;

    this->particles.add(vec(START_X, START_Y), vec(0.0, 0.0), GRAVITY_MASS_BIG);
}

Gravity::~Gravity()
{
    for (auto &n : this->nodes_pool.pool)
        delete n;
}

void Gravity::Add(int count)
{
    for (int i = 0; i < count; ++i)
    {
        double x = Random(2.0 * START_SIZE) - START_SIZE;
        double y = Random(2.0 * START_SIZE) - START_SIZE;
        double r = sqrt(x * x + y * y);

        this->particles.add(vec(START_X + x, START_Y + y), vec(y * 0.02 / sqrt(r), -x * 0.02 / sqrt(r)), GRAVITY_MASS);
    }
}

void Gravity::DoGravity(double dt)
{
    if (this->integration == INTEGRATION_EULER1)
        this->calc_gravity_euler1(dt);
    else if (this->integration == INTEGRATION_EULER2)
        this->calc_gravity_euler2(dt);
    else if (this->integration == INTEGRATION_RUNGE_KUTTA4)
        this->calc_gravity_runge_kutta4(dt);
}

void Gravity::DoCollisions()
{
    for (int i = 0; i < this->particles.count; ++i)
        for (int j = i + 1; j < this->particles.count; ++j)
        {
            vec dist = this->particles.pos[j] - this->particles.pos[i];
            double min_r = this->particles.info[i].radius + this->particles.info[j].radius;

            if (glm::dot(dist, dist) > min_r * min_r)
                continue;

            this->particles.pos[i] = this->particles.pos[i] * this->particles.info[i].mass
                    + this->particles.pos[j] * this->particles.info[j].mass;
            this->particles.vel[i] = this->particles.vel[i] * this->particles.info[i].mass
                    + this->particles.vel[j] * this->particles.info[j].mass;

            this->particles.info[i].mass += this->particles.info[j].mass;
            this->particles.info[i].radius = 0.007 * sqrt(sqrt(this->particles.info[i].mass));
            this->particles.circles[i].setRadius(this->particles.info[i].radius);

            this->particles.pos[i] /= this->particles.info[i].mass;
            this->particles.vel[i] /= this->particles.info[i].mass;

            this->particles.del(j);
            j--;
        }
}

void Gravity::DoClean()
{
    this->particles.center = vec(0.0, 0.0);
    this->particles.mass_all = 0.0;

    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.center += this->particles.pos[i] * this->particles.info[i].mass;
        this->particles.mass_all += this->particles.info[i].mass;
    }

    this->particles.center /= this->particles.mass_all;

    for (int i = 0; i < this->particles.count; ++i)
    {
        vec dist = this->particles.pos[i] - this->particles.center;
        if (glm::dot(dist, dist) > MAX_DISTANCE * MAX_DISTANCE)
        {
            this->particles.del(i);
            i--;
        }
    }
}

void Gravity::Draw(sf::RenderWindow &window)
{
    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.circles[i].setPosition(this->particles.pos[i].x - this->particles.info[i].radius,
                                               this->particles.pos[i].y - this->particles.info[i].radius);
        window.draw(this->particles.circles[i]);
    }
}

void Gravity::GetCenter(float &x, float &y)
{
    x = this->particles.center.x;
    y = this->particles.center.y;
}

int Gravity::GetCount()
{
    return this->particles.count;
}

void Gravity::calc_gravity(std::vector<vec> &pos, std::vector<vec> &accel)
{
    if (this->optimisation == OPTIMISATION_NONE)
    {
        for (int i = 0; i < this->particles.count; ++i)
            this->calc_gravity_all(i, pos, accel);
    }
    else if (this->optimisation == OPTIMISATION_TREE)
    {
        this->calc_tree();
        for (int i = 0; i < this->particles.count; ++i)
            this->calc_gravity_tree(this->root, i, pos, accel);
    }
}

void Gravity::calc_gravity_all(int particle_idx, std::vector<vec> &pos, std::vector<vec> &accel)
{
    for (int j = particle_idx + 1; j < this->particles.count; ++j)
    {
        vec dist = pos[j] - pos[particle_idx];
        double r2 = glm::dot(dist, dist);
        double k = 1.0 / (r2 * sqrt(r2) + GRAVITY_E);

        accel[particle_idx] += this->particles.info[j].mass * k * dist;
        accel[j] -= this->particles.info[particle_idx].mass * k * dist;
    }
    accel[particle_idx] *= GRAVITY_G;
}

void Gravity::calc_gravity_tree(TreeNode *node, int particle_idx, std::vector<vec> &pos, std::vector<vec> &accel)
{
    if (node->particle_idx == particle_idx)
        return;

    vec dist = node->mass_center - pos[particle_idx];
    double r = glm::length(dist);

    int bbox_size = node->bbox_size.x;
    if (bbox_size < node->bbox_size.y)
        bbox_size = node->bbox_size.y;

    if (node->particle_idx < 0
            && (r < 2.0 * this->particles.info[particle_idx].radius
                    || bbox_size / (r - this->particles.info[particle_idx].radius) > THETA))
    {
        for (int i = 0; i < 4; ++i)
            if (node->nodes[i] != nullptr)
                this->calc_gravity_tree(node->nodes[i], particle_idx, pos, accel);
        return;
    }

    double k = 1.0 / (r * r * r + GRAVITY_E);

    accel[particle_idx] += GRAVITY_G * node->mass * k * dist;
}

void Gravity::calc_gravity_euler1(double dt)
{
    std::vector<vec> accel(this->particles.count, vec(0.0, 0.0));

    this->calc_gravity(this->particles.pos, accel);

    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.vel[i] += accel[i] * dt;
        this->particles.pos[i] += this->particles.vel[i] * dt;
    }
}

void Gravity::calc_gravity_euler2(double dt)
{
    std::vector<vec> accel0(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(this->particles.pos, accel0);

    std::vector<vec> pos1 = this->particles.pos;
    std::vector<vec> vel1 = this->particles.vel;
    for (int i = 0; i < this->particles.count; ++i)
    {
        pos1[i] += dt * this->particles.vel[i];
        vel1[i] += dt * accel0[i];
    }

    std::vector<vec> accel1(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(pos1, accel1);

    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.pos[i] += dt * 0.5 * (vel1[i] + this->particles.vel[i]);
        this->particles.vel[i] += dt * 0.5 * (accel0[i] + accel1[i]);
    }
}

void Gravity::calc_gravity_runge_kutta4(double dt)
{
    std::vector<vec> accel(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(this->particles.pos, accel);

    std::vector<vec> pos(this->particles.count);
    std::vector<vec> vel(this->particles.count);

    std::vector<vec> p1(this->particles.count);
    std::vector<vec> v1(this->particles.count);
    for (int i = 0; i < this->particles.count; ++i)
    {
        p1[i] = dt * this->particles.vel[i];
        v1[i] = dt * accel[i];
        pos[i] = this->particles.pos[i] + 0.5 * p1[i];
        vel[i] = this->particles.vel[i] + 0.5 * v1[i];
    }

    accel.resize(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(pos, accel);

    std::vector<vec> p2(this->particles.count);
    std::vector<vec> v2(this->particles.count);
    for (int i = 0; i < this->particles.count; ++i)
    {
        p2[i] = dt * vel[i];
        v2[i] = dt * accel[i];
        pos[i] = this->particles.pos[i] + 0.5 * p2[i];
        vel[i] = this->particles.vel[i] + 0.5 * v2[i];
    }

    accel.resize(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(pos, accel);

    std::vector<vec> p3(this->particles.count);
    std::vector<vec> v3(this->particles.count);
    for (int i = 0; i < this->particles.count; ++i)
    {
        p3[i] = dt * vel[i];
        v3[i] = dt * accel[i];
        pos[i] = this->particles.pos[i] + p3[i];
        vel[i] = this->particles.vel[i] + v3[i];
    }

    accel.resize(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(pos, accel);

    std::vector<vec> p4(this->particles.count);
    std::vector<vec> v4(this->particles.count);
    for (int i = 0; i < this->particles.count; ++i)
    {
        p4[i] = dt * vel[i];
        v4[i] = dt * accel[i];
    }

    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.pos[i] += (p1[i] + 2.0 * p2[i] + 2.0 * p3[i] + p4[i]) / 6.0;
        this->particles.vel[i] += (v1[i] + 2.0 * v2[i] + 2.0 * v3[i] + v4[i]) / 6.0;
    }
}

void Gravity::calc_tree()
{
    this->nodes_pool.reset();
    this->root = this->nodes_pool.get_node(0);

    this->calc_bbox();

    for (int i = 1; i < this->particles.count; ++i)
        this->add_node1(this->root, i);

    this->fill_tree(this->root);
}

void Gravity::calc_bbox()
{
    vec min = this->particles.pos[0];
    vec max = min;

    for (int i = 1; i < this->particles.count; ++i)
    {
        if (this->particles.pos[i].x < min.x)
            min.x = this->particles.pos[i].x;
        if (this->particles.pos[i].x > max.x)
            max.x = this->particles.pos[i].x;

        if (this->particles.pos[i].y < min.y)
            min.y = this->particles.pos[i].y;
        if (this->particles.pos[i].y > max.y)
            max.y = this->particles.pos[i].y;
    }

    this->root->bbox_size = 0.5 * (max - min);
    this->root->bbox = min + this->root->bbox_size;
}

void Gravity::add_node1(TreeNode *node, int partice_idx)
{
    if (node->particle_idx >= 0)
    {
        this->add_node2(node, node->particle_idx);
        node->particle_idx = -1;
    }
    this->add_node2(node, partice_idx);
}

void Gravity::add_node2(TreeNode *node, int partice_idx)
{
    int idx = 0;
    if (this->particles.pos[partice_idx].x > node->bbox.x)
        idx += 1;
    if (this->particles.pos[partice_idx].y > node->bbox.y)
        idx += 2;
    if (node->nodes[idx] != nullptr)
    {
        this->add_node1(node->nodes[idx], partice_idx);
        return;
    }
    node->nodes[idx] = this->nodes_pool.get_node(partice_idx);
    node->nodes[idx]->bbox_size = 0.5 * node->bbox_size;
    node->nodes[idx]->bbox.x = node->bbox.x
            + ((this->particles.pos[partice_idx].x > node->bbox.x) ?
                    node->nodes[idx]->bbox_size.x : -node->nodes[idx]->bbox_size.x);
    node->nodes[idx]->bbox.y = node->bbox.y
            + ((this->particles.pos[partice_idx].y > node->bbox.y) ?
                    node->nodes[idx]->bbox_size.y : -node->nodes[idx]->bbox_size.y);
}

void Gravity::fill_tree(TreeNode *node)
{
    if (node->particle_idx >= 0)
    {
        node->mass_center = this->particles.pos[node->particle_idx];
        node->mass = this->particles.info[node->particle_idx].mass;
        return;
    }

    node->mass_center = vec(0.0, 0.0);
    node->mass = 0.0;

    for (int i = 0; i < 4; ++i)
    {
        if (node->nodes[i] == nullptr)
            continue;
        this->fill_tree(node->nodes[i]);
        node->mass_center += node->nodes[i]->mass_center * node->nodes[i]->mass;
        node->mass += node->nodes[i]->mass;
    }

    node->mass_center /= node->mass;
}

void Gravity::Particles::add(vec pos, vec vel, double mass)
{
    this->count++;
    this->pos.push_back(pos);
    this->vel.push_back(vel);
    this->info.push_back( { mass, 0.007 * sqrt(sqrt(mass)) });
    this->circles.push_back(sf::CircleShape(1, 16));
    this->circles.back().setRadius(this->info.back().radius);
}

void Gravity::Particles::del(int idx)
{
    this->count--;

    this->pos[idx] = this->pos[this->count];
    this->vel[idx] = this->vel[this->count];
    this->info[idx] = this->info[this->count];
    this->circles[idx] = this->circles[this->count];

    this->pos.pop_back();
    this->vel.pop_back();
    this->info.pop_back();
    this->circles.pop_back();
}

Gravity::TreeNode* Gravity::NodesPool::get_node(int particle_idx)
{
    if (this->idx == this->pool.size())
        this->pool.push_back(new TreeNode);

    TreeNode *node = this->pool[this->idx++];
    node->particle_idx = particle_idx;
    for (int i = 0; i < 4; ++i)
        node->nodes[i] = nullptr;

    return node;
}

void Gravity::NodesPool::reset()
{
    this->idx = 0;
}

double Gravity::Random(double q)
{
    return q * static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) + 1.0);
}

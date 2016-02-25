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

Gravity::Gravity()
{
    srand((unsigned) time(nullptr));

    this->particles.count = 0;
    this->particles.center = vec(START_X, START_Y);
    this->particles.mass_all = 0.0;

    this->particles.add(vec(START_X, START_Y), vec(0.0, 0.0), GRAVITY_MASS_BIG);
}

Gravity::~Gravity()
{
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
//    this->calc_gravity_euler1(dt);
    this->calc_gravity_euler2(dt);
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
    sf::CircleShape circle(1, 16);
    for (int i = 0; i < this->particles.count; ++i)
    {
        circle.setPosition(this->particles.pos[i].x - this->particles.info[i].radius,
                           this->particles.pos[i].y - this->particles.info[i].radius);
        circle.setRadius(this->particles.info[i].radius);
        window.draw(circle);
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
    for (int i = 0; i < this->particles.count; ++i)
    {
        for (int j = i + 1; j < this->particles.count; ++j)
        {
            vec dist = pos[j] - pos[i];
            double r2 = glm::dot(dist, dist);
            double k = 1.0 / (r2 * sqrt(r2) + GRAVITY_E);

            accel[i] += this->particles.info[j].mass * k * dist;
            accel[j] -= this->particles.info[i].mass * k * dist;
        }
        accel[i] *= GRAVITY_G;
    }
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
    std::vector<vec> accel(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(this->particles.pos, accel);

    std::vector<vec> pos1 = this->particles.pos;
    std::vector<vec> vel1 = this->particles.vel;
    for (int i = 0; i < this->particles.count; ++i)
    {
        pos1[i] += dt * this->particles.vel[i];
        vel1[i] += dt * accel[i];
    }

    std::vector<vec> accel1(this->particles.count, vec(0.0, 0.0));
    this->calc_gravity(pos1, accel1);

    for (int i = 0; i < this->particles.count; ++i)
    {
        this->particles.pos[i] += dt * 0.5 * (vel1[i] + this->particles.vel[i]);
        this->particles.vel[i] += dt * 0.5 * (accel[i] + accel1[i]);
    }
}

void Gravity::Particles::add(vec pos, vec vel, double mass)
{
    this->count++;
    this->pos.push_back(pos);
    this->vel.push_back(vel);
    this->info.push_back( { mass, 0.007 * sqrt(sqrt(mass)) });
}

void Gravity::Particles::del(int idx)
{
    this->count--;

    this->pos[idx] = this->pos[this->count];
    this->vel[idx] = this->vel[this->count];
    this->info[idx] = this->info[this->count];

    this->pos.pop_back();
    this->vel.pop_back();
    this->info.pop_back();
}

double Gravity::Random(double q)
{
    return q * static_cast<double>(rand()) / (static_cast<double>(RAND_MAX) + 1.0);
}
